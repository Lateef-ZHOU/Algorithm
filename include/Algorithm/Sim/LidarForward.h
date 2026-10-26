#pragma once
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <numbers>
#include <random>
#include <vector>

#include "Algorithm/Geometry/GimbalModel.h" // WorldModel（光束的世界朝向）
#include "Algorithm/Wind/VerticalWindProfile.h" // 等高间隔 ENU 风场
#include "Algorithm/Wind/RangeGate.h"       // 门 → 斜距，全仓一处

// ============================================================
// LidarForward —— 虚拟探测器的正向模型：风场 + 光束朝向 → 一帧功率谱
//
// 反演链（GaussFit → RadialWind → WindProfile）怎么读谱，本模型就怎么写谱，两边同一套口径：
//   · 逐门：第 k 门的斜距 R = RangeGate::rangeOf(k, 门距)，离地高度 z = R·sin(高度角)，光束单位矢量
//     d = (sin az·cos el, cos az·cos el, sin el)（ENU，与 WindProfile 设计矩阵那一行同式），
//     径向风 v_r = wind(z)·d，正 = 远离雷达。
//   · 频点：x₀ = (v_r − intercept) / slope，与 RadialWind 的 ws = x₀·slope + intercept 互逆。
//   · 谱形：噪声底 d 上一个高斯峰 A·exp(−½((i−x₀)/w)²)。CNR 按 RadialWind 的定义反着写——它算
//     CNR = 10·lg(√(2π)·|A·w| / (|d|·cols))，于是给定目标 CNR 就能解出 A。
//   · 累积 = 逐帧求和（用户定）：N 帧之和 = N × 单帧均值，每个频点的相对起伏按 1/√N 给（高斯乘性）；量化成整数，超过
//     计数上限顶满（16 位累加器饱和，真卡同理：累积越久近场越先顶满）。CNR 是比值，不随 N 变。
// CNR 随距离：1 km 处的值，减几何衰减 20·lg(R / 1 km)，减双程消光，加近场重叠函数 10·lg(1 − e^{−(R/R_ov)²})，
// 减环境折减。重叠函数是望远镜收发视场在近场逐渐重合的标准写法：没有它，1/R² 从几十米到几十公里跨 50 多 dB，
// 固定噪声底之上的谱峰在近场必定顶到计数上限（探针实测顶满）。口径、激光功率这些探测器属性全折进
// "1 km 处的 CNR"这一个数——一周交付，先把链路跑通（用户定）。
//
// **环境折减是声明了假设的模型，不是实测**（用户定）：温度在舒适区之外每度折减若干 dB，数值由页面给，
// 文档注明来源是经验或文献。它体现为 CNR 下降，远端门先失效。
//
// 只依赖算法层与标准库，不认识任何设备；输出行主序 rows × cols 的计数，与 RadialWind 读谱的
// buf + row·cols 一致。
// ============================================================
namespace LidarForward
{
    // 探测器属性：形状与标定来自要模仿的那张卡，灵敏度与噪声来自页面。
    struct Detector
    {
        int rows;                 // 距离门数
        int cols;                 // 频点数
        double gateSpacingM;      // 门距（m）
        double freq2WsSlope;      // 频点 → m/s
        double freq2WsIntercept;  // 第 0 频点对应的 m/s
        double cnrAt1kmDb;        // 1 km 处的 CNR（dB）
        double extinctionDbPerKm; // 大气消光（双程，dB/km）
        double overlapRangeM;     // 近场重叠函数的特征距离（m）
        double floorCounts;       // 单帧噪声底（计数）；累积 N 帧后的底 = N·d
        double peakWidthBins;     // 谱峰 σ（频点）
        int accumFrames;          // 累积帧数 N：计数 ∝ N（求和），相对噪声 1/√N
        double maxCounts;         // 量化上限（计数）
    };

    // 环境：先只有温度（甲方点名的那一项）；要加别的因素就加字段、加一条折减。
    struct Environment
    {
        double temperatureC;
    };

    // 环境折减（dB）：温度在 [lowC, highC] 之外每度折减的 dB，两侧各一个斜率。
    struct Degradation
    {
        double lowC;
        double highC;
        double dbPerDegBelow;
        double dbPerDegAbove;

        [[nodiscard]] double penaltyDb(const Environment &env) const
        {
            if (env.temperatureC > highC)
                return (env.temperatureC - highC) * dbPerDegAbove;
            if (env.temperatureC < lowC)
                return (lowC - env.temperatureC) * dbPerDegBelow;
            return 0.0;
        }
    };

    // 光束单位矢量（ENU）。
    [[nodiscard]] inline LinAlg::Vec3 beamVector(const WorldModel &beam)
    {
        constexpr double D2R = std::numbers::pi / 180.0;
        const double h = beam.azimuth() * D2R;
        const double v = beam.elevation() * D2R;
        const double cv = std::cos(v);
        return {std::sin(h) * cv, std::cos(h) * cv, std::sin(v)};
    }

    // 第 k 门的真值径向风（m/s）。高度取 R·sin(el)，与风廓线的高度口径一致。
    [[nodiscard]] inline double radialWind(const VerticalWindProfileUnit &field, const WorldModel &beam, const double gateSpacingM,
                                           const int gate)
    {
        const LinAlg::Vec3 d = beamVector(beam);
        const double range = RangeGate::rangeOf(gate, gateSpacingM);
        return LinAlg::Vec3::dot(field.windAt(range * d.z()), d);
    }

    // 该斜距处的 CNR（dB）。
    [[nodiscard]] inline double cnrDb(const Detector &det, const double rangeM, const double penaltyDb)
    {
        const double km = rangeM / 1000.0;
        const double q = rangeM / det.overlapRangeM;
        const double overlap = 1.0 - std::exp(-q * q);
        return det.cnrAt1kmDb - 20.0 * std::log10(std::max(km, 1e-6)) - det.extinctionDbPerKm * km +
               10.0 * std::log10(std::max(overlap, 1e-12)) - penaltyDb;
    }

    // 一帧谱（累积 N 帧之和）：rows × cols，行主序。rng 由调用方持有（固定种子即可复现）。
    [[nodiscard]] inline std::vector<std::uint64_t> synthesize(const Detector &det, const double penaltyDb,
                                                               const VerticalWindProfileUnit &field, const WorldModel &beam,
                                                               std::mt19937 &rng)
    {
        std::vector<std::uint64_t> out(static_cast<std::size_t>(det.rows) * static_cast<std::size_t>(det.cols));
        std::normal_distribution<double> noise(0.0, 1.0);
        const double n = static_cast<double>(std::max(det.accumFrames, 1));
        const double relSigma = 1.0 / std::sqrt(n);
        const double w = det.peakWidthBins;
        const double d = det.floorCounts;
        for (int k = 0; k < det.rows; ++k)
        {
            const double vr = radialWind(field, beam, det.gateSpacingM, k);
            // 频点位置。超出频带的速度钳在带边：真机是混叠，这里先不模拟。
            const double x0 = std::clamp((vr - det.freq2WsIntercept) / det.freq2WsSlope, 0.0, det.cols - 1.0);
            const double cnr =
                std::pow(10.0, cnrDb(det, RangeGate::rangeOf(k, det.gateSpacingM), penaltyDb) / 10.0);
            const double A = cnr * d * det.cols / (std::sqrt(2.0 * std::numbers::pi) * w);
            std::uint64_t *row = out.data() + static_cast<std::size_t>(k) * static_cast<std::size_t>(det.cols);
            for (int i = 0; i < det.cols; ++i)
            {
                const double dx = (i - x0) / w;
                const double mean = n * (d + A * std::exp(-0.5 * dx * dx)); // N 帧之和的均值
                const double y = std::clamp(mean * (1.0 + relSigma * noise(rng)), 0.0, det.maxCounts);
                row[i] = static_cast<std::uint64_t>(std::llround(y));
            }
        }
        return out;
    }
} // namespace LidarForward
