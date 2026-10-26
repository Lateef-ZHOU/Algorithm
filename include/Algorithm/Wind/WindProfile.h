#pragma once
#include <cmath>
#include <numbers>
#include <optional>
#include <span>
#include <QVector>

#include "Algorithm/Wind/RangeGate.h"
#include "Algorithm/Wind/VerticalWindProfile.h"
#include "Algorithm/Geometry/GimbalModel.h"

// ============================================================
// WindProfile — 多波束径向风反演风廓线
//
// 一轮 ≥4 个不同方位的波束（各自世界系方位角/高度角 + 一列径向风速），最小二乘解出
// 每个高度门的 ENU 向量 (u,v,w)，与起始高度、层间距一起交出。
//
// Input.ws 裸值自持（QVector 是 COW 值类，不套指针壳）。各 worker 把最近 kBeams 束
// 存在定长环里，Calc 以 span **只看不拷**——每束带着一份 optional<WorldModel>（内含堆向量），
// 先拷成一份 QVector<Input> 再算就是白搬。
// ============================================================
namespace WindProfile
{
    // **一组反演用几束**——全栈唯一事实源。这个数同时被三处使用，三处必须一致，否则风廓线
    // 要么永远凑不齐、要么拿半组去反演：
    //   · 本文件 Calc 的下限（少于它，方位矩阵解不出三个风分量）；
    //   · worker 的风廓线**环长**（HfnlWorker / LidarQWorker）；
    //   · 四轴计划的**束数**（CycleTask / SnakeTask 的固定束数），协程的同步位间隔正建立在两者相等上。
    // 它既是算法的下限，也是系统的选择（正好取下限）。要改成 6 束扫描，只改这一行。
    inline constexpr int kBeams = 4;

    struct Input
    {
        /// 该波束的**世界朝向**（只读方位与高度；翻滚与光束无关）。**空 = 朝向未知**。
        /// **参考系由类型钉住**：只收 `WorldModel`，由编排层经云台族口 `worldOf` 换好再交——
        /// 收一袋裸数的话，喂进来的若是云台角，基座一歪就指错且不响。
        std::optional<WorldModel> dir;
        QVector<double> ws; ///< 该波束径向风速（每门一个；COW 值类，裸值自持）
    };

    // 反演直接生成带到料时刻的只读整帧；层数组属于该帧，失败返回空份额。
    [[nodiscard]] inline VerticalWindProfilePtr
    Calc(const std::span<const Input> beams, const double gateSpacingMeters, const std::chrono::system_clock::time_point at)
    {
        const int N = static_cast<int>(beams.size());
        if (N < kBeams) return {};
        for (const auto &b : beams)
            if (!b.dir) return {}; // 任一束朝向未知（无云台 / 到位斗上还没有实测位形）→ 整组不反演

        constexpr double DEG_TO_RAD = std::numbers::pi / 180.0;

        // 设计矩阵每行 = 波束世界指向的单位矢量（ENU 球面→笛卡尔，方位/高度，本就是风反演的内核）。
        QVector<double> A(N * 3);
        for (int i = 0; i < N; ++i)
        {
            const double h = beams[i].dir->azimuth() * DEG_TO_RAD;
            const double v = beams[i].dir->elevation() * DEG_TO_RAD;
            const double cv = std::cos(v);
            A[i * 3 + 0] = std::sin(h) * cv;
            A[i * 3 + 1] = std::cos(h) * cv;
            A[i * 3 + 2] = std::sin(v);
        }

        LinAlg::Mat3 AtA = LinAlg::Mat3::Zeros();
        for (int r = 0; r < 3; ++r)
            for (int c = 0; c < 3; ++c)
                for (int i = 0; i < N; ++i)
                    AtA(r, c) += A[i * 3 + r] * A[i * 3 + c];

        const LinAlg::Mat3 AtA_inv = AtA.inverse();
        if (!AtA_inv) return {}; // 伪逆 (AᵀA)⁻¹ 奇异 → 几何退化

        QVector<double> Pinv(3 * N, 0.0);
        for (int r = 0; r < 3; ++r)
            for (int i = 0; i < N; ++i)
            {
                double val = 0.0;
                for (int c = 0; c < 3; ++c)
                    val += AtA_inv(r, c) * A[i * 3 + c];
                Pinv[r * N + i] = val;
            }

        const int nGates = beams[0].ws.size(); // ws 是裸值、必在；空列即 0 门
        if (nGates == 0) return {};
        for (int i = 1; i < N; ++i)
        {
            if (beams[i].ws.size() != nGates) return {};
        }

        // 高度 = 斜距 × **sin**(高度角)。高度角自水平面起、抬头为正（WorldModel 的 el），
        // 与上面设计矩阵的竖直分量 sin(v) 是同一个口径。
        // ⚠ 别写成 cos：高度 = 斜距 × cos 是**天顶角**口径的式子，而喂进来的是高度角。写反了不会
        // 报错——60° 高度角下算出的高度只有正确值的 58%（cos60° / sin60°），还随风廓线一起落盘；
        // 只有 45° 时两者碰巧相等。
        double meanEl = 0.0;
        for (int i = 0; i < N; ++i) meanEl += beams[i].dir->elevation() * DEG_TO_RAD;
        meanEl /= N;
        const double sinEl = std::sin(meanEl);
        auto result = std::make_shared<VerticalWindProfileUnit>(at,
            RangeGate::rangeOf(0, gateSpacingMeters) * sinEl,
            gateSpacingMeters * sinEl,
            std::vector<LinAlg::Vec3>(nGates));
        const double *pinvU = Pinv.constData() + 0 * N;
        const double *pinvV = Pinv.constData() + 1 * N;
        const double *pinvW = Pinv.constData() + 2 * N;

#pragma omp parallel for
        for (int k = 0; k < nGates; ++k)
        {
            double uk = 0.0, vk = 0.0, wk = 0.0;
            bool valid = true;
            for (int i = 0; i < N; ++i)
            {
                const double vr = beams[i].ws[k]; // const 元素直读
                if (!std::isfinite(vr)) { valid = false; break; }
                uk += pinvU[i] * vr;
                vk += pinvV[i] * vr;
                wk += pinvW[i] * vr;
            }
            result->winds[k] = valid ? LinAlg::Vec3(uk, vk, wk) : LinAlg::Vec3::Empty();
        }
        return result;
    }
} // namespace WindProfile
