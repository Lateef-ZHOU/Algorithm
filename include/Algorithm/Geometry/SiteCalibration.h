#pragma once
#include <cmath>
#include <numbers>
#include <QList>
#include <QVector>

#include "Algorithm/Geometry/SO3.h"

// ============================================================
// SiteCalibration — 站点标定（手眼标定 AX=XB，分离地基与手机偏差）
//
// 物理模型（全程旋转、roll 全在场）：手机贴云台顶，第 i 点已知云台姿态 gimbal_i、手机实测世界系
// 全朝向 measured_i，满足
//     R_meas_i = R_base · R_gimbal_i · R_offset
//   R_base  ：地基安装姿态（= 零位三元组，待求）
//   R_offset：手机相对云台顶的固定偏差（待求）
//
// 解法（两步，皆为旋转，不假设 pan 等分、不丢 roll）：
//   ① 手眼：相邻对的相对运动 A=R_meas_iᵀR_meas_j、B=R_gimbal_iᵀR_gimbal_j 满足 B=R_off·A·R_offᵀ，
//      即旋转向量 log(B)=R_off·log(A)——问题变成「找一个旋转，把一组向量 a 转到 b」（Wahba 问题）：
//      累加 M=Σ log(B)·log(A)ᵀ，离它最近的旋转即 R_offset。
//      **需 ≥2 个不平行的相对转轴** → 标定时云台除 pan 还须 tilt（俯仰激励）；转轴全共线（只转了
//      水平）时 R_offset 与 R_base 绕 pan 轴的分量不可分（规范自由度），此时**退化为 R_offset=I**、
//      只解 R_base——手机绕水平轴贴歪的那几度会原样算进地基的方位。
//      **补秩**：两个不平行的方向本就足以定出一个旋转，但那时 M 的秩只有 2，求最近旋转的迭代（要对
//      M 求逆）起不了头。「四点、隔点抬俯仰」恰是这种情形：第 1→2 段与第 3→4 段的相对转动相同，
//      独立方向只有两个。故先把 M 加上自己的余因子阵（按 ‖M‖ 归一）：
//          cof(Σ b·aᵀ) = Σ_{i<j} (b_i×b_j)·(a_i×a_j)ᵀ          （Cauchy–Binet）
//      即把每两对观测的**叉积**也当成一对观测补进去——旋转保叉积（b=R·a ⇒ b_i×b_j=R·(a_i×a_j)），
//      第三个方向就有了。写成带符号的奇异值分解 M=U·diag(d₁,d₂,d₃)·Vᵀ（U、V 取正常旋转，d₃ 可负），
//      cof(M)=U·diag(d₂d₃,d₁d₃,d₁d₂)·Vᵀ 与 M 同一组 U、V：加上去只是把三个 d 都抬成正的，最近旋转
//      就是 U·Vᵀ——Wahba 问题在 SO(3) 里的最优解。M 本来满秩且自洽时结果不变；秩 1（转轴全共线）时
//      余因子阵为零、补不出来，照旧退化。
//   ② 回代：R_base = 平均_i( R_meas_i · (R_gimbal_i·R_offset)ᵀ )，求和后投影回 SO(3)。
//
// ⚠ **本类不判「激励够不够」**：俯仰只抬几度时第二个方向淹在读数噪声里，R_offset 绕主轴那一圈由
// 噪声说了算，而残差照样很小。够不够由流程保证（云台族基 GMBDeviceBase::kCalibTiltExcite* 限住激励角）。
// 量级（合成数据：手机贴歪约 6°，读数噪声方位 0.5°、倾角 0.1°，四点隔点抬俯仰）——地基姿态的误差
// 均值：抬 40° 约 0.26°、抬 15° 约 0.45°、抬 5° 约 1.2°；只转水平约 5°（即手机绕水平轴贴歪的量）。
// 方位读数的噪声到 3°（手机罗盘的实际水平）时，抬 40° 四点约 1.3°——瓶颈在罗盘，不在算法。
//
// 样本与结果一律用**朝向 SO3**（旋转矩阵本体），本类自始至终不碰欧拉角、不知装配约定——角 ↔ 矩阵
// 的转换由调用方按云台模型（PanTiltModel / WorldModel 的 toOrientation）先做好再喂进来。
// **measured 必须是全朝向**（方位 / 高度 / 翻滚三个角都要）：算法吃的是相邻两点之间的相对转动，
// 只有指向的两个角拼不出来。样本按运动顺序传入（相邻相对旋转避免接近 180°，否则旋转向量轴向
// 不可靠）。声明 / 定义分开，inline 类外（header-only，单文件）。
// **目前全仓没有调用点**：标定页尚未接线，本类是备好的算法。
// ============================================================
class SiteCalibration
{
public:
    struct Sample
    {
        SO3 gimbal;   ///< 已知云台姿态（朝向）
        SO3 measured; ///< 手机实测世界系朝向
    };

    struct Estimate
    {
        SO3 base;   ///< 地基零位姿态 R_base
        SO3 offset; ///< 手机相对云台顶偏差 R_offset；数据仅单轴时退化为单位阵
        bool offsetResolved = false; ///< true=俯仰激励足够、R_offset 被显式解出；false=退化为 I
        double residualDeg = 0.0;  ///< 拟合 RMS 残差 (°)，越小越可信
        bool ok = false;
    };

    // N<3 / 旋转求解退化 → ok=false。
    [[nodiscard]] static Estimate estimateFromSamples(const QList<Sample> &samples);
};

// ============================================================
// SiteCalibration — inline 实现（类外，header-only）
// ============================================================

inline SiteCalibration::Estimate SiteCalibration::estimateFromSamples(const QList<Sample> &samples)
{
    using namespace LinAlg;
    Estimate res;
    const int N = samples.size();
    if (N < 3) return res;

    QVector<Mat3> Rm(N), Rg(N); // 各点：手机实测旋转 / 云台已知旋转（SO3 切到 Mat3 做线代）
    for (int i = 0; i < N; ++i)
    {
        Rm[i] = samples[i].measured;
        Rg[i] = samples[i].gimbal;
    }

    // ① 手眼解 R_offset：相邻对的相对旋转向量 b = R_off·a，累加 Wahba 阵 M = Σ b·aᵀ。
    Mat3 M = Mat3::Zeros();
    for (int i = 0; i + 1 < N; ++i)
    {
        const Vec3 a = SO3::logRotation(Rm[i].transpose() * Rm[i + 1]);
        const Vec3 b = SO3::logRotation(Rg[i].transpose() * Rg[i + 1]);
        const double bv[3] = {b.x(), b.y(), b.z()};
        const double av[3] = {a.x(), a.y(), a.z()};
        for (int r = 0; r < 3; ++r) // M += b·aᵀ（外积累加）
            for (int c = 0; c < 3; ++c)
                M(r, c) += bv[r] * av[c];
    }
    // 补秩（道理见文件头 ①）：M += cof(M) / ‖M‖。余因子阵的第 i 行 = M 另外两行的叉积（轮换序）。
    {
        const Vec3 r0{M(0, 0), M(0, 1), M(0, 2)}, r1{M(1, 0), M(1, 1), M(1, 2)}, r2{M(2, 0), M(2, 1), M(2, 2)};
        const Vec3 cof[3] = {Vec3::cross(r1, r2), Vec3::cross(r2, r0), Vec3::cross(r0, r1)};
        const double scale = std::sqrt(r0.norm2() + r1.norm2() + r2.norm2()); // ‖M‖（Frobenius 范数）
        if (scale > 0.0)
            for (int r = 0; r < 3; ++r)
            {
                M(r, 0) += cof[r].x() / scale;
                M(r, 1) += cof[r].y() / scale;
                M(r, 2) += cof[r].z() / scale;
            }
    }
    Mat3 Roff = Mat3::Identity();
    const Mat3 cand = SO3::nearestRotation(M);
    if (cand.determinant() > 0.5) // 至少两个不平行的转轴 → 采纳；转轴全共线（秩 1，补不出来）→ 退回 I
    {
        Roff = cand;
        res.offsetResolved = true;
    }

    // ② 回代旋转平均求 R_base：S = Σ R_meas_i·(R_gimbal_i·R_off)ᵀ，投影回 SO(3)。
    Mat3 S = Mat3::Zeros();
    for (int i = 0; i < N; ++i)
    {
        const Mat3 Rbi = Rm[i] * (Rg[i] * Roff).transpose();
        for (int r = 0; r < 3; ++r)
            for (int c = 0; c < 3; ++c)
                S(r, c) += Rbi(r, c);
    }
    const Mat3 Rbase = SO3::nearestRotation(S);
    if (Rbase.determinant() < 0.5) return res; // 求和退化（朝向自相矛盾）→ 失败

    // ③ 残差：predicted = R_base·R_gimbal_i·R_off，与实测的失配旋转取角，RMS。
    double sse = 0.0;
    for (int i = 0; i < N; ++i)
    {
        const Mat3 err = (Rbase * Rg[i] * Roff).transpose() * Rm[i];
        const double ang = SO3::logRotation(err).norm(); // 失配角（弧度）
        sse += ang * ang;
    }

    res.base = SO3::fromApprox(Rbase);   // Rbase/Roff 已是 nearestRotation 输出(正交)，
    res.offset = SO3::fromApprox(Roff);   // fromApprox 再投影一次幂等、安全(不走裸构造)
    res.residualDeg = std::sqrt(sse / N) * (180.0 / std::numbers::pi);
    res.ok = true;
    return res;
}
