#pragma once
#include <cmath>
#include <numbers>
#include <utility>
#include <QVector>
#include <QtGlobal>

#include "Algorithm/Wind/GaussFit.h"

// ============================================================
// RadialWind — 功率谱 → 逐距离门的径向风速 + CNR
//
// 一帧谱是 rows × cols 的行主序矩阵（**row = 距离门、col = 频率 bin**）：每一行是该门上的
// 一条频谱，高斯拟合出峰位 x0（bin）即多普勒频移，线性换算成径向风速；峰的能量与基线之比
// 即 CNR。
//   信号功率 = sqrt(2π)·|A·w|；噪声 = |d|·cols；CNR = 10·log10(signal/noise)。
//
// **为什么住 Algorithm/**：它一个型号字节都不碰——谱的行列数与频率 → 风速的斜率截距全部
// 由调用方传入，与哪张卡、几个通道无关。两张卡的 worker（HfnlWorker / LidarQWorker）用的是
// 同一份，故与 GaussFit / WindProfile 并列住算法层，而不是在每个型号夹里各放一份数值代码。
//
// 谱统一是 double 列表；算法不依赖图形库，也不把小数高度转回计数。
// ============================================================
namespace RadialWind
{
// 一帧的两列产物（裸值对：QVector 是 COW 值类，不套指针壳）。
struct Result
{
    QVector<double> ws;  ///< 每门径向风速（长度 = rows）
    QVector<double> cnr; ///< 每门 CNR（dB，长度 = rows）
};

// spec 行主序 rows × cols；slope/intercept 为频率 bin → 风速（m/s）的线性换算。
[[nodiscard]] inline Result
Calc(const QVector<double> &spec,
     const int rows, const int cols, const double slope, const double intercept)
{
    const double *const buf = spec.constData();
    QVector<double> ws(rows);
    QVector<double> cnr(rows);
    double *const wsOut = ws.data();   // 循环前取写指针：detach 检查一次在此，omp 体内零 COW 判
    double *const cnrOut = cnr.data(); // （各卡提取谱时也只在循环前取一次 data）
#pragma omp parallel for
    for (int row = 0; row < rows; ++row)
    {
        const double *rowBuf = buf + row * cols;
        const GaussFit::Params p = GaussFit(rowBuf, cols).params; // 构造即拟合，取参数
        const double signal = std::sqrt(2.0 * std::numbers::pi) * std::abs(p.A * p.w);
        const double noise = std::abs(p.d) * cols;
        wsOut[row] = p.x0 * slope + intercept;
        cnrOut[row] = noise > 1e-30 ? 10.0 * std::log10(signal / noise) : 0.0;
    }
    return {std::move(ws), std::move(cnr)}; // 拥有型交付：move 归对（非 const 局部，move 合法）
}
} // namespace RadialWind
