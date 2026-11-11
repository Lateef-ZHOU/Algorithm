#pragma once
#include <cmath>
#include <numbers>
#include <utility>
#include <vector>
#include <QtGlobal>

#include "GaussFit.h"
#include "LinearArray.h"

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
// 谱为二维线性数组（行 = 距离门、列 = 频点）；算法不依赖图形库，也不把小数高度转回计数。
// ============================================================
namespace RadialWind
{
// 相干多普勒测风：径向风速 v = λ·f_D / 2，每 MHz 频移对应 λ/2 × 10⁶ m/s（λ 以米计；1550 nm → 0.775）。
[[nodiscard]] constexpr double metersPerSecondPerMHz(const double wavelengthM) noexcept
{
    return wavelengthM / 2.0 * 1e6;
}

// 一帧的两列产物：轴取功率谱的行轴（斜距）。
struct Result
{
    LinearArray<double> ws;  ///< 每门径向风速（m/s）
    LinearArray<double> cnr; ///< 每门 CNR（dB）
};

// spectrum 行 = 距离门、列 = 频点（列轴单位 MHz）。峰位按列轴换成频率，减零风速参考后乘 λ/2 系数得风速。
// 网格不完整或轴无效时返回空数组。
[[nodiscard]] inline Result
Calc(const LinearArray2D<double> &spectrum, const double metersPerSecondPerMHz, const double zeroWindMHz)
{
    if (!spectrum.complete() || !spectrum.hasValidAxes() || spectrum.columns < 2)
        return {};
    const int rows = spectrum.rows, cols = spectrum.columns;
    const double slope = spectrum.columnStep * metersPerSecondPerMHz;
    const double intercept = (spectrum.columnStart - zeroWindMHz) * metersPerSecondPerMHz;
    Result result{LinearArray<double>(spectrum.rowStart, spectrum.rowStep, std::vector<double>(rows)),
                  LinearArray<double>(spectrum.rowStart, spectrum.rowStep, std::vector<double>(rows))};
    const double *const buf = spectrum.data();
    double *const wsOut = result.ws.data(); // 循环前取写指针，omp 体内只按下标写
    double *const cnrOut = result.cnr.data();
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
    return result;
}
} // namespace RadialWind
