#pragma once

// ============================================================
// RangeGate.h —— **距离门下标 → 沿波束的斜距**，全仓只写这一处。
//
// 门下标就是数组下标，从 0 起；第 k 门的斜距 = 起点 + k·门距。起点由各卡的通道规格显式给出
// （用户定：现有各卡都以每门起点为坐标，第 0 门在 0 m）。
//
// **为什么单独成叶子头**：这条约定同时被 Algorithm/WindProfile（高度 = 斜距 × sin 高度角）、
// LidarForward（正向模型）与 ChartCore 的剖面图纵轴用到，而那几层互不依赖。
// ⚠ 别在使用处各写各的：那样迟早写出两个口径，同一个门在两张图上差一个门距，而两边都"看着合理"。
// ============================================================
namespace RangeGate
{
// 第 gate 门（数组下标，从 0 起）的斜距，米。
[[nodiscard]] constexpr double rangeOf(const int gate, const double startMeters, const double spacingMeters) noexcept
{
    return startMeters + gate * spacingMeters;
}
} // namespace RangeGate
