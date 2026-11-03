#pragma once

// ============================================================
// RangeGate.h —— **距离门下标 → 沿波束的斜距**，全仓只写这一处。
//
// 门下标就是数组下标，**从 0 起**；第 k 门的数据代表斜距 (k+1)·门距 处。
// （用户依物理测量定：下标从 0 起 ⇒ (k+1)·门距；若从 1 起才是 k·门距——本仓库各卡的逐门
//   数组都是前者。）
//
// **为什么单独成叶子头**：这条约定同时被 Algorithm/WindProfile（高度 = 斜距 × sin 高度角）
// 与 UI/Chart/Wind2DItem（剖面图纵轴的刻度文字）用到，而那两层互不依赖。
// ⚠ 别在使用处各写各的：那样迟早写出两个口径（一处 k+1、一处 k），同一个门在两张图上差一个
// 门距，而两边都"看着合理"、谁也发现不了。**换算只有一份，才谈得上统一。**
// ============================================================
namespace RangeGate
{
// 第 gate 门（数组下标，从 0 起）的斜距，米。对 gate 是线性的：rangeOf(k) = rangeOf(0) + k·门距。
[[nodiscard]] constexpr double rangeOf(const int gate, const double spacingMeters) noexcept
{
    return (gate + 1) * spacingMeters;
}
} // namespace RangeGate
