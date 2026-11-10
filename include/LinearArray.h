#pragma once
#include <cmath>
#include <cstddef>
#include <utility>
#include <vector>

// 等间隔一维数组：第 i 个元素位于 start + i·step（单位由使用处定，如斜距 m、高度 m）。
// 位置只由 start/step 定，不另存位置数组；元素就是 std::vector 本身。
template <class T>
struct LinearArray : std::vector<T>
{
    double start = 0.0;
    double step = 0.0;

    LinearArray() = default;
    LinearArray(const double start, const double step, std::vector<T> values = {})
        : std::vector<T>(std::move(values)), start(start), step(step) {}

    [[nodiscard]] double positionAt(const std::size_t index) const
    {
        return start + static_cast<double>(index) * step;
    }
    // 轴有效：起点、步长有限且步长为正（不检查元素）。
    [[nodiscard]] bool hasValidAxis() const
    {
        return std::isfinite(start) && std::isfinite(step) && step > 0.0;
    }
    // 同一网格：起点、步长与个数都相同（元素类型可不同）。
    template <class U>
    [[nodiscard]] bool sameAxis(const LinearArray<U> &other) const
    {
        return start == other.start && step == other.step && this->size() == other.size();
    }
};
