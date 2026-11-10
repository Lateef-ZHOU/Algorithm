#pragma once
#include <cmath>
#include <cstddef>
#include <span>
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

// 行主序二维等间隔数组（长扁平数组，长度 = rows × columns）：第 r 行 c 列位于
// (rowPosition(r), columnPosition(c)) = (rowStart + r·rowStep, columnStart + c·columnStep)。
// 元素为空、只带行列数与两条轴时，可作只描述网格的轴对象。
template <class T>
struct LinearArray2D : std::vector<T>
{
    int rows = 0;
    int columns = 0;
    double rowStart = 0.0;
    double rowStep = 0.0;
    double columnStart = 0.0;
    double columnStep = 0.0;

    LinearArray2D() = default;
    LinearArray2D(const int rows, const int columns, const double rowStart, const double rowStep, const double columnStart,
                  const double columnStep, std::vector<T> values = {})
        : std::vector<T>(std::move(values)), rows(rows), columns(columns), rowStart(rowStart), rowStep(rowStep),
          columnStart(columnStart), columnStep(columnStep) {}

    [[nodiscard]] double rowPosition(const std::size_t row) const { return rowStart + static_cast<double>(row) * rowStep; }
    [[nodiscard]] double columnPosition(const std::size_t column) const
    {
        return columnStart + static_cast<double>(column) * columnStep;
    }
    // 两条轴有效：行列数为正，起点、步长有限且步长为正（不检查元素）。
    [[nodiscard]] bool hasValidAxes() const
    {
        return rows > 0 && columns > 0 && std::isfinite(rowStart) && std::isfinite(rowStep) && rowStep > 0.0
               && std::isfinite(columnStart) && std::isfinite(columnStep) && columnStep > 0.0;
    }
    // 同一网格：行列数与两条轴的起点、步长都相同（元素类型可不同、不比较元素）。
    template <class U>
    [[nodiscard]] bool sameAxes(const LinearArray2D<U> &other) const
    {
        return rows == other.rows && columns == other.columns && rowStart == other.rowStart && rowStep == other.rowStep
               && columnStart == other.columnStart && columnStep == other.columnStep;
    }
    // 只取网格（行列数与两条轴），不带元素。
    [[nodiscard]] LinearArray2D axes() const { return LinearArray2D(rows, columns, rowStart, rowStep, columnStart, columnStep); }
    // 元素恰好填满 rows × columns。
    [[nodiscard]] bool complete() const
    {
        return rows > 0 && columns > 0 && this->size() == static_cast<std::size_t>(rows) * static_cast<std::size_t>(columns);
    }
    [[nodiscard]] std::span<const T> row(const std::size_t index) const
    {
        return std::span<const T>(this->data() + index * static_cast<std::size_t>(columns), static_cast<std::size_t>(columns));
    }
};
