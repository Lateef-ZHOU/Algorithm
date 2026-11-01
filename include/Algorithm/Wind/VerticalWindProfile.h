#pragma once
#include <cmath>
#include <cstddef>
#include <memory>
#include <stdexcept>
#include <utility>
#include <vector>
#include <chrono>

#include "Algorithm/Geometry/Matrix.h"

// 一个时刻的等高间隔风廓线。起始高度和正层间距以米为单位，ENU 为东、北、上（米/秒）。
// 数组按值属于本帧；发布后通过 const 份额共享，不另存高度数组、风速或风向。
struct VerticalWindProfileUnit
{
    std::chrono::system_clock::time_point at;
    double startHeightM;
    double heightStepM;
    std::vector<LinAlg::Vec3> winds;

    [[nodiscard]] bool hasValidGrid() const
    {
        return !winds.empty() && std::isfinite(startHeightM) && std::isfinite(heightStepM)
               && heightStepM > 0.0 && std::isfinite(heightAt(winds.size() - 1));
    }

    // 调用者给出有效下标；单层仍保存采样间距，计算不依赖层数，也不做除法。
    [[nodiscard]] double heightAt(const std::size_t index) const
    {
        return startHeightM + static_cast<double>(index) * heightStepM;
    }

    // 按均匀网格直接定位相邻样本，以 ENU 分量线性插值；覆盖外取端值。
    // 空/非法网格或非有限查询返回无效向量；无效邻点由 Vec3 运算传播。
    [[nodiscard]] LinAlg::Vec3 windAt(const double heightM) const
    {
        if (!hasValidGrid() || !std::isfinite(heightM)) return LinAlg::Vec3::Empty();
        if (heightM <= startHeightM || winds.size() == 1) return winds.front();
        if (heightM >= heightAt(winds.size() - 1)) return winds.back();
        const double position = (heightM - startHeightM) / heightStepM;
        const auto low = static_cast<std::size_t>(position);
        if (low >= winds.size() - 1) return winds.back(); // 端点附近可能舍入到末下标。
        // 用同一高度映射识别节点，避免小数步长的反算舍入把精确落点误当作插值。
        if (heightM == heightAt(low)) return winds[low];
        if (heightM == heightAt(low + 1)) return winds[low + 1];
        const double fraction = position - static_cast<double>(low);
        return winds[low] * (1.0 - fraction) + winds[low + 1] * fraction;
    }
};

using VerticalWindProfilePtr = std::shared_ptr<const VerticalWindProfileUnit>;

// 非均匀整帧：高度与风矢量分数组连续保存。同下标对应同一层。
// 消费入口只在接收时检查等长、有限且递增的高度；发布后生产者不得修改任一数组。
struct NonuniformVerticalWindProfileUnit
{
    std::chrono::system_clock::time_point at;
    std::vector<double> heights;
    std::vector<LinAlg::Vec3> winds;
};
using NonuniformVerticalWindProfilePtr = std::shared_ptr<const NonuniformVerticalWindProfileUnit>;

// 定长历史：只保存完整帧的只读份额，满后覆盖最旧槽；同步由所属线程/使用方负责。
template <class Profile>
class WindProfileHistory
{
public:
    explicit WindProfileHistory(std::size_t capacity) : slots_(capacity)
    {
        if (capacity == 0) throw std::invalid_argument("Wind profile history capacity must be positive.");
    }
    [[nodiscard]] std::size_t size() const noexcept { return count_; }
    [[nodiscard]] std::size_t capacity() const noexcept { return slots_.size(); }
    [[nodiscard]] std::size_t slotAt(std::size_t index) const noexcept { return (head_ + index) % slots_.size(); }
    [[nodiscard]] const std::shared_ptr<const Profile> &at(std::size_t index) const { return slots_[slotAt(index)]; }
    [[nodiscard]] std::shared_ptr<const Profile> latest() const { return count_ ? at(count_ - 1) : nullptr; }
    std::size_t append(std::shared_ptr<const Profile> frame)
    {
        if (!frame) throw std::invalid_argument("Wind profile history requires a frame.");
        const auto slot = slotAt(count_);
        slots_[slot] = std::move(frame);
        if (count_ < slots_.size()) ++count_;
        else head_ = (head_ + 1) % slots_.size();
        return slot;
    }

private:
    std::vector<std::shared_ptr<const Profile>> slots_;
    std::size_t head_ = 0;
    std::size_t count_ = 0;
};


using VerticalWindProfileHistory = WindProfileHistory<VerticalWindProfileUnit>;
using NonuniformVerticalWindProfileHistory = WindProfileHistory<NonuniformVerticalWindProfileUnit>;
