#pragma once
#include <algorithm>
#include <cmath>
#include <numbers>
#include "Algorithm/Wind/VerticalWindProfile.h"

// 虚拟风场的输入工厂。人工输入采用水平风速/气象来向，进入帧时转为 ENU。
// 场为水平均匀、仅随高度变化；插值和覆盖外取端值统一由风廓线帧实现。
namespace WindField
{
    [[nodiscard]] inline LinAlg::Vec3 velocity(const double speed, const double dirDeg, const double vertical)
    {
        const double rad = dirDeg * std::numbers::pi / 180.0;
        return {-speed * std::sin(rad), -speed * std::cos(rad), vertical};
    }

    [[nodiscard]] inline VerticalWindProfileUnit uniform(VerticalWindProfileUnit profile, const double speed,
                                                         const double dirDeg, const double vertical = 0.0)
    {
        std::fill(profile.winds.begin(), profile.winds.end(), velocity(speed, dirDeg, vertical));
        return profile;
    }

    // 两端风由地面参数和每千米增量确定，按目标网格填入 ENU；覆盖外取端值。
    [[nodiscard]] inline VerticalWindProfileUnit linearShear(VerticalWindProfileUnit profile, const double speedAtGround,
        const double speedGainPerKm, const double dirAtGroundDeg, const double dirGainPerKm,
        const double topM, const double vertical = 0.0)
    {
        if (!(topM > 0.0) || !std::isfinite(topM)) {
            std::fill(profile.winds.begin(), profile.winds.end(), LinAlg::Vec3::Empty());
            return profile;
        }
        const double km = topM / 1000.0;
        const auto low = velocity(speedAtGround, dirAtGroundDeg, vertical);
        const auto high = velocity(speedAtGround + speedGainPerKm * km, dirAtGroundDeg + dirGainPerKm * km, vertical);
        for (std::size_t i = 0; i < profile.winds.size(); ++i) {
            const double fraction = std::clamp(profile.heightAt(i) / topM, 0.0, 1.0);
            profile.winds[i] = low * (1.0 - fraction) + high * fraction;
        }
        return profile;
    }
}
