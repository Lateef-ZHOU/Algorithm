#include <array>
#include <chrono>
#include <cmath>
#include <memory>
#include <type_traits>
#include <vector>
#include <QCoreApplication>
#include "WindProfile.h"
#include "LidarForward.h"
using RxTime = std::chrono::system_clock::time_point;
int main(int argc, char **argv)
{
    QCoreApplication app(argc, argv);
    const auto now = std::chrono::system_clock::now();
    VerticalWindProfileUnit frame{now, 100.0, 25.0, {{1.0, 2.0, 3.0}}};
    if (!frame.hasValidGrid() || frame.heightAt(0) != 100.0 || frame.at != now) return 1;
    const WorldModel world(35.0, 60.0);
    const WorldPoseHistory history{{now, world}};
    if (history.front().first != now) return 2;
    VerticalWindProfileHistory history2(2);
    auto column = std::make_shared<const VerticalWindProfileUnit>(frame);
    std::weak_ptr<const VerticalWindProfileUnit> old = column;
    history2.append(column);
    if (history2.at(0).get() != column.get()) return 3;
    column.reset();
    history2.append(std::make_shared<const VerticalWindProfileUnit>(now + std::chrono::seconds(1), 40.0, 10.0,
        std::vector<LinAlg::Vec3>{{1, 0, 0}}));
    if (history2.size() != 2 || old.expired()) return 4;
    history2.append(std::make_shared<const VerticalWindProfileUnit>(now + std::chrono::seconds(2), 50.0, 10.0,
        std::vector<LinAlg::Vec3>{{2, 0, 0}}));
    if (history2.size() != 2 || !old.expired() || history2.at(0)->heightAt(0) != 40
        || history2.latest()->heightAt(0) != 50) return 5;
    if (RangeGate::rangeOf(0, 0, 15) != 0 || RangeGate::rangeOf(3, 10, 15) != 55) return 6;
    std::array<WindProfile::Input, 4> beams;
    for (int i = 0; i < 4; ++i)
        beams[i] = {WorldModel(90.0 * i, 30.0), {1.0, 2.0, 3.0}};
    // 30° 仰角：高度 = 斜距 / 2；起点 0、门距 15 m。
    const auto column3 = WindProfile::Calc(beams, 0.0, 15.0, now);
    if (!column3 || column3->at != now || column3->startHeightM != 0 || std::abs(column3->heightAt(2) - 15) > 1e-9
        || std::abs(column3->winds[2].z() - 6) > 1e-9) return 7;
    if (WindProfile::Calc(beams, 0.0, 0.0, now) || WindProfile::Calc(beams, 0.0, -15.0, now)) return 8;
    return 0;
}
