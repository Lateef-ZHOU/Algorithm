#include <chrono>
#include <memory>
#include <type_traits>
#include <vector>
#include <QCoreApplication>
#include "Algorithm/Wind/WindProfile.h"
#include "Algorithm/Sim/LidarForward.h"
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
    NonuniformVerticalWindProfileHistory irregular(2);
    auto column = std::make_shared<const NonuniformVerticalWindProfileUnit>(
        now, std::vector<double>{15, 30, 300},
        std::vector<LinAlg::Vec3>{{1, 2, 3}, {4, 5, 6}, {7, 8, 9}});
    std::weak_ptr<const NonuniformVerticalWindProfileUnit> old = column;
    irregular.append(column);
    if (irregular.at(0).get() != column.get() || irregular.at(0)->heights[2] != 300) return 3;
    column.reset();
    irregular.append(std::make_shared<const NonuniformVerticalWindProfileUnit>(
        now + std::chrono::seconds(1), std::vector<double>{40}, std::vector<LinAlg::Vec3>{{1, 0, 0}}));
    if (irregular.size() != 2 || old.expired()) return 4;
    irregular.append(std::make_shared<const NonuniformVerticalWindProfileUnit>(
        now + std::chrono::seconds(2), std::vector<double>{50}, std::vector<LinAlg::Vec3>{{2, 0, 0}}));
    if (irregular.size() != 2 || !old.expired() || irregular.at(0)->heights[0] != 40
        || irregular.latest()->heights[0] != 50) return 5;
    VerticalWindProfileHistory uniform(2);
    uniform.append(std::make_shared<const VerticalWindProfileUnit>(frame));
    if (uniform.latest()->heightAt(0) != 100) return 6;
    return 0;
}
