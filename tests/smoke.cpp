#include <array>
#include <chrono>
#include <cmath>
#include <memory>
#include <type_traits>
#include <vector>
#include <QCoreApplication>
#include "RadialWind.h"
#include "WindProfile.h"
#include "LidarForward.h"
using RxTime = std::chrono::system_clock::time_point;
int main(int argc, char **argv)
{
    QCoreApplication app(argc, argv);
    const auto now = std::chrono::system_clock::now();
    using Winds = LinearArray<LinAlg::Vec3>;
    VerticalWindProfileUnit frame{now, Winds(100.0, 25.0, {{1.0, 2.0, 3.0}})};
    if (!frame.hasValidGrid() || frame.winds.positionAt(0) != 100.0 || frame.at != now) return 1;
    const WorldModel world(35.0, 60.0);
    const WorldPoseHistory history{{now, world}};
    if (history.front().first != now) return 2;
    VerticalWindProfileHistory history2(2);
    auto column = std::make_shared<const VerticalWindProfileUnit>(frame);
    std::weak_ptr<const VerticalWindProfileUnit> old = column;
    history2.append(column);
    if (history2.at(0).get() != column.get()) return 3;
    column.reset();
    history2.append(std::make_shared<const VerticalWindProfileUnit>(
        VerticalWindProfileUnit{now + std::chrono::seconds(1), Winds(40.0, 10.0, {{1, 0, 0}})}));
    if (history2.size() != 2 || old.expired()) return 4;
    history2.append(std::make_shared<const VerticalWindProfileUnit>(
        VerticalWindProfileUnit{now + std::chrono::seconds(2), Winds(50.0, 10.0, {{2, 0, 0}})}));
    if (history2.size() != 2 || !old.expired() || history2.at(0)->winds.positionAt(0) != 40
        || history2.latest()->winds.positionAt(0) != 50) return 5;
    const LinearArray<double> axis(10, 15, {0, 0, 0, 0});
    if (axis.positionAt(3) != 55 || !axis.hasValidAxis() || LinearArray<double>(0, 0).hasValidAxis()
        || !axis.sameAxis(LinearArray<int>(10, 15, {1, 2, 3, 4})) || axis.sameAxis(LinearArray<double>(10, 15, {0}))) return 6;
    std::array<WindProfile::Input, 4> beams;
    for (int i = 0; i < 4; ++i)
        beams[i] = {WorldModel(90.0 * i, 30.0), LinearArray<double>(0.0, 15.0, {1.0, 2.0, 3.0})};
    // 30° 仰角：高度 = 斜距 / 2；起点 0、门距 15 m。
    const auto column3 = WindProfile::Calc(beams, now);
    if (!column3 || column3->at != now || column3->winds.start != 0 || std::abs(column3->winds.positionAt(2) - 15) > 1e-9
        || std::abs(column3->winds[2].z() - 6) > 1e-9) return 7;
    // 各束斜距网格不一致或门距非正：不反演。
    beams[3].ws.step = 30.0;
    if (WindProfile::Calc(beams, now)) return 8;
    for (auto &beam : beams) beam.ws.step = 0.0;
    if (WindProfile::Calc(beams, now)) return 9;
    // 二维谱 → 径向风：两门、九个频点，列轴 10 MHz 起、步长 2 MHz；峰在第 4 列（18 MHz）。
    std::vector<double> values;
    for (int gate = 0; gate < 2; ++gate)
        for (int bin = 0; bin < 9; ++bin)
            values.push_back(1.0 + 100.0 * std::exp(-(bin - 4.0) * (bin - 4.0) / 2.0));
    const LinearArray2D<double> spectrum(2, 9, 30.0, 15.0, 10.0, 2.0, values);
    if (!spectrum.complete() || spectrum.columnPosition(8) != 26.0 || spectrum.row(1).size() != 9) return 10;
    const auto calm = RadialWind::Calc(spectrum, 0.775, 18.0);
    if (calm.ws.size() != 2 || calm.ws.start != 30.0 || calm.ws.step != 15.0 || !calm.cnr.sameAxis(calm.ws)
        || std::abs(calm.ws[0]) > 0.01) return 11;
    const auto shifted = RadialWind::Calc(spectrum, 0.775, 16.0);
    if (std::abs(shifted.ws[1] - 2.0 * 0.775) > 0.01) return 12;
    const auto outOfBand = RadialWind::Calc(spectrum, 0.775, 80.0);
    if (RadialWind::inBand(spectrum, 80.0) || !std::isnan(outOfBand.ws[0]) || outOfBand.cnr[0] != calm.cnr[0]) return 13;
    if (!RadialWind::Calc(LinearArray2D<double>(2, 9, 0, 15, 10, 2, {1, 2}), 0.775, 18).ws.empty()) return 14;
    return 0;
}
