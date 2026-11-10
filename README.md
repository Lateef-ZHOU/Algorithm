# Algorithm

几何换算、径向风与风廓线反演、风场和谱合成的数学工具。Windows / Linux，C++20，CMake ≥ 3.28。

## 边界与依赖

- Qt ≥ 6.8.3；Qt 由使用环境提供，不在仓库中打包。
- 上游：无其他项目模块。已有同名 CMake 目标时复用；否则 FetchContent 拉取固定完整提交。
- 本仓库不依赖 LiDAR_DataACQ 源码、应用目录或设备型号；不包含第三方二进制。
- 只依赖 Qt Core 容器/数值工具和 OpenMP；不依赖设备、QML 或图表。时间使用标准 system_clock::time_point。

## 独立构建和验证

安装对应 Qt 开发包，然后执行（将路径替换为本机 Qt 前缀）：

```sh
cmake -S . -B build -DCMAKE_PREFIX_PATH=/path/to/Qt/6.8.3/platform -DCMAKE_BUILD_TYPE=Release
cmake --build build --config Release --parallel 2
ctest --test-dir build -C Release --output-on-failure
```

Windows 在安装 MSVC 的环境运行以上命令，Qt 选择对应 MSVC 构建；Linux 使用 GCC 与 OpenMP 开发环境。
测试逐个独立编译公开头、拒绝意外加载 GUI/QML，并运行契约检查。GitHub Actions 在 Windows/Linux 使用 Qt 6.11.2 验证；最低版本声明与最新验证版本是不同概念。

## 接入应用

```cmake
include(FetchContent)
FetchContent_Declare(lidar_algorithm
    GIT_REPOSITORY https://github.com/Lateef-ZHOU/Algorithm.git
    GIT_TAG <经过验证的完整提交号>)
FetchContent_MakeAvailable(lidar_algorithm)
target_link_libraries(my_app PRIVATE Lidar::Algorithm)
```

从仓库提交列表选择固定版本；示例占位值必须替换，不使用浮动分支保证可复现。多个模块共用已有上游目标，应用若需覆盖上游版本应先统一声明并创建它们。标准 `FETCHCONTENT_SOURCE_DIR_LIDAR_模块名`（模块名大写无分隔）可用于本地开发，发布验证应取消覆盖并实际拉取远端。

## 线性数组与风廓线帧

`LinearArray<T>`（`LinearArray.h`）继承 `std::vector<T>`，另带 `start`、`step`：第 i 个元素位于 `positionAt(i) = start + i·step`；`hasValidAxis()` 要求两者有限且步长为正，`sameAxis()` 比较起点、步长与个数。凡是传达"等间隔位置 + 逐点数值"的签名都用它表达。

- `VerticalWindProfileUnit{at, LinearArray<Vec3> winds}`：第 i 层高度 = `winds.positionAt(i)`；`VerticalWindProfileHistory` 为定长整帧环。
- `WindProfile::Input{dir, LinearArray<double> ws}`：`ws` 的轴是斜距。`Calc(beams, at)` 要求各束同一斜距网格，高度起点与层距 = 斜距起点与门距 × sin(平均仰角)。
- `RadialWind::Result` 为逐门 `std::vector<double>`，距离轴由调用方的通道规格给出。
- 斜距约定：起点由各卡通道规格显式给出（现有各卡第 0 门在 0 m）；`LidarForward::Detector` 只有门起点与门距两个数，就地按 起点 + k·门距 计算。`RangeGate.h` 已删除。

2026-11-10 验证：Windows / MSVC / Qt 6.11.2，build-1110 Release /W4 /WX 零诊断、公开头独立编译、CTest 1/1。

## 来源与许可

从 [LiDAR_DataACQ](https://github.com/Lateef-ZHOU/LiDAR_DataACQ) 的模块分离提交 `94ad6f8` 提取；初始分仓保持运行逻辑，独立整理构建、验证与发布入口。

自有代码采用 [MIT](LICENSE)，Copyright (c) 2026 周涛 (ZHOU TAO)。第三方依赖不随本仓库改为 MIT，参见 [依赖说明](THIRD_PARTY_NOTICES.md)。

## 源码目录（2026-11-03）

公开文件直接放在 `include/`，引用时使用文件名，不再携带原单仓库的业务目录前缀。CMake 目标名称、C++ 类型和 QML 模块名保持不变；使用方须同步更新头文件路径及固定依赖版本。

目录整理验证（2026-11-03）：Windows Qt 6.11.2 / MSVC：公开头独立编译、CTest 1/1 通过。本轮通过本地源码覆盖验证，固定依赖同步更新到本地提交；尚未推送，不代表已验证远端拉取。
