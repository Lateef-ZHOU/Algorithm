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

## 风廓线帧与斜距约定

`VerticalWindProfileUnit`：时刻、起始高度、正层间距、ENU 向量数组；`VerticalWindProfileHistory` 为定长整帧环，满后覆盖最旧帧，只保存只读共享份额。

斜距约定只写在 `RangeGate::rangeOf(gate, start, spacing) = start + gate·spacing`：起点由各卡通道规格显式给出（2026-11-10 用户定：现有各卡都以每门起点为坐标，第 0 门在 0 m）。`WindProfile::Calc(beams, startM, stepM, at)` 复用四束 ENU 求解，高度起点 = 起点 × sin(平均仰角)，层间距 = 门距 × sin(平均仰角)；门距非正或仰角不在水平面之上返回空。`LidarForward::Detector` 同样带门起点 `gateStartM`。

2026-11-10 起删除非均匀风廓线帧 `NonuniformVerticalWindProfileUnit`、其历史环与按逐门斜距的 `Calc(ranges)`：所有采集卡的距离与频率都是等间隔。验证：Windows / MSVC / Qt 6.11.2，build-1110 Release 严格警告零诊断、公开头独立编译、CTest 1/1。

2026-11-01 本地验证：Windows / MSVC / Qt 6.11.2，Release 严格警告构建、公开头独立编译和 CTest 1/1 通过，覆盖两类帧及整帧环的容量内填充、覆盖释放、数组对应关系。证据 `build/nonuniform-build.log` 与 `build/Testing/Temporary/LastTest.log`。本轮未重跑 Linux；未推送。

## 来源与许可

从 [LiDAR_DataACQ](https://github.com/Lateef-ZHOU/LiDAR_DataACQ) 的模块分离提交 `94ad6f8` 提取；初始分仓保持运行逻辑，独立整理构建、验证与发布入口。

自有代码采用 [MIT](LICENSE)，Copyright (c) 2026 周涛 (ZHOU TAO)。第三方依赖不随本仓库改为 MIT，参见 [依赖说明](THIRD_PARTY_NOTICES.md)。

## 源码目录（2026-11-03）

公开文件直接放在 `include/`，引用时使用文件名，不再携带原单仓库的业务目录前缀。CMake 目标名称、C++ 类型和 QML 模块名保持不变；使用方须同步更新头文件路径及固定依赖版本。

目录整理验证（2026-11-03）：Windows Qt 6.11.2 / MSVC：公开头独立编译、CTest 1/1 通过。本轮通过本地源码覆盖验证，固定依赖同步更新到本地提交；尚未推送，不代表已验证远端拉取。
