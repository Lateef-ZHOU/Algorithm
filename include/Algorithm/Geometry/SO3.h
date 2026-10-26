#pragma once
#include <cmath>

#include "Algorithm/Geometry/Matrix.h" // LinAlg::Vec3 / Mat3 + 线性代数底座

// ════════════════════════════════════════════════════════════
// 【单位契约】全工程角量：**跨边界一律「角度(°)」，弧度只活在 SO(3) 内核**
//
//   · canonical = 角度：凡跨边界 / 被存储 / 被传递的角量——各模型的关节当前值、规格参数、
//     Joint 作动范围、WorldModel 三轴范围、落盘文件、UI 显示——全部是「度」。为什么选度而非弧度：
//     UI / 文件 / 电机协议三个外边界天生是度，只有最内核的 trig / exp 需要弧度；选度 → 换算
//     只发生在内核门口一处，且内存里存的值人能读。
//   · 弧度 = SO(3) 内核的实现细节，**绝不外逃**：仅 expmap / logRotation 及局部 trig 吃 / 吐弧度，
//     它们在自己首 / 末行换算（`D2R` = deg→rad 进门、`R2D` = rad→deg 出门）。`D2R` / `R2D` /
//     `DEG_TO_RAD` / `numbers::pi` / `qDegreesToRadians` 就是全工程的「弧度前沿」——grep 它们即可看全部弧度触点
//     （RadialWind 的 √2π 是高斯面积，不是角）：
//     GimbalModel.h 的 toOrientation / intrinsicXYZ / PanTiltModel::fromOrientation、WindProfile、
//     SiteCalibration 的残差输出、风羽图画箭头处、罗盘画束位处（BeamCompassItem）。各是一次门口换算，无泄漏。
//   · 铁律：任何函数若需弧度，就地换、就地用、就地换回，**不得把弧度写进关节值或跨模块传递**。
// ════════════════════════════════════════════════════════════

// ============================================================
// SO3 — 朝向对象（继承 LinAlg::Mat3，本体即 SO(3) 群元素 = 旋转矩阵）+ 通用几何工具门面
//
// 线代是通用设施（LinAlg），旋转是其应用，独立成类、脱离 LinAlg。SO3 本体 IS-A 旋转矩阵
// （朝向矩阵即内部状态），把"复合 / 求逆"做成方法；Wrap360 / SO(3) 流形工具
// （expmap / nearestRotation / logRotation）作静态方法。**"作动量↔朝向"的换算不在这里——它依赖装配，
// 交给 GimbalModel**（见 GimbalModel.h，本类只管"朝向本身"，与具体装配无关）。
//
// 声明与定义分开：类内只留声明 + 一行式平凡成员，复杂 inline 定义集中在类外（仍 header-only，单文件）。
// ============================================================
class SO3 : public LinAlg::Mat3
{
public:
    SO3() : LinAlg::Mat3(LinAlg::Mat3::Identity()) {} // 单位朝向

    // 构造只走「保证正交」的路:轴角指数映射 / 把近似阵投影回最近旋转 / 拷贝既有 SO3。**不开放裸
    // Mat3 构造**——随手给个 Mat3 未必正交,而 inverse()=转置、各种复合都依赖"本体恒正交"这条不变量
    // (否则转置≠逆、复合脱离 SO(3))。裸 Mat3→SO3 仅作私有,供已知正交的内部复合包装(见末尾)。
    [[nodiscard]] static SO3 exp(const LinAlg::Vec3 &rotvec) { return SO3(expmap(rotvec)); }       // 轴角(旋量/四元数同源)→旋转
    [[nodiscard]] static SO3 fromApprox(const LinAlg::Mat3 &m) { return SO3(nearestRotation(m)); } // 近似/带噪阵→最近旋转(Higham)；不查奇异 / 反射，调用方先自检（SiteCalibration 以 det 把关）

    [[nodiscard]] SO3 inverse() const { return SO3(transpose()); } // 本体恒正交 → 转置即逆
    // 复合保型（左乘）：旋转不可交换，顺序固定。
    [[nodiscard]] SO3 operator*(const SO3 &o) const { return SO3(LinAlg::Mat3::operator*(o)); }

    // 角度回绕到 [0, 360)（关节进窗、反解收尾用）。
    [[nodiscard]] static double Wrap360(double deg)
    {
        deg = std::fmod(deg, 360.0);
        return deg < 0.0 ? deg + 360.0 : deg;
    }

    // ── SO(3) 流形工具（手眼标定 / 装配模型用；复杂实现见类外）──
    [[nodiscard]] static LinAlg::Mat3 expmap(const LinAlg::Vec3 &rotvec);   // 指数映射（Rodrigues，弧度）
    [[nodiscard]] static LinAlg::Mat3 nearestRotation(const LinAlg::Mat3 &M); // 最近旋转（Higham 极分解）
    [[nodiscard]] static LinAlg::Vec3 logRotation(const LinAlg::Mat3 &R);     // 旋转对数（SO(3)→so(3)，弧度）

private:
    // 裸包装:仅供已知正交的内部路径(exp/fromApprox/operator*/inverse)用，不对外——外部进不来非正交阵。
    explicit SO3(const LinAlg::Mat3 &m) : LinAlg::Mat3(m) {}
};

// ── SO3 inline 实现（类外，header-only）────────────────────────
// 指数映射（Rodrigues，logRotation 的逆）：旋转向量 axis·angle（弧度）→ 旋转矩阵。PoE 正向的积木。
//   R = cosθ·I + sinθ·[û]× + (1−cosθ)·ûûᵀ，θ=‖rotvec‖。近零 → 单位阵。
inline LinAlg::Mat3 SO3::expmap(const LinAlg::Vec3 &rotvec)
{
    const double th = rotvec.norm();
    if (th < 1e-12) return LinAlg::Mat3::Identity();
    const LinAlg::Vec3 u = (1.0 / th) * rotvec;
    const double c = std::cos(th), s = std::sin(th), k = 1.0 - c;
    const double x = u.x(), y = u.y(), z = u.z();
    return {c + k * x * x, k * x * y - s * z, k * x * z + s * y,
            k * y * x + s * z, c + k * y * y, k * y * z - s * x,
            k * z * x - s * y, k * z * y + s * x, c + k * z * z};
}

// 最近旋转（Higham 正交极分解）：把带噪/求和后的矩阵投影回最近正交阵。det 可能 ±1，调用方自检。
inline LinAlg::Mat3 SO3::nearestRotation(const LinAlg::Mat3 &M)
{
    LinAlg::Mat3 R = M;
    for (int it = 0; it < 32; ++it)
    {
        const LinAlg::Mat3 Rinv = R.inverse();
        if (!Rinv)
            break; // M 奇异 → 提前退出，det 自检会判退化
        const LinAlg::Mat3 Rit = Rinv.transpose();
        double diff = 0.0;
        LinAlg::Mat3 next;
        for (int r = 0; r < 3; ++r)
            for (int c = 0; c < 3; ++c)
            {
                next(r, c) = 0.5 * (R(r, c) + Rit(r, c));
                diff += std::abs(next(r, c) - R(r, c));
            }
        R = next;
        if (diff < 1e-12)
            break;
    }
    return R;
}

// 旋转对数（SO(3) → so(3)）：axis·angle 向量（弧度）。角 0 → 零向量；角 = π 时退化为零向量（轴向与模长一起丢），近 π 轴向不可靠。
inline LinAlg::Vec3 SO3::logRotation(const LinAlg::Mat3 &R)
{
    double cosA = (R(0, 0) + R(1, 1) + R(2, 2) - 1.0) * 0.5;
    cosA = cosA < -1.0 ? -1.0 : (cosA > 1.0 ? 1.0 : cosA);
    const double angle = std::acos(cosA);
    const LinAlg::Vec3 w{R(2, 1) - R(1, 2), R(0, 2) - R(2, 0), R(1, 0) - R(0, 1)};
    const double s = w.norm();
    if (s < 1e-12)
        return LinAlg::Vec3::Zeros(); // 近单位 → 零旋转
    return (angle / s) * w;           // = angle · 单位轴
}
