#pragma once
#include <algorithm>
#include <array>
#include <cmath>
#include <chrono>
#include <numbers>
#include <utility>
#include <vector>

#include "Algorithm/Geometry/SO3.h" // 朝向对象 + 流形工具；角量单位契约也写在那里

// GimbalModel.h — 云台的**位形模型**一族：Joint（关节）→ GimbalModel（关节链 + 正解）→
// WorldModel{az,el,roll} / PanTiltModel{pan,tilt}（各自的反解）。
// SO3 是"朝向本身"，本文件是"哪种装配、作动到哪儿会得到那个朝向"——两层，上层包含下层。
// 模型实例即一个位形（谈位形不谈指向）；角量一律是度。

// ============================================================
// Joint — 一根关节：带符号单位轴向（IS-A LinAlg::Vec3）+ **窗的偏置** + **当前作动值**（°）
//
// 把"轴向"、"这根轴的角从哪儿数起"、"它此刻转到哪"绑在一起：轴向继续当 Vec3 用（PoE 连乘按基类切片，正向不变）。
// **关节没有范围、只有窗**（用户定）：一根转轴的作动值本来就是整整一圈，窗只定这一圈从哪儿数起（偏置 offset）——方位 / 水平
// 是罗盘方向，窗 [0°, 360°)；高度 / 翻滚 / 俯仰以水平为零，窗 [−180°, 180°)。收进来的数进门时折进窗，不拒、不坍缩
// ——(130°, 120°, 0°) 与 (310°, 60°, 180°) 是同一个朝向的两个位形，都合法；此后读到的就是窗里的数，界面用同一张窗，
// 中间没有换算。设备行程是另一回事，归官邸按协议定、造斗 / 组包时查。
//
// **关节自带当前值**，于是**模型实例就是一个位形**——`WorldModel` 的实例是一个世界朝向，
// `PanTiltModel` 的实例是这台云台的一个位形。参考系由类型说，长度由构造口保证；作动量不另装
// 一袋裸 double 满世界传（那样"这袋数是谁的、长度对不对"无人知晓，传反了编译器也不响）。
// 构造后不可变：位形是值，要换就整台换。
// ============================================================
class Joint : public LinAlg::Vec3
{
public:
    // 轴向 + 窗的偏置 offset（整圈从 offset 数起，窗 = [offset, offset+360°)；它不是限位，只定从哪儿数）+ 作动值
    //（任意实数；已在窗里的原样收下，否则折进窗）。非 explicit：便于 GimbalModel 列表构造 {{x,y,z}, offset, deg}。
    Joint(const LinAlg::Vec3 &axis, const double offset, const double deg)
        : LinAlg::Vec3(axis),
          value_(deg >= offset && deg < offset + 360.0 ? deg : offset + SO3::Wrap360(deg - offset))
    {
    }

    // 本轴此刻的作动值，落在本轴的窗里。
    [[nodiscard]] double value() const { return value_; }

private:
    double value_;
};

// ============================================================
// GimbalModel — 一台云台的**一个位形**（关节列表，每根关节自带当前值）：位形 ↔ 朝向
//
// 定朝向需「轴序 / 内外旋 / 符号」约定——本质是把每根关节的（带符号）单位轴向按复合顺序排成一串。
// 关节列表**就是模型的核心数据**，轴数即表长、不写死；每个 Joint 自带作动范围。指数积(PoE)/旋量：
//     R(θ₁,…,θ_N) = ∏ exp(θ_i[û_i]×)   （内旋=动轴=正序连乘）
//
// 正向 toOrientation 通用且简单，收在基类。**反向（朝向 → 位形）是真逆运动学**，通用解复杂度极高、
// 可能多解 / 无解（存在性看雅可比的秩、唯一性看冗余自由度），故基类不提供，由每个具体模型自写
// 简单闭式的静态口 `fromOrientation`：两个正交轴模型都复用 protected intrinsicXYZ（双轴截前两个分量）；
// 将来非正交轴子类另写 Paden–Kahan。
// **基类不可实例化**（构造 / 析构皆 protected）：一袋任意关节没有已知反解，不算一台模型。
//
// **关节表私有持有，不继承 `std::vector<Joint>`**。"关节列表就是模型的核心数据"听上去该是 IS-A，
// 但 IS-A vector 会把两样东西一并继承下来：
//   · **写入面**：任何拿着模型的人都能 `clear()` / `push_back()` / `resize()`，而清空关节表会让
//     `orthoBasis()` 去读空表的第 0 个元素。关节表是构造期就定下的事实，写入口一概不开。
//   · **`std::vector::operator==`**：⚠ 实测踩过——带 Q_OBJECT 的类一旦经由模型 IS-A vector，Qt 的
//     `QMetaType::fromType<T>` 会探测"此型别可否相等比较"，它顺着基类找到这个声明、认定可行、
//     再实例化，然后在 `std::equal` 里硬错（`Joint` 没有 `operator==`）。典型的 SFINAE 逃逸：
//     探测看得见声明、看不见体，一旦承诺就无从退回。
// 私有继承本可封门，但本仓库不用私有继承；私有持有一举两清。代价是自己开一扇**只读**门
// （axis）。
//
// 「两个模型间换算」（世界 ↔ 云台）**不做成成员**——它是两台模型 + 一个安装朝向的关系，不属于
// 任一模型；就是几行 SO3 复合（`toOrientation` / `inverse` / `fromOrientation`），住在持有安装标定的
// 云台族基（GMBDeviceBase::gimbalOf / worldOf）。
// ============================================================
class GimbalModel
{
public:

    // 第 i 根关节（**只读**）：各模型自己的读口（azimuth() / pan() 等）与几何换算走它。关节没有范围（整圈），
    // 写入口一概不开。
    [[nodiscard]] const Joint &axis(const std::size_t i) const { return joints_[i]; }

    // 两个位形之间的**距离**（°）：逐轴取圆周差，取最大者——量的是"每根关节离目标值多远"，不是两个朝向差多远
    //（后者是 SO(3) 上的测地角；判到位不用它，云台的每根关节都得到位，用户定）。
    //   · **圆周差**：每根关节都是一圈，0° 与 360° 是同一处，直接相减在跨零时会得出近一整圈的假距离——实机：
    //     命令 0°、停在 359.97°，差 359.97 实为 0.03。算法（用户提出）：两角差的绝对值超过半周，就取 360 − |差|。
    //   · **取最大者**：保持"每轴都在容差内"的含义，族级收集器的容差照旧按单轴度数给。
    // 模板只为把"两台同型"钉在编译期（WorldModel 与 PanTiltModel 比不了）；调用面写具体型（PanTiltModel::distance）。
    template <class Model>
    [[nodiscard]] static double distance(const Model &a, const Model &b)
    {
        double worst = 0.0;
        for (std::size_t i = 0; i < a.joints_.size(); ++i)
        {
            const double d = std::fmod(std::abs(a.joints_[i].value() - b.joints_[i].value()), 360.0);
            worst = std::max(worst, d > 180.0 ? 360.0 - d : d);
        }
        return worst;
    }

    // **本位形**的朝向（PoE 正向，通用）：R = ∏ exp(θ_i·û_i)，θ_i 取自各关节自己的当前值。
    // **无参**：位形就在自己身上。两个位形就是两个对象，换算写成
    //     mounting.toOrientation().inverse() * worldTarget.toOrientation()
    [[nodiscard]] SO3 toOrientation() const;

    // ── 反解（朝向 → 一台摆到那个位形的模型）**不在本类声明**，由各具体模型提供一个静态
    //    具名口：`static Self fromOrientation(const SO3 &R)`。**按值返回、恒有解**（用户定）：都走 intrinsicXYZ，
    //    输出落在典范窗（方位整圈、高度 ±90°、翻滚整圈），两极多解收成一个代表（第三轴折进第一轴）；双轴模型
    //    对翻滚无计可施，目标翻滚直接忽略、截前两个分量。
    //
    // **为什么是静态而不是虚**：它要**返自身的具体型**，而协变返回不认值；
    // 做成虚口就只能返基类指针，调用方再 `static_cast` 下行一次——那次转型恒真却写在每个调用点上，
    // 还白搭一次堆分配。而调用点（云台族基的 gimbalOf / worldOf）全都静态知道自己要哪型，没有
    // 一个真走擦除。副作用是好的：本类没有虚函数，位形是被到处按值造的东西，省下的是每个实例
    // 一个 vptr。

protected:
    explicit GimbalModel(std::vector<Joint> joints) : joints_(std::move(joints)) {}
    // **非虚 + protected**：本类一个虚函数都没有，而位形从不经基类指针删除（它们是被按值造、
    // 按值传的东西）。protected 顺带封死 `delete (GimbalModel *)p`——同 Collector / GMBDeviceBase。
    ~GimbalModel() = default;

    // 极点带宽：|cos θ₂| ≤ kPoleCos 按"光指着极点"收（intrinsicXYZ 用；两模型的反解都建在它上面）。
    // **带宽取 ±0.1°，不是"正好 90°"**（用户定）：严格的 90° 谁也到不了——浮点会丢精度，设备更有误差。
    // 落在带内就按极点收，换来的是读数**稳定**：带外 0.05° 处俯仰动一点方位就甩几十度。
    // 代价明码：带内报出的三角是**代表**而非精确分解，重新合成会与原朝向差一点（量级 ≈ |θ₃|·sin(带宽)；
    // 指向误差不超过 2·sin(带宽) ≈ 0.2°，且只发生在"光已经指着天顶"的时候）。
    // **带宽该盖住"水平电机误差投到光束上的那个量"**（实机定标）：铭辰下单天顶时水平差 0.11°，
    // 光束离云台自身极点 28°，投影后 0.11°·sin(28°) ≈ 0.05° —— 光就停在距天顶 0.05° 处。带宽取 0.05°
    // 恰好卡在边界上、时灵时不灵，故取 0.1° 留一倍余量。换设备（水平精度不同）就该重算这个数。
    static constexpr double kPoleCos = 1.7453284e-3; // = sin(0.1°) = cos(89.9°)

    // 三轴正交内旋提取（两模型的 fromOrientation 都用它：WorldModel 三个分量全要，PanTiltModel 截前两个；
    // 非正交轴子类不用它）：把关节轴补成正交 U∈SO(3)（缺的第 3 列 = 前两轴叉积），E=UᵀRU 还原标准内旋 XYZ，
    // 返回 {θ₁,θ₂,θ₃}(°)，θ₁ ∈ [0°,360°)、θ₂ ∈ [−90°,90°]、θ₃ ∈ (−180°,180°]——这就是反解的典范窗（进关节后按窗 [−180°,180°) 收，180° 读成 −180°）。
    [[nodiscard]] std::array<double, 3> intrinsicXYZ(const SO3 &R) const;

private:
    // 关节基下的朝向 E = UᵀRU（U = orthoBasis）：在 E 里第 i 根关节就是第 i 根坐标轴。intrinsicXYZ 建在它上面。
    [[nodiscard]] LinAlg::Mat3 jointBasis(const SO3 &R) const;

    // 把关节轴补成正交 3×3：第 3 列缺则用前两轴叉积补(双轴→û₁×û₂)。要求各轴正交、N≥2(任何 pan-tilt 皆是)。
    [[nodiscard]] LinAlg::Mat3 orthoBasis() const;

    // 关节表：轴向 + 当前值，按复合顺序排列。**私有持有**（见类注）。
    std::vector<Joint> joints_;
};

// ── GimbalModel inline 实现（类外，header-only）─────────────────
inline SO3 GimbalModel::toOrientation() const
{
    constexpr double D2R = std::numbers::pi / 180.0;
    SO3 R; // 单位
    for (const Joint &j : joints_)
        R = R * SO3::exp(j.value() * D2R * j); // 逐轴轴角指数(各保证正交)再复合
    return R;
}

inline LinAlg::Mat3 GimbalModel::jointBasis(const SO3 &R) const
{
    const LinAlg::Mat3 U = orthoBasis();
    return U.transpose() * R * U;
}

// ── 退化（奇异）的完整分类：**奇异 ⟺ 第一根轴与第三根轴共线** ──────────────────
// 三轴链 E = R₁(θ₁)·R₂(θ₂)·R₃(θ₃)，中间那根轴把首尾两根拧到共线时，首尾就只剩一个自由度：
//   · **Tait–Bryan（三轴互异，本模型 x-y-z）**：θ₂ = ±90° 时首尾共线。+90° 只定 θ₁+θ₃，−90° 只定 θ₃−θ₁。
//   · **Euler（首尾同轴，如 z-x-z）**：θ₂ = 0° 或 180° 时首尾共线。0° 定 θ₁+θ₃，180° 定 θ₁−θ₃。
//   本仓两个模型都是 Tait–Bryan，故**只有中间轴 ±90° 这一族**；将来真加了 Euler 型子类，
//   按同一条规矩办即可：**把第三轴的量折进与它共线的那根（三轴链里就是第一根），第三轴取 0**，
//   同向相加、反向相减。（双轴的 PanTiltModel 正解是单射、位形本身从不含糊；它的反解截取本函数的前两个
//   分量，两极的折法随之同一套。）
//
// **两极是多解不是无解**（用户定）：θ₂ = ±90° 时 cos θ₂ = 0，常规两式用到的
// E(1,2)/E(2,2)/E(0,1)/E(0,0) 同时趋零，atan2 退化成 atan2(0,0)——数值上是噪声，而数学上是
// **一族解**：θ₂ = +90° 只定 θ₁+θ₃，θ₂ = −90° 只定 θ₃−θ₁。既然是多解，就该在算法层收成**一个
// 确定的代表**，而不是把烂摊子丢给显示层去判"解不出"：取 **θ₃ = 0，把它折进 θ₁**，
// 两支合成一式 θ₁ = atan2(sinθ₂·E(1,0), E(1,1))。
// 于是同一个朝向恒得同一组三角，正反解在两极也闭合。
// **注意它治的是"正好在极点"**：离极点还有零点零几度时公式并不退化、解是唯一的，但方位角对指向的
// 条件数会炸（高度 89.95° 时俯仰动 0.05°、方位甩 44°）。那是这套参数化的固有性质，不是缺陷。
inline std::array<double, 3> GimbalModel::intrinsicXYZ(const SO3 &R) const
{
    constexpr double R2D = 180.0 / std::numbers::pi;
    // 极点带宽 kPoleCos（基类 protected 常量，选值理由见它的声明处）。
    const LinAlg::Mat3 E = jointBasis(R);
    const double e02 = E(0, 2) < -1.0 ? -1.0 : (E(0, 2) > 1.0 ? 1.0 : E(0, 2)); // sin θ₂
    const double c2 = std::hypot(E(0, 0), E(0, 1));                             // |cos θ₂|
    if (c2 <= kPoleCos)
        return {SO3::Wrap360(std::atan2(e02 * E(1, 0), E(1, 1)) * R2D), // θ₁（θ₃ 已折进来）
                std::asin(e02) * R2D,                                   // θ₂ = ±90°
                0.0};                                                   // θ₃ 取 0 作代表
    return {SO3::Wrap360(std::atan2(-E(1, 2), E(2, 2)) * R2D), // θ₁
            std::asin(e02) * R2D,                              // θ₂
            std::atan2(-E(0, 1), E(0, 0)) * R2D};              // θ₃
}

inline LinAlg::Mat3 GimbalModel::orthoBasis() const
{
    const LinAlg::Vec3 a0 = joints_[0];
    const LinAlg::Vec3 a1 = joints_[1];
    const LinAlg::Vec3 a2 = joints_.size() > 2 ? LinAlg::Vec3(joints_[2]) : LinAlg::Vec3::cross(a0, a1);
    return {a0.x(), a1.x(), a2.x(),
            a0.y(), a1.y(), a2.y(),
            a0.z(), a1.z(), a2.z()};
}

// ============================================================
// WorldModel — 世界系 ENU，本质是一台「三正交轴、内旋」的虚拟云台
//
// 关键洞见：世界系不是独立概念，它就是个 GimbalModel——az/el/roll 即这台三轴虚拟云台的作动量。
// 三轴取 û₁=-ẑ(方位,真北顺时针)、û₂=x̂(高度,抬头为正)、û₃=û₁×û₂=-ŷ(翻滚)，正是 orthoBasis 补出的正交基。
// 反解：三正交轴直接 XYZ 提取。
//
// **一个实例就是一个世界朝向**：安装标定是一个，本次目标是另一个，两者同类不同值。
// 构造收三个角（任意实数，各折进自己的窗、不坍缩）；裸构造给出归零位，即水平正北。
// **翻滚有默认值 0**：指向是 S² 上的一点、朝向是 SO(3) 上的一点，同一指向的全部朝向构成陪集
// R·SO(2)，翻滚是那根纤维的坐标。对"我要指哪"而言它不关心，给个代表元即可；对"基座架成什么样"
// 而言它是实打实的物理，那一处必须显式给。
// ============================================================
class WorldModel : public GimbalModel
{
public:
    // 摆到这个世界朝向：三个角任意实数，各折进自己的窗、不坍缩。翻滚默认 0（见类注）。
    // 窗写在关节上：方位 [0°, 360°)（罗盘方向），高度 / 翻滚 [−180°, 180°)（以水平为零）；读口与界面用同一张窗。
    WorldModel(const double azimuth, const double elevation, const double roll = 0.0)
        : GimbalModel({{{0.0, 0.0, -1.0}, 0.0, azimuth},     // az 方位
                       {{1.0, 0.0, 0.0}, -180.0, elevation}, // el 高度
                       {{0.0, -1.0, 0.0}, -180.0, roll}})    // roll 翻滚
    {
    }
    WorldModel() : WorldModel(0.0, 0.0, 0.0) {} // 归零位：水平正北

    [[nodiscard]] double azimuth() const { return axis(0).value(); }
    [[nodiscard]] double elevation() const { return axis(1).value(); }
    [[nodiscard]] double roll() const { return axis(2).value(); }

    // 朝向 → 世界三角，恒有解、落在典范窗（三轴双解靠高度 ±90° 那一刀收成唯一，两极多解折成一个代表）。
    [[nodiscard]] static WorldModel fromOrientation(const SO3 &R)
    {
        const auto e = WorldModel().intrinsicXYZ(R);
        return {e[0], e[1], e[2]};
    }
};

// ============================================================
// PanTiltModel — 双轴 pan(绕 -ẑ) + tilt(绕 x̂)，内旋，ENU。R = Rz(-pan)·Rx(tilt)。
//
// **本类是位形，不是指向**（用户定；算法层只谈位形，拒绝拿"指向"混淆概念）：
//   · **位形空间** = 两根转轴各一圈 = 环面 T²。构造收任意实数、各折进自己的窗（水平 [0°,360°)、俯仰
//     [−180°,180°)）、不坍缩。它是双轴云台这一数学对象的全部位形，不是哪台设备的行程。
//   · **正解是单射**：Rz(−p)·Rx(t) = Rz(−p')·Rx(t') ⇒ Rz(p'−p) = Rx(t'−t) ⇒ 两边都是单位阵。环面上
//     每个位形各得一个 SO(3) 朝向，**两极也不例外**——"两极处水平角任意"说的是指向，不是位形。
//   · **像是 SO(3) 里的二维子集**：双轴表达不了一般姿态。它恰是世界模型（前两轴相同）把第三轴
//     （翻滚，−ŷ）钉在 0 的那部分。
//
// ⚠ **别把位形空间收窄**（俯仰收到 [0°,90°] 或 [−90°,90°]）：理由听上去是"(p,t) 与 (p+180°,180°−t) 指向相同、
// 二义"，但那是把指向当成了位形——两者**朝向不同**（翻滚差 180°），在位形空间里本就是两个点。收窄的
// 后果是半个环面装不进来（实测俯仰 92° 换不出世界角）。**反解取 |俯仰| ≤ 90° 那支是另一回事**：那是
// "目标翻滚被忽略"之后对指向的唯一代表，位形空间本身仍是整个环面。
//
// **设备行程不在这儿**：某台云台够不够得着由**官邸**按协议定、造斗 / 组包时查（McztVocab::TiltAngleMin/Max、
// SqafConfig::tiltMin/Max），与本类的范围互不抄数（用户定："小解耦"）。发给电机的最后一道闸仍是
// `McztVocab::*2Code` 里的钳位。
//
// **反解 `fromOrientation`：朝向 → 位形**：目标的翻滚直接忽略（用户定：双轴对翻滚无计可施，需求只要指向），
// 取 intrinsicXYZ 的前两个分量——水平 [0°,360°)、俯仰 [−90°,90°]，两极折进水平角，与 WorldModel 同一套。
// ============================================================
class PanTiltModel : public GimbalModel
{
public:
    // 摆到这个位形：两个角任意实数，各折进自己的窗、不坍缩；设备行程不在此判（归官邸）。
    // 窗写在关节上：水平 [0°, 360°)，俯仰 [−180°, 180°)（以水平为零）；读口与界面用同一张窗。
    PanTiltModel(const double pan, const double tilt)
        : GimbalModel({{{0.0, 0.0, -1.0}, 0.0, pan},     // pan：绕 −ẑ
                       {{1.0, 0.0, 0.0}, -180.0, tilt}}) // tilt：绕 x̂
    {
    }
    PanTiltModel() : PanTiltModel(0.0, 0.0) {} // 归零位，一个真实合法的位形

    [[nodiscard]] double pan() const { return axis(0).value(); }
    [[nodiscard]] double tilt() const { return axis(1).value(); }

    // 朝向 → 位形：目标翻滚忽略，取 intrinsicXYZ 的前两个分量（见类注）；设备行程由官邸造斗时查。
    [[nodiscard]] static PanTiltModel fromOrientation(const SO3 &R)
    {
        const auto e = PanTiltModel().intrinsicXYZ(R);
        return {e[0], e[1]};
    }
};

// 按到料时刻排列的世界指向；采样者按前值保持解释相邻回报之间的区间。
using WorldPoseHistory = std::vector<std::pair<std::chrono::system_clock::time_point, WorldModel>>;
