#pragma once
#include <cmath>

// ============================================================
// LinAlg — 3D 线性代数原语（最底座，零几何 / 坐标系 / 角度语义）
//
// 只管「数」本身：Vec3 / Mat3 及其加法、标量乘、点叉积、范数、矩阵积、转置、行列式、逆。
// 这是能搬去任何项目的**通用设施**，故命名空间叫 LinAlg（线性代数），不挂"几何"招牌——
// 几何 / 旋转是线代的一个应用，反过来拿应用命名底座是错位。凡是"知道矩阵是旋转"的东西
// （指数映射、最近旋转、旋转对数、角 ↔ 朝向的约定）都在 SO3.h：本文件管裸矩阵，那个管旋转。
//
// Vec3 / Mat3 是**封装类**：数据私有，对外只给构造函数 + 访问器（Vec3 的 x()/y()/z()、
// Mat3 的 operator()(r,c)）+ 运算。有效性用私有 bool valid_，经 explicit operator bool 暴露
// （if(M)/if(!M) 即判），Empty() 造无效实例（inverse 奇异时返回，此后经运算逐级传播）；不抢占
// double 的 NaN 位、对 /fp:fast 免疫。一步构建结果靠 9 参（Mat3）/ 3 参（Vec3）构造函数——
// `return {e0, …}` 仍是列表初始化命中该构造函数，写起来如数学原样：R = A * B * C。
// 「左操作数即本类型」的运算做成成员（A*B、a+b、v*s、A.transpose()/inverse()/determinant()）；
// 对称的 dot / cross 做成**静态双参**（挂单参成员会制造"谁点谁"的假不对称）；s*v 是唯一的
// 自由函数（左操作数是 double，当不了成员）。
// ============================================================
namespace LinAlg
{
    class Vec3
    {
    public:
        Vec3() = default;
        Vec3(const double x, const double y, const double z) : v_{x, y, z} {}

        [[nodiscard]] double x() const { return v_[0]; }
        [[nodiscard]] double y() const { return v_[1]; }
        [[nodiscard]] double z() const { return v_[2]; }

        // 有效性：私有 bool，经 operator bool 暴露（if(v)/if(!v)）。
        [[nodiscard]] explicit operator bool() const { return valid_; }

        // 常用实例：一律静态工厂。无典范"单位向量"（单位向量不唯一），故不给 Identity。Zeros()/Empty()
        // 复用私有 bool 构造（本体恒零向量，bool 只定有效性）。
        [[nodiscard]] static Vec3 Zeros() { return Vec3(false); }
        [[nodiscard]] static Vec3 Empty() { return Vec3(true); }

        // ---- 代数：能以本类作 lhs 的，一律单参成员（lhs = *this 隐式传入）。任一操作数无效 → 失效
        // 传播：向量结果返回 Empty()，标量结果（dot）返回 NaN。成员能读同类实例私有，直接碰 b.v_。
        // 加法：风场按高度插值要把两层的风按权相加（LidarForward 那一侧）。
        [[nodiscard]] Vec3 operator+(const Vec3 &b) const
        {
            if (!valid_ || !b.valid_) return Empty();
            return {v_[0] + b.v_[0], v_[1] + b.v_[1], v_[2] + b.v_[2]};
        }
        // 标量乘严谨写法是 v*s（被缩放对象在左、标量在右）：s*v 的本义是"s 个 v 相加"，标量在左是计数；
        // 缩放向量时数学定义就是 v*s。成员 lhs 必为本类，恰好只表达 v*s——与数学定义同构。
        [[nodiscard]] Vec3 operator*(const double s) const
        {
            if (!valid_) return Empty();
            return {v_[0] * s, v_[1] * s, v_[2] * s};
        }

        // dot/cross 都不是运算符（无中缀好处），单参成员会读成"a 点/叉 b"的假不对称——故都做成**静态
        // 双参**对称形式。dot 失效返回 NaN（double 无 Empty；NaN 令后续比较失败，安全）。norm2/norm 建在
        // dot/自身上，失效经 NaN 自动传播（sqrt(NaN)=NaN），无需各自再判。
        [[nodiscard]] static double dot(const Vec3 &a, const Vec3 &b)
        {
            if (!a.valid_ || !b.valid_) return NAN;
            return a.v_[0] * b.v_[0] + a.v_[1] * b.v_[1] + a.v_[2] * b.v_[2];
        }
        [[nodiscard]] static Vec3 cross(const Vec3 &a, const Vec3 &b)
        {
            if (!a.valid_ || !b.valid_) return Empty();
            return {a.v_[1] * b.v_[2] - a.v_[2] * b.v_[1], a.v_[2] * b.v_[0] - a.v_[0] * b.v_[2], a.v_[0] * b.v_[1] - a.v_[1] * b.v_[0]};
        }
        // 模方 = ‖v‖²（v·v）：比大小 / 判退化用它，省一次 sqrt（数学上 ‖a‖<‖b‖ ⇔ ‖a‖²<‖b‖²，sqrt
        // 吃算力，能不开方就不开方）。模长 ‖v‖：真要长度（归一化 / 取夹角弧度）才用，走 norm2 再开方。
        [[nodiscard]] double norm2() const { return dot(*this, *this); }
        [[nodiscard]] double norm() const { return std::sqrt(norm2()); }

    private:
        // 私有便捷构造：本体恒为零向量，bool 只定有效性（invalid=true → valid_=false）。Zeros()=
        // Vec3(false)、Empty()=Vec3(true) 复用它。私有 + explicit 挡隐式转换。与 Mat3(bool) 同形。
        explicit Vec3(const bool invalid) : valid_(!invalid) {}

        double v_[3]{};
        bool valid_ = true;
    };

    // s*v：唯一当不了成员的（lhs 是 double），双参自由函数转发到严谨的 v*s。交换律的"坏习惯"在此兜底，
    // 让调用方写 s*v 也能用，但实现仍归结到 v*s 那一份（只 forward，不碰私有，故非 friend、非成员）。
    [[nodiscard]] inline Vec3 operator*(const double s, const Vec3 &v) { return v * s; }

    class Mat3
    {
    public:
        Mat3() = default;
        // 行优先 9 元素构造：让上层把矩阵式一步写出，如 transpose()/inverse() 及 expmap 的 return。
        Mat3(const double a00, const double a01, const double a02,
             const double a10, const double a11, const double a12,
             const double a20, const double a21, const double a22)
            : m_{{a00, a01, a02}, {a10, a11, a12}, {a20, a21, a22}} {}

        // 元素访问器：const 读 / 非 const 读写（支持 M(r,c) += …）。原始访问器，不查 valid_
        // （非 const 版须返回真实左值，无从返回 Empty）；有效性在运算层 / 经 operator bool 判。
        [[nodiscard]] double operator()(const int r, const int c) const { return m_[r][c]; }
        [[nodiscard]] double &operator()(const int r, const int c) { return m_[r][c]; }

        // ---- 左操作数即 Mat3 的运算：成员（任一操作数无效 → 失效传播，返回 Empty）----
        [[nodiscard]] Mat3 operator*(const Mat3 &B) const
        {
            if (!valid_ || !B.valid_)
                return Empty();
            return {m_[0][0] * B.m_[0][0] + m_[0][1] * B.m_[1][0] + m_[0][2] * B.m_[2][0],
                    m_[0][0] * B.m_[0][1] + m_[0][1] * B.m_[1][1] + m_[0][2] * B.m_[2][1],
                    m_[0][0] * B.m_[0][2] + m_[0][1] * B.m_[1][2] + m_[0][2] * B.m_[2][2],
                    m_[1][0] * B.m_[0][0] + m_[1][1] * B.m_[1][0] + m_[1][2] * B.m_[2][0],
                    m_[1][0] * B.m_[0][1] + m_[1][1] * B.m_[1][1] + m_[1][2] * B.m_[2][1],
                    m_[1][0] * B.m_[0][2] + m_[1][1] * B.m_[1][2] + m_[1][2] * B.m_[2][2],
                    m_[2][0] * B.m_[0][0] + m_[2][1] * B.m_[1][0] + m_[2][2] * B.m_[2][0],
                    m_[2][0] * B.m_[0][1] + m_[2][1] * B.m_[1][1] + m_[2][2] * B.m_[2][1],
                    m_[2][0] * B.m_[0][2] + m_[2][1] * B.m_[1][2] + m_[2][2] * B.m_[2][2]};
        }

        // 转置（正交阵的转置即逆）。
        [[nodiscard]] Mat3 transpose() const
        {
            if (!valid_)
                return Empty();
            return {m_[0][0], m_[1][0], m_[2][0],
                    m_[0][1], m_[1][1], m_[2][1],
                    m_[0][2], m_[1][2], m_[2][2]};
        }

        // 标量返回的运算（determinant/dot/norm…）无 Empty 可返，失效则返回 NaN：double 没有"空"状态，
        // NaN 是 IEEE 的"此实数不存在"，且任何 <eps 比较对 NaN 均为 false，失效不会误过阈值。这是临时
        // 计算信号；对象自身的有效位不借 NaN 存（见文件头）。
        [[nodiscard]] double determinant() const
        {
            if (!valid_)
                return NAN;
            return m_[0][0] * (m_[1][1] * m_[2][2] - m_[1][2] * m_[2][1]) -
                   m_[0][1] * (m_[1][0] * m_[2][2] - m_[1][2] * m_[2][0]) +
                   m_[0][2] * (m_[1][0] * m_[2][1] - m_[1][1] * m_[2][0]);
        }

        // 逆：先算行列式判奇异（|det|<1e-15 → 返回 Empty()，调用方 if(!inv) 取舍）；非奇异才值得
        // 算伴随阵（= 余子式转置）。返回阵每项 = 对应代数余子式 / det，一步构造成形。
        [[nodiscard]] Mat3 inverse() const
        {
            if (!valid_)
                return Empty();
            const double det = determinant();
            if (std::abs(det) < 1e-15)
                return Empty();
            const double s = 1.0 / det;
            return {(m_[1][1] * m_[2][2] - m_[1][2] * m_[2][1]) * s,
                    (m_[0][2] * m_[2][1] - m_[0][1] * m_[2][2]) * s,
                    (m_[0][1] * m_[1][2] - m_[0][2] * m_[1][1]) * s,
                    (m_[1][2] * m_[2][0] - m_[1][0] * m_[2][2]) * s,
                    (m_[0][0] * m_[2][2] - m_[0][2] * m_[2][0]) * s,
                    (m_[0][2] * m_[1][0] - m_[0][0] * m_[1][2]) * s,
                    (m_[1][0] * m_[2][1] - m_[1][1] * m_[2][0]) * s,
                    (m_[0][1] * m_[2][0] - m_[0][0] * m_[2][1]) * s,
                    (m_[0][0] * m_[1][1] - m_[0][1] * m_[1][0]) * s};
        }

        [[nodiscard]] explicit operator bool() const { return valid_; }

        // 常用实例：一律静态工厂（统一出口：不让 Identity 是常量、Empty 却是函数）。
        // Identity()/Empty() 复用下面的私有 bool 构造；Zeros() 走公开 9 参构造（valid_ 默认真）。
        [[nodiscard]] static Mat3 Identity() { return Mat3(false); }
        [[nodiscard]] static Mat3 Empty() { return Mat3(true); }
        [[nodiscard]] static Mat3 Zeros() { return {0, 0, 0, 0, 0, 0, 0, 0, 0}; }

    private:
        // 私有便捷构造：本体恒为单位阵，bool 只定有效性（invalid=true → valid_=false，此时本体长啥样
        // 已无所谓）。Identity()=Mat3(false)、Empty()=Mat3(true) 复用它。私有 + explicit 挡隐式转换。
        explicit Mat3(const bool invalid)
            : m_{{1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, {0.0, 0.0, 1.0}}, valid_(!invalid) {}

        double m_[3][3]{};
        bool valid_ = true;
    };
    // 注：SO(3) 流形工具（expmap / nearestRotation / logRotation）在 SO3.h——它们"知道矩阵是旋转"，
    // 属旋转层；本头是零几何语义的纯线代。绕基轴的初等旋转不单列：它只是 expmap 的特例
    // （绕 x̂ 转 a ≡ expmap(a·x̂)）。
} // namespace LinAlg
