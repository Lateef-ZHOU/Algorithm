#pragma once
#include <cmath>
#include <utility>
#include <QtGlobal>

// ============================================================
// GaussFit — 单峰高斯拟合（谱线 → 峰位/峰宽/幅值/基底）
//
// 「构造即拟合」：传入一段谱线（频点索引空间，与世界坐标无关），ctor 内部完成初值估计 +
// Levenberg-Marquardt 迭代。每一步算法都是 **private static** 纯函数——本类对外只有"喂一段谱、
// 读拟合结果"这一条机制化用途，不暴露任何中间步骤。
//
// 结果（params）作 **公开 const 数据成员**铺开，构造期定死、不可变：
//   GaussFit fit(rowBuf, cols);
//   fit.params.x0 ...
// 实现：公开 ctor 在初始化列表里调静态 fit() 直接初始化 const 成员（const 成员只能初始化列表
// 赋值，静态函数入列表合法）。由 RadialWind::Calc 逐行（每个距离门一条谱）调用。
//
// 声明与定义分开：类内只留声明，inline 定义集中在类外（仍 header-only，单文件）。
// ============================================================
class GaussFit
{
public:
    struct Params
    {
        double A = 1.0;  ///< 峰幅值
        double x0 = 0.0; ///< 峰位（bin）
        double w = 1.0;  ///< 高斯宽（σ，bin）
        double d = 0.0;  ///< 基底
    };

    // 构造即拟合（n<4 → 仅以初值落定）。
    explicit GaussFit(const double *y, int n, int maxIter = 50, double tol = 1e-6);

    const Params params; ///< 拟合参数，构造期定死

private:
    // ── 私有算法步骤（全部静态纯函数）──
    [[nodiscard]] static double gauss(double x, double A, double x0, double w, double d);
    [[nodiscard]] static double sumSquaredError(const double *y, int n, const Params &p);
    [[nodiscard]] static Params estimateInitial(const double *y, int n);
    [[nodiscard]] static Params fit(const double *y, int n, int maxIter, double tol);
};

// ============================================================
// GaussFit — inline 实现（类外，header-only）
// ============================================================

// 构造即拟合：静态 fit 的产物直接初始化 const 成员。
inline GaussFit::GaussFit(const double *y, const int n, const int maxIter, const double tol)
    : params(fit(y, n, maxIter, tol))
{
}

inline double GaussFit::gauss(const double x, const double A, const double x0,
                              const double w, const double d)
{
    const double dx = (x - x0) / w;
    return A * std::exp(-0.5 * dx * dx) + d;
}

// 残差平方和（SSE）：模型 gauss(i; p) 与样本 y 的逐点差平方累加，拟合优度度量。
inline double GaussFit::sumSquaredError(const double *y, const int n, const Params &p)
{
    double sse = 0.0;
    for (int i = 0; i < n; ++i)
    {
        const double r = y[i] - gauss(i, p.A, p.x0, p.w, p.d);
        sse += r * r;
    }
    return sse;
}

// 初值估计：基底取两端均值，峰位/幅值取最大点，宽由半高全宽 / 2.3548 估。
inline GaussFit::Params GaussFit::estimateInitial(const double *y, const int n)
{
    Params p;
    if (n <= 0) return p;

    p.d = (y[0] + y[n - 1]) * 0.5;

    int iMax = 0;
    double yMax = y[0];
    for (int i = 1; i < n; ++i)
    {
        if (y[i] > yMax)
        {
            yMax = y[i];
            iMax = i;
        }
    }
    p.A = yMax - p.d;
    p.x0 = static_cast<double>(iMax);

    const double halfMax = p.d + p.A * 0.5;
    int left = iMax, right = iMax;
    while (left > 0 && y[left] > halfMax) --left;
    while (right < n - 1 && y[right] > halfMax) ++right;
    const double fwhm = right - left;
    p.w = fwhm > 0.0 ? fwhm / 2.3548 : 1.0;

    if (p.A <= 0.0) p.A = 1.0;
    if (p.w <= 0.0) p.w = 1.0;

    return p;
}

// 拟合主体（Levenberg-Marquardt）：estimateInitial → 迭代，返回参数。
// n<4（参数个数）→ 仅以初值落定。
inline GaussFit::Params GaussFit::fit(const double *y, const int n, const int maxIter, const double tol)
{
    constexpr int NP = 4;

    Params p0 = estimateInitial(y, n);
    if (n < NP) return p0;

    double A = p0.A, x0 = p0.x0, w = p0.w, d = p0.d;

    double lambda = 1e-3;
    double sse = sumSquaredError(y, n, {A, x0, w, d});

    for (int iter = 0; iter < maxIter; ++iter)
    {
        double JtJ[NP][NP] = {};
        double Jtr[NP] = {};

        for (int i = 0; i < n; ++i)
        {
            const double xi = i;
            const double dx = (xi - x0) / w;
            const double expv = std::exp(-0.5 * dx * dx);
            const double ri = y[i] - (A * expv + d);

            double J[NP];
            J[0] = expv;
            J[1] = A * expv * dx / w;
            J[2] = A * expv * dx * dx / w;
            J[3] = 1.0;

            for (int r = 0; r < NP; ++r)
            {
                Jtr[r] += J[r] * ri;
                for (int c = 0; c < NP; ++c)
                    JtJ[r][c] += J[r] * J[c];
            }
        }

        double H[NP][NP];
        for (int r = 0; r < NP; ++r)
        {
            for (int c = 0; c < NP; ++c)
                H[r][c] = JtJ[r][c];
            H[r][r] *= 1.0 + lambda;
        }

        double aug[NP][NP + 1];
        for (int r = 0; r < NP; ++r)
        {
            for (int c = 0; c < NP; ++c)
                aug[r][c] = H[r][c];
            aug[r][NP] = Jtr[r];
        }

        bool singular = false;
        for (int col = 0; col < NP; ++col)
        {
            int pivotRow = col;
            for (int r = col + 1; r < NP; ++r)
                if (std::abs(aug[r][col]) > std::abs(aug[pivotRow][col]))
                    pivotRow = r;
            if (pivotRow != col)
                for (int c = 0; c <= NP; ++c)
                    std::swap(aug[col][c], aug[pivotRow][c]);
            if (std::abs(aug[col][col]) < 1e-15)
            {
                singular = true;
                break;
            }
            for (int r = col + 1; r < NP; ++r)
            {
                const double factor = aug[r][col] / aug[col][col];
                for (int c = col; c <= NP; ++c)
                    aug[r][c] -= factor * aug[col][c];
            }
        }

        if (singular)
        {
            lambda *= 10.0;
            continue;
        }

        double dp[NP];
        for (int r = NP - 1; r >= 0; --r)
        {
            dp[r] = aug[r][NP];
            for (int c = r + 1; c < NP; ++c)
                dp[r] -= aug[r][c] * dp[c];
            dp[r] /= aug[r][r];
        }

        const double A_new = A + dp[0];
        const double x0_new = x0 + dp[1];
        double w_new = w + dp[2];
        const double d_new = d + dp[3];

        if (w_new <= 0.0) w_new = std::abs(w_new) + 1e-6;

        if (const double sse_new = sumSquaredError(y, n, {A_new, x0_new, w_new, d_new}); sse_new < sse)
        {
            A = A_new;
            x0 = x0_new;
            w = w_new;
            d = d_new;

            const double improvement = (sse - sse_new) / (sse + 1e-30);
            sse = sse_new;
            lambda *= 0.1;

            if (improvement < tol)
                break;
        }
        else
        {
            lambda *= 10.0;
        }
    }

    return {A, x0, std::abs(w), d};
}
