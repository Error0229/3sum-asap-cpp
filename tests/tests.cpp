// Randomized cross-checks of every layer against brute force.

#include <algorithm>
#include <cstdint>
#include <cstdio>
#include <limits>
#include <optional>
#include <random>
#include <stdexcept>
#include <vector>

#include "asap/exact_triangle.hpp"
#include "asap/schonhage.hpp"
#include "asap/thin_product.hpp"
#include "asap/three_sum.hpp"

namespace {

int failures = 0;
int checks = 0;

#define CHECK(cond)                                                                \
    do {                                                                           \
        ++checks;                                                                  \
        if (!(cond)) {                                                             \
            ++failures;                                                            \
            std::fprintf(stderr, "%s:%d: CHECK failed: %s\n", __FILE__, __LINE__, #cond); \
        }                                                                          \
    } while (0)

using asap::Matrix;
using asap::Position;

std::mt19937_64 rng(20261008);

std::int64_t uniform(std::int64_t lo, std::int64_t hi) {
    return std::uniform_int_distribution<std::int64_t>(lo, hi)(rng);
}

Matrix random_matrix(std::size_t r, std::size_t c, std::int64_t lo, std::int64_t hi) {
    Matrix M(r, c);
    for (auto& v : M.data) v = uniform(lo, hi);
    return M;
}

std::uint64_t binom(int n, int k) {
    if (k < 0 || k > n) return 0;
    std::uint64_t r = 1;
    for (int i = 0; i < k; ++i) r = r * static_cast<std::uint64_t>(n - i) / static_cast<std::uint64_t>(i + 1);
    return r;
}

std::uint64_t pow_u(std::uint64_t b, int e) {
    std::uint64_t r = 1;
    while (e-- > 0) r *= b;
    return r;
}

// Lemma 2.1: the coefficient of every monomial s t z in sum_λ φ_λ ψ_λ χ_λ equals its coefficient
// in G + E, with p̂, q̂, G and E written out independently of the library's tables.
void test_schonhage_identity() {
    namespace sch = asap::schonhage;
    auto hat_p = [](int i, int j, int s) {  // coefficient of left variable s in p̂_ij
        if (s < 3 || j == 2) return 0;
        const int pi = (s - 3) / 2, pj = (s - 3) % 2;
        if (i < 2) return (pi == i && pj == j) ? 1 : 0;
        return pj == j ? -1 : 0;
    };
    auto hat_q = [](int i, int j, int t) {  // coefficient of right variable t in q̂_ij
        if (t < 3 || i == 2) return 0;
        const int qi = (t - 3) / 2, qj = (t - 3) % 2;
        if (j < 2) return (qi == i && qj == j) ? 1 : 0;
        return qi == i ? -1 : 0;
    };
    for (int s = 0; s < sch::kVars; ++s) {
        for (int t = 0; t < sch::kVars; ++t) {
            for (int z = 0; z < sch::kTerms; ++z) {
                int lhs = 0;
                for (int term = 0; term < sch::kTerms; ++term) {
                    if (sch::contributes(term, z)) lhs += sch::kTables.phi[term][s] * sch::kTables.psi[term][t];
                }
                int rhs;
                if (z == sch::kZero) {
                    rhs = (s >= 3 && s == t) ? 1 : 0;  // p_ij q_ij z_0
                } else {
                    const int i = z / 3, j = z % 3;
                    const int x = s == i, y = t == j;
                    rhs = x * y + x * hat_q(i, j, t) + hat_p(i, j, s) * y + hat_p(i, j, s) * hat_q(i, j, t);
                }
                CHECK(lhs == rhs);
            }
        }
    }
}

std::int64_t naive_entry(const Matrix& X, const Matrix& Y, std::size_t i, std::size_t j) {
    std::int64_t sum = 0;
    for (std::size_t k = 0; k < X.cols; ++k) sum += X(i, k) * Y(k, j);
    return sum;
}

// Theorem 2.1, checked against inner products for many shapes, recursion depths and densities.
void test_thin_product() {
    struct Case {
        std::size_t n_rows, inner, n_cols;
        int L;
        double density;
        std::int64_t bound;
    };
    const std::vector<Case> cases = {
        {5, 4, 7, 1, 1.0, 9},       {30, 4, 30, 2, 0.5, 9},     {40, 3, 25, 3, 0.3, 9},
        {50, 4, 50, 4, 0.2, 9},     {20, 16, 20, 2, 1.0, 9},    {45, 16, 45, 4, 0.3, 9},
        {60, 10, 60, 5, 0.1, 9},    {30, 16, 30, 5, 0.25, 1},   {10, 64, 12, 3, 0.5, 9},
        {35, 16, 35, 4, 0.5, 1 << 28},
    };
    for (const Case& cs : cases) {
        const Matrix X = random_matrix(cs.n_rows, cs.inner, -cs.bound, cs.bound);
        const Matrix Y = random_matrix(cs.inner, cs.n_cols, -cs.bound, cs.bound);
        std::vector<Position> W;
        for (std::size_t i = 0; i < cs.n_rows; ++i) {
            for (std::size_t j = 0; j < cs.n_cols; ++j) {
                if (std::uniform_real_distribution<double>(0, 1)(rng) < cs.density) W.push_back({i, j});
            }
        }
        if (!W.empty()) W.push_back(W.front());  // a duplicate position
        asap::ThinProductOptions opt;
        opt.L = cs.L;
        asap::ThinProductStats st;
        const auto got = asap::thin_product_entries(X, Y, W, opt, &st);
        bool all = true;
        for (std::size_t k = 0; k < W.size(); ++k) all &= got[k] == naive_entry(X, Y, W[k].row, W[k].col);
        CHECK(all);

        // Lemma 2.6: a tile with wanted set U visits at most sum_d min(|U| α_d, β_d) leaves,
        // where α_d = C(m,d) 9^d and β_d = C(L, m-d) 9^(L-m+d).
        for (const auto& [u, leaves] : st.per_tile) {
            std::uint64_t bound = 0;
            for (int d = 0; d <= st.m; ++d) {
                const std::uint64_t alpha = binom(st.m, d) * pow_u(9, d);
                const std::uint64_t beta = binom(st.L, st.m - d) * pow_u(9, st.L - st.m + d);
                bound += std::min(u * alpha, beta);
            }
            CHECK(leaves <= bound);
        }
    }

    // A single wanted entry uses exactly its 10^m contributing leaves.
    {
        const Matrix X = random_matrix(40, 16, -3, 3);
        const Matrix Y = random_matrix(16, 40, -3, 3);
        asap::ThinProductOptions opt;
        opt.L = 4;
        asap::ThinProductStats st;
        const auto got = asap::thin_product_entries(X, Y, {{17, 33}}, opt, &st);
        CHECK(got[0] == naive_entry(X, Y, 17, 33));
        CHECK(st.leaves == 100);
    }

    // The paper's parameters (L = 19m) need 10^38 encoded numbers and are refused.
    {
        const Matrix X = random_matrix(4, 16, 0, 1);
        const Matrix Y = random_matrix(16, 4, 0, 1);
        bool threw = false;
        try {
            asap::thin_product_entries(X, Y, {{0, 0}});
        } catch (const std::length_error&) {
            threw = true;
        }
        CHECK(threw);
    }

    // An empty request needs no encoding, so even the paper's parameters are fine.
    {
        const Matrix X = random_matrix(1, 1, 0, 1);
        const Matrix Y = random_matrix(1, 1, 0, 1);
        bool threw = false;
        std::vector<std::int64_t> got;
        try {
            got = asap::thin_product_entries(X, Y, {});
        } catch (const std::exception&) {
            threw = true;
        }
        CHECK(!threw && got.empty());
    }
}

asap::ExactTriangleOptions small_triangle_options(int L) {
    asap::ExactTriangleOptions opt;
    opt.D = 16;
    opt.g = 1;
    opt.thin.L = L;
    return opt;
}

std::int64_t weight(const asap::ExactTriangleInstance& in, std::size_t a, std::size_t b, std::size_t c) {
    return in.w_ab(a, b) + in.w_bc(b, c) + in.w_ac(a, c);
}

asap::ExactTriangleInstance random_triangle_instance(std::size_t nA, std::size_t nB, std::size_t nC, std::int64_t bound) {
    return {random_matrix(nA, nB, -bound, bound), random_matrix(nB, nC, -bound, bound),
            random_matrix(nA, nC, -bound, bound)};
}

// Theorems 3.4 and 3.5, checked against exhaustive search.
void test_exact_triangle() {
    struct Config {
        std::uint64_t D, g;
        int L;
    };
    const std::vector<Config> configs = {{16, 1, 2}, {16, 2, 4}, {16, 4, 3}, {64, 1, 4}, {64, 3, 3}, {256, 2, 4}};
    for (int round = 0; round < 120; ++round) {
        const Config cfg = configs[round % configs.size()];
        const std::size_t nA = uniform(1, 24), nB = uniform(1, 24), nC = uniform(1, 24);
        const std::int64_t bound = round % 3 == 0 ? 1000 : uniform(2, 30);
        asap::ExactTriangleInstance in = random_triangle_instance(nA, nB, nC, bound);
        if (round % 4 == 0) {  // plant a zero triangle
            const std::size_t a = uniform(0, nA - 1), b = uniform(0, nB - 1), c = uniform(0, nC - 1);
            in.w_ab(a, b) = -(in.w_bc(b, c) + in.w_ac(a, c));
        }
        asap::ExactTriangleOptions opt;
        opt.D = cfg.D;
        opt.g = cfg.g;
        opt.thin.L = cfg.L;
        asap::ExactTriangleStats st;
        const auto got = asap::exact_triangle(in, opt, &st);
        const auto want = asap::exact_triangle_brute_force(in);
        CHECK(got.has_value() == want.has_value());
        if (got) CHECK(weight(in, got->a, got->b, got->c) == 0);
        CHECK(!st.brute_force);

        // The chosen prime has the fewest triples with S ≡ 0 (mod p) in [sqrt(D)/2, sqrt(D)), and
        // every failed scan meets a distinct false positive of it.
        auto count_mod = [&](std::int64_t p, bool false_only) {
            std::uint64_t count = 0;
            for (std::size_t a = 0; a < nA; ++a)
                for (std::size_t b = 0; b < nB; ++b)
                    for (std::size_t c = 0; c < nC; ++c) {
                        const std::int64_t s = weight(in, a, b, c);
                        if (s % p == 0 && (!false_only || s != 0)) ++count;
                    }
            return count;
        };
        std::uint64_t sqrt_d = 1;
        while (sqrt_d * sqrt_d < cfg.D) ++sqrt_d;
        const std::uint64_t chosen = count_mod(static_cast<std::int64_t>(st.prime), false);
        for (std::uint64_t q = sqrt_d / 2; q < sqrt_d; ++q) {
            bool prime = q >= 2;
            for (std::uint64_t d = 2; d * d <= q; ++d) prime &= q % d != 0;
            if (prime) CHECK(chosen <= count_mod(static_cast<std::int64_t>(q), false));
        }
        CHECK(st.failed_scans <= count_mod(static_cast<std::int64_t>(st.prime), true));
    }

    // Weights at the limit ±2^61 are handled exactly; weights beyond it, such as INT64_MIN, are rejected.
    {
        constexpr std::int64_t big = asap::kMaxAbsWeight;
        asap::ExactTriangleInstance in = random_triangle_instance(12, 12, 12, big);
        in.w_ab(3, 5) = big;
        in.w_bc(5, 7) = -big;
        in.w_ac(3, 7) = 0;
        const auto got = asap::exact_triangle(in, small_triangle_options(4));
        CHECK(got.has_value());
        if (got) CHECK(weight(in, got->a, got->b, got->c) == 0);
        CHECK(asap::exact_triangle_brute_force(in).has_value());

        for (const std::int64_t bad : {big + 1, std::numeric_limits<std::int64_t>::min()}) {
            in.w_bc(0, 0) = bad;
            bool threw = false;
            try {
                asap::exact_triangle(in, small_triangle_options(4));
            } catch (const std::out_of_range&) {
                threw = true;
            }
            CHECK(threw);
            threw = false;
            try {
                asap::exact_triangle_brute_force(in);
            } catch (const std::out_of_range&) {
                threw = true;
            }
            CHECK(threw);
        }
    }

    // With the paper's parameters, D = 16 needs n >= 16^18, so every feasible instance is a base case.
    {
        const asap::ExactTriangleInstance in = random_triangle_instance(10, 10, 10, 5);
        asap::ExactTriangleStats st;
        const auto got = asap::exact_triangle(in, {}, &st);
        CHECK(st.brute_force);
        CHECK(got.has_value() == asap::exact_triangle_brute_force(in).has_value());
    }
}

void test_convolution_three_sum() {
    for (int round = 0; round < 60; ++round) {
        const std::size_t n = uniform(1, 70);
        const std::int64_t bound = uniform(3, 40);
        asap::ConvolutionThreeSumInstance in;
        for (auto* v : {&in.x, &in.y, &in.z}) {
            for (std::size_t i = 0; i < n; ++i) {
                if (uniform(0, 9) < 7) {
                    v->push_back(uniform(-bound, bound));
                } else {
                    v->push_back(std::nullopt);
                }
            }
        }
        bool want = false;
        for (std::size_t i = 0; i < n; ++i)
            for (std::size_t j = 0; i + j < n; ++j)
                want |= in.x[i] && in.y[j] && in.z[i + j] && *in.x[i] + *in.y[j] == *in.z[i + j];
        const auto got = asap::convolution_three_sum(in, small_triangle_options(round % 2 ? 3 : 4));
        CHECK(got.has_value() == want);
        if (got) {
            CHECK(got->i + got->j < n);
            CHECK(in.x[got->i] && in.y[got->j] && in.z[got->i + got->j]);
            CHECK(*in.x[got->i] + *in.y[got->j] == *in.z[got->i + got->j]);
        }
    }
}

std::vector<std::int64_t> random_set(std::size_t n, std::int64_t bound) {
    std::vector<std::int64_t> v(n);
    for (auto& x : v) x = uniform(-bound, bound);
    return v;
}

bool contains(const std::vector<std::int64_t>& v, std::int64_t x) { return std::find(v.begin(), v.end(), x) != v.end(); }

void test_three_sum() {
    for (int round = 0; round < 60; ++round) {
        const std::size_t n = uniform(1, 45);
        const std::int64_t bound = round % 3 == 0 ? 1'000'000'000'000LL : uniform(n, 6 * n);
        std::vector<std::int64_t> A = random_set(n, bound), B = random_set(uniform(1, 45), bound),
                                  C = random_set(uniform(1, 45), bound);
        if (round % 3 == 0 && round % 2 == 0) C.push_back(A[uniform(0, A.size() - 1)] + B[uniform(0, B.size() - 1)]);
        asap::ThreeSumOptions opt;
        opt.triangle = small_triangle_options(3);
        // 0 is automatic; 64 exceeds every bucket, so every solution must come through
        // Convolution 3SUM, Exact Triangle and the thin matrix product.
        opt.bucket_cap = round % 5 == 4 ? 64 : static_cast<std::uint64_t>(round % 4);
        asap::ThreeSumStats st;
        const auto got = asap::three_sum(A, B, C, opt, &st);
        const auto want = asap::three_sum_quadratic(A, B, C);
        if (opt.bucket_cap == 64) CHECK(st.heavy == 0);
        CHECK(got.has_value() == want.has_value());
        if (got) {
            CHECK(got->a + got->b == got->c);
            CHECK(contains(A, got->a) && contains(B, got->b) && contains(C, got->c));
        }
    }

    // Distinct-index, single-list form.
    for (int round = 0; round < 40; ++round) {
        const std::size_t n = uniform(0, 30);
        const std::vector<std::int64_t> nums = random_set(n, uniform(1, 3 * n + 1));
        bool want = false;
        for (std::size_t i = 0; i < n; ++i)
            for (std::size_t j = i + 1; j < n; ++j)
                for (std::size_t k = j + 1; k < n; ++k) want |= nums[i] + nums[j] + nums[k] == 0;
        asap::ThreeSumOptions opt;
        opt.triangle = small_triangle_options(3);
        const auto got = asap::three_sum_zero(nums, opt);
        CHECK(got.has_value() == want);
        if (got) {
            const auto [i, j, k] = *got;
            CHECK(i < j && j < k && k < n);
            CHECK(nums[i] + nums[j] + nums[k] == 0);
        }
    }
}

}  // namespace

int main() {
    test_schonhage_identity();
    test_thin_product();
    test_exact_triangle();
    test_convolution_three_sum();
    test_three_sum();
    std::printf("%d checks, %d failures\n", checks, failures);
    return failures == 0 ? 0 : 1;
}
