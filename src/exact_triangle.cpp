#include "asap/exact_triangle.hpp"

#include <algorithm>
#include <limits>
#include <stdexcept>
#include <vector>

namespace asap {
namespace {

using u64 = std::uint64_t;

std::int64_t mod(std::int64_t x, std::int64_t p) {
    const std::int64_t r = x % p;
    return r < 0 ? r + p : r;
}

bool is_prime(u64 x) {
    if (x < 2) return false;
    for (u64 d = 2; d * d <= x; ++d) {
        if (x % d == 0) return false;
    }
    return true;
}

// x^e, saturating at the largest u64.
u64 saturating_pow(u64 x, int e) {
    u64 r = 1;
    for (int i = 0; i < e; ++i) {
        if (x != 0 && r > std::numeric_limits<u64>::max() / x) return std::numeric_limits<u64>::max();
        r *= x;
    }
    return r;
}

// Square matrices over the ring Z[x]/(x^p - 1): entry (i, j) is its vector of p coefficients.
struct PolyMatrix {
    std::size_t n = 0;
    std::size_t p = 0;
    std::vector<std::int64_t> c;

    PolyMatrix(std::size_t n_, std::size_t p_) : n(n_), p(p_), c(n_ * n_ * p_, 0) {}
    std::int64_t* at(std::size_t i, std::size_t j) { return c.data() + (i * n + j) * p; }
    const std::int64_t* at(std::size_t i, std::size_t j) const { return c.data() + (i * n + j) * p; }
};

PolyMatrix quadrant(const PolyMatrix& M, std::size_t qi, std::size_t qj) {
    const std::size_t h = M.n / 2;
    PolyMatrix Q(h, M.p);
    for (std::size_t i = 0; i < h; ++i) {
        std::copy_n(M.at(qi * h + i, qj * h), h * M.p, Q.at(i, 0));
    }
    return Q;
}

void place(PolyMatrix& M, const PolyMatrix& Q, std::size_t qi, std::size_t qj) {
    for (std::size_t i = 0; i < Q.n; ++i) {
        std::copy_n(Q.at(i, 0), Q.n * Q.p, M.at(qi * Q.n + i, qj * Q.n));
    }
}

PolyMatrix combine(const PolyMatrix& X, const PolyMatrix& Y, int sign) {
    PolyMatrix Z = X;
    for (std::size_t k = 0; k < Z.c.size(); ++k) Z.c[k] += sign * Y.c[k];
    return Z;
}

PolyMatrix operator+(const PolyMatrix& X, const PolyMatrix& Y) { return combine(X, Y, 1); }
PolyMatrix operator-(const PolyMatrix& X, const PolyMatrix& Y) { return combine(X, Y, -1); }

// Strassen's algorithm (Section 2.1) over Z[x]/(x^p - 1). n must be a power of two.
PolyMatrix multiply(const PolyMatrix& A, const PolyMatrix& B) {
    const std::size_t n = A.n, p = A.p;
    PolyMatrix C(n, p);
    if (n <= 32) {
        for (std::size_t i = 0; i < n; ++i) {
            for (std::size_t k = 0; k < n; ++k) {
                const std::int64_t* a = A.at(i, k);
                for (std::size_t j = 0; j < n; ++j) {
                    const std::int64_t* b = B.at(k, j);
                    std::int64_t* out = C.at(i, j);
                    for (std::size_t u = 0; u < p; ++u) {
                        if (a[u] == 0) continue;
                        for (std::size_t v = 0; v < p; ++v) out[(u + v) % p] += a[u] * b[v];
                    }
                }
            }
        }
        return C;
    }
    const PolyMatrix A11 = quadrant(A, 0, 0), A12 = quadrant(A, 0, 1), A21 = quadrant(A, 1, 0), A22 = quadrant(A, 1, 1);
    const PolyMatrix B11 = quadrant(B, 0, 0), B12 = quadrant(B, 0, 1), B21 = quadrant(B, 1, 0), B22 = quadrant(B, 1, 1);
    const PolyMatrix S1 = multiply(A11 + A22, B11 + B22);
    const PolyMatrix S2 = multiply(A21 + A22, B11);
    const PolyMatrix S3 = multiply(A11, B12 - B22);
    const PolyMatrix S4 = multiply(A22, B21 - B11);
    const PolyMatrix S5 = multiply(A11 + A12, B22);
    const PolyMatrix S6 = multiply(A21 - A11, B11 + B12);
    const PolyMatrix S7 = multiply(A12 - A22, B21 + B22);
    place(C, S1 + S4 - S5 + S7, 0, 0);
    place(C, S3 + S5, 0, 1);
    place(C, S2 + S4, 1, 0);
    place(C, S1 - S2 + S3 + S6, 1, 1);
    return C;
}

std::int64_t triangle_weight(const ExactTriangleInstance& in, std::size_t a, std::size_t b, std::size_t c) {
    return in.w_ab(a, b) + in.w_bc(b, c) + in.w_ac(a, c);
}

// The number of triples (a, b, c) with S(a, b, c) ≡ 0 (mod p), that is F(p) + Z_0, by the
// small-weights trick of the proof of Theorem 3.4: with P[a, c] = x^(w(a,c) mod p) and
// Q[c, b] = x^(w(b,c) mod p), the coefficient of x^r in (PQ)[a, b] counts the c with
// w(a,c) + w(b,c) ≡ r, and we read it at r = -w(a,b) mod p.
u64 count_zero_mod_p(const ExactTriangleInstance& in, std::int64_t p) {
    const std::size_t nA = in.w_ab.rows, nB = in.w_ab.cols, nC = in.w_ac.cols;
    std::size_t size = 1;
    while (size < std::max({nA, nB, nC})) size *= 2;
    PolyMatrix P(size, static_cast<std::size_t>(p)), Q(size, static_cast<std::size_t>(p));
    for (std::size_t a = 0; a < nA; ++a) {
        for (std::size_t c = 0; c < nC; ++c) P.at(a, c)[mod(in.w_ac(a, c), p)] = 1;
    }
    for (std::size_t c = 0; c < nC; ++c) {
        for (std::size_t b = 0; b < nB; ++b) Q.at(c, b)[mod(in.w_bc(b, c), p)] = 1;
    }
    const PolyMatrix R = multiply(P, Q);
    u64 count = 0;
    for (std::size_t a = 0; a < nA; ++a) {
        for (std::size_t b = 0; b < nB; ++b) count += static_cast<u64>(R.at(a, b)[mod(-in.w_ab(a, b), p)]);
    }
    return count;
}

}  // namespace

std::optional<Triangle> exact_triangle_brute_force(const ExactTriangleInstance& in) {
    for (std::size_t a = 0; a < in.w_ab.rows; ++a) {
        for (std::size_t b = 0; b < in.w_ab.cols; ++b) {
            for (std::size_t c = 0; c < in.w_ac.cols; ++c) {
                if (triangle_weight(in, a, b, c) == 0) return Triangle{a, b, c};
            }
        }
    }
    return std::nullopt;
}

std::optional<Triangle> exact_triangle(const ExactTriangleInstance& in, const ExactTriangleOptions& options,
                                       ExactTriangleStats* stats) {
    const std::size_t nA = in.w_ab.rows, nB = in.w_ab.cols, nC = in.w_ac.cols;
    if (in.w_ac.rows != nA || in.w_bc.rows != nB || in.w_bc.cols != nC) {
        throw std::invalid_argument("exact_triangle: inconsistent part sizes");
    }
    ExactTriangleStats local;
    ExactTriangleStats& st = stats ? *stats : local;
    st = ExactTriangleStats{};
    if (nA == 0 || nB == 0 || nC == 0) return std::nullopt;

    u64 D = options.D;
    if (D == 0) {
        // The largest power of four with D^18 <= n. Smaller instances (D < 16) are solved by brute force.
        const u64 n = std::max({nA, nB, nC});
        D = 1;
        while (saturating_pow(4 * D, 18) <= n) D *= 4;
        if (D < 16) {
            st.brute_force = true;
            return exact_triangle_brute_force(in);
        }
    }
    int m = 0;
    for (u64 d = D; d > 1; d /= 4) {
        if (d % 4 != 0) throw std::invalid_argument("exact_triangle: D must be a power of four");
        ++m;
    }
    if (m < 2) throw std::invalid_argument("exact_triangle: D must be at least 16");
    const u64 sqrt_d = u64{1} << m;
    u64 g = options.g;
    if (g == 0) {
        g = 1;
        while (saturating_pow(g, 36) < D) ++g;  // ceil(D^(1/36))
    }
    if (g < 1 || g > sqrt_d) throw std::invalid_argument("exact_triangle: g must be in [1, sqrt(D)]");
    st.D = D;
    st.g = g;

    // Hashing modulo a prime p in [sqrt(D)/2, sqrt(D)), the one with the fewest false positives.
    std::int64_t p = 0;
    u64 best = std::numeric_limits<u64>::max();
    for (u64 q = sqrt_d / 2; q < sqrt_d; ++q) {
        if (!is_prime(q)) continue;
        const u64 count = count_zero_mod_p(in, static_cast<std::int64_t>(q));
        if (count < best) {
            best = count;
            p = static_cast<std::int64_t>(q);
        }
    }
    st.prime = static_cast<u64>(p);

    // The instances. C is cut into pieces of at most ceil(s/g) vertices, s = floor(sqrt(D)), and
    // W_ρ = {(a, b) : w(a, b) ≡ ρ} into chunks of at most |A||B|/sqrt(D) query pairs. The
    // instance of a chunk of W_ρ and a piece C_k has the middle part C_k x Z_p, with
    //   a ~ (c, σ)  iff  σ ≡ w(a, c) + ρ,    (c, σ) ~ b  iff  σ ≡ -w(b, c)   (mod p).
    const u64 piece = (sqrt_d + g - 1) / g;
    const u64 chunk = std::max<u64>(1, static_cast<u64>(nA) * nB / sqrt_d);
    std::vector<std::vector<Position>> by_residue(static_cast<std::size_t>(p));
    for (std::size_t a = 0; a < nA; ++a) {
        for (std::size_t b = 0; b < nB; ++b) by_residue[mod(in.w_ab(a, b), p)].push_back({a, b});
    }
    std::vector<Matrix> right;  // the M x B side depends only on the piece
    for (std::size_t c0 = 0; c0 < nC; c0 += piece) {
        Matrix Y(D, nB);
        for (std::size_t c = c0; c < std::min<std::size_t>(nC, c0 + piece); ++c) {
            for (std::size_t b = 0; b < nB; ++b) Y((c - c0) * p + mod(-in.w_bc(b, c), p), b) = 1;
        }
        right.push_back(std::move(Y));
    }

    for (std::int64_t rho = 0; rho < p; ++rho) {
        const std::vector<Position>& pairs = by_residue[rho];
        if (pairs.empty()) continue;
        for (std::size_t k = 0; k < right.size(); ++k) {
            const std::size_t c0 = k * piece, c1 = std::min<std::size_t>(nC, c0 + piece);
            Matrix X(nA, D);
            for (std::size_t a = 0; a < nA; ++a) {
                for (std::size_t c = c0; c < c1; ++c) X(a, (c - c0) * p + mod(in.w_ac(a, c) + rho, p)) = 1;
            }
            for (std::size_t first = 0; first < pairs.size(); first += chunk) {
                const std::vector<Position> W(pairs.begin() + first,
                                              pairs.begin() + std::min<std::size_t>(pairs.size(), first + chunk));
                ThinProductStats thin_stats;
                const std::vector<std::int64_t> common = count_lopsided_sparse_triangles(X, right[k], W, options.thin, &thin_stats);
                ++st.instances;
                st.query_pairs += W.size();
                st.leaves += thin_stats.leaves;

                // Witnesses: scan the piece for every accepted query pair.
                for (std::size_t i = 0; i < W.size(); ++i) {
                    if (common[i] == 0) continue;
                    ++st.accepted;
                    for (std::size_t c = c0; c < c1; ++c) {
                        if (triangle_weight(in, W[i].row, W[i].col, c) == 0) return Triangle{W[i].row, W[i].col, c};
                    }
                    ++st.failed_scans;
                }
            }
        }
    }
    return std::nullopt;
}

}  // namespace asap
