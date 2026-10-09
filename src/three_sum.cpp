#include "asap/three_sum.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <map>
#include <stdexcept>

namespace asap {
namespace {

using u64 = std::uint64_t;

std::int64_t mod(std::int64_t x, std::int64_t m) {
    const std::int64_t r = x % m;
    return r < 0 ? r + m : r;
}

bool is_prime(u64 x) {
    if (x < 2) return false;
    for (u64 d = 2; d * d <= x; ++d) {
        if (x % d == 0) return false;
    }
    return true;
}

void check_range(const std::vector<std::int64_t>& v) {
    for (std::int64_t x : v) {
        if (x > kMaxAbsValue || x < -kMaxAbsValue) throw std::out_of_range("three_sum: |value| exceeds kMaxAbsValue");
    }
}

std::vector<std::int64_t> sorted_unique(std::vector<std::int64_t> v) {
    std::sort(v.begin(), v.end());
    v.erase(std::unique(v.begin(), v.end()), v.end());
    return v;
}

// The number of pairs (x, x') in X^2 with x ≡ x' (mod m), summed over the sets.
u64 collisions(const std::vector<const std::vector<std::int64_t>*>& sets, std::int64_t m) {
    u64 total = 0;
    std::vector<std::int64_t> residues;
    for (const auto* set : sets) {
        residues.clear();
        for (std::int64_t x : *set) residues.push_back(mod(x, m));
        std::sort(residues.begin(), residues.end());
        for (std::size_t i = 0; i < residues.size();) {
            std::size_t j = i;
            while (j < residues.size() && residues[j] == residues[i]) ++j;
            total += static_cast<u64>((j - i) * (j - i));
            i = j;
        }
    }
    return total;
}

// Deterministic additive hashing with few collisions [Chan-Lewenstein, Chan-He], as stated in
// Fischer, Kaliciak, Polak (Lemma 3.1): m starts at 1 and is multiplied, one at a time, by the
// prime of size about n^δ that minimises the collisions of m*p, until m >= n^(1-δ)/2; then m is
// multiplied by ceil(n/m), which puts it in [n, 2n).
std::int64_t choose_modulus(const std::vector<const std::vector<std::int64_t>*>& sets, u64 n, double delta) {
    const double size = std::pow(static_cast<double>(n), delta);
    u64 lo = std::max<u64>(2, static_cast<u64>(std::ceil(size)));
    u64 hi = std::max<u64>(lo + 1, static_cast<u64>(std::floor(2 * size)));
    std::vector<u64> primes;
    auto extend = [&] {  // make sure an unused prime is available
        while (primes.empty()) {
            for (u64 q = lo; q < hi; ++q) {
                if (is_prime(q)) primes.push_back(q);
            }
            lo = hi;
            hi *= 2;
        }
    };
    const double stop = 0.5 * std::pow(static_cast<double>(n), 1 - delta);
    std::int64_t m = 1;
    while (true) {
        extend();
        std::size_t best = 0;
        u64 best_count = std::numeric_limits<u64>::max();
        for (std::size_t i = 0; i < primes.size(); ++i) {
            const u64 count = collisions(sets, m * static_cast<std::int64_t>(primes[i]));
            if (count < best_count) {
                best_count = count;
                best = i;
            }
        }
        m *= static_cast<std::int64_t>(primes[best]);
        primes.erase(primes.begin() + static_cast<std::ptrdiff_t>(best));
        if (static_cast<double>(m) >= stop) break;
    }
    const std::int64_t target = static_cast<std::int64_t>(std::max<u64>(n, 1));
    return m * ((target + m - 1) / m);
}

// Values of a sorted set grouped by residue modulo m.
std::vector<std::vector<std::int64_t>> buckets(const std::vector<std::int64_t>& set, std::int64_t m) {
    std::vector<std::vector<std::int64_t>> out(static_cast<std::size_t>(m));
    for (std::int64_t x : set) out[mod(x, m)].push_back(x);
    return out;
}

// b in B, c in C with c - b = a, for sorted B and C.
std::optional<ThreeSumWitness> find_difference(std::int64_t a, const std::vector<std::int64_t>& B,
                                               const std::vector<std::int64_t>& C) {
    std::size_t i = 0, j = 0;
    while (i < B.size() && j < C.size()) {
        const std::int64_t diff = C[j] - B[i];
        if (diff == a) return ThreeSumWitness{a, B[i], C[j]};
        if (diff < a) {
            ++j;
        } else {
            ++i;
        }
    }
    return std::nullopt;
}

// a in A, b in B with a + b = c, for sorted A and B.
std::optional<ThreeSumWitness> find_sum(std::int64_t c, const std::vector<std::int64_t>& A,
                                        const std::vector<std::int64_t>& B) {
    if (A.empty() || B.empty()) return std::nullopt;
    std::size_t i = 0, j = B.size();
    while (i < A.size() && j > 0) {
        const std::int64_t sum = A[i] + B[j - 1];
        if (sum == c) return ThreeSumWitness{A[i], B[j - 1], c};
        if (sum < c) {
            ++i;
        } else {
            --j;
        }
    }
    return std::nullopt;
}

}  // namespace

std::optional<ConvolutionWitness> convolution_three_sum(const ConvolutionThreeSumInstance& in,
                                                        const ExactTriangleOptions& options, ThreeSumStats* stats) {
    const std::size_t n = in.x.size();
    if (in.y.size() != n || in.z.size() != n) {
        throw std::invalid_argument("convolution_three_sum: x, y, z must have equal lengths");
    }
    std::int64_t max_abs = 0;
    for (const auto* v : {&in.x, &in.y, &in.z}) {
        for (const auto& e : *v) {
            if (!e) continue;
            if (*e > kMaxAbsValue || *e < -kMaxAbsValue) {
                throw std::out_of_range("convolution_three_sum: |value| exceeds kMaxAbsValue");
            }
            max_abs = std::max(max_abs, *e < 0 ? -*e : *e);
        }
    }
    if (n == 0) return std::nullopt;

    // Write i = α s + β and j = γ s + (t - β) with t = (i mod s) + (j mod s), so that
    // i + j = (α + γ) s + t. For each t in [0, 2s - 1), the graph on α in A, γ in B, β in C with
    //   w(α, β) = x[α s + β],   w(γ, β) = y[γ s + t - β],   w(α, γ) = -z[(α + γ) s + t]
    // has a zero triangle exactly when some i, j with these digits have x[i] + y[j] = z[i + j].
    // Absent entries get a weight that lies in no zero triangle.
    const std::size_t s = static_cast<std::size_t>(std::ceil(std::sqrt(static_cast<double>(n))));
    const std::size_t h = (n + s - 1) / s;
    const std::int64_t missing = 2 * max_abs + 1;
    for (std::size_t t = 0; t + 1 < 2 * s; ++t) {
        ExactTriangleInstance tri{Matrix(h, h, missing), Matrix(h, s, missing), Matrix(h, s, missing)};
        bool any_ab = false, any_bc = false, any_ac = false;
        for (std::size_t alpha = 0; alpha < h; ++alpha) {
            for (std::size_t gamma = 0; gamma < h; ++gamma) {
                const std::size_t k = (alpha + gamma) * s + t;
                if (k < n && in.z[k]) {
                    tri.w_ab(alpha, gamma) = -*in.z[k];
                    any_ab = true;
                }
            }
        }
        for (std::size_t beta = 0; beta < s; ++beta) {
            if (beta > t || t - beta >= s) continue;
            for (std::size_t gamma = 0; gamma < h; ++gamma) {
                const std::size_t j = gamma * s + (t - beta);
                if (j < n && in.y[j]) {
                    tri.w_bc(gamma, beta) = *in.y[j];
                    any_bc = true;
                }
            }
        }
        for (std::size_t alpha = 0; alpha < h; ++alpha) {
            for (std::size_t beta = 0; beta < s; ++beta) {
                const std::size_t i = alpha * s + beta;
                if (i < n && in.x[i]) {
                    tri.w_ac(alpha, beta) = *in.x[i];
                    any_ac = true;
                }
            }
        }
        if (!any_ab || !any_bc || !any_ac) continue;

        ExactTriangleStats tri_stats;
        const std::optional<Triangle> found = exact_triangle(tri, options, &tri_stats);
        if (stats) {
            ++stats->triangle_instances;
            stats->triangle_brute_force += tri_stats.brute_force ? 1 : 0;
            stats->sparse_triangle_instances += tri_stats.instances;
            stats->leaves += tri_stats.leaves;
        }
        if (found) return ConvolutionWitness{found->a * s + found->c, found->b * s + (t - found->c)};
    }
    return std::nullopt;
}

std::optional<ThreeSumWitness> three_sum(const std::vector<std::int64_t>& A_in, const std::vector<std::int64_t>& B_in,
                                         const std::vector<std::int64_t>& C_in, const ThreeSumOptions& options,
                                         ThreeSumStats* stats) {
    check_range(A_in);
    check_range(B_in);
    check_range(C_in);
    ThreeSumStats local;
    ThreeSumStats& st = stats ? *stats : local;
    st = ThreeSumStats{};
    const std::vector<std::int64_t> A = sorted_unique(A_in), B = sorted_unique(B_in), C = sorted_unique(C_in);
    if (A.empty() || B.empty() || C.empty()) return std::nullopt;
    const u64 n = std::max({A.size(), B.size(), C.size()});

    // Hash x -> x mod m. Since (a mod m) + (b mod m) is (a + b) mod m or that plus m, a solution
    // a + b = c puts a at position a mod m, b at b mod m, and c at c mod m or c mod m + m of
    // arrays of length 2m.
    const std::int64_t m = choose_modulus({&A, &B, &C}, n, options.hash_delta);
    const u64 R = options.bucket_cap != 0
                      ? options.bucket_cap
                      : std::max<u64>(1, static_cast<u64>(std::llround(std::pow(static_cast<double>(n), options.bucket_exponent))));
    st.modulus = static_cast<u64>(m);
    st.bucket_cap = R;
    const auto bucket_a = buckets(A, m), bucket_b = buckets(B, m), bucket_c = buckets(C, m);

    // Elements of buckets with more than R elements are few, since m has few collisions; check
    // each of them against the whole input in linear time.
    for (const auto& bucket : bucket_a) {
        if (bucket.size() <= R) continue;
        for (std::int64_t a : bucket) {
            ++st.heavy;
            if (auto w = find_difference(a, B, C)) return w;
        }
    }
    for (const auto& bucket : bucket_b) {
        if (bucket.size() <= R) continue;
        for (std::int64_t b : bucket) {
            ++st.heavy;
            if (auto w = find_difference(b, A, C)) return ThreeSumWitness{w->b, b, w->c};
        }
    }
    for (const auto& bucket : bucket_c) {
        if (bucket.size() <= R) continue;
        for (std::int64_t c : bucket) {
            ++st.heavy;
            if (auto w = find_sum(c, A, B)) return w;
        }
    }

    // The light elements: one Convolution 3SUM instance for every (i, j, k) in [R]^3, holding the
    // i-th element of every light bucket of A, the j-th of B, and the k-th of C.
    auto fill = [&](const std::vector<std::vector<std::int64_t>>& bkts, u64 index, bool twice,
                    std::vector<std::optional<std::int64_t>>& out) {
        out.assign(static_cast<std::size_t>(2 * m), std::nullopt);
        bool any = false;
        for (std::size_t r = 0; r < bkts.size(); ++r) {
            if (bkts[r].size() > R || index >= bkts[r].size()) continue;
            out[r] = bkts[r][index];
            if (twice) out[r + static_cast<std::size_t>(m)] = bkts[r][index];
            any = true;
        }
        return any;
    };
    ConvolutionThreeSumInstance conv;
    for (u64 i = 0; i < R; ++i) {
        if (!fill(bucket_a, i, false, conv.x)) continue;
        for (u64 j = 0; j < R; ++j) {
            if (!fill(bucket_b, j, false, conv.y)) continue;
            for (u64 k = 0; k < R; ++k) {
                if (!fill(bucket_c, k, true, conv.z)) continue;
                ++st.convolution_instances;
                if (auto w = convolution_three_sum(conv, options.triangle, &st)) {
                    return ThreeSumWitness{*conv.x[w->i], *conv.y[w->j], *conv.z[w->i + w->j]};
                }
            }
        }
    }
    return std::nullopt;
}

std::optional<ThreeSumWitness> three_sum_quadratic(const std::vector<std::int64_t>& A_in,
                                                   const std::vector<std::int64_t>& B_in,
                                                   const std::vector<std::int64_t>& C_in) {
    check_range(A_in);
    check_range(B_in);
    check_range(C_in);
    const std::vector<std::int64_t> A = sorted_unique(A_in), B = sorted_unique(B_in), C = sorted_unique(C_in);
    for (std::int64_t c : C) {
        if (auto w = find_sum(c, A, B)) return w;
    }
    return std::nullopt;
}

std::optional<std::array<std::size_t, 3>> three_sum_zero(const std::vector<std::int64_t>& nums,
                                                         const ThreeSumOptions& options, ThreeSumStats* stats) {
    check_range(nums);
    if (stats) *stats = ThreeSumStats{};
    const std::size_t n = nums.size();
    if (n < 3) return std::nullopt;
    int bits = 0;
    while ((std::size_t{1} << bits) < n) ++bits;

    for (int b1 = 0; b1 < bits; ++b1) {
        for (int b2 = 0; b2 < bits; ++b2) {
            if (b2 == b1) continue;
            for (std::size_t v = 0; v < 2; ++v) {
                // A: indices with bit b1 = v; B and C: the others, split by bit b2. The search is
                // for a + b = c' with c' = -nums[k].
                std::map<std::int64_t, std::size_t> index_a, index_b, index_c;
                std::vector<std::int64_t> A, B, C;
                for (std::size_t i = 0; i < n; ++i) {
                    if (((i >> b1) & 1) == v) {
                        index_a.emplace(nums[i], i);
                        A.push_back(nums[i]);
                    } else if (((i >> b2) & 1) == 0) {
                        index_b.emplace(nums[i], i);
                        B.push_back(nums[i]);
                    } else {
                        index_c.emplace(-nums[i], i);
                        C.push_back(-nums[i]);
                    }
                }
                if (A.empty() || B.empty() || C.empty()) continue;
                ThreeSumStats round;
                const auto w = three_sum(A, B, C, options, &round);
                if (stats) {
                    stats->modulus = std::max(stats->modulus, round.modulus);
                    stats->bucket_cap = std::max(stats->bucket_cap, round.bucket_cap);
                    stats->heavy += round.heavy;
                    stats->convolution_instances += round.convolution_instances;
                    stats->triangle_instances += round.triangle_instances;
                    stats->triangle_brute_force += round.triangle_brute_force;
                    stats->sparse_triangle_instances += round.sparse_triangle_instances;
                    stats->leaves += round.leaves;
                }
                if (w) {
                    std::array<std::size_t, 3> idx{index_a.at(w->a), index_b.at(w->b), index_c.at(w->c)};
                    std::sort(idx.begin(), idx.end());
                    return idx;
                }
            }
        }
    }
    return std::nullopt;
}

}  // namespace asap
