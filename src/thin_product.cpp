#include "asap/thin_product.hpp"

#include <algorithm>
#include <array>
#include <cstring>
#include <limits>
#include <map>
#include <numeric>
#include <stdexcept>
#include <string>

#include "asap/schonhage.hpp"

namespace asap {
namespace {

using u64 = std::uint64_t;
namespace sch = schonhage;

u64 checked_pow(u64 base, int exp) {
    u64 r = 1;
    for (int i = 0; i < exp; ++i) {
        if (r > std::numeric_limits<u64>::max() / base) {
            throw std::length_error("thin_product_entries: 10^L overflows; L is too large");
        }
        r *= base;
    }
    return r;
}

u64 binomial(int n, int k) {
    u64 r = 1;
    for (int i = 0; i < k; ++i) r = r * static_cast<u64>(n - i) / static_cast<u64>(i + 1);
    return r;
}

u64 isqrt(u64 x) {
    u64 r = 0;
    while ((r + 1) * (r + 1) <= x) ++r;
    return r;
}

// Strings of variables are indexed with level 1 as the most significant digit: a left or right
// string u is the base-7 number sum_l u_l 7^(L-l), and a leaf τ the base-10 number
// sum_l τ_l 10^(L-l). Output strings are kept as byte strings over the digits 0..9.

// The shape of one tile (Sections 2.3.3 and 2.3.4). The block product in grid cell (r, c) of a
// tile uses the subset Q number r*K0 + c, the same in every tile.
struct Layout {
    int m = 0;
    int L = 0;
    u64 D = 0;
    u64 N0 = 0;
    u64 K = 0;
    u64 K0 = 0;
    u64 band = 0;          // K0 * N0 rows (or columns) per band
    u64 left_size = 0;     // 7^L
    std::vector<std::vector<int>> inner_levels;  // the levels of Q, increasing
    std::vector<std::vector<int>> outer_levels;  // the other levels, increasing
    // Index contribution, in a 7^L array, of the outer part of a string (by its row of X_Q, or
    // column of Y_Q, in [0, N0)) and of its inner part (by its column of X_Q, or row of Y_Q,
    // in [0, D)). Left and right strings share these, since x_i, y_j are 0..2 and p_ij, q_ij 3..6.
    std::vector<std::vector<u64>> outer_offset;
    std::vector<std::vector<u64>> inner_offset;
};

Layout make_layout(int m, int L) {
    Layout lay;
    lay.m = m;
    lay.L = L;
    lay.D = checked_pow(4, m);
    lay.N0 = checked_pow(3, L - m);
    lay.K = binomial(L, m);
    lay.K0 = isqrt(lay.K);
    lay.band = lay.K0 * lay.N0;
    lay.left_size = checked_pow(7, L);

    std::vector<u64> weight(L);
    for (int level = 0; level < L; ++level) weight[level] = checked_pow(7, L - 1 - level);

    // The first K0^2 subsets of size m, in lexicographic order.
    std::vector<int> comb(m);
    std::iota(comb.begin(), comb.end(), 0);
    for (u64 q = 0; q < lay.K0 * lay.K0; ++q) {
        std::vector<int> outer;
        for (int level = 0, k = 0; level < L; ++level) {
            if (k < m && comb[k] == level) {
                ++k;
            } else {
                outer.push_back(level);
            }
        }
        std::vector<u64> outer_off(lay.N0), inner_off(lay.D);
        for (u64 row = 0; row < lay.N0; ++row) {
            u64 off = 0, rest = row;
            for (int k = L - m - 1; k >= 0; --k, rest /= 3) off += (rest % 3) * weight[outer[k]];
            outer_off[row] = off;
        }
        for (u64 col = 0; col < lay.D; ++col) {
            u64 off = 0, rest = col;
            for (int k = m - 1; k >= 0; --k, rest /= 4) off += (sch::kOuterVars + rest % 4) * weight[comb[k]];
            inner_off[col] = off;
        }
        lay.inner_levels.push_back(comb);
        lay.outer_levels.push_back(std::move(outer));
        lay.outer_offset.push_back(std::move(outer_off));
        lay.inner_offset.push_back(std::move(inner_off));

        int i = m - 1;
        while (i >= 0 && comb[i] == L - m + i) --i;
        if (i < 0) break;
        ++comb[i];
        for (int k = i + 1; k < m; ++k) comb[k] = comb[k - 1] + 1;
    }
    return lay;
}

// The output string of the entry in row `row` and column `col` of the block product X_Q Y_Q
// (Section 2.3.3): z_0 at the levels of Q, and z_ij at the others, where the x_i form the row
// and the y_j form the column.
std::string output_string(const Layout& lay, u64 q, u64 row, u64 col) {
    std::string s(lay.L, static_cast<char>(sch::kZero));
    const auto& outer = lay.outer_levels[q];
    for (int k = lay.L - lay.m - 1; k >= 0; --k, row /= 3, col /= 3) {
        s[outer[k]] = static_cast<char>(3 * (row % 3) + col % 3);
    }
    return s;
}

// The left input array of a row band: X_Q is the row block r of the band for the subset Q of
// grid cell (r, c), and the array is zero at every other string (Sections 2.3.3 and 2.3.4).
std::vector<u64> left_array(const Matrix& X, const Layout& lay, u64 band) {
    std::vector<u64> a(lay.left_size, 0);
    for (u64 r = 0; r < lay.K0; ++r) {
        for (u64 c = 0; c < lay.K0; ++c) {
            const u64 q = r * lay.K0 + c;
            for (u64 row = 0; row < lay.N0; ++row) {
                const u64 I = band * lay.band + r * lay.N0 + row;
                if (I >= X.rows) break;
                const u64 base = lay.outer_offset[q][row];
                for (u64 col = 0; col < X.cols; ++col) {
                    a[base + lay.inner_offset[q][col]] = static_cast<u64>(X(I, col));
                }
            }
        }
    }
    return a;
}

// The right input array of a column band: Y_Q is the column block c of the band.
std::vector<u64> right_array(const Matrix& Y, const Layout& lay, u64 band) {
    std::vector<u64> b(lay.left_size, 0);
    for (u64 r = 0; r < lay.K0; ++r) {
        for (u64 c = 0; c < lay.K0; ++c) {
            const u64 q = r * lay.K0 + c;
            for (u64 col = 0; col < lay.N0; ++col) {
                const u64 J = band * lay.band + c * lay.N0 + col;
                if (J >= Y.cols) break;
                const u64 base = lay.outer_offset[q][col];
                for (u64 row = 0; row < Y.rows; ++row) {
                    b[base + lay.inner_offset[q][row]] = static_cast<u64>(Y(row, J));
                }
            }
        }
    }
    return b;
}

// The encoding of an input array (Section 2.4.1): the 10^L numbers Φ_τ(a) (or Ψ_τ(b)), computed
// level by level as the encoding steps of the recursion do (Yates' algorithm). Level l turns an
// array of shape 10^(l-1) x 7 x 7^(L-l) into one of shape 10^(l-1) x 10 x 7^(L-l), forming
// A_λ = sum_s φ_λ(s) a_s for every term λ.
std::vector<u64> encode(std::vector<u64> cur, int L, const std::array<std::array<int, sch::kVars>, sch::kTerms>& coef) {
    u64 head = 1;
    u64 tail = cur.size();
    for (int level = 0; level < L; ++level) {
        tail /= sch::kVars;
        std::vector<u64> next(head * sch::kTerms * tail, 0);
        for (u64 h = 0; h < head; ++h) {
            for (int term = 0; term < sch::kTerms; ++term) {
                u64* out = next.data() + (h * sch::kTerms + term) * tail;
                for (int s = 0; s < sch::kVars; ++s) {
                    const int c = coef[term][s];
                    if (c == 0) continue;
                    const u64* in = cur.data() + (h * sch::kVars + s) * tail;
                    if (c > 0) {
                        for (u64 t = 0; t < tail; ++t) out[t] += in[t];
                    } else {
                        for (u64 t = 0; t < tail; ++t) out[t] -= in[t];
                    }
                }
            }
        }
        cur.swap(next);
        head *= sch::kTerms;
    }
    return cur;
}

// The pruned recursion Pruned (Section 2.4.2). A call at vertex τ_1...τ_k receives the sorted set
// S of the output-string suffixes wanted from it, and returns the array of its outputs on S.
class PrunedRecursion {
public:
    PrunedRecursion(const u64* left_encoding, const u64* right_encoding)
        : ea_(left_encoding), eb_(right_encoding) {}

    // S holds `count` distinct strings of length `len`, sorted, stored back to back.
    std::vector<u64> run(u64 vertex, int len, const std::vector<std::uint8_t>& S, std::size_t count) {
        ++calls;
        // (1) At a leaf, multiply the two encoded numbers.
        if (len == 0) {
            ++leaves;
            return {ea_[vertex] * eb_[vertex]};
        }
        const int sub = len - 1;
        // The slices S_z are contiguous, since S is sorted: S_z is [first[z], first[z+1]).
        std::array<std::size_t, sch::kTerms + 1> first{};
        std::size_t scan = 0;
        for (int z = 0; z <= sch::kTerms; ++z) {
            while (scan < count && S[scan * len] < z) ++scan;
            first[z] = scan;
        }
        const std::size_t z0_begin = first[sch::kZero];
        const std::size_t z0_count = count - z0_begin;
        auto suffix = [&](std::size_t i) { return S.data() + i * len + 1; };

        std::array<std::vector<u64>, sch::kTerms> C;
        // pos_z0[λ][k]: position, in S_λ, of the suffix of the k-th string of the slice S_{z_0}.
        std::array<std::vector<std::size_t>, sch::kTerms> pos_z0;
        // pos_own[i]: for a string z_ij w' of S, the position of w' in S_{P_ij}.
        std::vector<std::size_t> pos_own(z0_begin);

        for (int term = 0; term < sch::kTerms; ++term) {
            // (2) S_{P_ij} = S_{z_ij} ∪ S_{z_0} and S_{P_0} = S_{z_0}, merged without duplicates.
            std::vector<std::uint8_t> T;
            std::size_t t_count = 0;
            pos_z0[term].resize(z0_count);
            std::size_t i = term == sch::kZero ? first[sch::kZero] : first[term];
            const std::size_t i_end = term == sch::kZero ? first[sch::kZero] : first[term + 1];
            std::size_t k = 0;
            while (i < i_end || k < z0_count) {
                int cmp;
                if (i == i_end) {
                    cmp = 1;
                } else if (k == z0_count) {
                    cmp = -1;
                } else {
                    cmp = std::memcmp(suffix(i), suffix(z0_begin + k), sub);
                }
                const std::uint8_t* src = cmp <= 0 ? suffix(i) : suffix(z0_begin + k);
                T.insert(T.end(), src, src + sub);
                if (cmp <= 0) pos_own[i++] = t_count;
                if (cmp >= 0) pos_z0[term][k++] = t_count;
                ++t_count;
            }
            // (3) Recurse only into the children from which some output is wanted.
            if (t_count > 0) C[term] = run(vertex * sch::kTerms + term, sub, T, t_count);
        }

        // (4) Decode: c_{z_ij} = C_{P_ij}, and c_{z_0} = sum over all terms of C_λ.
        std::vector<u64> out(count);
        for (std::size_t i = 0; i < z0_begin; ++i) out[i] = C[S[i * len]][pos_own[i]];
        for (std::size_t k = 0; k < z0_count; ++k) {
            u64 sum = 0;
            for (int term = 0; term < sch::kTerms; ++term) sum += C[term][pos_z0[term][k]];
            out[z0_begin + k] = sum;
        }
        return out;
    }

    u64 leaves = 0;
    u64 calls = 0;

private:
    const u64* ea_;
    const u64* eb_;
};

}  // namespace

std::vector<std::int64_t> thin_product_entries(const Matrix& X, const Matrix& Y, const std::vector<Position>& W,
                                               const ThinProductOptions& options, ThinProductStats* stats) {
    if (X.cols != Y.rows) throw std::invalid_argument("thin_product_entries: X.cols must equal Y.rows");
    for (const Position& w : W) {
        if (w.row >= X.rows || w.col >= Y.cols) throw std::out_of_range("thin_product_entries: position outside XY");
    }

    int m = 1;
    u64 D = 4;
    while (D < X.cols) {
        D *= 4;
        ++m;
    }
    const int L = options.L == 0 ? 19 * m : options.L;
    if (L < m) throw std::invalid_argument("thin_product_entries: L must be at least m = log4(D)");
    ThinProductStats local;
    ThinProductStats& st = stats ? *stats : local;
    st = ThinProductStats{};
    st.m = m;
    st.L = L;
    st.D = D;
    // Nothing is wanted, so nothing is encoded, however large 10^L is.
    if (W.empty()) return {};

    const u64 leaves_per_encoding = checked_pow(10, L);
    if (leaves_per_encoding > options.max_encoded_values) {
        throw std::length_error("thin_product_entries: an encoding of 10^L numbers exceeds max_encoded_values");
    }

    const Layout lay = make_layout(m, L);
    st.N0 = lay.N0;
    st.K = lay.K;
    st.K0 = lay.K0;

    std::vector<std::int64_t> result(W.size());

    // Locate the tile and the output string of every wanted position (Section 2.4.4).
    struct Wanted {
        u64 row_band;
        u64 col_band;
        std::string str;
        std::size_t index;
    };
    std::vector<Wanted> wanted;
    wanted.reserve(W.size());
    for (std::size_t i = 0; i < W.size(); ++i) {
        const u64 I = W[i].row, J = W[i].col;
        const u64 r = (I % lay.band) / lay.N0, c = (J % lay.band) / lay.N0;
        wanted.push_back({I / lay.band, J / lay.band,
                          output_string(lay, r * lay.K0 + c, I % lay.N0, J % lay.N0), i});
    }
    std::sort(wanted.begin(), wanted.end(), [](const Wanted& x, const Wanted& y) {
        if (x.row_band != y.row_band) return x.row_band < y.row_band;
        if (x.col_band != y.col_band) return x.col_band < y.col_band;
        return x.str < y.str;
    });

    // Every band's encoding is computed once and shared by all its tiles (Section 2.4.1). Only
    // the bands that contain wanted positions are encoded; the others would never be read.
    std::map<u64, std::vector<u64>> col_encodings;
    for (const Wanted& w : wanted) col_encodings.emplace(w.col_band, std::vector<u64>{});
    if ((col_encodings.size() + 1) > options.max_encoded_values / leaves_per_encoding) {
        throw std::length_error("thin_product_entries: the shared encodings exceed max_encoded_values");
    }
    for (auto& [band, enc] : col_encodings) {
        enc = encode(right_array(Y, lay, band), L, sch::kTables.psi);
        ++st.encodings;
    }

    std::vector<u64> row_encoding;
    u64 encoded_row_band = std::numeric_limits<u64>::max();
    for (std::size_t begin = 0; begin < wanted.size();) {
        std::size_t end = begin;
        while (end < wanted.size() && wanted[end].row_band == wanted[begin].row_band &&
               wanted[end].col_band == wanted[begin].col_band) {
            ++end;
        }
        if (wanted[begin].row_band != encoded_row_band) {
            encoded_row_band = wanted[begin].row_band;
            row_encoding = encode(left_array(X, lay, encoded_row_band), L, sch::kTables.phi);
            ++st.encodings;
        }

        // The set W_T of the tile's wanted output strings, sorted and without duplicates.
        std::vector<std::uint8_t> S;
        std::vector<std::size_t> slot(end - begin);
        std::size_t count = 0;
        for (std::size_t i = begin; i < end; ++i) {
            if (i == begin || wanted[i].str != wanted[i - 1].str) {
                S.insert(S.end(), wanted[i].str.begin(), wanted[i].str.end());
                ++count;
            }
            slot[i - begin] = count - 1;
        }

        PrunedRecursion pruned(row_encoding.data(), col_encodings.at(wanted[begin].col_band).data());
        const std::vector<u64> values = pruned.run(0, L, S, count);
        for (std::size_t i = begin; i < end; ++i) {
            result[wanted[i].index] = static_cast<std::int64_t>(values[slot[i - begin]]);
        }

        ++st.tiles;
        st.leaves += pruned.leaves;
        st.calls += pruned.calls;
        st.per_tile.emplace_back(count, pruned.leaves);
        begin = end;
    }
    return result;
}

}  // namespace asap
