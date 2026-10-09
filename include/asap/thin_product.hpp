#pragma once

// Section 2 / Theorem 2.1: computing a given set W of entries of a thin matrix product XY,
// X in Z^{N x D} and Y in Z^{D x N}, by Schönhage's identity applied recursively (Section 2.3),
// tiled over the product (Section 2.3.4), with the encodings shared among tiles (Section 2.4.1)
// and the decoding recursion pruned to the wanted entries (Section 2.4.2).
//
// The paper sets L = 19m for D = 4^m and assumes N >= D^18 for its running-time bound. The
// algorithm is correct for every L >= m, so L is a parameter here; small L makes it runnable.

#include <cstdint>
#include <utility>
#include <vector>

#include "asap/matrix.hpp"

namespace asap {

struct ThinProductOptions {
    // Recursion depth L (number of levels of Schönhage's identity). 0 selects the paper's
    // L = 19m. Must satisfy L >= m, where D = 4^m is the inner dimension padded to a power of four.
    int L = 0;
    // Refuse to run when the shared encodings would hold more than this many numbers
    // (each encoding has 10^L of them). The paper's parameters need about 10^(19m).
    std::uint64_t max_encoded_values = std::uint64_t{1} << 27;
};

struct ThinProductStats {
    int m = 0;                     // D = 4^m
    int L = 0;                     // recursion depth
    std::uint64_t D = 0;           // padded inner dimension
    std::uint64_t N0 = 0;          // 3^(L-m): side of one block product
    std::uint64_t K = 0;           // C(L, m): products computed by one run of the recursion
    std::uint64_t K0 = 0;          // floor(sqrt(K)): a tile is a K0 x K0 grid of block products
    std::uint64_t encodings = 0;   // band encodings computed (each 10^L numbers)
    std::uint64_t tiles = 0;       // tiles with at least one wanted entry
    std::uint64_t leaves = 0;      // leaves visited by the pruned recursion, summed over tiles
    std::uint64_t calls = 0;       // calls of the pruned recursion, summed over tiles
    // For each tile: (|U|, |Leaves(U)|), the number of wanted output strings and leaves visited.
    std::vector<std::pair<std::uint64_t, std::uint64_t>> per_tile;
};

// Returns (XY)[I, J] for every (I, J) in W, in the order of W.
//
// Arithmetic is done modulo 2^64; every step is a ring operation, so the results are exact
// whenever the true entries fit in an int64_t. X.cols must equal Y.rows; it is padded with
// zeros to a power of four D >= 4.
std::vector<std::int64_t> thin_product_entries(const Matrix& X, const Matrix& Y,
                                               const std::vector<Position>& W,
                                               const ThinProductOptions& options = {},
                                               ThinProductStats* stats = nullptr);

// Definitions 3.1 and 3.2 (Lopsided All-Edges Sparse Triangle, and its counting version).
// X is the 0/1 biadjacency matrix between A and the middle part M, and Y between M and B.
// For every query pair (a, b) in W, returns the number of common neighbours of a and b in M.
inline std::vector<std::int64_t> count_lopsided_sparse_triangles(const Matrix& X, const Matrix& Y,
                                                                 const std::vector<Position>& W,
                                                                 const ThinProductOptions& options = {},
                                                                 ThinProductStats* stats = nullptr) {
    return thin_product_entries(X, Y, W, options, stats);
}

}  // namespace asap
