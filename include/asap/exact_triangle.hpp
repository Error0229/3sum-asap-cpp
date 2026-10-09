#pragma once

// Section 3.2 / Theorem 3.4: the deterministic reduction from Exact Triangle to Lopsided
// All-Edges Sparse Triangle, and Section 3.3 / Theorem 3.5: Exact Triangle solved with it and
// the thin matrix product algorithm of Theorem 2.1.

#include <cstdint>
#include <optional>

#include "asap/matrix.hpp"
#include "asap/thin_product.hpp"

namespace asap {

// A tripartite graph on parts A, B, C with an integer weight on every edge.
// w_ab is |A| x |B|, w_bc is |B| x |C|, and w_ac is |A| x |C|. A zero triangle is (a, b, c) with
// S(a, b, c) = w_ab(a, b) + w_bc(b, c) + w_ac(a, c) = 0. A missing edge can be given a weight
// larger than twice the largest absolute weight of the other edges, so it lies in no zero triangle.
struct ExactTriangleInstance {
    Matrix w_ab;
    Matrix w_bc;
    Matrix w_ac;
};

struct Triangle {
    std::size_t a = 0;
    std::size_t b = 0;
    std::size_t c = 0;
};

struct ExactTriangleOptions {
    // Bound D on the middle part of the Lopsided Sparse Triangle instances; a power of four,
    // at least 16. 0 selects the paper's choice (proof of Theorem 3.5, by Theorem 2.1): the
    // largest power of four with D <= n^(1/18), where instances with D < 16, i.e. with
    // n < 16^18, are solved by brute force.
    std::uint64_t D = 0;
    // Trade-off g between the number of instances and the witness scans, 1 <= g <= sqrt(D).
    // 0 selects the paper's g = ceil(D^(1/36)).
    std::uint64_t g = 0;
    // Options of the thin matrix product that solves each instance (L = 0 is the paper's 19m).
    ThinProductOptions thin;
};

struct ExactTriangleStats {
    bool brute_force = false;     // solved by brute force instead of the reduction
    std::uint64_t D = 0;
    std::uint64_t g = 0;
    std::uint64_t prime = 0;      // the prime p chosen by counting false positives
    std::uint64_t instances = 0;  // Lopsided Sparse Triangle instances solved
    std::uint64_t query_pairs = 0;
    std::uint64_t accepted = 0;   // query pairs with a common neighbour, each scanned for a witness
    std::uint64_t failed_scans = 0;  // scans that met only false positives of p
    std::uint64_t leaves = 0;     // leaves visited by the thin matrix products
};

// Finds a zero triangle, or returns nullopt if there is none.
std::optional<Triangle> exact_triangle(const ExactTriangleInstance& instance, const ExactTriangleOptions& options = {},
                                       ExactTriangleStats* stats = nullptr);

// The O(|A| |B| |C|) exhaustive search.
std::optional<Triangle> exact_triangle_brute_force(const ExactTriangleInstance& instance);

}  // namespace asap
