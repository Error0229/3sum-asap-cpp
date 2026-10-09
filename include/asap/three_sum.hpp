#pragma once

// Section 3.4 / Theorem 3.7: 3SUM through the known deterministic reductions to Exact Triangle
// (Theorem 3.6(a)), solved by Theorem 3.5:
//
//   3SUM --[Chan-He: hashing x mod m, m a product of primes chosen greedily for few collisions]-->
//   Convolution 3SUM --[Vassilevska Williams-Williams: split indices into sqrt(n) digits]-->
//   Exact Triangle --[Theorem 3.4]--> Lopsided All-Edges Sparse Triangle --[Theorem 2.1]

#include <array>
#include <cstdint>
#include <optional>
#include <vector>

#include "asap/exact_triangle.hpp"

namespace asap {

// Input numbers must have absolute value at most this, so that no sum below overflows.
inline constexpr std::int64_t kMaxAbsValue = std::int64_t{1} << 59;

struct ThreeSumOptions {
    // Options of the Exact Triangle algorithm (and of the thin matrix products below it).
    ExactTriangleOptions triangle;
    // The hashing modulus m in [n, 2n) is a product of primes of size about n^hash_delta,
    // chosen one at a time to minimise the number of collisions.
    double hash_delta = 0.25;
    // Buckets with more than R elements are "heavy"; their elements are checked one by one in
    // O(n) time each, and the light ones go into R^3 Convolution 3SUM instances. 0 selects
    // R = max(1, round(n^bucket_exponent)).
    std::uint64_t bucket_cap = 0;
    // The R^3 instances must cost less than what the algorithm saves on each of them. By
    // Theorem 2.1, Exact Triangle saves n^(1/648), so Convolution 3SUM (sqrt(n) instances on
    // sqrt(n) vertices) saves n^(1/1296) ≈ n^0.00077, and R^3 = n^0.00075 stays below it.
    double bucket_exponent = 0.00025;
};

struct ThreeSumStats {
    std::uint64_t modulus = 0;
    std::uint64_t bucket_cap = 0;
    std::uint64_t heavy = 0;                      // elements checked one by one
    std::uint64_t convolution_instances = 0;
    std::uint64_t triangle_instances = 0;
    std::uint64_t triangle_brute_force = 0;       // of which solved by brute force (too small)
    std::uint64_t sparse_triangle_instances = 0;  // Lopsided Sparse Triangle instances
    std::uint64_t leaves = 0;                     // leaves visited by the thin matrix products
};

// a in A, b in B, c in C with a + b = c.
struct ThreeSumWitness {
    std::int64_t a = 0;
    std::int64_t b = 0;
    std::int64_t c = 0;
};

// Convolution 3SUM: given x, y, z of equal length n, whose entries may be absent, find i, j
// with i + j < n and x[i] + y[j] = z[i + j].
struct ConvolutionThreeSumInstance {
    std::vector<std::optional<std::int64_t>> x;
    std::vector<std::optional<std::int64_t>> y;
    std::vector<std::optional<std::int64_t>> z;
};

struct ConvolutionWitness {
    std::size_t i = 0;
    std::size_t j = 0;
};

// 3SUM: finds a in A, b in B, c in C with a + b = c, or returns nullopt.
std::optional<ThreeSumWitness> three_sum(const std::vector<std::int64_t>& A, const std::vector<std::int64_t>& B,
                                         const std::vector<std::int64_t>& C, const ThreeSumOptions& options = {},
                                         ThreeSumStats* stats = nullptr);

// The textbook O(n^2) algorithm (sort, then two pointers for every c).
std::optional<ThreeSumWitness> three_sum_quadratic(const std::vector<std::int64_t>& A,
                                                   const std::vector<std::int64_t>& B,
                                                   const std::vector<std::int64_t>& C);

// The single-list form: distinct indices i < j < k with nums[i] + nums[j] + nums[k] = 0.
// Reduces to O(log^2 n) calls of three_sum on disjoint index classes A, B, C: for distinct
// i, j, k some bit separates one of them from the other two, and another bit separates those two.
std::optional<std::array<std::size_t, 3>> three_sum_zero(const std::vector<std::int64_t>& nums,
                                                         const ThreeSumOptions& options = {},
                                                         ThreeSumStats* stats = nullptr);

// Convolution 3SUM through Exact Triangle. Adds to *stats rather than resetting it.
std::optional<ConvolutionWitness> convolution_three_sum(const ConvolutionThreeSumInstance& instance,
                                                        const ExactTriangleOptions& options = {},
                                                        ThreeSumStats* stats = nullptr);

}  // namespace asap
