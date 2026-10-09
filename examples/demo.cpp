// Runs 3SUM three ways on a random instance: the textbook O(n^2) algorithm, the paper's
// pipeline with the paper's parameters, and the same pipeline with small parameters that make
// every layer (down to the thin matrix product) actually run.
//
//   demo [n] [seed] [planted]
//
// All numbers are odd, so a + b is even and never equals c: the search has to rule out every
// triple. With "planted", one solution is added to C.

#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <random>
#include <string>

#include "asap/three_sum.hpp"

namespace {

template <class F>
double seconds(F&& f) {
    const auto start = std::chrono::steady_clock::now();
    f();
    return std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
}

void print(const char* name, const std::optional<asap::ThreeSumWitness>& w, double secs) {
    if (w) {
        std::printf("%-22s %8.3fs  found %lld + %lld = %lld\n", name, secs, static_cast<long long>(w->a),
                    static_cast<long long>(w->b), static_cast<long long>(w->c));
    } else {
        std::printf("%-22s %8.3fs  no solution\n", name, secs);
    }
}

void print_stats(const asap::ThreeSumStats& st) {
    std::printf("    modulus m=%llu, bucket cap R=%llu, heavy elements=%llu, Convolution 3SUM instances=%llu\n",
                static_cast<unsigned long long>(st.modulus), static_cast<unsigned long long>(st.bucket_cap),
                static_cast<unsigned long long>(st.heavy), static_cast<unsigned long long>(st.convolution_instances));
    std::printf("    Exact Triangle instances=%llu (brute-force base cases: %llu), sparse triangle instances=%llu, leaves=%llu\n",
                static_cast<unsigned long long>(st.triangle_instances),
                static_cast<unsigned long long>(st.triangle_brute_force),
                static_cast<unsigned long long>(st.sparse_triangle_instances),
                static_cast<unsigned long long>(st.leaves));
}

}  // namespace

int main(int argc, char** argv) {
    const std::size_t n = argc > 1 ? std::stoul(argv[1]) : 100;
    const unsigned seed = argc > 2 ? static_cast<unsigned>(std::stoul(argv[2])) : 1;
    const bool planted = argc > 3 && std::string(argv[3]) == "planted";
    std::mt19937_64 rng(seed);
    const auto bound = static_cast<std::int64_t>(n) * static_cast<std::int64_t>(n);
    std::uniform_int_distribution<std::int64_t> value(-bound, bound);
    std::vector<std::int64_t> A(n), B(n), C(n);
    for (auto* v : {&A, &B, &C}) {
        for (auto& x : *v) x = 2 * value(rng) + 1;
    }
    if (planted && n > 0) C.push_back(A[rng() % n] + B[rng() % n]);
    std::printf("3SUM on |A| = |B| = %zu, |C| = %zu random odd integers in [-%lld, %lld]%s\n\n", n, C.size(),
                static_cast<long long>(2 * bound + 1), static_cast<long long>(2 * bound + 1),
                planted ? ", one solution planted" : ", no solution");

    std::optional<asap::ThreeSumWitness> w;
    print("textbook O(n^2)", w, seconds([&] { w = asap::three_sum_quadratic(A, B, C); }));

    asap::ThreeSumStats st;
    print("paper parameters", w, seconds([&] { w = asap::three_sum(A, B, C, {}, &st); }));
    print_stats(st);
    std::printf("    (the paper's D <= n^(1/18) with D >= 16 needs Exact Triangle instances with n >= 16^18,\n"
                "     so every instance of feasible size is a brute-force base case)\n");

    asap::ThreeSumOptions small;
    small.triangle.D = 16;
    small.triangle.g = 1;
    small.triangle.thin.L = 4;
    small.bucket_cap = 3;
    print("small parameters", w, seconds([&] { w = asap::three_sum(A, B, C, small, &st); }));
    print_stats(st);
    std::printf("    (D = 16, g = 1, L = 4, R = 3: every layer runs, but the savings need L >= 19 log4(D))\n");
    return 0;
}
