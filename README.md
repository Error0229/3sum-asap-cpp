<img src="assets/icon.svg" width="64" height="64" alt="">

# 3sum-asap-cpp

A C++20 implementation, using only the standard library, of the deterministic truly subquadratic
3SUM algorithm of Josh Alman and Virginia Vassilevska Williams,
[*Truly Subquadratic 3SUM and Truly Subcubic APSP via Triangles in Sparse Lopsided Graphs*](https://arxiv.org/abs/2610.06783)
(arXiv:2610.06783).

## The pipeline

The 3SUM algorithm is a chain of reductions that ends in the paper's new algorithm for computing a
sparse set of entries of a thin matrix product:

| Step | Paper | Code |
| --- | --- | --- |
| 3SUM → Convolution 3SUM: hash `x mod m`, with `m` a product of primes chosen greedily for few collisions; heavy buckets are checked one element at a time | Thm 3.6(a), cited from Chan–He [CH20] | `src/three_sum.cpp`: `three_sum` |
| Convolution 3SUM → Exact Triangle: write indices as two digits of size about √n | Thm 3.6(a), cited from [VW13] | `src/three_sum.cpp`: `convolution_three_sum` |
| Exact Triangle → Lopsided All-Edges Sparse Triangle: hash the weights mod a prime `p` (chosen by counting false positives with Strassen over Z[x]/(x^p − 1)), build the instances, scan for witnesses | Thm 3.4, Thm 3.5 | `src/exact_triangle.cpp` |
| Lopsided Sparse Triangle → wanted entries of a thin product `XY` | Def 3.1–3.2, Cor 3.3 | `count_lopsided_sparse_triangles` |
| Wanted entries of `XY`: Schönhage's identity applied recursively, Coppersmith-style tiling, shared encodings, pruned decoding | Section 2, Thm 2.1 | `include/asap/schonhage.hpp`, `src/thin_product.cpp` |

## Build

```sh
cmake -S . -B build -G Ninja
cmake --build build
ctest --test-dir build          # randomized cross-checks of every layer against brute force
./build/asap3sum_demo 100       # add "planted" as the third argument to plant a solution
```

Any C++20 compiler works, for example: `g++ -std=c++20 -O2 -Iinclude src/*.cpp examples/demo.cpp`.

## Usage

```cpp
#include "asap/three_sum.hpp"

// a in A, b in B, c in C with a + b = c
std::optional<asap::ThreeSumWitness> w = asap::three_sum(A, B, C);

// distinct indices i < j < k with nums[i] + nums[j] + nums[k] = 0
std::optional<std::array<std::size_t, 3>> idx = asap::three_sum_zero(nums);
```

The lower layers are public too: `exact_triangle`, `convolution_three_sum`, and
`thin_product_entries` (Theorem 2.1 on its own).

## Parameters, and what to expect in practice

The algorithm's savings are asymptotic, and with the paper's constants they appear only at sizes
no computer can reach:

- **Paper parameters (the defaults).** Theorem 3.5 takes `D` to be the largest power of four with
  `D ≤ n^(1/18)`, and needs `D ≥ 16`. So an Exact Triangle instance needs `n ≥ 16^18 ≈ 4.7·10^21`
  vertices per part, which means a 3SUM input of about `10^43` numbers. Smaller instances are solved
  by brute force, as in the paper's proof. Even past that point, `L = 19m` makes each encoding hold
  `10^38` numbers. With the defaults, then, every feasible input runs the reductions down to the
  brute-force base case.
- **Small parameters** (`options.triangle.D = 16`, `.g = 1`, `.thin.L = 3..5`) make every layer
  run, down to the pruned Schönhage recursion. The demo does this. The results are correct but slower
  than the textbook `O(n²)` algorithm. Theorem 2.1 saves work only when the number of leaves of
  order `d` decays, i.e. when `L` is well above `10m` (Section 2.4.3). With `D ≥ 4`, that already
  means encodings of more than `10^11` numbers per band.

`ThinProductOptions::max_encoded_values` (default `2^27`) makes the library refuse encodings it
cannot hold, rather than try to allocate them.

## Notes on the implementation

- Arithmetic in the thin product is done modulo `2^64`. Every step is a ring operation, so the
  results are exact whenever the true entries fit in an `int64_t` (the paper uses `O(log N)`-bit words).
- Only bands that contain wanted entries are encoded. The paper encodes all of them, but the
  others are never read.
- Only Theorem 2.1 is implemented, with saving `D^(1/18)`, which gives the paper's
  `n^(2 − 1/1296 + o(1))` bound for 3SUM. Section 4's data structure, with saving `D^0.063`, is not.
- The paper cites the 3SUM → Exact Triangle reduction without restating it. This code follows
  Chan–He's deterministic hashing as described by Fischer, Kaliciak and Polak
  ([arXiv:2310.12913](https://arxiv.org/abs/2310.12913), Lemma 3.1), with Pătraşcu's bucket
  construction and the standard digit-splitting reduction to Exact Triangle. The details may differ
  from the original papers.
- Input numbers must satisfy `|x| ≤ 2^59` (`asap::kMaxAbsValue`).
