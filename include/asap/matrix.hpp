#pragma once

#include <cstddef>
#include <cstdint>
#include <vector>

namespace asap {

// Dense row-major integer matrix.
struct Matrix {
    std::size_t rows = 0;
    std::size_t cols = 0;
    std::vector<std::int64_t> data;

    Matrix() = default;
    Matrix(std::size_t r, std::size_t c, std::int64_t fill = 0) : rows(r), cols(c), data(r * c, fill) {}

    std::int64_t& operator()(std::size_t i, std::size_t j) { return data[i * cols + j]; }
    std::int64_t operator()(std::size_t i, std::size_t j) const { return data[i * cols + j]; }
};

// A position (row, column) of a matrix product.
struct Position {
    std::size_t row = 0;
    std::size_t col = 0;
};

}  // namespace asap
