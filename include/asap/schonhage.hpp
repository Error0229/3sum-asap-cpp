#pragma once

// Schönhage's identity for an inner product and an outer product (Section 2.2, Lemma 2.1).
//
//   G = sum_{i,j<=3} x_i y_j z_ij + (sum_{i,j<=2} p_ij q_ij) z_0
//
// is computed, up to the harmless error E, by the ten terms
//
//   P_ij : (x_i + p̂_ij)(y_j + q̂_ij)(z_ij + z_0)        for 1 <= i,j <= 3
//   P_0  : -(x_1 + x_2 + x_3)(y_1 + y_2 + y_3) z_0
//
// where p̂ is the 2x2 matrix of p's extended so that its columns sum to zero, and q̂ the
// 2x2 matrix of q's extended so that its rows sum to zero.

#include <array>

namespace asap::schonhage {

// Variable numbering. Left variables x1,x2,x3,p11,p12,p21,p22 are 0..6, and right variables
// y1,y2,y3,q11,q12,q21,q22 are 0..6 in the same pattern. Indices 0..2 are the outer
// variables, 3..6 the inner ones (p_ij and q_ij are 3 + 2(i-1) + (j-1)).
inline constexpr int kVars = 7;
inline constexpr int kOuterVars = 3;
inline constexpr int kInnerVars = 4;

// Terms P_ij are 3(i-1) + (j-1), and P_0 is 9. The output variables z_ij and z_0 use the
// same numbering: every term contributes to z_0, and only P_ij contributes to z_ij.
inline constexpr int kTerms = 10;
inline constexpr int kZero = 9;  // P_0, and the output variable z_0

struct Tables {
    // phi[λ][s] and psi[λ][t]: the coefficients φ_λ(s), ψ_λ(t) in {-1, 0, 1}.
    std::array<std::array<int, kVars>, kTerms> phi{};
    std::array<std::array<int, kVars>, kTerms> psi{};
};

constexpr Tables make_tables() {
    Tables t;
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            const int term = 3 * i + j;
            // φ_{P_ij} = x_i + p̂_ij
            t.phi[term][i] += 1;
            if (j < 2) {
                if (i < 2) {
                    t.phi[term][3 + 2 * i + j] += 1;  // p_ij
                } else {
                    t.phi[term][3 + j] -= 1;      // -p_1j
                    t.phi[term][3 + 2 + j] -= 1;  // -p_2j
                }
            }
            // ψ_{P_ij} = y_j + q̂_ij
            t.psi[term][j] += 1;
            if (i < 2) {
                if (j < 2) {
                    t.psi[term][3 + 2 * i + j] += 1;  // q_ij
                } else {
                    t.psi[term][3 + 2 * i] -= 1;      // -q_i1
                    t.psi[term][3 + 2 * i + 1] -= 1;  // -q_i2
                }
            }
        }
    }
    // φ_{P_0} = -(x1 + x2 + x3), ψ_{P_0} = y1 + y2 + y3
    for (int s = 0; s < kOuterVars; ++s) {
        t.phi[kZero][s] = -1;
        t.psi[kZero][s] = 1;
    }
    return t;
}

inline constexpr Tables kTables = make_tables();

// Whether term λ contributes to output variable z, i.e. z appears in χ_λ.
constexpr bool contributes(int term, int z) { return term == z || z == kZero; }

}  // namespace asap::schonhage
