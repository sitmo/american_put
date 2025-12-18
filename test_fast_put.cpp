// MIT License
//
// Copyright (c) Thijs van den Berg
//
// This code implements the method described in:
// Leif Andersen, Mark Lake, and Dimitri Offengenden,
// "High-performance American option pricing",
// Journal of Computational Finance (2015).
// https://ideas.repec.org/a/rsk/journ0/2464632.html
//
// -----------------------------------------------------------------------------
// Parameter mapping to the paper (very important)
//
// The paper reports configurations as (l, m, n) and p:
//
//   - n : number of Chebyshev nodes used to represent the early-exercise boundary
//         (time discretization / Chebyshev expansion degree).
//
//   - l : number of tanh–sinh quadrature nodes used in the boundary fixed-point
//         integrals (the integrals in the fixed-point operator, e.g. eqs. 21–24).
//
//   - m : number of outer fixed-point iterations used to converge the boundary.
//
//   - p : number of tanh–sinh quadrature nodes used in the price integral
//         (for the early-exercise premium term, e.g. Vee(T)).
//
// In THIS codebase these correspond to:
//
//   Paper n  -> CHEBYSHEV_NODES
//   Paper l  -> QUAD_NODES_FP
//   Paper m  -> FIXED_POINT_ITER
//   Paper p  -> QUAD_NODES_PRICE
//
// Example (from paper): FP-A (l,m,n) = (15,3,7), p = 41
//   => QUAD_NODES_FP=15, FIXED_POINT_ITER=3, CHEBYSHEV_NODES=7, QUAD_NODES_PRICE=41
// -----------------------------------------------------------------------------

#include "fast_put.hpp"
#include <cstdlib>
#include <cstring>

// -----------------------------------------------------------------------------
// main
// -----------------------------------------------------------------------------
int main(int argc, char *argv[]) {
    RealType S = 100;
    RealType K = 100;
    RealType T = 1.0;
    RealType sigma = 0.25;
    RealType r = 0.05;
    RealType q = 0.05;
    int fixed_point_iterations = FIXED_POINT_ITER;  // Default from constants

    int arg_pairs = (argc - 1) / 2;
    for (int i = 0; i < arg_pairs; ++i) {
        const char* key = argv[2 * i + 1];
        const char* val = argv[2 * i + 2];

        if (std::strcmp(key, "S") == 0)     S = std::strtod(val, nullptr);
        if (std::strcmp(key, "K") == 0)     K = std::strtod(val, nullptr);
        if (std::strcmp(key, "T") == 0)     T = std::strtod(val, nullptr);
        if (std::strcmp(key, "sigma") == 0) sigma = std::strtod(val, nullptr);
        if (std::strcmp(key, "r") == 0)     r = std::strtod(val, nullptr);
        if (std::strcmp(key, "q") == 0)     q = std::strtod(val, nullptr);
        if (std::strcmp(key, "m") == 0)    fixed_point_iterations = std::atoi(val);
    }

    FastPut<RealType, CHEBYSHEV_NODES, QUAD_NODES_FP, QUAD_NODES_PRICE> 
        put(S, K, sigma, r, q, T, fixed_point_iterations);
    RealType ans = put.calc();
    return 0;
}

