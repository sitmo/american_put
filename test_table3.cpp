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
// Replicates Table 3 from the paper: Estimated 3-year put option price values
// We replicate the FP-B column for spot prices S = 80, 100, 120

#include "fast_put.hpp"
#include <iostream>
#include <iomanip>
#include <chrono>
#include <cmath>
#include <vector>

// Helper function to compute American put price without JSON output
template <typename Real, size_t NumChebyshevNodes, size_t NumQuadNodesFP, size_t NumQuadNodesPrice>
std::pair<Real, double> compute_price_and_time(
    Real S, Real K, Real sigma, Real r, Real q, Real T, int fixed_point_iterations) {
    
    FastPut<Real, NumChebyshevNodes, NumQuadNodesFP, NumQuadNodesPrice> 
        put(S, K, sigma, r, q, T, fixed_point_iterations);
    
    auto t1 = std::chrono::high_resolution_clock::now();
    
    put.initalGuessB();
    put.initH();
    put.initChebyshevInterpolation();
    
    Real veur = compute_european_put_price(S, K, sigma, r, q, T);
    
    for (int ww = 0; ww < fixed_point_iterations; ++ww) {
        for (size_t i = 0; i + 1 < put.NC; ++i) {
            put.Bv[i] = put.FixedPoint_B_StepOrdinary(put.tv[i]);
        }
        put.initH();
        put.initChebyshevInterpolation();
    }
    
    Real vee = put.Vee(T);
    Real american_price = veur + vee;
    
    auto t2 = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = t2 - t1;
    
    return {american_price, elapsed.count()};
}

int main() {
    // Model parameters from Table 3 caption
    const RealType K = 100.0;
    const RealType T = 3.0;
    const RealType sigma = 0.2;
    const RealType r = 0.04;
    const RealType q = 0.04;
    
    // FP-B parameters: (l, m, n) = (10, 3, 7), p = 25
    const int l = 10;  // QUAD_NODES_FP
    const int m = 3;   // fixed_point_iterations
    const int n = 7;   // CHEBYSHEV_NODES
    const int p = 25;  // QUAD_NODES_PRICE
    
    // True prices (from FP-A with (l,m,n) = (131,16,64))
    struct SpotData {
        RealType S;
        RealType true_price;
    };
    
    const SpotData spots[] = {
        {80.0, 23.22834},
        {100.0, 12.60521},
        {120.0, 6.482425}
    };
    
    // Print table header
    std::cout << "Table 3: Estimated 3-year put option price values for K = 100 and S as listed in the table.\n";
    std::cout << "Model settings: r = q = 4%, σ = 0.2\n";
    std::cout << "FP-B parameters: (l, m, n) = (" << l << ", " << m << ", " << n << "), p = " << p << "\n";
    std::cout << "True prices computed with FP-A (l,m,n) = (131,16,64)\n\n";
    
    // Print header rows
    std::cout << std::left << std::setw(10) << "Spot S"
              << std::setw(12) << "True Price"
              << std::right << std::setw(12) << "FP-B"
              << "\n";
    std::cout << std::left << std::setw(10) << ""
              << std::setw(12) << ""
              << std::right << std::setw(12) << ("(" + std::to_string(l) + "," + std::to_string(m) + "," + std::to_string(n) + ")")
              << "\n";
    std::cout << std::string(34, '-') << "\n";
    
    // Compute and print each row
    for (const auto& spot : spots) {
        auto [price, cpu_time] = compute_price_and_time<RealType, n, l, p>(
            spot.S, K, sigma, r, q, T, m);
        
        RealType absolute_error = std::abs(price - spot.true_price);
        
        // Format output: price on first line, error on second line (in italics style)
        std::cout << std::left << std::fixed << std::setprecision(0) << std::setw(10) << spot.S
                  << std::setprecision(5) << std::setw(12) << spot.true_price
                  << std::right << std::setw(12) << price
                  << "\n";
        std::cout << std::left << std::setw(10) << ""
                  << std::setw(12) << ""
                  << std::right << std::scientific << std::setprecision(2) << std::setw(12) << absolute_error
                  << "\n";
    }
    
    return 0;
}

