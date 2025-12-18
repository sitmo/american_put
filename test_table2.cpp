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
// Replicates Table 2 from the paper: Estimated 1-year American premium for K = S = 100

#include "fast_put.hpp"
#include <iostream>
#include <iomanip>
#include <chrono>
#include <cmath>

// Helper function to compute American premium without JSON output
template <typename Real, size_t NumChebyshevNodes, size_t NumQuadNodesFP, size_t NumQuadNodesPrice>
std::pair<Real, double> compute_premium_and_time(
    Real S, Real K, Real sigma, Real r, Real q, Real T, int fixed_point_iterations) {
    
    FastPut<Real, NumChebyshevNodes, NumQuadNodesFP, NumQuadNodesPrice> 
        put(S, K, sigma, r, q, T, fixed_point_iterations);
    
    auto t1 = std::chrono::high_resolution_clock::now();
    
    put.initalGuessB();
    put.initH();
    put.initChebyshevInterpolation();
    
    // European price not needed for premium calculation, but computation path is same
    (void)compute_european_put_price(S, K, sigma, r, q, T);
    
    for (int ww = 0; ww < fixed_point_iterations; ++ww) {
        for (size_t i = 0; i + 1 < put.NC; ++i) {
            put.Bv[i] = put.FixedPoint_B_StepOrdinary(put.tv[i]);
        }
        put.initH();
        put.initChebyshevInterpolation();
    }
    
    Real vee = put.Vee(T);
    
    auto t2 = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = t2 - t1;
    
    return {vee, elapsed.count()};
}

int main() {
    // Model parameters from Table 2 caption
    const RealType S = 100.0;
    const RealType K = 100.0;
    const RealType T = 1.0;
    const RealType sigma = 0.25;
    const RealType r = 0.05;
    const RealType q = 0.05;
    
    // Reference premium computed with (l,m,n) = (201, 16, 64) and p = 201
    const RealType reference_premium = 0.106952702747;
    
    // Table 2 data: (l, m, n) and p values
    struct TableRow {
        int l, m, n, p;
    };
    
    const TableRow rows[] = {
        {5, 1, 4, 15},
        {7, 2, 5, 20},
        {11, 2, 5, 31},
        {15, 2, 6, 41},
        {15, 3, 7, 41},
        {25, 4, 9, 51},
        {25, 5, 12, 61},
        {25, 6, 15, 61},
        {35, 8, 16, 81},
        {51, 8, 24, 101},
        {65, 8, 32, 101}
    };
    
    // Print table header
    std::cout << "Table 2: Estimated 1-year American premium for K = S = 100\n";
    std::cout << "Model: r = q = 5%, σ = 0.25\n";
    std::cout << "Reference premium: " << std::fixed << std::setprecision(12) << reference_premium 
              << " (computed with (l,m,n) = (201,16,64), p = 201)\n\n";
    
    std::cout << std::left << std::setw(12) << "(l,m,n)"
              << std::setw(6) << "p"
              << std::right << std::setw(20) << "American Premium"
              << std::setw(18) << "Relative Error"
              << std::setw(15) << "CPU Seconds"
              << "\n";
    std::cout << std::string(71, '-') << "\n";
    
    // Compute and print each row
    // Note: Template parameters must be compile-time constants, so we handle each case
    auto compute_row = [&](int l, int m, int n, int p) {
        RealType premium;
        double cpu_time;
        
        // Dispatch to correct template instantiation based on parameters
        if (n == 4 && l == 5 && p == 15) {
            std::tie(premium, cpu_time) = compute_premium_and_time<RealType, 4, 5, 15>(S, K, sigma, r, q, T, m);
        } else if (n == 5 && l == 7 && p == 20) {
            std::tie(premium, cpu_time) = compute_premium_and_time<RealType, 5, 7, 20>(S, K, sigma, r, q, T, m);
        } else if (n == 5 && l == 11 && p == 31) {
            std::tie(premium, cpu_time) = compute_premium_and_time<RealType, 5, 11, 31>(S, K, sigma, r, q, T, m);
        } else if (n == 6 && l == 15 && p == 41) {
            std::tie(premium, cpu_time) = compute_premium_and_time<RealType, 6, 15, 41>(S, K, sigma, r, q, T, m);
        } else if (n == 7 && l == 15 && p == 41) {
            std::tie(premium, cpu_time) = compute_premium_and_time<RealType, 7, 15, 41>(S, K, sigma, r, q, T, m);
        } else if (n == 9 && l == 25 && p == 51) {
            std::tie(premium, cpu_time) = compute_premium_and_time<RealType, 9, 25, 51>(S, K, sigma, r, q, T, m);
        } else if (n == 12 && l == 25 && p == 61) {
            std::tie(premium, cpu_time) = compute_premium_and_time<RealType, 12, 25, 61>(S, K, sigma, r, q, T, m);
        } else if (n == 15 && l == 25 && p == 61) {
            std::tie(premium, cpu_time) = compute_premium_and_time<RealType, 15, 25, 61>(S, K, sigma, r, q, T, m);
        } else if (n == 16 && l == 35 && p == 81) {
            std::tie(premium, cpu_time) = compute_premium_and_time<RealType, 16, 35, 81>(S, K, sigma, r, q, T, m);
        } else if (n == 24 && l == 51 && p == 101) {
            std::tie(premium, cpu_time) = compute_premium_and_time<RealType, 24, 51, 101>(S, K, sigma, r, q, T, m);
        } else if (n == 32 && l == 65 && p == 101) {
            std::tie(premium, cpu_time) = compute_premium_and_time<RealType, 32, 65, 101>(S, K, sigma, r, q, T, m);
        } else {
            std::cerr << "Error: Unsupported parameter combination (" << l << "," << m << "," << n << "), p=" << p << "\n";
            return;
        }
        
        RealType relative_error = std::abs((premium - reference_premium) / reference_premium);
        
        // Format output with proper alignment
        std::cout << std::left << std::setw(12) << ("(" + std::to_string(l) + "," + std::to_string(m) + "," + std::to_string(n) + ")")
                  << std::setw(6) << p
                  << std::right << std::fixed << std::setprecision(12) << std::setw(20) << premium
                  << std::scientific << std::setprecision(2) << std::setw(18) << relative_error
                  << std::setw(15) << cpu_time
                  << "\n";
    };
    
    for (const auto& row : rows) {
        compute_row(row.l, row.m, row.n, row.p);
    }
    
    return 0;
}

