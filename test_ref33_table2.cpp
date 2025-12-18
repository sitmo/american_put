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
// Replicates Table 2 from reference [33]: Prices of American Put Options
// Compares FastPut method against true prices

#include "fast_put.hpp"
#include <iostream>
#include <iomanip>
#include <chrono>
#include <cmath>
#include <vector>
#include <iomanip>
#include <sstream>

// Helper function to compute American put price
template <typename Real, size_t NumChebyshevNodes, size_t NumQuadNodesFP, size_t NumQuadNodesPrice>
Real compute_price(Real S, Real K, Real sigma, Real r, Real q, Real T, int fixed_point_iterations) {
    FastPut<Real, NumChebyshevNodes, NumQuadNodesFP, NumQuadNodesPrice> 
        put(S, K, sigma, r, q, T, fixed_point_iterations);
    
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
    return veur + vee;
}

// Compute true price using high-precision FastPut (equivalent to BIN-BS 15000)
template <typename Real>
Real compute_true_price(Real S, Real K, Real sigma, Real r, Real q, Real T) {
    // Use high-precision parameters: (l,m,n) = (65,8,32), p = 101
    return compute_price<Real, 32, 65, 101>(S, K, sigma, r, q, T, 8);
}

int main() {
    // Model parameters from table title: K = 100, δ = q = 0.04, r = 0.04
    const RealType K = 100.0;
    const RealType r = 0.04;
    const RealType q = 0.04;
    
    // FastPut parameters - using (l,m,n) = (15,3,7), p = 41 (from Table 2 paper)
    const int l = 15;  // QUAD_NODES_FP
    const int m = 3;   // fixed_point_iterations
    const int n = 7;   // CHEBYSHEV_NODES
    const int p = 41;  // QUAD_NODES_PRICE
    
    // Table structure: (T, sigma, Spot, True Price)
    struct TableRow {
        RealType T;
        RealType sigma;
        RealType S;
        RealType true_price;
    };
    
    // True prices from reference [33] Table 2
    // Note: These are the "True" values from the table (computed with BIN-BS 15000)
    const TableRow rows[] = {
        // T = 0.5, σ = 0.2
        {0.5, 0.2, 80.0, 20.14372},
        {0.5, 0.2, 100.0, 5.54634},
        {0.5, 0.2, 120.0, 0.70724},
        
        // T = 0.5, σ = 0.5
        {0.5, 0.5, 80.0, 24.67736},
        {0.5, 0.5, 100.0, 13.80581},
        {0.5, 0.5, 120.0, 7.28738},
        
        // T = 3, σ = 0.2
        {3.0, 0.2, 80.0, 23.22834},
        {3.0, 0.2, 100.0, 12.60521},
        {3.0, 0.2, 120.0, 6.482425},
        
        // T = 3, σ = 0.5
        {3.0, 0.5, 80.0, 37.97483},
        {3.0, 0.5, 100.0, 30.74247},
        {3.0, 0.5, 120.0, 25.21333}
    };
    
    // Print table header
    std::cout << "Table 2 – Prices of American Put Options (K = 100, δ = 0.04, r = 0.04)\n";
    std::cout << "FastPut method with (l, m, n) = (" << l << ", " << m << ", " << n << "), p = " << p << "\n\n";
    
    std::cout << std::left << std::setw(6) << "T"
              << std::setw(8) << "σ"
              << std::setw(8) << "Spot"
              << std::setw(12) << "True"
              << std::setw(15) << "FastPut"
              << std::setw(15) << "Error"
              << "\n";
    std::cout << std::string(64, '-') << "\n";
    
    std::vector<RealType> errors;
    std::vector<RealType> relative_errors;
    
    // Process each row
    for (const auto& row : rows) {
        RealType computed_price = compute_price<RealType, n, l, p>(
            row.S, K, row.sigma, r, q, row.T, m);
        
        RealType absolute_error = std::abs(computed_price - row.true_price);
        RealType relative_error = absolute_error / row.true_price;
        
        errors.push_back(absolute_error);
        relative_errors.push_back(relative_error);
        
        // Format output
        std::cout << std::fixed << std::setprecision(2) << std::setw(6) << row.T
                  << std::setw(8) << row.sigma
                  << std::setprecision(0) << std::setw(8) << row.S
                  << std::setprecision(5) << std::setw(12) << row.true_price
                  << std::setw(15) << computed_price
                  << std::scientific << std::setprecision(2) << std::setw(15) << absolute_error
                  << "\n";
    }
    
    // Compute RMSE and RRMSE
    RealType rmse = 0.0;
    RealType rrmse = 0.0;
    for (size_t i = 0; i < errors.size(); ++i) {
        rmse += errors[i] * errors[i];
        rrmse += relative_errors[i] * relative_errors[i];
    }
    rmse = std::sqrt(rmse / errors.size());
    rrmse = std::sqrt(rrmse / relative_errors.size());
    
    std::cout << std::string(64, '-') << "\n";
    std::cout << std::left << std::setw(22) << "RMSE:"
              << std::scientific << std::setprecision(3) << rmse << "\n";
    std::cout << std::setw(22) << "RRMSE:"
              << std::scientific << std::setprecision(3) << rrmse << "\n";
    
    return 0;
}

