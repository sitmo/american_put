// MIT License
//
// Copyright (c) Thijs van den Berg
//
// This code implements the method described in:
// Leif Andersen, Mark Lake, and Dimitri Offengenden,
// "High-performance American option pricing",
// Journal of Computational Finance (2015).
// https://ideas.repec.org/a/rsk/journ0/2464632.html

#pragma once

#include "constants.hpp"
#include "black_scholes.hpp"
#include "math_utils.hpp"
#include <cmath>
#include <cstdint>

// -----------------------------------------------------------------------------
// QD+ Boundary Initialization
//
// Computes an initial estimate of the early exercise boundary for an
// American put option using the QD+ (quadratic dividend) method.
// -----------------------------------------------------------------------------

/// Computes initial early exercise boundary estimate using QD+ method.
///
/// **Parameters:**
/// - `K`: Strike price
/// - `sigma`: Volatility
/// - `r`: Risk-free rate
/// - `q`: Dividend yield
/// - `t`: Time to expiration
/// - `tol`: Convergence tolerance (default: 1E-6)
/// - `max_steps`: Maximum number of iterations (default: 10)
///
/// **Returns:** Estimated early exercise boundary
template<typename Real>
Real compute_initial_early_exercise_boundary(Real K, Real sigma, Real r, Real q, Real t,
                                            Real tol = 1E-6, unsigned max_steps = 10) {
    if (t == Real(0)) {
        if (q == Real(0)) return K;
        return K * std::min(Real(1), r / q);
    }

    // Precompute common values
    const Real e_rt = std::exp(-r * t);
    const Real h = Real(1) - e_rt;
    const Real s2 = sigma * sigma;
    const Real s2h = s2 * h;
    const Real w_1 = Real(2) * (r - q) / s2 - Real(1);
    const Real s_ = std::sqrt(w_1 * w_1 + Real(8) * r / s2h);
    const Real lambda = Real(0.5) * (-w_1 - s_);
    const Real lambda_ = Real(2) * r / (s2h * h * s_);

    const Real c0a = Real(2) * lambda + w_1;
    const Real c0b = -(Real(1) - h) * Real(2) * r / s2 / c0a;
    const Real c0c = Real(1) / h + lambda_ / c0a;

    Real put_price, put_theta;
    Real yu, ym, yd;
    Real B = K;

    uint32_t step = 0;
    while (step < max_steps) {
        ++step;

        const Real dB = B * Real(1e-6);
        const Real sqrt_t = std::sqrt(t);
        const Real sigma_sqrt_t = sigma * sqrt_t;

        // MID
        Real S_mid = B;
        compute_european_put_price_and_theta(S_mid, K, sigma, r, q, t, put_price, put_theta);
        
        const Real d1_mid = (std::log(S_mid / K) + (r - q + Real(0.5) * sigma * sigma) * t) / sigma_sqrt_t;
        const Real e_qt = std::exp(-q * t);
        const Real N_d1_mid = normal_cdf(-d1_mid);
        
        ym = -e_qt * N_d1_mid
            + c0b * (-put_theta / (r * e_rt) + c0c * (K - S_mid - put_price)) / S_mid
            + lambda * (K - S_mid - put_price) / S_mid
            + Real(1);

        if (std::fabs(ym) < tol) break;

        // UP
        Real S_up = B + dB;
        compute_european_put_price_and_theta(S_up, K, sigma, r, q, t, put_price, put_theta);
        
        const Real d1_up = (std::log(S_up / K) + (r - q + Real(0.5) * sigma * sigma) * t) / sigma_sqrt_t;
        const Real N_d1_up = normal_cdf(-d1_up);
        
        yu = -e_qt * N_d1_up
            + c0b * (-put_theta / (r * e_rt) + c0c * (K - S_up - put_price)) / S_up
            + lambda * (K - S_up - put_price) / S_up
            + Real(1);

        // DOWN
        Real S_down = B - dB;
        compute_european_put_price_and_theta(S_down, K, sigma, r, q, t, put_price, put_theta);
        
        const Real d1_down = (std::log(S_down / K) + (r - q + Real(0.5) * sigma * sigma) * t) / sigma_sqrt_t;
        const Real N_d1_down = normal_cdf(-d1_down);
        
        yd = -e_qt * N_d1_down
            + c0b * (-put_theta / (r * e_rt) + c0c * (K - S_down - put_price)) / S_down
            + lambda * (K - S_down - put_price) / S_down
            + Real(1);

        const Real dy = Real(0.5) * (yu - yd) / dB;
        const Real ddy = ((yu - ym) - (ym - yd)) / (dB * dB);

        const Real denom = (Real(2) * dy * dy - ym * ddy);
        if (denom == Real(0)) break;

        Real Bnew = B - (Real(2) * ym * dy) / denom;
        if (Bnew <= Real(0)) Bnew = K * Real(0.5);

        B = Bnew;
    }

    return B;
}
