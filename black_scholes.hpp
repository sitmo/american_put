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
#include "math_utils.hpp"
#include <cmath>

// -----------------------------------------------------------------------------
// Black-Scholes European Put Option Pricing
// -----------------------------------------------------------------------------

/// Computes the price of a European put option using the Black-Scholes formula
/// @param S Current stock price
/// @param K Strike price
/// @param sigma Volatility
/// @param r Risk-free rate
/// @param q Dividend yield
/// @param t Time to expiration
/// @return European put option price
template <typename Real>
inline Real compute_european_put_price(Real S, Real K, Real sigma, Real r, Real q, Real t) {
    if (t <= Real(0)) {
        return std::max(Real(0), K - S);
    }

    const Real sqrt_t = std::sqrt(t);
    const Real sigma_sqrt_t = sigma * sqrt_t;
    const Real d1 = (std::log(S / K) + (r - q + Real(0.5) * sigma * sigma) * t) / sigma_sqrt_t;
    const Real d2 = d1 - sigma_sqrt_t;
    
    const Real N_d1 = normal_cdf(-d1);
    const Real N_d2 = normal_cdf(-d2);
    const Real e_rt = std::exp(-r * t);
    const Real e_qt = std::exp(-q * t);
    
    return N_d2 * K * e_rt - N_d1 * S * e_qt;
}

/// Computes the price and theta (time derivative) of a European put option
/// @param S Current stock price
/// @param K Strike price
/// @param sigma Volatility
/// @param r Risk-free rate
/// @param q Dividend yield
/// @param t Time to expiration
/// @param[out] price European put option price
/// @param[out] theta Time derivative of the option price (dP/dt)
template <typename Real>
inline void compute_european_put_price_and_theta(Real S, Real K, Real sigma, Real r, Real q, Real t,
                                                 Real& price, Real& theta) {
    if (t <= Real(0)) {
        price = std::max(Real(0), K - S);
        theta = Real(0);
        return;
    }

    const Real sqrt_t = std::sqrt(t);
    const Real sigma_sqrt_t = sigma * sqrt_t;
    const Real d1 = (std::log(S / K) + (r - q + Real(0.5) * sigma * sigma) * t) / sigma_sqrt_t;
    const Real d2 = d1 - sigma_sqrt_t;
    
    const Real N_d1 = normal_cdf(-d1);
    const Real N_d2 = normal_cdf(-d2);
    const Real e_rt = std::exp(-r * t);
    const Real e_qt = std::exp(-q * t);
    
    const Real term1 = N_d2 * K * e_rt;
    const Real term2 = N_d1 * S * e_qt;
    
    price = term1 - term2;
    theta = r * term1 - q * term2 - sigma * S / (Real(2) * sqrt_t) * e_qt * normal_pdf(d1);
}
