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

// -------------------------------
// Black-Scholes building blocks
// -------------------------------
template <typename Real>
struct Put {
    Real S, K, sigma, r, q, t;
};

template <typename Real>
struct Put_1 {
    Real e_rt, e_qt, sqrt_t, sigma_sqrt_t;

    Put_1(Put<Real>& p) { operator()(p); }

    void operator()(Put<Real>& p) {
        sqrt_t = std::sqrt(p.t);
        sigma_sqrt_t = p.sigma * sqrt_t;
        e_rt = std::exp(-p.r * p.t);
        e_qt = std::exp(-p.q * p.t);
    }
};

template <typename Real>
struct Put_1_1 {
    Real d1, d2, ncdf_d1, ncdf_d2;

    Put_1_1() {}
    Put_1_1(Put<Real>& p, Put_1<Real>& p1) { operator()(p, p1); }

    void operator()(Put<Real>& p, Put_1<Real>& p1) {
        d1 = (std::log(p.S / p.K) + (p.r - p.q + Real(0.5) * p.sigma * p.sigma) * p.t) / p1.sigma_sqrt_t;
        d2 = d1 - p1.sigma_sqrt_t;
        ncdf_d1 = normal_cdf(-d1);
        ncdf_d2 = normal_cdf(-d2);
    }
};

template <typename Real>
static inline void bs_european_put_price_theta(Put<Real>& p, Put_1<Real>& p1, Put_1_1<Real>& p2,
                                              Real& price, Real& theta) {
    Real term1 = p2.ncdf_d2 * p.K * p1.e_rt;
    Real term2 = p2.ncdf_d1 * p.S * p1.e_qt;
    price = term1 - term2;
    theta = p.r * term1 - p.q * term2
          - p.sigma * p.S / (Real(2) * p1.sqrt_t) * p1.e_qt * normal_pdf(p2.d1);
}

template <typename Real>
static inline Real bs_european_put_price(Put<Real>& p, Put_1<Real>& p1, Put_1_1<Real>& p2) {
    return p2.ncdf_d2 * p.K * p1.e_rt - p2.ncdf_d1 * p.S * p1.e_qt;
}

