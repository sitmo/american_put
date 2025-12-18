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
#include <cmath>
#include <cstdint>

// -------------------------------
// QD+ boundary init
// -------------------------------
template<typename Real>
Real gbm_american_put_oeb_qd_plus(Put<Real>& p_, Real tol = 1E-6, unsigned max_steps = 10) {

    Put<Real> p(p_);
    Put_1<Real> p1(p);
    Put_1_1<Real> p2;

    if (p.t == Real(0)) {
        if (p.q == Real(0)) return p.K;
        return p.K * std::min(Real(1), p.r / p.q);
    }

    Real h = Real(1) - p1.e_rt;
    Real s2 = p.sigma * p.sigma;
    Real s2h = s2 * h;
    Real w_1 = Real(2) * (p.r - p.q) / s2 - Real(1);
    Real s_ = std::sqrt(w_1 * w_1 + Real(8) * p.r / s2h);
    Real lambda  = Real(0.5) * (-w_1 - s_);
    Real lambda_ = Real(2) * p.r / (s2h * h * s_);

    Real c0a = Real(2) * lambda + w_1;
    Real c0b = -(Real(1) - h) * Real(2) * p.r / s2 / c0a;
    Real c0c = Real(1) / h + lambda_ / c0a;

    Real put_price, put_theta;
    Real yu, ym, yd;
    Real B = p.K;

    uint32_t step = 0;
    while (step < max_steps) {
        ++step;

        Real dB = B * Real(1e-6);

        // MID
        p.S = B;
        p2(p, p1);
        bs_european_put_price_theta(p, p1, p2, put_price, put_theta);
        ym =
            -p1.e_qt * p2.ncdf_d1
            + c0b * (
                -put_theta / (p.r * p1.e_rt)
                + c0c * (p.K - p.S - put_price)
            ) / p.S
            + lambda * (p.K - p.S - put_price) / p.S
            + Real(1);

        if (std::fabs(ym) < tol) break;

        // UP
        p.S = B + dB;
        p2(p, p1);
        bs_european_put_price_theta(p, p1, p2, put_price, put_theta);
        yu =
            -p1.e_qt * p2.ncdf_d1
            + c0b * (
                -put_theta / (p.r * p1.e_rt)
                + c0c * (p.K - p.S - put_price)
            ) / p.S
            + lambda * (p.K - p.S - put_price) / p.S
            + Real(1);

        // DOWN
        p.S = B - dB;
        p2(p, p1);
        bs_european_put_price_theta(p, p1, p2, put_price, put_theta);
        yd =
            -p1.e_qt * p2.ncdf_d1
            + c0b * (
                -put_theta / (p.r * p1.e_rt)
                + c0c * (p.K - p.S - put_price)
            ) / p.S
            + lambda * (p.K - p.S - put_price) / p.S
            + Real(1);

        Real dy  = Real(0.5) * (yu - yd) / dB;
        Real ddy = ((yu - ym) - (ym - yd)) / (dB * dB);

        Real denom = (Real(2) * dy * dy - ym * ddy);
        if (denom == Real(0)) break;

        Real Bnew = B - (Real(2) * ym * dy) / denom;
        if (Bnew <= Real(0)) Bnew = p.K * Real(0.5);

        B = Bnew;
    }

    return B;
}

