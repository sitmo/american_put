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
#include "quadrature.hpp"
#include "math_utils.hpp"
#include "black_scholes.hpp"
#include "boundary_init.hpp"
#include <array>
#include <cmath>
#include <limits>
#include <iostream>
#include <iomanip>
#include <chrono>
#include <cstdint>

// -------------------------------
// Main pricer (Chebyshev/boundary arrays are now compile-time std::array)
// -------------------------------
template <typename Real>
class AloPut {
public:
    static constexpr size_t NC = size_t(CHEBYSHEV_NODES) + 1;

    Real S, K, sigma, r, q, T;

    std::array<Real, NC> zv{}; // Chebyshev nodes [-1,1]
    std::array<Real, NC> xv{}; // rescaled [0,sqrt(T)]
    std::array<Real, NC> tv{}; // nodes [0,T]
    std::array<Real, NC> Bv{}; // boundary values
    std::array<Real, NC> Hv{}; // H(sqrt(t))
    std::array<Real, NC> av{}; // Chebyshev coeffs
    std::array<Real, NC> bv{}; // Clenshaw temp

    // cached constants
    Real sqrtT, inv_sqrtT;
    Real logK, log_max, log_min;
    Real sig2;

    AloPut(Real S, Real K, Real sigma, Real r, Real q, Real T)
        : S(S), K(K), sigma(sigma), r(r), q(q), T(T)
    {
        sqrtT     = std::sqrt(T);
        inv_sqrtT = Real(1) / sqrtT;

        logK    = std::log(K);
        log_max = std::log(std::numeric_limits<Real>::max());
        log_min = std::log(std::numeric_limits<Real>::min());

        sig2 = sigma * sigma;

        initChebyshev();
    }

    void initChebyshev() {
        const Real x_factor = Real(0.5) * sqrtT;
        const Real pi_over_n = Real(PI) / Real(CHEBYSHEV_NODES);

        for (size_t i = 0; i < NC; ++i) {
            zv[i] = std::cos(pi_over_n * Real(i));
            xv[i] = x_factor * (Real(1) + zv[i]);
            tv[i] = xv[i] * xv[i];
        }
        tv[0] = T;
        tv[NC - 1] = Real(0);
    }

    inline void d_plus_minus_sqrt(Real t, Real sqrt_t, Real z, Real& d_plus, Real& d_minus) {
        const Real eps = Real(1e-12);
        if (t <= eps) {
            Real lz = std::log(z);
            if (std::fabs(lz) < Real(1e-30)) {
                d_plus = Real(0);
                d_minus = Real(0);
                return;
            }
            Real big = Real(1e30);
            d_plus  = big * lz;
            d_minus = big * lz;
            return;
        }

        Real inv = Real(1) / (sigma * sqrt_t);
        Real left  = std::log(z) + (r - q) * t;
        Real right = Real(0.5) * sig2 * t;
        d_plus  = inv * (left + right);
        d_minus = inv * (left - right);
    }

    void initalGuessB() {
        const Real tol = Real(1e-4);
        const uint32_t max_iter = 5;

        for (size_t i = 0; i < NC; ++i) {
            Put<Real> p = {S, K, sigma, r, q, tv[i]};
            Bv[i] = gbm_american_put_oeb_qd_plus<Real>(p, tol, max_iter);
        }
        if (q == Real(0)) Bv[NC - 1] = K;
        else Bv[NC - 1] = K * std::min(Real(1), r / q);
    }

    inline void initH() {
        for (size_t i = 0; i < NC; ++i) {
            Real g = std::log(Bv[i]) - logK;
            Hv[i] = g * std::fabs(g);
        }
    }

    inline void initChebyshevInterpolation() {
        constexpr uint32_t n = uint32_t(CHEBYSHEV_NODES);
        const Real c = Real(2) / Real(n);
        const uint32_t two_n = 2 * n;

        for (uint32_t k = 0; k <= n; ++k) {
            av[k] = Real(0.5) * Hv[0];
            const Real sgnk = (k & 1) ? Real(-1) : Real(1);
            av[k] += Real(0.5) * Hv[n] * sgnk;

            for (uint32_t i = 1; i < n; ++i) {
                uint32_t ki = (k * i) % two_n;
                if (ki > n) ki = two_n - ki;
                av[k] += Hv[i] * zv[ki];
            }
            av[k] *= c;
        }
    }

    inline Real q_c(Real z) {
        constexpr uint32_t n = uint32_t(CHEBYSHEV_NODES);
        bv[n] = Real(0.5) * av[n];
        bv[n - 1] = av[n - 1] + Real(2) * z * bv[n];
        for (int32_t k = int32_t(n) - 2; k >= 0; --k) {
            bv[size_t(k)] = av[size_t(k)] + Real(2) * z * bv[size_t(k) + 1] - bv[size_t(k) + 2];
        }
        return Real(0.5) * (bv[0] - bv[2]);
    }

    inline Real Bt_from_sqrt(Real sqrt_ti) {
        Real z = Real(2) * sqrt_ti * inv_sqrtT - Real(1);
        Real H_x = q_c(z);

        Real expo = std::sqrt(std::fabs(H_x)) * Real(sgn(H_x));
        if (logK + expo > log_max) expo = log_max - logK;
        if (logK + expo < log_min) expo = log_min - logK;

        return K * std::exp(expo);
    }

    Real FixedPoint_B_StepOrdinary(Real t) {
        Real d_plus, d_minus;

        const Real sqrt_t = std::sqrt(t);
        const Real B1 = Bt_from_sqrt(sqrt_t);

        d_plus_minus_sqrt(t, sqrt_t, B1 / K, d_plus, d_minus);
        Real N_tB = normal_cdf(d_minus);
        Real D_tB = normal_cdf(d_plus);

        Real Nint = Real(0);
        Real Dint = Real(0);

        for (size_t i = 0; i < QUAD_NODES_FP; ++i) {
            const Real u = fixed_point_time_coeffs.alpha[i] * t;
            const Real sqrt_u = fixed_point_time_coeffs.sqrt_alpha[i] * sqrt_t;

            const Real Bu = Bt_from_sqrt(sqrt_u);

            Real dt = fixed_point_time_coeffs.beta[i] * t;
            if (dt < Real(1e-18)) dt = Real(1e-18);
            const Real sqrt_dt = std::sqrt(dt);

            d_plus_minus_sqrt(dt, sqrt_dt, B1 / Bu, d_plus, d_minus);

            Nint += fixed_point_quadrature_rule.weights[i] * std::exp(r * u) * normal_cdf(d_minus);
            Dint += fixed_point_quadrature_rule.weights[i] * std::exp(q * u) * normal_cdf(d_plus);
        }

        Nint *= t * Real(0.5);
        Dint *= t * Real(0.5);

        N_tB += r * Nint;
        D_tB += q * Dint;

        const Real K_star = K * std::exp(-(r - q) * t);
        return K_star * N_tB / D_tB;
    }

    Real Vee(Real t) {
        Real d_plus, d_minus;

        const Real sqrt_t = std::sqrt(t);
        Real ans = Real(0);

        for (size_t i = 0; i < QUAD_NODES_PRICE; ++i) {
            const Real u = price_time_coeffs.alpha[i] * t;
            const Real sqrt_u = price_time_coeffs.sqrt_alpha[i] * sqrt_t;

            const Real Bu = Bt_from_sqrt(sqrt_u);

            Real dt = price_time_coeffs.beta[i] * t;
            if (dt < Real(1e-18)) dt = Real(1e-18);
            const Real sqrt_dt = std::sqrt(dt);

            d_plus_minus_sqrt(dt, sqrt_dt, S / Bu, d_plus, d_minus);

            const Real y =
                r * K * std::exp(-r * (t - u)) * normal_cdf(-d_minus)
              - q * S * std::exp(-q * (t - u)) * normal_cdf(-d_plus);

            ans += price_quadrature_rule.weights[i] * y;
        }

        return ans * t * Real(0.5);
    }

    Real Main(int m) {
        auto t1 = std::chrono::high_resolution_clock::now();

        initalGuessB();
        initH();
        initChebyshevInterpolation();

        Put<Real> put = {S, K, sigma, r, q, T};
        Put_1<Real> p1(put);
        Put_1_1<Real> p2(put, p1);
        Real veur = bs_european_put_price(put, p1, p2);

        for (int ww = 0; ww < m; ++ww) {
            for (size_t i = 0; i + 1 < NC; ++i) {
                Bv[i] = FixedPoint_B_StepOrdinary(tv[i]);
            }
            initH();
            initChebyshevInterpolation();
        }

        Real vee = Vee(T);

        auto t2 = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double, std::milli> fp_ms = t2 - t1;

        std::cout.precision(2);
        std::cout << std::fixed;

        std::cout << "{\n"
            << "\t\"l_FixedPointQuadNodes\": " << QUAD_NODES_FP << ",\n"
            << "\t\"m_FixedPointIter\": " << m << ",\n"
            << "\t\"n_ChebyshevNodes\": " << CHEBYSHEV_NODES << ",\n"
            << "\t\"p_PriceQuadNodes\": " << QUAD_NODES_PRICE << ",\n"
            << "\t\"S\": " << S << ",\n"
            << "\t\"K\": " << K << ",\n"
            << "\t\"T\": " << T << ",\n";

        std::cout.precision(4);
        std::cout << std::fixed;
        std::cout
            << "\t\"sigma\": " << sigma << ",\n"
            << "\t\"r\": " << r << ",\n"
            << "\t\"q\": " << q << ",\n"
            << "\t\"time_ms\": " << fp_ms.count() << ",\n";

        std::cout.precision(12);
        std::cout << std::scientific;
        std::cout
            << "\t\"european_put\": " << veur << ",\n"
            << "\t\"american_put\": " << (veur + vee) << ",\n"
            << "\t\"ee_premium\": " << vee << "\n"
            << "\t\"err\": " << vee  - 0.106952702747 << "\n"
            << "}\n";

        return (veur + vee);
    }
};

