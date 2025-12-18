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
#include <cmath>

// -------------------------------
// math helpers
// -------------------------------
template <typename T>
static inline int sgn(T val) {
    return (T(0) < val) - (val < T(0));
}

template <typename Real>
static inline Real normal_pdf(Real x) {
    return std::exp(Real(-0.5) * x * x) * Real(sqrt2piinv);
}

template <typename Real>
static inline Real normal_cdf(Real x) {
    return Real(0.5) + Real(0.5) * std::erf(x * Real(sqrt2inv));
}

