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
#include <array>
#include <cmath>

// -----------------------------------------------------------------------------
// Tanh-Sinh Quadrature
//
// Generates quadrature nodes and weights for numerical integration using
// the tanh-sinh (double exponential) method. The nodes are in [-1, 1].
// -----------------------------------------------------------------------------

/// Computes a single tanh-sinh quadrature node and weight pair.
///
/// **Parameters:**
/// - `step_size`: Step size in the transformed domain
/// - `transformed_index`: Index in the transformed domain
///
/// **Returns:** A pair containing (node_location, weight)
template<typename Real>
inline std::pair<Real, Real> compute_tanh_sinh_node_and_weight(Real step_size, Real transformed_index) {
    const Real half_pi = Real(HALF_PI);
    Real z = half_pi * std::sinh(transformed_index);
    Real cosh_z = std::cosh(z);

    Real node_location = std::tanh(z);
    Real weight = half_pi * std::cosh(transformed_index) * step_size / (cosh_z * cosh_z);
    return {node_location, weight};
}

/// Container for quadrature nodes and weights.
///
/// Stores the node locations and corresponding weights for numerical integration.
template<size_t N, typename Real = RealType>
struct QuadratureRule {
    std::array<Real, N> nodes{};      // Quadrature node locations in [-1, 1]
    std::array<Real, N> weights{};    // Corresponding weights
};

/// Generates a tanh-sinh quadrature rule with N nodes.
///
/// **Parameters:**
/// - `truncation_parameter`: Controls the range of integration (typically 3.0-3.5)
///
/// **Returns:** Quadrature rule with nodes in [-1, 1] and corresponding weights
template<size_t N, typename Real = RealType>
inline QuadratureRule<N, Real> generate_tanh_sinh_quadrature(Real truncation_parameter) {
    QuadratureRule<N, Real> rule{};
    Real step_size = Real(2) * truncation_parameter / Real(N - 1);

    // Generate symmetric nodes (only need to compute half)
    for (size_t i = 0; i <= N / 2; ++i) {
        Real transformed_index = -truncation_parameter + Real(i) * step_size;
        auto [node, weight] = compute_tanh_sinh_node_and_weight<Real>(step_size, transformed_index);

        rule.nodes[i] = node;
        rule.nodes[N - 1 - i] = -node;  // Symmetric
        rule.weights[i] = weight;
        rule.weights[N - 1 - i] = weight;  // Symmetric
    }
    return rule;
}

// -----------------------------------------------------------------------------
// Time Transformation Precomputation
//
// For time-domain integrals, we need to transform quadrature nodes from [-1, 1]
// to [0, t]. The transformation is:
//   u = alpha[i] * t        (time point in [0, t])
//   dt = beta[i] * t        (remaining time = t - u)
// where alpha = 0.5 * (1 + node) maps [-1, 1] to [0, 1]
// -----------------------------------------------------------------------------

/// Precomputed time transformation coefficients for efficient integration.
///
/// Stores coefficients that transform quadrature nodes from [-1, 1] to [0, t] time points.
template<size_t N, typename Real>
struct TimeTransformationCoefficients {
    std::array<Real, N> alpha{};           // Maps node to [0, 1]: u = alpha[i] * t
    std::array<Real, N> beta{};            // Remaining time: dt = beta[i] * t = (1 - alpha[i]) * t
    std::array<Real, N> sqrt_alpha{};      // Precomputed sqrt(alpha[i])
    std::array<Real, N> sqrt_beta{};       // Precomputed sqrt(beta[i])
};

/// Precomputes time transformation coefficients from a quadrature rule.
///
/// This allows efficient transformation from [-1, 1] nodes to [0, t] time points.
///
/// **Parameters:**
/// - `rule`: The quadrature rule to transform
///
/// **Returns:** Time transformation coefficients with precomputed alpha, beta, and their square roots
template<size_t N, typename Real>
inline TimeTransformationCoefficients<N, Real> 
precompute_time_transformation(const QuadratureRule<N, Real>& rule) {
    TimeTransformationCoefficients<N, Real> coeffs{};
    for (size_t i = 0; i < N; ++i) {
        Real alpha = Real(0.5) * (Real(1) + rule.nodes[i]);
        coeffs.alpha[i] = alpha;
        coeffs.beta[i] = Real(1) - alpha;
        coeffs.sqrt_alpha[i] = std::sqrt(alpha);
        coeffs.sqrt_beta[i] = std::sqrt(coeffs.beta[i]);
    }
    return coeffs;
}

// -----------------------------------------------------------------------------
// Compile-time Generated Quadrature Rules
//
// These are generated at compile-time for the fixed-point iteration and
// price computation. The truncation parameters (3.3 and 3.5) are chosen
// to balance accuracy and numerical stability.
// -----------------------------------------------------------------------------

// Quadrature rule for fixed-point boundary iteration
const QuadratureRule<QUAD_NODES_FP, RealType> fixed_point_quadrature_rule =
    generate_tanh_sinh_quadrature<QUAD_NODES_FP, RealType>(RealType(3.3));

// Quadrature rule for price computation
const QuadratureRule<QUAD_NODES_PRICE, RealType> price_quadrature_rule =
    generate_tanh_sinh_quadrature<QUAD_NODES_PRICE, RealType>(RealType(3.5));

// Precomputed time transformation coefficients
const TimeTransformationCoefficients<QUAD_NODES_FP, RealType> fixed_point_time_coeffs =
    precompute_time_transformation(fixed_point_quadrature_rule);

const TimeTransformationCoefficients<QUAD_NODES_PRICE, RealType> price_time_coeffs =
    precompute_time_transformation(price_quadrature_rule);
