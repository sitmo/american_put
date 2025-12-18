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

// Type definition
typedef double RealType;

// Mathematical constants
#define PI         3.141592653589793238462643383279502884197169L
#define sqrt2inv   0.707106781186547524400844362104849039284836L
#define sqrt2piinv 0.398942280401432677939946059934381868475859L

// Compile-time configurable node counts
#ifndef QUAD_NODES_FP
#define QUAD_NODES_FP 9
#endif
#ifndef FIXED_POINT_ITER
#define FIXED_POINT_ITER 2
#endif
#ifndef CHEBYSHEV_NODES
#define CHEBYSHEV_NODES 5
#endif
#ifndef QUAD_NODES_PRICE
#define QUAD_NODES_PRICE 21
#endif

static_assert(CHEBYSHEV_NODES >= 2, "CHEBYSHEV_NODES must be >= 2 (Clenshaw uses bv[2]).");

