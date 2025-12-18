# High-Performance American Put Option Pricing

A C++ implementation of the fast American put option pricing method described in:

> Leif Andersen, Mark Lake, and Dimitri Offengenden,  
> "High-performance American option pricing",  
> Journal of Computational Finance (2015).  
> https://ideas.repec.org/a/rsk/journ0/2464632.html

## Overview

This repository implements a high-performance numerical method for pricing American put options using:
- **Chebyshev interpolation** for the early-exercise boundary representation
- **Tanh-sinh quadrature** for numerical integration
- **Fixed-point iteration** to converge the boundary

The method achieves high accuracy (relative errors typically < 1e-6) with very fast computation times (microseconds to milliseconds per option).

## Performance

### Accuracy
- **Table 2 replication**: American premium convergence study shows relative errors decreasing from ~2.6e-02 to ~1.3e-07 as parameters increase
- **Table 3 replication**: Comparison against reference prices shows RMSE of ~1.1e-04 and RRMSE of ~7.3e-06 across various market conditions

### Speed
- Typical computation time: **10-400 microseconds** per option (depending on parameter settings)
- Example: With (l,m,n) = (15,3,7), p = 41, pricing takes ~60-100 microseconds

### Example Results

**Table 2** (1-year American premium, K = S = 100, r = q = 5%, σ = 0.25):
- With (l,m,n) = (15,3,7), p = 41: Premium = 0.10695120, Relative Error = 1.40e-05, Time = 60μs
- With (l,m,n) = (35,8,16), p = 81: Premium = 0.10695272, Relative Error = 1.85e-07, Time = 346μs

**Table 3** (3-year put options, various spot prices):
- RMSE: 1.095e-04
- RRMSE: 7.261e-06
- All test cases show errors < 2e-04

## Code Structure

The codebase is organized into clean, modular headers:

- **`constants.hpp`**: Type definitions, mathematical constants, and compile-time configuration
- **`math_utils.hpp`**: Mathematical utility functions (sign, normal PDF/CDF)
- **`quadrature.hpp`**: Tanh-sinh quadrature generation and time transformation
- **`black_scholes.hpp`**: Black-Scholes European put option pricing
- **`boundary_init.hpp`**: QD+ method for initial early-exercise boundary estimation
- **`fast_put.hpp`**: Main `FastPut` class implementing the American put pricer

### Template Parameters

The `FastPut` class is templated for flexibility:

```cpp
template <typename Real, size_t NumChebyshevNodes, size_t NumQuadNodesFP, size_t NumQuadNodesPrice>
class FastPut {
    // ...
    int fixed_point_iterations;  // Runtime parameter
};
```

**Parameters:**
- `Real`: Floating point type (e.g., `double`, `float`)
- `NumChebyshevNodes`: Number of Chebyshev nodes (n in paper)
- `NumQuadNodesFP`: Number of quadrature nodes for fixed-point iteration (l in paper)
- `NumQuadNodesPrice`: Number of quadrature nodes for price computation (p in paper)
- `fixed_point_iterations`: Number of fixed-point iterations (m in paper) - **runtime parameter**

## Building

```bash
make
```

This builds all executables:
- `test_fast_put`: Basic test program
- `test_table2`: Replicates Table 2 from the paper
- `test_table3`: Replicates Table 3 from the paper (FP-B column)
- `test_ref33_table2`: Comparison against reference [33] Table 2
- `test_ref33_table3`: Comparison against reference [33] Table 3
- `reference_tree_put`: Binomial tree reference implementation

## Usage

### Basic Usage

```cpp
#include "fast_put.hpp"

// Model parameters
double S = 100.0;      // Spot price
double K = 100.0;      // Strike price
double T = 1.0;        // Time to expiration
double sigma = 0.25;   // Volatility
double r = 0.05;       // Risk-free rate
double q = 0.05;       // Dividend yield
int m = 3;             // Fixed-point iterations

// Create pricer with template parameters: (n, l, p) = (7, 15, 41)
FastPut<double, 7, 15, 41> put(S, K, sigma, r, q, T, m);

// Compute American put price
double price = put.calc();
```

### Command-Line Interface

```bash
# Basic test with default parameters
./test_fast_put

# With custom parameters
./test_fast_put S 110 K 100 T 1.0 sigma 0.3 r 0.05 q 0.05 m 4

# Run table replications
./test_table2
./test_table3
./test_ref33_table2
./test_ref33_table3
```

## Parameter Mapping

The paper uses notation `(l, m, n)` and `p`:

| Paper Notation | Code Parameter | Description |
|----------------|----------------|-------------|
| `n` | `NumChebyshevNodes` | Chebyshev nodes for boundary representation |
| `l` | `NumQuadNodesFP` | Quadrature nodes for fixed-point iteration |
| `m` | `fixed_point_iterations` | Number of fixed-point iterations (runtime) |
| `p` | `NumQuadNodesPrice` | Quadrature nodes for price computation |

**Example**: Paper configuration `(l,m,n) = (15,3,7), p = 41` corresponds to:
```cpp
FastPut<double, 7, 15, 41> put(S, K, sigma, r, q, T, 3);
```

## Test Programs

- **`test_fast_put.cpp`**: Basic functionality test with configurable parameters
- **`test_table2.cpp`**: Replicates Table 2 showing convergence of American premium
- **`test_table3.cpp`**: Replicates Table 3 FP-B column for 3-year options
- **`test_ref33_table2.cpp`**: Comparison against reference [33] Table 2 (12 test cases)
- **`test_ref33_table3.cpp`**: Comparison against reference [33] Table 3 (12 test cases)
- **`reference_tree_put.cpp`**: Binomial tree implementation for comparison

## Dependencies

- C++17 compiler (tested with GCC and Clang)
- Standard library only (no external dependencies)

## License

MIT License - see LICENSE file for details.

## References

1. Andersen, L., Lake, M., & Offengenden, D. (2015). High-performance American option pricing. *Journal of Computational Finance*.
2. Reference [33] from the paper (for Table 2 and Table 3 comparisons)

