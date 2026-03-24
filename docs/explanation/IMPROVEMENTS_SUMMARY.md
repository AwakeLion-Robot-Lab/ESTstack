# ESKFOM Improvements Summary

## Overview

This document summarizes the improvements made to the ESKFOM (Error-State Kalman Filter on Manifolds) implementation based on the reference implementation from the manif library's SE2 localization example.

**Reference:** https://github.com/artivis/manif/blob/devel/examples/se2_localization.cpp

## Key Improvements

### 1. Configurable Covariance Update Method (eskfom.hpp:242-260)

Added support for two covariance update strategies:

#### a) Joseph Form (Default - More Stable)
```cpp
P = (I - KH) * P * (I - KH)^T + K * Hw * R * Hw^T * K^T
```
- Guarantees positive semi-definite covariance
- Better numerical stability
- Prevents filter divergence
- Recommended for production systems

#### b) Standard Form (manif-compatible - Simpler)
```cpp
S = H * P * H^T + Hw * R * Hw^T
P = P - K * S * K^T
```
- Simpler and more efficient
- Matches canonical textbook formulation
- Used in manif se2_localization.cpp
- Sufficient for well-conditioned systems

**Configuration:**
```cpp
eststack::solution::ESKFOM<StateT>::Options options;
options.use_joseph_form = false;  // Use standard form
```

### 2. Optional Reset Jacobian (eskfom.hpp:263-281)

The reset Jacobian correction can now be disabled:

**With Reset (Default):**
```cpp
G = I ± 0.5 * ad_dx
P = G * P * G^T
```

**Without Reset (manif-compatible):**
- No correction applied
- Acceptable for small error states
- Used in manif se2_localization.cpp

**Configuration:**
```cpp
options.apply_reset_jacobian = false;  // Disable for simple cases
```

### 3. Fixed dt Parameter Handling (eskfom.hpp:116-119, 187-197)

**Previous bug:** Model methods were called without the required `dt` parameter

**Fixed:**
- `predictImpl` now includes `dt` parameter (default: 1.0)
- `updateImpl` now includes `dt` parameter (default: 1.0)
- All model method calls properly pass `dt`

**Impact:** This fixes compatibility with models that depend on time step (like CT model)

### 4. Innovation Tracking (eskfom.hpp:88-103, 199-215)

Added methods to extract filter statistics for debugging/analysis:

```cpp
// Store innovation and innovation covariance during update
last_innovation_ = y;               // y = z - h(x)
last_innovation_cov_ = S;           // S = H * P * H^T + R

// Retrieve for analysis
auto innovation = eskf.getInnovation<MeasurementT>();
auto S = eskf.getInnovationCovariance<2>();
```

**Use cases:**
- Consistency checking (Mahalanobis distance)
- Outlier detection
- Filter health monitoring
- Debugging estimation issues

This pattern follows the manif example where E and Z are computed separately for analysis.

### 5. Optimized Kalman Gain Computation (eskfom.hpp:217-229)

**Improvement:** Compute innovation covariance S once and reuse:

```cpp
// Compute S = H * P * H^T + Hw * R * Hw^T
const auto S = (H * this->P_ * H.transpose() + Hw * R * Hw.transpose()).eval();
last_innovation_cov_ = S;  // Store for analysis

// Use pre-computed S in gain calculation
const auto K = this->P_ * H.transpose() * S.inverse();
```

**Benefits:**
- Eliminates redundant computation
- More efficient for standard form update
- Cleaner code structure

## New Files Created

### 1. SE2 Motion Model
**File:** `include/eststack/model/motion/se2.hpp`

Simple SE2 motion model:
```cpp
X(k+1) = X(k) * exp(u) = X(k) ⊞ u
```

Uses manif's `rplus()` operation with automatic Jacobian computation.

### 2. SE2 Landmark Measurement Model
**File:** `include/eststack/model/measurement/landmark_se2.hpp`

Landmark measurement model:
```cpp
y = X^{-1} * b
```

**Key features:**
- Chain rule Jacobian: H = J_h_xi * J_xi_x
- Uses manif's `inverse()` and `act()` with automatic Jacobians
- Follows the exact pattern from manif se2_localization.cpp

### 3. SE2 Localization Example
**File:** `examples/se2_localization_eskfom.cpp`

Complete working example demonstrating:
- ESKFOM with SE2 state
- Multiple landmark measurements
- Comparison: filtered vs unfiltered vs ground truth
- Configuration of Joseph form and reset Jacobian

### 4. Comprehensive Documentation
**Files:**
- `docs/explanation/eskfom.md` - Updated with configuration guide
- `docs/explanation/eskfom_improvements.md` - Detailed improvement explanation
- `examples/README.md` - Example usage guide

## Comparison: Before vs After

| Feature | Before | After |
|---------|--------|-------|
| Covariance update | Joseph form only | Joseph or standard (configurable) |
| Reset Jacobian | Always applied | Optional (configurable) |
| dt parameter | Missing (bug) | Properly passed with defaults |
| Innovation tracking | Not available | Stored and retrievable |
| Innovation cov. | Computed twice | Computed once, reused |
| SE2 support | No example | Complete example |
| Configuration | Hardcoded | Runtime configurable |

## Usage Examples

### Example 1: manif-compatible mode (simple, fast)
```cpp
// Configure to match manif se2_localization.cpp
eststack::solution::ESKFOM<manif::SE2d>::Options options;
options.use_joseph_form = false;
options.apply_reset_jacobian = false;

eststack::solution::ESKFOM<manif::SE2d> filter(x0, P0, options);
```

### Example 2: Production mode (stable, robust)
```cpp
// Use defaults for maximum stability
eststack::solution::ESKFOM<manif::SE2d> filter(x0, P0);  // Both options true by default
```

### Example 3: Innovation monitoring
```cpp
filter.update(meas_model, z, R);

// Check innovation for outlier detection
auto innovation = filter.getInnovation<Eigen::Vector2d>();
auto S = filter.getInnovationCovariance<2>();
double mahalanobis = innovation.transpose() * S.inverse() * innovation;

if (mahalanobis > threshold) {
    // Reject measurement as outlier
}
```

## Technical Details

### Chain Rule Jacobian Pattern (from manif)

The measurement model demonstrates proper chain rule application:

```cpp
// h(x, b) = x^{-1} * b
// Chain rule: dh/dx = (dh/d(x^{-1})) * (d(x^{-1})/dx)

manif::SE2d::Jacobian J_xi_x;  // J_xi_x = d(x^{-1})/dx
manif::SE2d x_inv = x.inverse(J_xi_x);

Eigen::Matrix<double, 2, 3> J_h_xi;  // J_h_xi = dh/d(x^{-1})
x_inv.act(landmark_position_, J_h_xi);

H = J_h_xi * J_xi_x;  // Chain rule
```

This pattern is directly from manif se2_localization.cpp lines 151-153.

### Innovation Covariance Computation

The manif example computes innovation covariance in two steps (lines 149-153):
```cpp
E = H * P * H.transpose();  // Projected state covariance
Z = E + R;                   // Innovation covariance
```

Our implementation combines this as:
```cpp
S = H * P * H^T + Hw * R * Hw^T  // Full innovation covariance
```

Both are equivalent when Hw = I (additive noise).

## Build and Test

### Building the example:
```bash
xmake config --mode=debug
xmake build ESTstack-example-se2_localization_eskfom
```

### Running:
```bash
xmake run ESTstack-example-se2_localization_eskfom
```

**Expected output:** 10 time steps showing simulated, estimated, and unfiltered trajectories. The ESKFOM estimate should track the ground truth better than the unfiltered estimate.

## References

1. **manif library SE2 localization:** https://github.com/artivis/manif/blob/devel/examples/se2_localization.cpp
2. **Solà et al. (2018):** "A micro Lie theory for state estimation in robotics" - arXiv:1812.01537
3. **Barfoot (2024):** "State Estimation for Robotics" (2nd ed.), Chapter 2.2.15 (Joseph form)
4. **Deray & Solà (2020):** "manif: A small header-only library for Lie theory" - JOSS

## Summary

These improvements make ESKFOM:
1. **More flexible** - configurable to match different reference implementations
2. **More correct** - fixes dt parameter bug
3. **More efficient** - eliminates redundant computation
4. **More observable** - provides innovation statistics
5. **Better documented** - clear examples and configuration guide

The implementation now properly balances:
- **Numerical stability** (Joseph form, reset Jacobian) for production
- **Simplicity** (standard form, no reset) for prototyping/learning
- **Performance** (SMW lemma, optimized computation) for efficiency
