# ESKFOM Improvements Based on manif Reference

This document describes improvements made to the ESKFOM implementation based on the reference implementation in the manif library ([se2_localization.cpp](https://github.com/artivis/manif/blob/devel/examples/se2_localization.cpp)).

## Summary of Improvements

### 1. Configurable Covariance Update Method

**Location:** `include/eststack/solution/eskfom.hpp:205-229`

The ESKFOM implementation now supports two covariance update methods:

#### Joseph Form (Default, more stable)
```cpp
P = (I - KH) * P * (I - KH)^T + K * R * K^T
```

**Advantages:**
- Guarantees positive semi-definite covariance
- Better numerical stability
- Prevents filter divergence in ill-conditioned systems

#### Standard Form (as in manif example)
```cpp
S = H * P * H^T + R
P = P - K * S * K^T
```

**Advantages:**
- Simpler and more efficient
- Matches canonical Kalman filter equations
- Sufficient for well-conditioned systems

**Configuration:**
```cpp
eststack::solution::ESKFOM<StateT>::Options options;
options.use_joseph_form = false;  // Use standard form like manif
```

### 2. Optional Reset Jacobian

**Location:** `include/eststack/solution/eskfom.hpp:226-249`

The reset Jacobian correction can now be disabled for simple cases:

**With Reset Jacobian (Default):**
```cpp
G = I ± 0.5 * ad_dx
P = G * P * G^T
```

**Without Reset Jacobian (as in manif example):**
```cpp
// No correction applied
```

The manif se2_localization.cpp example doesn't apply the reset Jacobian correction. This is acceptable for:
- Simple state spaces (SE2, SO3)
- Small error state increments
- Short-duration filtering

**Configuration:**
```cpp
eststack::solution::ESKFOM<StateT>::Options options;
options.apply_reset_jacobian = false;  // Disable reset Jacobian
```

### 3. SE2 Motion Model

**New file:** `include/eststack/model/motion/se2.hpp`

A clean SE2 motion model implementation:
- State: SE(2) pose (x, y, theta)
- Control: Twist in se(2) tangent space
- Motion equation: X(k+1) = X(k) ⊞ u = X(k) * exp(u)

This model uses manif's built-in Jacobian computation through the `rplus()` operation.

### 4. SE2 Landmark Measurement Model

**New file:** `include/eststack/model/measurement/landmark_se2.hpp`

Implementation of landmark measurements in robot frame:
- Measurement equation: y = X^{-1} * b
- Uses chain rule for Jacobian: H = J_h_xi * J_xi_x
- Leverages manif's `inverse()` and `act()` operations with automatic Jacobian computation

**Key insight from manif example:** The chain rule pattern:
```cpp
manif::SE2d::Jacobian J_xi_x;
manif::SE2d x_inv = x.inverse(J_xi_x);

Eigen::Matrix<double, 2, 3> J_h_xi;
x_inv.act(landmark_position_, J_h_xi);

H = J_h_xi * J_xi_x;  // Chain rule
```

### 5. Complete SE2 Localization Example

**New file:** `examples/se2_localization_eskfom.cpp`

A complete working example that:
- Simulates a robot moving in 2D with control noise
- Simulates landmark measurements with measurement noise
- Runs ESKFOM to estimate the robot pose
- Compares filtered vs unfiltered vs ground truth trajectories
- Demonstrates the configuration options

This example directly corresponds to the manif se2_localization.cpp reference.

## Comparison: ESKFOM vs manif Example

| Feature | manif se2_localization.cpp | ESTstack ESKFOM |
|---------|---------------------------|-----------------|
| State representation | SE(2) | Generic (any Lie group) |
| Covariance update | Standard form | Both Joseph & standard |
| Reset Jacobian | Not applied | Optional (configurable) |
| Perturbation | Local (implicit) | Both global & local |
| SMW lemma | Not used | Automatic for high-dim measurements |
| Type system | manif types | Generic concepts + templates |

## Implementation Notes

### Why Joseph Form is Better

From Barfoot (2024), Chapter 2.2.15:

The Joseph form guarantees that the updated covariance remains positive semi-definite even in the presence of:
- Numerical round-off errors
- Ill-conditioned state covariance matrices
- Model linearization errors

The standard form `P = P - K*S*K^T` can produce slightly negative eigenvalues due to numerical errors, which can cause filter instability over time.

### When to Disable Reset Jacobian

The reset Jacobian accounts for second-order effects in the error state covariance redistribution. It can be safely disabled when:

1. Error state increments are small (||dx|| << 1)
2. State space is simple (SE2, SO3 without large perturbations)
3. Short-duration filtering (a few updates)

For production systems with long-running filters or large state perturbations, keeping the reset Jacobian enabled is recommended.

## Testing

Build and run the example:
```bash
xmake config --mode=debug
xmake build ESTstack-example-se2_localization_eskfom
xmake run ESTstack-example-se2_localization_eskfom
```

Expected output: The filtered estimate should track the simulated trajectory better than the unfiltered estimate.

## References

1. Solà, J., Deray, J., & Atchuthan, D. (2018). "A micro Lie theory for state estimation in robotics." arXiv:1812.01537
2. Barfoot, T. D. (2024). "State Estimation for Robotics" (2nd ed.)
3. Deray, J., & Solà, J. (2020). "manif: A small header-only library for Lie theory." JOSS
4. manif library SE2 example: https://github.com/artivis/manif/blob/devel/examples/se2_localization.cpp
