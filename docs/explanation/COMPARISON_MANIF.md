# Detailed Comparison: ESKFOM vs manif SE2 Localization

This document provides a line-by-line comparison between the ESTstack ESKFOM implementation and the manif library's se2_localization.cpp reference.

## Algorithm Flow Comparison

### Prediction Step

| Step | manif se2_localization.cpp | ESTstack ESKFOM |
|------|---------------------------|-----------------|
| 1. Compute Jacobians | `X.plus(u, J_x, J_u)` | `model.computeStateJacobian(x, u, dt)` |
| | | `model.computeNoiseJacobian(x, u, dt)` |
| 2. Propagate state | `X = X.plus(u, J_x, J_u)` | `x = model.compute(x, u, dt)` |
| 3. Perturbation conv. | Implicit (local) | Explicit (global/local configurable) |
| 4. Propagate cov. | `P = J_x * P * J_x^T + J_u * U * J_u^T` | `P = Fx * P * Fx^T + Fw * Q * Fw^T` |

**Key difference:** ESKFOM separates Jacobian computation from state propagation, allowing generic models.

### Update Step

| Step | manif se2_localization.cpp | ESTstack ESKFOM |
|------|---------------------------|-----------------|
| 1. Expectation | `e = X.inverse(J_xi_x).act(b, J_e_xi)` | `e = model.compute(x, dt)` |
| 2. Jacobian | `H = J_e_xi * J_xi_x` (chain rule) | `H = model.computeMeasJacobian(x, dt)` |
| 3. Innovation cov. | `E = H * P * H^T; Z = E + R` | `S = H * P * H^T + Hw * R * Hw^T` |
| 4. Innovation | `z = y - e` | `y = z - e` |
| 5. Kalman gain | `K = P * H^T * Z.inverse()` | `K = P * H^T * S.inverse()` or SMW |
| 6. Error state | `dx = K * z` | `dx = K * y` |
| 7. State update | `X = X + dx` (rplus) | `x = x.rplus(dx)` or `x.lplus(dx)` |
| 8. Cov. update | `P = P - K * Z * K^T` | **Configurable:** |
| | | Joseph: `P = (I-KH)*P*(I-KH)^T + K*R*K^T` |
| | | Standard: `P = P - K * S * K^T` |
| 9. Reset | **Not applied** | **Optional:** `P = G * P * G^T` |

## Code Comparison

### manif: Prediction (lines 95-98)
```cpp
X = X.plus(u_est, J_x, J_u);
P = J_x * P * J_x.transpose() + J_u * U * J_u.transpose();
```

### ESKFOM: Prediction (eskfom.hpp:116-141)
```cpp
const auto Fx_raw = model.computeStateJacobian(this->x_, u, dt);
const auto Fw_raw = model.computeNoiseJacobian(this->x_, u, dt);
const State priori_x = model.compute(this->x_, u, dt);

// Perturbation convention conversion
const auto Fx = (P == Global) ? priori_x.adj() * Fx_raw * x.inv().adj() : Fx_raw;
const auto Fw = (P == Global) ? priori_x.adj() * Fw_raw : Fw_raw;

this->x_ = priori_x;
this->P_ = Fx * this->P_ * Fx.transpose() + Fw * Q * Fw.transpose();
```

**Difference:** ESKFOM supports generic models and both perturbation conventions.

---

### manif: Update (lines 105-122)
```cpp
e = X.inverse(J_xi_x).act(b, J_e_xi);
H = J_e_xi * J_xi_x;
E = H * P * H.transpose();

z = y - e;
Z = E + R;

K = P * H.transpose() * Z.inverse();

dx = K * z;
X = X + dx;
P = P - K * Z * K.transpose();
```

### ESKFOM: Update (eskfom.hpp:187-260)
```cpp
const auto H_raw = model.computeMeasJacobian(this->x_, dt);
const auto Hw = model.computeNoiseJacobian(this->x_, dt);

// Perturbation convention conversion
const auto H = (P == Global) ? H_raw * x.inverse().adj() : H_raw;

const auto y = (z - model.compute(this->x_, dt)).eval();
const auto S = (H * P * H^T + Hw * R * Hw^T).eval();

// SMW lemma for high-dimensional measurements
const auto K = (MeasDim > 6*StateDim)
    ? (P^{-1} + H^T * R^{-1} * H)^{-1} * H^T * R^{-1}
    : P * H^T * S^{-1};

const auto dx = K * y;
x = (P == Global) ? x.lplus(dx) : x.rplus(dx);

// Configurable covariance update
if (use_joseph_form)
    P = (I - K*H) * P * (I - K*H)^T + K * R * K^T;
else
    P = P - K * S * K^T;

// Optional reset Jacobian
if (apply_reset_jacobian)
    P = G * P * G^T;  // G = I ± 0.5 * ad_dx
```

**Differences:**
1. Generic models instead of hardcoded operations
2. Perturbation convention support
3. SMW lemma optimization
4. Configurable covariance update
5. Optional reset Jacobian

## Performance Characteristics

| Aspect | manif (lines of code) | ESKFOM (lines of code) |
|--------|----------------------|------------------------|
| Prediction | ~3 lines | ~30 lines (+ genericity) |
| Update | ~12 lines | ~70 lines (+ options) |
| Type safety | Runtime | Compile-time (concepts) |
| Flexibility | Fixed SE2 | Any Lie group |
| Optimization | None | SMW lemma |
| Stability | Standard | Configurable (Joseph) |

## Recommendations

**Use manif-compatible mode when:**
- Prototyping or learning
- Simple state spaces (SE2, SO3)
- Short-duration filtering
- Reference checking against manif

**Use production mode when:**
- Long-running filters
- Complex state spaces (SE_2(3), high-dim Bundles)
- Ill-conditioned systems
- Stability is critical

## Configuration Guide

```cpp
// manif-compatible (matches se2_localization.cpp exactly)
ESKFOM<SE2d>::Options manif_mode;
manif_mode.use_joseph_form = false;
manif_mode.apply_reset_jacobian = false;

// Production mode (maximum stability)
ESKFOM<SE2d>::Options prod_mode;  // defaults are fine
prod_mode.use_joseph_form = true;
prod_mode.apply_reset_jacobian = true;

// Balanced mode (stable update, skip reset for speed)
ESKFOM<SE2d>::Options balanced_mode;
balanced_mode.use_joseph_form = true;
balanced_mode.apply_reset_jacobian = false;
```

## Conclusion

The improved ESKFOM implementation:
1. ✅ Maintains compatibility with manif reference (via configuration)
2. ✅ Adds production-grade stability features (Joseph form)
3. ✅ Fixes critical bugs (dt parameter passing)
4. ✅ Provides observability (innovation tracking)
5. ✅ Optimizes performance (SMW lemma, efficient computation)
6. ✅ Supports generic state spaces (template-based)
7. ✅ Includes complete example and documentation

All improvements are based on established literature (Solà, Barfoot) and reference implementations (manif).
