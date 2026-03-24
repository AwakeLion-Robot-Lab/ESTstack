# ESKFOM Explanation

## Reset Jacobian Derivation

we need some tricks about lie bracket $[\cdot,\cdot]$ to derive how small adjoint related to wedge mapping, here provides derivation in $\boldsymbol{SO}(3)$:

![reset jacobian derivation](./assets/reset_jacobian_derivation.jpg)

**NOTE:** if manifold is $\boldsymbol{SE}(3)$ or $\boldsymbol{SE}_2(3)$, the skew-symmetric matrix would not equal to small adjoint matrix since $(\cdot)^{\wedge}: \mathbb{R}^{n} \to \mathbb{R}^{m\times m}$, but $ \boldsymbol{ad}_{(\cdot)}: \mathbb{R}^n \to \mathbb{R}^{n\times n}$, here we define $m$ is dimension of lie group (`manif::Dim`) and $n$ is dimension of lie algebra (`manif::DoF`).

## Covariance Update Methods

The ESKFOM implementation supports two covariance update methods:

### 1. Joseph Form (Default)
**Formula:** $P = (I - KH) P (I - KH)^T + K R K^T$

**Advantages:**
- Guarantees positive semi-definite covariance matrix
- More numerically stable for ill-conditioned systems
- Recommended for production systems

**Use case:** Complex state spaces, long-running filters, or when numerical stability is critical

### 2. Standard Form (as in manif se2_localization.cpp)
**Formula:** $P = P - K S K^T$, where $S = H P H^T + R$

**Advantages:**
- Simpler and more computationally efficient
- Matches canonical Kalman filter formulation
- Used in many reference implementations (e.g., manif examples)

**Use case:** Simple state spaces, short-duration filters, or when performance is critical

### Configuration

```cpp
eststack::solution::ESKFOM<StateT>::Options options;
options.use_joseph_form = true;      // Use Joseph form (default)
options.apply_reset_jacobian = true; // Apply reset Jacobian (default)

// Or use standard form as in manif examples:
options.use_joseph_form = false;
options.apply_reset_jacobian = false;

eststack::solution::ESKFOM<StateT> filter(x0, P0, options);
```

## Reset Jacobian Application

The reset Jacobian correction $G$ accounts for the redistribution of error state covariance after state injection:

**Global perturbation:** $G = I + \frac{1}{2} \boldsymbol{ad}_{\delta x}$

**Local perturbation:** $G = I - \frac{1}{2} \boldsymbol{ad}_{\delta x}$

This correction is theoretically important for maintaining consistency, but can be disabled for simple cases where the error state increments are small (as done in the manif se2_localization.cpp example).

## References

- See `manif` library SE2 localization example: [se2_localization.cpp](https://github.com/artivis/manif/blob/devel/examples/se2_localization.cpp)
- Solà et al., "A micro Lie theory for state estimation in robotics" (2018)
- Barfoot, "State Estimation for Robotics" (2024), Chapter 2.2.15 for Joseph form
