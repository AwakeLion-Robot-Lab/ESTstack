# ESTstack Examples

This directory contains example applications demonstrating the use of ESTstack algorithms.

## SE2 Localization Example

**File:** `se2_localization_eskfom.cpp`

**Description:** 2D robot localization using ESKFOM with fixed landmarks. This example is based on the reference implementation from the manif library ([se2_localization.cpp](https://github.com/artivis/manif/blob/devel/examples/se2_localization.cpp)).

**What it demonstrates:**
- ESKFOM prediction step with SE2 motion model
- ESKFOM update step with landmark measurements
- Comparison between filtered, unfiltered, and ground truth trajectories
- Configuration options for covariance update (Joseph form vs standard form)
- Configuration options for reset Jacobian application

**Models used:**
- **Motion model:** SE2 motion with velocity twist control (see `include/eststack/model/motion/se2.hpp`)
- **Measurement model:** Landmark position in robot frame (see `include/eststack/model/measurement/landmark_se2.hpp`)

**Building:**
```bash
xmake config --mode=debug
xmake build ESTstack-example-se2_localization_eskfom
```

**Running:**
```bash
xmake run ESTstack-example-se2_localization_eskfom
```

## Key Differences from manif Example

The ESTstack ESKFOM implementation provides additional features compared to the basic manif example:

1. **Configurable covariance update**: Choose between Joseph form (more stable) or standard form (simpler)
2. **Optional reset Jacobian**: Can enable/disable the reset Jacobian correction
3. **Perturbation convention**: Support for both global and local perturbation
4. **SMW lemma optimization**: Automatic use of matrix inversion lemma for high-dimensional measurements
5. **Generic state types**: Works with any Lie group state (SE2, SE3, Bundle, etc.)

## References

- Solà, J., et al. "A micro Lie theory for state estimation in robotics." arXiv:1812.01537 (2018)
- Deray, J., & Solà, J. "manif: A small header-only library for Lie theory." JOSS (2020)
- Barfoot, T. D. "State Estimation for Robotics" (2024)
