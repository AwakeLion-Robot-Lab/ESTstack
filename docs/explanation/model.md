# Model Explanation

**All models here are based on tangent space, since ONLY can do linear operation on tangent space but not on manifold directly.**

## Auto-Computable

We use `plus(t, J_mout_m, J_mout_t)` in [manif](https://github.com/artivis/manif/blob/devel/include/manif/impl/lie_group_base.h#L184) to auto compute state jacobian $\mathbf{F}_x$ and noise jacobian $\mathbf{F}_i$, but this can only work in normal lie group, in fact, extended lie groups like $\boldsymbol{SE}_2(3)$ and $\textbf{Bundle}$ do not meet this API.

That's why I define `AutoComputable` concept to set constraint for you guys about `plus()`. 

### Transition Model

A general state transition on manifold can be written as:

$$\mathbf{x}^- = f(\mathbf{x}, \mathbf{u}, \mathbf{i}) = \mathbf{x} \oplus \boldsymbol{\tau}(\mathbf{x}, \mathbf{u}, \mathbf{i}) = \mathbf{x} \circ \textbf{exp}(\mathbf{\tau})$$

where $\boldsymbol{\tau}\in \mathbf{T}_{I}\mathcal{M}$ is a tangent vector means **small increment on manifold**, and $\delta \mathbf{x}\in \mathbf{T}_{I}\mathcal{M}$ is also a tangent but with different meaning which is **error-state(right-perturbation) denoted by $\mathbf{x}_{\text{true}} = \mathbf{x} \oplus \delta \mathbf{x}$**.

Calling `x.plus(tau, J_x, J_u)` gives two DoF $\times$ DoF Jacobians:

- $J_\mathbf{x} = \frac{\partial(\mathbf{x} \oplus \boldsymbol{\tau})}{\partial \mathbf{x}}$
- $J_\mathbf{u} = \frac{\partial(\mathbf{x} \oplus \boldsymbol{\tau})}{\partial \boldsymbol{\tau}}$

The ESKF Jacobians are obtained via **chain rule**:

$$\mathbf{F}_x = \frac{\partial f}{\partial \delta \mathbf{x}} = J_\mathbf{x} + J_\mathbf{u} \cdot \frac{\partial \boldsymbol{\tau}}{\partial \delta\mathbf{x}}, \qquad \mathbf{F}_i =\frac{\partial f}{\partial \mathbf{i}} = J_\mathbf{u} \cdot \frac{\partial \boldsymbol{\tau}}{\partial \mathbf{i}}$$

where:

| Symbol                                                         | Size                    | Meaning                                        |
| -------------------------------------------------------------- | ----------------------- | ---------------------------------------------- |
| $J_\mathbf{x}, J_\mathbf{u}$                                   | DoF $\times$ DoF        | from `plus()`                                  |
| $\frac{\partial \boldsymbol{\tau}}{\partial \delta\mathbf{x}}$ | DoF $\times$ DoF        | how $\boldsymbol{\tau}$ depends on error state |
| $\frac{\partial \boldsymbol{\tau}}{\partial \mathbf{i}}$       | DoF $\times$ dim(noise) | noise-to-tangent mapping                       |

> [!NOTE]
> In manif, Dim = dim(Manifold) and DoF = dim(Tangent) = **dim(Error-State)**.

#### Auto-Computable Case

When **both** of the following hold:

1. `ControlInput = State::Tangent` (control lives in full tangent space)
2. $\boldsymbol{\tau} = \mathbf{u}$ (tangent vector is the control, independent of state)

the chain rule simplifies to:

$$\frac{\partial \boldsymbol{\tau}}{\partial \delta\mathbf{x}} = \mathbf{0} \implies \mathbf{F}_x = J_\mathbf{x}, \qquad \frac{\partial \boldsymbol{\tau}}{\partial \mathbf{i}} = \mathbf{I} \implies \mathbf{F}_i = J_\mathbf{u}$$

So `x.plus(u, Fx, Fw)` directly gives both Jacobians.

##### Example

$\boldsymbol{SE}(3)$ with velocity control $\mathbf{u} = [\mathbf{v}, \boldsymbol{\omega}] \in \mathbb{R}^6$:

$$
\boldsymbol{\tau} = \mathbf{u} \cdot dt, \quad \frac{\partial\boldsymbol{\tau}}{\partial\delta\mathbf{x}} = \mathbf{0}, \quad \frac{\partial\boldsymbol{\tau}}{\partial\mathbf{i}} = \mathbf{I}, \mathbf{F}_x = J_\mathbf{x}, \mathbf{F}_u = J_\mathbf{u}
$$

```cpp
// Fx = J_x, Fi = J_u directly from manif
priori_x = x.plus(u * dt, Fx, Fi);
```

#### Non-Auto-Computable Case

When `dim(ControlInput) != dim(Error-State)` or $\boldsymbol{\tau}$ depends on state, **you must hand-derive** $\frac{\partial \boldsymbol{\tau}}{\partial \delta\mathbf{x}}$ and $\frac{\partial \boldsymbol{\tau}}{\partial \mathbf{i}}$, then compose with $J_\mathbf{u}$ from `plus()`. Or simply hand-write $\mathbf{F}_x$ and $\mathbf{F}_i$ directly (note to discretize).

##### Example

$\boldsymbol{SE}_2(3)$ with IMU input $\mathbf{u} = [\boldsymbol{\alpha}, \boldsymbol{\omega}] \in \mathbb{R}^6$:

$$
\boldsymbol{\tau} = \begin{bmatrix} \frac{1}{2}\mathbf{a}dt^2 + \mathbf{R}^\top\mathbf{v}dt \\ \boldsymbol{\omega}dt \\ \mathbf{a}dt \end{bmatrix}, \quad \mathbf{a} = \boldsymbol{\alpha} - \mathbf{R}^\top\mathbf{g}
$$

$$
\frac{\partial\boldsymbol{\tau}}{\partial\delta\mathbf{x}} = \begin{bmatrix} \mathbf{0} & (\frac{1}{2}\mathbf{a}dt^2 + \mathbf{R}^\top\mathbf{v}dt)^{\wedge} & \mathbf{R}^\top dt \\ \mathbf{0} & \mathbf{0} & \mathbf{0} \\ \mathbf{0} & -(\mathbf{R}^\top\mathbf{g})^{\wedge}dt & \mathbf{0} \end{bmatrix}
$$

$$
\frac{\partial\boldsymbol{\tau}}{\partial\mathbf{i}} = \begin{bmatrix} \mathbf{0} & \frac{1}{2}dt^2\mathbf{I} \\ dt\mathbf{I} & \mathbf{0} \\ \mathbf{0} & dt\mathbf{I} \end{bmatrix}
$$

```cpp
/* chain rule */
Fx = J_x + J_u * J_tau_dx;
Fw = J_u * J_tau_i;
```

### Measurement Model

A general measurement model on manifold:

$$
\mathbf{y} = h(\mathbf{x}) + \mathbf{v}
$$

where $\mathbf{v}$ is measurement noise vector.

The measurement Jacobian is obtained via **chain rule**:

$$
\mathbf{H} = \frac{\partial h}{\partial \delta\mathbf{x}}\bigg|_{\mathbf{x}} = \underbrace{\frac{\partial h}{\partial \mathbf{x}_t}\bigg|_{\mathbf{x}}}_{\mathbf{H}_\mathbf{x}} \cdot \underbrace{\frac{\partial \mathbf{x}_t}{\partial \delta\mathbf{x}}\bigg|_{\mathbf{x}}}_{\mathbf{X}_{\delta\mathbf{x}}} = \mathbf{H}_\mathbf{x} \cdot \mathbf{X}_{\delta\mathbf{x}}
$$

$\mathbf{x}_t$ means true state vector, the condition is $\mathbf{x}_t = \mathbf{x}$ here, where:

| Symbol                          | Size                              | Meaning                                                  |
| ------------------------------- | --------------------------------- | -------------------------------------------------------- |
| $\mathbf{H}_\mathbf{x}$         | dim(measurement) $\times$ **Dim** | Standard EKF Jacobian                                    |
| $\mathbf{X}_{\delta\mathbf{x}}$ | Dim $\times$ DoF                  | ESKF-specific, depends on state (like $\mathbf{\theta}$) |

Given that measurement model is hard to "integrate" into a black box since their is too many kinds of sensors nowadays, so no `autoCompute()` function for `MeasurementModel` class. **Don't forget to compute** $\mathbf{X}_{\delta\mathbf{x}}$ **in your own model!**