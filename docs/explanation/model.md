# Model Explanation

## auto-computable

The `autoComputable` concept identifies transition models where manif's `plus()` built-in Jacobian outputs can be **directly** used as $F_x$ and $F_i$ in ESKF covariance prediction:

$$\mathbf{P} \leftarrow \mathbf{F}_x \mathbf{P} \mathbf{F}_x^\top + \mathbf{F}_i \mathbf{Q}_i \mathbf{F}_i^\top$$

### General case

A general state transition on manifold can be written as:

$$\mathbf{x}^{-} = f(\mathbf{x}, \delta \mathbf{x}, \mathbf{u}_m, \mathbf{i}) = \mathbf{x} \oplus \boldsymbol{\tau}$$

where $\boldsymbol{\tau}\in \mathbf{T}_{I}\mathcal{M}$ is a tangent vector that may depend on state $\mathbf{x}$, control input $\mathbf{u}_m$ and noise $\mathbf{i}$.

Calling `x.plus(tau, J_x, J_u)` gives two DoF $\times$ DoF Jacobians:

- $J_x = \frac{\partial(\mathbf{x} \oplus \boldsymbol{\tau})}{\partial \mathbf{x}}$
- $J_u = \frac{\partial(\mathbf{x} \oplus \boldsymbol{\tau})}{\partial \boldsymbol{\tau}}$

The ESKF Jacobians are obtained via **chain rule**:

$$\mathbf{F}_x = \frac{\partial f}{\partial \delta \mathbf{x}} = J_x + J_u \cdot \frac{\partial \boldsymbol{\tau}}{\partial \delta\mathbf{x}}, \qquad \mathbf{F}_i =\frac{\partial f}{\partial \mathbf{i}} = J_u \cdot \frac{\partial \boldsymbol{\tau}}{\partial \mathbf{i}}$$

where:

| Symbol                                                                | Size                    | Meaning                                        |
| --------------------------------------------------------------------- | ----------------------- | ---------------------------------------------- |
| $J_x, J_u$                                                            | DoF $\times$ DoF        | from `plus()`                                  |
| $\frac{\partial \boldsymbol{\tau}}{\partial \delta\mathbf{x}}$        | DoF $\times$ DoF        | how $\boldsymbol{\tau}$ depends on error state |
| $\mathbf{G} = \frac{\partial \boldsymbol{\tau}}{\partial \mathbf{i}}$ | DoF $\times$ dim(noise) | noise-to-tangent mapping                       |

> [!NOTE]
> In manif, Dim = dim(Manifold) and DoF = dim(Tangent).

### auto-computable case

When **both** of the following hold:

1. `ControlInput = State::Tangent` (control lives in full tangent space)
2. $\boldsymbol{\tau} = \mathbf{u}$ (tangent vector is the control, independent of state)

the chain rule simplifies to:

$$\frac{\partial \boldsymbol{\tau}}{\partial \delta\mathbf{x}} = \mathbf{0} \implies \mathbf{F}_x = J_x, \qquad \mathbf{G} = \mathbf{I} \implies \mathbf{F}_i = J_u$$

So `x.plus(u, Fx, Fw)` directly gives both Jacobians. No manual computation needed.

**Typical examples:** models on $\boldsymbol{SO}(2)$, $\boldsymbol{SO}(3)$, $\boldsymbol{SE}(2)$, $\boldsymbol{SE}(3)$ with tangent-space control.

### Non-auto-computable case

When `dim(ControlInput) != DoF` or $\boldsymbol{\tau}$ depends on state, **you must hand-derive** $\frac{\partial \boldsymbol{\tau}}{\partial \delta\mathbf{x}}$ and $\mathbf{G}$, then compose with $J_u$ from `plus()`. Or simply hand-write $\mathbf{F}_x$ and $\mathbf{F}_i$ directly (often easier).

**Typical examples:**

| Model        | State                           | DoF | Control dim | Noise dim | Why not autoComputable                                                           |
| ------------ | ------------------------------- | --- | ----------- | --------- | -------------------------------------------------------------------------------- |
| IMU ESKF [4] | $\boldsymbol{SE}_2(3)$          | 9   | 6           | 12        | dim(control) $\neq$ DoF, $\boldsymbol{\tau}$ depends on $\mathbf{R}, \mathbf{b}$ |
| CT           | Bundle$\langle$SO2, R5$\rangle$ | 6   | 3           | 3         | dim(control) $\neq$ DoF, $\boldsymbol{\tau}$ depends on state                    |
| SingleArmor  | Bundle$\langle$SO2, R8$\rangle$ | 9   | 4           | 5         | dim(control) $\neq$ DoF                                                          |
