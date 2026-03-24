# 关于 Fx, Fw, H, Hw 的链式法则应用说明

## 问题回答

**问题：是否需要对于 Fx Fw 和 H Hw 应用链式法则？除了 ESKF 应该都不需要考虑扰动雅可比吧？**

**简短回答：**
- **Fx, Fw**: 模型实现中**通常需要**应用链式法则（除非是 auto-computable 的简单情况）
- **H, Hw**: 测量模型实现中**需要**应用链式法则（当观测方程是复合函数时）
- **扰动相关的 Adjoint 变换**: 这个**只在 ESKFOM 内部**处理，模型不需要关心

## 详细解释

### 1. 模型层（Motion/Measurement Model）的职责

#### 1.1 运动模型需要计算的 Jacobian

运动模型需要提供：
- **Fx = ∂f/∂δx**: 状态转移关于**误差状态**的 Jacobian
- **Fw = ∂f/∂w**: 状态转移关于**过程噪声**的 Jacobian

**关键点：** 模型返回的 Jacobian 应该是关于**误差状态 δx**（切空间）的，而不是关于名义状态 x（流形）的。

#### 案例分析

##### a) auto-computable 模型（简单情况）

如 SE2 运动模型 (`model/motion/se2.hpp`)：

```cpp
// 状态转移: X(k+1) = X(k) ⊞ u
// 其中 u ∈ se(2) 是切空间控制量

State computeImpl(const State& x, const ControlInput& u, const double& dt) const
{
    return x.rplus(u);  // X ⊞ u
}

StateJacobian computeStateJacobianImpl(const State& x, const ControlInput& u, const double& dt) const
{
    StateJacobian J_x;
    x.rplus(u, J_x);  // manif 自动计算 ∂(X ⊞ u)/∂X
    return J_x;       // 这就是 Fx，无需链式法则！
}

NoiseJacobian computeNoiseJacobianImpl(const State& x, const ControlInput& u, const double& dt) const
{
    StateJacobian J_u;
    manif::SE2d::Jacobian J_x_dummy;
    x.rplus(u, J_x_dummy, J_u);  // manif 自动计算 ∂(X ⊞ u)/∂u
    return J_u;                   // 这就是 Fw，无需链式法则！
}
```

**为什么不需要链式法则？**
- 控制量 u 直接是切空间元素
- τ = u（tangent vector 就是控制量本身）
- ∂τ/∂δx = 0（τ 不依赖状态）
- ∂τ/∂w = I（噪声直接作用在切空间）
- 因此 Fx = J_x, Fw = J_u，manif 的 `plus()` 输出可直接使用

##### b) 非 auto-computable 模型（需要链式法则）

如 CT 模型 (`model/motion/ct.hpp`)：

**问题在于：**
- 状态: Bundle<SO2, R5>（6 DoF）
- 控制: Vector3d（3 维）→ dim(control) ≠ DoF
- 实际的切空间扰动 τ 依赖于当前状态（vx, vy, omega）

**运动方程：**
```
theta(k+1) = theta(k) + omega * dt + dtheta
x(k+1) = x(k) + (vx * sin(ω*dt) - vy * (1-cos(ω*dt))) / ω + dx
y(k+1) = y(k) + (vy * sin(ω*dt) + vx * (1-cos(ω*dt))) / ω + dy
vx(k+1) = vx * cos(ω*dt) - vy * sin(ω*dt)
vy(k+1) = vx * sin(ω*dt) + vy * cos(ω*dt)
omega(k+1) = omega
```

**关键：** 位置和速度的更新依赖于**当前状态** (vx, vy, omega)，而不仅仅是控制输入！

因此，`computeStateJacobianImpl()` 中手工计算的 Fx 实际上**已经隐含应用了链式法则**：

```cpp
StateJacobian Fx = StateJacobian::Identity();

// 这些偏导数体现了链式法则：
Fx(1, 3) = sin_wt * inv_w;        // ∂x(k+1)/∂vx
Fx(1, 4) = -(1.0 - cos_wt) * inv_w;  // ∂x(k+1)/∂vy
Fx(1, 5) = ...;                    // ∂x(k+1)/∂omega (复杂表达式)
```

每一项都是通过对复合函数求导得到的。

#### 1.2 测量模型需要计算的 Jacobian

测量模型需要提供：
- **H = ∂h/∂δx**: 观测方程关于**误差状态**的 Jacobian
- **Hw = ∂h/∂v**: 观测方程关于**测量噪声**的 Jacobian

##### 案例：Landmark SE2 测量模型

观测方程：y = X^{-1} * b（地标在机器人坐标系下的位置）

这是一个**复合函数**：
1. 先求逆：ξ = X^{-1}
2. 再作用：y = ξ * b

**必须使用链式法则：**

```cpp
MeasJacobian computeMeasJacobianImpl(const State& x, const double& dt) const
{
    // 步骤1: 计算 X^{-1} 及其 Jacobian
    manif::SE2d::Jacobian J_xi_x;  // ∂(X^{-1})/∂X
    manif::SE2d x_inv = x.inverse(J_xi_x);

    // 步骤2: 计算 act 及其 Jacobian
    Eigen::Matrix<double, 2, 3> J_h_xi;  // ∂(ξ * b)/∂ξ
    x_inv.act(landmark_position_, J_h_xi);

    // 步骤3: 链式法则
    return J_h_xi * J_xi_x;  // H = ∂h/∂ξ · ∂ξ/∂X
}
```

这是 manif se2_localization.cpp 中的标准模式（第 151-153 行）。

### 2. ESKFOM 层的职责

ESKFOM 收到模型提供的 Fx, Fw, H, Hw 后，还需要进行**扰动约定转换**（如果使用 Global perturbation）。

#### 2.1 Predict 阶段的 Adjoint 变换

**代码位置：** `eskfom.hpp:140-154`

```cpp
// 模型提供的是局部扰动下的 Jacobian: F_l
const auto Fx_raw = model.computeStateJacobian(this->x_, u, dt);
const auto Fw_raw = model.computeNoiseJacobian(this->x_, u, dt);

// 如果使用全局扰动，需要转换
if constexpr (P == eststack::Perturbation::Global)
{
    Fx = priori_x.adj() * Fx_raw * this->x_.inverse().adj();
    Fw = priori_x.adj() * Fw_raw;
}
else
{
    Fx = Fx_raw;  // 局部扰动，直接使用
    Fw = Fw_raw;
}
```

**这不是链式法则，而是坐标变换：**
- 局部扰动：δx_l = Adj(X)^{-1} · δx_g
- 将局部扰动下的 F_l 转换为全局扰动下的 F_g

#### 2.2 Update 阶段的 Adjoint 变换

**代码位置：** `eskfom.hpp:189-195`

```cpp
const auto H_raw = model.computeMeasJacobian(this->x_, dt);

if constexpr (P == eststack::Perturbation::Global)
{
    H = H_raw * this->x_.inverse().adj();  // H_g = H_l * Adj(X)^{-1}
}
else
{
    H = H_raw;  // 局部扰动，直接使用
}
```

### 3. 总结：谁负责什么？

| 组件 | 需要应用链式法则？ | 需要考虑扰动约定？ | 说明 |
|------|-------------------|-------------------|------|
| **运动模型 Fx** | ✅ **需要**（非 auto-computable 情况） | ❌ **不需要** | 手工推导时已包含链式法则；auto-computable 时 manif 自动处理 |
| **运动模型 Fw** | ✅ **需要**（如果噪声映射复杂） | ❌ **不需要** | 简单情况下 Fw 常为单位阵或简单映射 |
| **测量模型 H** | ✅ **需要**（复合函数情况） | ❌ **不需要** | 如 SE2 landmark 需要链式法则 |
| **测量模型 Hw** | ⚠️ **通常不需要** | ❌ **不需要** | 多数情况下 Hw = I（加性噪声） |
| **ESKFOM Adjoint 变换** | ❌ **不是链式法则** | ✅ **需要** | 这是扰动约定的坐标变换，只在 ESKFOM 内部 |

### 4. 具体示例对比

#### 示例 1: SE2 模型（auto-computable）

```cpp
// 不需要手动链式法则，manif 的 plus() 自动处理
StateJacobian Fx;
x.rplus(u, Fx);  // Fx = ∂(X ⊞ u)/∂X，已经是关于误差状态的！
return Fx;
```

**原因：** τ = u（控制直接是切空间），manif 的 plus() 输出的 J_x 就是 ∂(x⊕u)/∂δx。

#### 示例 2: CT 模型（非 auto-computable）

```cpp
// 需要手动推导，因为位置更新依赖速度状态
Fx(1, 3) = sin(ω*dt) / ω;           // ∂x(k+1)/∂vx
Fx(1, 5) = vx*(dt*cos - sin/ω)/ω²   // ∂x(k+1)/∂ω
```

**原因：** 这些偏导数是通过对 x(k+1) = x(k) + f(vx, vy, ω, dt) 求导得到，本质上已经应用了链式法则。

#### 示例 3: Landmark SE2 测量（需要显式链式法则）

```cpp
// 观测方程: y = X^{-1} * b（复合函数）

// 步骤1: X^{-1} 的 Jacobian
manif::SE2d x_inv = x.inverse(J_xi_x);  // J_xi_x = ∂ξ/∂X

// 步骤2: act 的 Jacobian
x_inv.act(b, J_h_xi);  // J_h_xi = ∂y/∂ξ

// 步骤3: 链式法则
H = J_h_xi * J_xi_x;  // H = (∂y/∂ξ) · (∂ξ/∂X)
```

**原因：** 观测是两步操作的复合，必须显式应用链式法则组合 Jacobian。

### 5. 关于扰动约定的说明

**重要：** 模型实现时**不需要**考虑扰动约定（Global vs Local）！

扰动约定的处理**完全在 ESKFOM 内部**：

```cpp
// eskfom.hpp:140-154 (predict)
// eskfom.hpp:189-195 (update)

// ESKFOM 收到模型的 F_raw (局部扰动) 后：
if (Global perturbation) {
    Fx = Adj(x+) * Fx_raw * Adj(x)^{-1}  // 坐标变换
    Fw = Adj(x+) * Fw_raw
}
```

**这不是链式法则，而是误差状态的坐标变换：**
- δx_l = Adj(X)^{-1} · δx_g
- 将局部切空间的 Jacobian 转换到全局切空间

**模型开发者只需：**
1. 假设使用局部扰动（x+ = x ⊞ δx）
2. 计算 Fx = ∂f/∂δx 和 Fw = ∂f/∂w
3. 返回这些 Jacobian
4. ESKFOM 会根据配置自动转换

### 6. 实际建议

#### 开发运动模型时：

**情况 A: 控制量直接是切空间元素，且独立于状态**
- ✅ 使用 manif 的 `plus()` 自动计算
- ✅ 不需要手动链式法则
- 例如：SE2, SE3, SO3 的简单运动

**情况 B: 控制维度 ≠ DoF，或状态转移依赖当前状态**
- ⚠️ 必须手动推导 Fx 和 Fw
- ⚠️ 推导过程中会自然包含链式法则
- 例如：CT, SingleArmor, IMU 模型

#### 开发测量模型时：

**情况 A: 观测方程是简单函数**
- 可能不需要显式链式法则
- 例如：直接测量位置 h(x) = x.position()

**情况 B: 观测方程是复合函数**
- ✅ 必须应用链式法则
- ✅ 使用 manif 提供的中间 Jacobian
- 例如：Landmark SE2 (inverse + act)

### 7. 代码验证

让我们看看现有代码是否正确处理了这些情况：

#### CT 模型检查

```cpp
// ct.hpp:119-176
// computeStateJacobianImpl 手动计算了所有偏导数
Fx(1, 3) = sin_wt * inv_w;        // ∂x/∂vx
Fx(1, 5) = vx*(dt*cos - ...) ...  // ∂x/∂ω
```

✅ **正确** - 这些偏导数已经通过对复合函数求导得到，隐含应用了链式法则。

#### Landmark SE2 测量检查

```cpp
// landmark_se2.hpp:83-94
manif::SE2d x_inv = x.inverse(J_xi_x);     // 第一步
x_inv.act(landmark_position_, J_h_xi);     // 第二步
return J_h_xi * J_xi_x;                    // 链式法则组合
```

✅ **正确** - 显式应用了链式法则，与 manif 参考实现一致。

#### ESKFOM 扰动转换检查

```cpp
// eskfom.hpp:140-154
if constexpr (P == eststack::Perturbation::Global)
    return (priori_x.adj() * Fx_raw * this->x_.inverse().adj()).eval();
```

✅ **正确** - 这是扰动约定的 Adjoint 变换，不是链式法则。模型不需要关心这个。

## 最终结论

1. **Fx, Fw 的链式法则：**
   - Auto-computable 模型：manif 自动处理，模型代码简单
   - 非 auto-computable 模型：手动推导时自然包含，但需要正确推导所有偏导数

2. **H, Hw 的链式法则：**
   - 复合函数观测（如 inverse + act）：必须显式应用
   - 简单观测：可能不需要

3. **扰动约定（Adjoint 变换）：**
   - 只在 ESKFOM 内部处理
   - 模型实现时**完全不需要**考虑
   - 模型只需假设局部扰动，返回对应的 Jacobian

4. **实践建议：**
   - 开发模型时专注于正确的物理/几何推导
   - 不要在模型层混入扰动约定的考虑
   - 让 ESKFOM 负责扰动约定的转换
   - 使用 manif 提供的自动 Jacobian 功能尽可能简化代码

## 参考文献

- `docs/explanation/model.md` - autoComputable 概念详细说明
- Solà et al. (2018) - "A micro Lie theory for state estimation in robotics"
- manif 库 se2_localization.cpp - 链式法则标准实践
