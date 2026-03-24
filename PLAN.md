# ESTstack 实现规划

## 总体策略

遵循**自底向上**的实现顺序：Core → Problem → Solution

当前执行重点（按当前进度）：**ESKFOM 已完成并改进，当前实现 RM 运动模型 + ESKF Jacobian 转换 + 最小闭环测试**。

本迭代采用里程碑优先：为尽快形成可运行滤波闭环，暂时按 **Solution-first** 推进，再回补 Core 通用算法。

## 最新更新 (2026-03-24)

**ESKFOM 改进（基于 manif se2_localization.cpp 参考实现）：**
- ✅ 添加可配置的协方差更新方法（Joseph form vs 标准形式）
- ✅ 添加可选的 reset Jacobian 应用（可针对简单情况禁用）
- ✅ 修复 dt 参数传递 bug
- ✅ 添加 innovation 和 innovation covariance 跟踪
- ✅ 创建 SE2 运动模型 (`model/motion/se2.hpp`)
- ✅ 创建 SE2 地标测量模型 (`model/measurement/landmark_se2.hpp`)
- ✅ 创建 SE2 定位示例 (`examples/se2_localization_eskfom.cpp`)
- ✅ 更新文档（`docs/explanation/eskfom.md`, `eskfom_improvements.md`, `COMPARISON_MANIF.md`）

详见：`docs/explanation/eskfom_improvements.md` 和 `docs/explanation/COMPARISON_MANIF.md`

当前代码进度快照（基于最近提交与文件内容）：

**已完成：**
- `types.hpp` — concepts（LieGroupDoF / EigenRow / DimAtCompileTime）、Jacobian / Covariance 类型别名、Perturbation 枚举、getDimAtCompileTime consteval
- `model/base_model.hpp` — LieGroupState / TransitionModel / MeasurementModel concepts + BaseTransitionModel / BaseMeasurementModel CRTP 基类
- `solution/base_kf.hpp` — CRTP 基类，提供 predict / update 分发、state / covariance 存取
- `solution/eskfom.hpp` — **已改进** 完整 predict（Fx / Fw 协方差传播）+ update（SMW lemma 大维度优化、可配置 Joseph form / 标准协方差更新、可选 reset Jacobian via smallAdj、innovation 跟踪）
- `model/motion/se2.hpp` — **新增** SE2 运动模型（用于示例和测试）
- `model/motion/ct.hpp` — **已完成** CT 模型完整实现（Jacobians 已细化）
- `model/motion/single_armor.hpp` — **已完成** 单装甲板模型完整实现
- `model/measurement/landmark_se2.hpp` — **新增** SE2 地标测量模型（链式法则 Jacobian）
- `problem/base_problem.hpp` — EstProblem concept + BaseProblem CRTP 基类（isInitialized / setSolution 已实现）
- `docs/explanation/eskfom.md` — Reset Jacobian 推导说明 + 协方差更新方法配置指南
- `docs/explanation/eskfom_improvements.md` — **新增** ESKFOM 改进详细说明
- `docs/explanation/COMPARISON_MANIF.md` — **新增** ESKFOM 与 manif 参考实现详细对比
- `docs/explanation/references.md` — 15 条 BibTeX 参考文献
- `examples/se2_localization_eskfom.cpp` — **新增** SE2 定位完整示例
- `examples/README.md` — **新增** 示例说明文档

**骨架（类声明完成，compute 方法仅声明无实现）：**
- `solution/sukfom.hpp` — SUKFOM 类壳，仅类型别名，无 predict / update

**空文件 / 未开始：**
- `model/motion/ctrv.hpp`、`model/imu/base_imu_model.hpp`
- `problem/motion_tracking.hpp`、`problem/pose_estimation.hpp`
- `solution/btc_tcf.hpp`、`solution/imm.hpp`、`solution/icp.hpp`、`solution/inekf.hpp`、`solution/ndt.hpp`、`solution/epnp.hpp`
- `core/kabsch.hpp`、`core/bar_itzhack.hpp`、`core/quasar.hpp`

实现原则：
1. 每个模块独立测试后再进行上层开发
2. 优先实现最常用、最基础的功能
3. 使用 Eigen 作为核心线性代数库
4. 充分利用 Manif 处理流形上的操作
5. 使用 autodiff 库进行自动微分（forward mode dual numbers）

---

## Phase 1: Core Layer - 基础数学工具 (优先级: 最高)

### 1.1 `core/kabsch.hpp` - Kabsch 算法
**目标**: 3D点云配准的SVD解法

**实现内容**:
```cpp
namespace eststack::core {
  // 计算两组对应点云之间的最优旋转和平移
  // 输入: source points (3xN), target points (3xN)
  // 输出: R (3x3), t (3x1)
  struct KabschResult {
    Eigen::Matrix3d R;
    Eigen::Vector3d t;
    double rmse;
  };

  KabschResult kabsch(
    const Eigen::Matrix3Xd& source,
    const Eigen::Matrix3Xd& target,
    const Eigen::VectorXd& weights = Eigen::VectorXd()
  );
}
```

**依赖**: Eigen (SVD)

**用途**: ICP算法的核心组件

---

### 1.2 `core/bar_itzhack.hpp` - DCM转四元数
**目标**: 稳健的DCM到四元数转换

**实现内容**:
```cpp
namespace eststack::core {
  // Bar-Itzhack方法: 通过特征值分解稳健地从DCM提取四元数
  Eigen::Quaterniond dcm_to_quaternion(const Eigen::Matrix3d& dcm);

  // 批量转换
  std::vector<Eigen::Quaterniond> dcm_to_quaternion_batch(
    const std::vector<Eigen::Matrix3d>& dcms
  );
}
```

**优点**: 比直接提取更数值稳定

**用途**: IMU滤波、姿态估计

---

### 1.3 `core/quasar.hpp` - QUASAR算法
**目标**: 鲁棒的Wahba问题求解器（姿态估计）

**实现内容**:
```cpp
namespace eststack::core {
  // QUASAR: Quaternion-based Unscented Attitude Solver
  // 从向量观测估计姿态 (类似Wahba问题)
  struct QuasarResult {
    Eigen::Quaterniond q;
    Eigen::Matrix3d covariance;
    bool converged;
  };

  QuasarResult solve_wahba(
    const std::vector<Eigen::Vector3d>& reference_vectors,
    const std::vector<Eigen::Vector3d>& body_vectors,
    const std::vector<double>& weights
  );
}
```

**依赖**: Eigen

**用途**: 姿态初始化、星敏感器数据处理

---

## Phase 1.5: Model Layer - RM 运动模型 (优先级: 当前最高)

运动模型方向已明确为 **RoboMaster 枪口系（turret-frame）单板运动模型**，不再采用通用 CV/CT/CTRV。

### 状态定义

枪口系下装甲板状态 `[x, y, z, vx, vy, vz, theta, omega, r]`（9 维）：
- `(x, y, z)` — 装甲板位置
- `(vx, vy, vz)` — 装甲板速度
- `theta` — 偏航角
- `omega` — 偏航角速度
- `r` — 旋转半径

流形结构待确认（可能为 `Bundle<SO2d, Rn<double, 7>>` 或纯 `Rn<double, 9>`，取决于对 theta 周期性的处理需求）。

### 文件规划

| 文件 | 当前状态 | 说明 |
|------|----------|------|
| `model/motion/cv.hpp` | 骨架 | 改为 RM 匀速模型（omega=0 简化） |
| `model/motion/ct.hpp` | 骨架 | 改为 RM 匀转模型（omega != 0） |
| `model/motion/ctrv.hpp` | 空文件 | 可选：带变转速的扩展模型 |

每个模型需实现 `computeImpl` / `computeStateJacobianImpl` / `computeNoiseJacobianImpl`，继承自 `BaseTransitionModel<Derived>`。

**注意**：模型只需提供 df/dx（对 nominal state 的 Jacobian），ESKF 将在内部完成到 df/d(delta_x) 的转换（参见 3.2 节设计待办）。

---

## Phase 2: Problem Layer - 问题接口定义

### 2.1 `problem/base_problem.hpp` - 基础问题接口 ✅ 部分完成
**目标**: 统一的估计问题抽象（CRTP + concept 设计，非虚函数）

**实际实现**:
- `EstProblem` concept — 约束 `isInitialized()` 和 `setSolution()`
- `BaseProblem<Derived, SolutionT>` CRTP 基类 — 持有 `unique_ptr<SolutionT>`，提供 `isInitialized` / `setSolution`

---

### 2.2 `problem/imu_filtering.hpp` - IMU滤波问题
**实现内容**:
```cpp
namespace eststack::problem {
  struct IMUState {
    manif::SE3d pose;           // 位姿 (on SE(3) manifold)
    Eigen::Vector3d velocity;
    Eigen::Vector3d accel_bias;
    Eigen::Vector3d gyro_bias;
  };

  struct IMUMeasurement {
    double timestamp;
    Eigen::Vector3d accel;
    Eigen::Vector3d gyro;
  };

  class IMUFiltering : public BaseProblem<IMUState, IMUMeasurement> {
    // 定义IMU运动模型
    // 定义观测模型
  };
}
```

**用途**: ESKFOM、SUKFM 的问题定义

---

### 2.3 `problem/pointcloud_registration.hpp` - 点云配准问题
```cpp
namespace eststack::problem {
  using PointCloud = Eigen::Matrix3Xd;  // 3xN

  struct RegistrationResult {
    Eigen::Matrix4d transformation;  // SE(3)
    double fitness_score;
    size_t num_inliers;
    bool converged;
  };

  class PointCloudRegistration : public BaseProblem<PointCloud, PointCloud> {
    // ICP, TCF等算法的统一接口
  };
}
```

---

### 2.4 `problem/object_tracking.hpp` - 目标跟踪问题
```cpp
namespace eststack::problem {
  struct TrackingState {
    Eigen::VectorXd position;
    Eigen::VectorXd velocity;
    // 可选：加速度、角速度等
  };

  struct Detection {
    double timestamp;
    Eigen::VectorXd measurement;
    Eigen::MatrixXd covariance;
  };

  class ObjectTracking : public BaseProblem<TrackingState, Detection> {
    // 定义运动模型（CV, CA, CTRV等）
  };
}
```

**用途**: IMM 多模型跟踪

---

### 2.5 `problem/pose_estimation.hpp` - 位姿估计问题
```cpp
namespace eststack::problem {
  struct PoseState {
    manif::SE3d pose;
  };

  struct VisualMeasurement {
    std::vector<Eigen::Vector2d> image_points;
    std::vector<Eigen::Vector3d> world_points;
  };

  class PoseEstimation : public BaseProblem<PoseState, VisualMeasurement> {
    // EPnP等算法的接口
  };
}
```

---

## Phase 3: Solution Layer - 滤波器实现

### 3.1 `solution/base_kf.hpp` - Kalman 滤波器基类 ✅
**目标**: 统一的 KF 接口（CRTP 设计，零虚函数开销）

**实际实现**:
- `BaseKF<Derived, StateT>` CRTP 基类
- `predict()` / `update()` 分发至 Derived 的 `predictImpl()` / `updateImpl()`
- `setState` / `getState` / `setStateCovariance` / `getStateCovariance`
- 使用 `eststack::Covariance<State>` 编译期固定大小协方差矩阵

---

### 3.2 `solution/eskfom.hpp` - Error-State Kalman Filter on Manifold ✅
**目标**: 基于流形的 ESKFOM 实现

**实际实现**:
- `ESKFOM<StateT, P>` 模板类，`P` 控制 Global / Local 扰动约定
- `predictImpl()` — 通过 TransitionModel 获取 Fx / Fw，传播协方差 `P = Fx P FxT + Fw Q FwT`
- `updateImpl()` — 完整更新流程：
  - 当 MeasDim > 6 * StateDim 时自动切换 SMW lemma（matrix inversion lemma）计算 Kalman gain
  - 误差状态注入：Global 用 `lplus`，Local 用 `rplus`
  - Joseph form 协方差更新保证正半定
  - Reset Jacobian：`G = I +/- 0.5 * dx.smallAdj()`

**设计待办 — ESKF Jacobian 转换**:

当前 ESKFOM 假设运动模型直接提供误差状态 Jacobian（df/d(delta_x) 和 df/dw），但更通用的做法是：
- 运动模型只提供对 nominal state 的 Jacobian：df/dx（对用户更自然）
- ESKF 内部负责链式法则转换：`df/d(delta_x) = df/dx * dx/d(delta_x)`

其中 `dx/d(delta_x)` 取决于流形结构和扰动约定：
- Local 扰动：`dx/d(delta_x) = Jr(delta_x)`（右 Jacobian）
- Global 扰动：`dx/d(delta_x) = Jl(delta_x)`（左 Jacobian）

此设计将在后续迭代中实现，使模型层与滤波器层进一步解耦。

**依赖**: Manif, Eigen

**适用问题**: IMU 滤波、RM 目标跟踪、位姿估计

---

## Phase 4: Solution Layer - 算法实现

### 4.1 `solution/sukfm.hpp` - Scaled UKF on Manifold
**目标**: 在流形空间上实现的Scaled Unscented Kalman Filter

**核心特性**:
- 在非欧几里得流形（SE(3), SO(3)等）上进行状态估计
- 使用 Scaled Unscented Transformation 生成 sigma 点
- 利用 Manif 库处理流形上的加法/减法运算
- 相比ESKFOM，不需要线性化，对非线性系统更精确

**实现要点**:
```cpp
namespace eststack::solution {
  template<typename ProblemType>
  class SUKFM : public BaseKF<ProblemType> {
  private:
    State state_;           // 流形状态 (如 SE3d)
    Eigen::MatrixXd P_;     // 切空间上的协方差

    // Scaled UT 参数
    double alpha_;  // 控制sigma点分布 (typically 1e-3)
    double beta_;   // 先验知识参数 (2.0 for Gaussian)
    double kappa_;  // 二次参数 (typically 0 or 3-n)

  public:
    // 在流形上生成 scaled sigma 点
    std::vector<State> generate_scaled_sigma_points();

    // 流形上的加权平均（使用Manif的log/exp映射）
    State weighted_mean_on_manifold(
      const std::vector<State>& sigma_points,
      const Eigen::VectorXd& weights
    );

    // 预测步骤：在流形上传播sigma点
    void predict(double dt) override;

    // 更新步骤：测量更新
    void update(const Measurement& z) override;
  };
}
```

**依赖**: Manif, Eigen

**用途**: 高精度姿态估计、IMU滤波、机器人导航

---

### 4.2 `solution/imm.hpp` - Interacting Multiple Model
**目标**: 多模型自适应跟踪

**实现内容**:
```cpp
namespace eststack::solution {
  class IMM {
  private:
    std::vector<std::unique_ptr<BaseKF>> models_;  // 多个滤波器
    Eigen::VectorXd model_probs_;  // 模型概率
    Eigen::MatrixXd transition_matrix_;  // 模型转移矩阵

  public:
    void add_model(std::unique_ptr<BaseKF> model);

    void predict(double dt);
    void update(const Measurement& z);

  private:
    void interaction();  // 模型交互
    void mode_probability_update();  // 模型概率更新
    void state_fusion();  // 状态融合
  };
}
```

**用途**: 机动目标跟踪（如无人机、车辆）

---

### 4.3 `solution/icp.hpp` - Iterative Closest Point
**目标**: 经典点云配准算法

**实现内容**:
```cpp
namespace eststack::solution {
  class ICP {
  public:
    struct Config {
      size_t max_iterations = 50;
      double tolerance = 1e-6;
      double max_correspondence_distance = 0.5;
      bool use_point_to_plane = false;
    };

    RegistrationResult align(
      const PointCloud& source,
      const PointCloud& target,
      const Eigen::Matrix4d& initial_guess = Eigen::Matrix4d::Identity()
    );

  private:
    // 使用 KDTree 查找最近邻
    std::vector<int> find_correspondences(
      const PointCloud& source,
      const PointCloud& target
    );

    // 调用 Kabsch 算法求解变换
    Eigen::Matrix4d estimate_transformation(/*...*/);
  };
}
```

**依赖**: core/kabsch.hpp, nanoflann (KDTree)

---

### 4.4 `solution/btc_tcf.hpp` - BTC-augmented TCF
**目标**: 使用BTC描述子的紧耦合框架

**实现内容**:
```cpp
namespace eststack::solution {
  // BTC: Binary Tree Coding descriptor
  class BTCDescriptor {
  public:
    Eigen::VectorXf compute(const PointCloud& local_patch);
    double match_score(const Eigen::VectorXf& desc1, const Eigen::VectorXf& desc2);
  };

  class BTC_TCF {
  public:
    // Tight Coupling Framework
    // 联合优化位姿和描述子匹配
    RegistrationResult align_with_optimization(
      const PointCloud& source,
      const PointCloud& target
    );

  private:
    // 使用 autodiff 或手写优化求解非线性问题
    void solve_joint_optimization();
  };
}
```

**依赖**: Eigen（Ceres Solver 暂未引入，如需非线性优化后续评估）

**特点**: 比ICP更鲁棒，适合退化环境

---

## Phase 5: 测试与示例

### 5.1 单元测试 (使用 Google Test)
```
tests/
├── core/
│   ├── test_kabsch.cpp
│   ├── test_bar_itzhack.cpp
│   └── test_quasar.cpp
├── solution/
│   ├── test_eskfom.cpp
│   ├── test_sukfm.cpp
│   ├── test_imm.cpp
│   ├── test_icp.cpp
│   └── test_btc_tcf.cpp
└── integration/
    └── test_vio_pipeline.cpp
```

### 5.2 示例程序
```
examples/
├── imu_eskfom_demo.cpp       # IMU滤波演示
├── icp_registration_demo.cpp  # 点云配准
├── imm_tracking_demo.cpp      # 多模型跟踪
├── btc_tcf_demo.cpp           # BTC-TCF 配准演示
└── vio_pipeline.cpp           # 完整VIO管道
```

---

## 实现优先级排序

### P0（当前迭代：RM 运动模型 + Jacobian 转换 + 闭环测试）
1. ✅ `types.hpp` - concepts、类型别名、Perturbation 枚举
2. ✅ `model/base_model.hpp` - 概念约束 + CRTP 基类
3. ✅ `solution/base_kf.hpp` - CRTP 基类
4. ✅ `solution/eskfom.hpp` - 完整 predict / update（SMW lemma、Joseph form、reset Jacobian）
5. ✅ `problem/base_problem.hpp` - EstProblem concept + BaseProblem CRTP 基类（部分完成）
6. `model/motion/cv.hpp` - RM 枪口系匀速模型实现（computeImpl + Jacobian）
7. `model/motion/ct.hpp` - RM 枪口系匀转模型实现（computeImpl + Jacobian）
8. `solution/eskfom.hpp` 增加 Jacobian 转换 — ESKF 内部 `df/d(delta_x) = df/dx * dx/d(delta_x)`
9. `test/test_eskfom_minimal.cpp` - 最小闭环测试（编译 + 一次 predict / update）

### P1（ESKFOM 后：BTC_TCF）
10. `solution/btc_tcf.hpp` - 建立 BTC 描述子 + TCF 接口骨架并接入优化入口
11. `test/test_btc_tcf.cpp` - 接口与基础数据流测试

### P2（回补 Core 与扩展算法）
12. `core/kabsch.hpp` - SVD 配准基础工具
13. `core/bar_itzhack.hpp` - 稳健 DCM→Quaternion
14. `solution/icp.hpp` - 依赖 kabsch 的经典配准
15. `core/quasar.hpp` / `solution/sukfm.hpp` / `solution/imm.hpp` - 后续增强

---

## CMake 构建配置建议

```cmake
# CMakeLists.txt 需要添加：

# 依赖查找
find_package(Eigen3 REQUIRED)
find_package(manif REQUIRED)
find_package(autodiff REQUIRED)

# 可选依赖
find_package(GTest)  # 测试
find_package(benchmark)  # 性能测试

# 创建接口库
add_library(eststack_core INTERFACE)
target_include_directories(eststack_core INTERFACE include)
target_link_libraries(eststack_core INTERFACE
  Eigen3::Eigen
  manif::manif
  autodiff::autodiff
)

# 如果有实现文件
add_library(eststack_impl
  src/core/kabsch.cpp
  src/solution/eskfom.cpp
  src/solution/btc_tcf.cpp
  # ...
)

# 示例和测试
if(BUILD_EXAMPLES)
  add_subdirectory(examples)
endif()

if(BUILD_TESTING)
  add_subdirectory(tests)
endif()
```

---

## 开发建议

### 编码规范
- 命名空间: `eststack::{core|problem|solution}`
- 使用 C++20 特性: concepts, ranges
- 模板编程: 尽量使用 CRTP 避免虚函数开销
- 常量表达式: 多用 `constexpr`

### 性能考虑
- 使用 Eigen 的表达式模板避免临时对象
- 滤波器状态用 `Eigen::Ref` 避免拷贝
- 关键路径考虑 SIMD 优化

### 文档
- 每个函数添加 Doxygen 注释
- 在 `docs/` 目录添加理论推导和使用教程

---

## 总结

这个框架的设计理念是：
- **分层解耦**: Core 不依赖 Problem / Solution
- **模板化**: 适应不同的状态和测量类型
- **流形支持**: 使用 Manif 处理李群
- **自动微分**: 使用 autodiff 库（forward mode dual numbers）
- **CRTP 优先**: 全链路编译期多态，零虚函数开销
