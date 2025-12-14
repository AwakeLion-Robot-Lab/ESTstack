# ESTstack 实现规划

## 总体策略

遵循**自底向上**的实现顺序：Core → Problem → Solution

实现原则：
1. 每个模块独立测试后再进行上层开发
2. 优先实现最常用、最基础的功能
3. 使用 Eigen 作为核心线性代数库
4. 充分利用 Manif 处理流形上的操作
5. 在需要优化求解时引入 Ceres Solver

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

## Phase 2: Problem Layer - 问题接口定义

### 2.1 `problem/base_problem.hpp` - 基础问题接口
**目标**: 统一的问题定义抽象

**实现内容**:
```cpp
namespace eststack::problem {
  template<typename StateType, typename MeasurementType>
  class BaseProblem {
  public:
    using State = StateType;
    using Measurement = MeasurementType;

    virtual ~BaseProblem() = default;

    // 状态维度
    virtual size_t state_dim() const = 0;
    virtual size_t measurement_dim() const = 0;

    // 问题描述
    virtual std::string name() const = 0;
    virtual std::string description() const = 0;
  };
}
```

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

**用途**: ESEKF、SUKFM 的问题定义

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

### 3.1 `solution/base_kf.hpp` - Kalman滤波器基类
**目标**: 统一的KF接口

**实现内容**:
```cpp
namespace eststack::solution {
  template<typename ProblemType>
  class BaseKF {
  public:
    using State = typename ProblemType::State;
    using Measurement = typename ProblemType::Measurement;

    virtual ~BaseKF() = default;

    // 核心接口
    virtual void predict(double dt) = 0;
    virtual void update(const Measurement& z) = 0;

    // 状态访问
    virtual State get_state() const = 0;
    virtual Eigen::MatrixXd get_covariance() const = 0;

    // 初始化
    virtual void initialize(const State& x0, const Eigen::MatrixXd& P0) = 0;
  };
}
```

---

### 3.2 `solution/esekf.hpp` - Error-State Extended Kalman Filter
**目标**: 基于流形的ESEKF实现

**核心特性**:
- 使用 Manif 库处理 SE(3) 流形
- Error-state formulation（误差状态形式）
- 适用于IMU + GPS/视觉融合

**实现要点**:
```cpp
namespace eststack::solution {
  template<typename ProblemType>
  class ESEKF : public BaseKF<ProblemType> {
  private:
    State nominal_state_;       // 名义状态
    Eigen::VectorXd error_state_;  // 误差状态
    Eigen::MatrixXd P_;         // 误差状态协方差

  public:
    void predict(double dt) override;
    void update(const Measurement& z) override;

  private:
    // 误差状态的传播
    Eigen::MatrixXd compute_state_transition(double dt);
    // 注入误差状态到名义状态（流形操作）
    void inject_error_state();
  };
}
```

**依赖**: Manif, Eigen

**适用问题**: IMUFiltering, PoseEstimation

---

## Phase 4: Solution Layer - 算法实现

### 4.1 `solution/sukfm.hpp` - Scaled UKF on Manifold
**目标**: 在流形空间上实现的Scaled Unscented Kalman Filter

**核心特性**:
- 在非欧几里得流形（SE(3), SO(3)等）上进行状态估计
- 使用 Scaled Unscented Transformation 生成 sigma 点
- 利用 Manif 库处理流形上的加法/减法运算
- 相比ESEKF，不需要线性化，对非线性系统更精确

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

### 4.4 `solution/tcf_btc.hpp` - TCF with BTC Descriptor
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

  class TCF_BTC {
  public:
    // Tight Coupling Framework
    // 联合优化位姿和描述子匹配
    RegistrationResult align_with_optimization(
      const PointCloud& source,
      const PointCloud& target
    );

  private:
    // 使用 Ceres 求解非线性优化问题
    void solve_joint_optimization();
  };
}
```

**依赖**: Ceres Solver, Eigen

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
│   ├── test_esekf.cpp
│   ├── test_sukfm.cpp
│   ├── test_imm.cpp
│   └── test_icp.cpp
└── integration/
    └── test_vio_pipeline.cpp
```

### 5.2 示例程序
```
examples/
├── imu_esekf_demo.cpp        # IMU滤波演示
├── icp_registration_demo.cpp  # 点云配准
├── imm_tracking_demo.cpp      # 多模型跟踪
└── vio_pipeline.cpp           # 完整VIO管道
```

---

## 实现优先级排序

### 高优先级 (立即开始)
1. ✅ `core/kabsch.hpp` - 最基础的工具
2. ✅ `core/bar_itzhack.hpp` - 姿态表示基础
3. ✅ `problem/base_problem.hpp` - 接口定义
4. ✅ `solution/base_kf.hpp` - 滤波器基类

### 中优先级 (第二阶段)
5. `problem/imu_filtering.hpp`
6. `solution/esekf.hpp` - 最常用的滤波器
7. `solution/icp.hpp` - 经典算法

### 低优先级 (后期完善)
8. `core/quasar.hpp` - 特定应用
9. `solution/sukfm.hpp` - 高级滤波器
10. `solution/imm.hpp` - 特定场景
11. `solution/tcf_btc.hpp` - 研究性算法

---

## CMake 构建配置建议

```cmake
# CMakeLists.txt 需要添加：

# 依赖查找
find_package(Eigen3 REQUIRED)
find_package(manif REQUIRED)
find_package(Ceres REQUIRED)

# 可选依赖
find_package(GTest)  # 测试
find_package(benchmark)  # 性能测试

# 创建接口库
add_library(eststack_core INTERFACE)
target_include_directories(eststack_core INTERFACE include)
target_link_libraries(eststack_core INTERFACE
  Eigen3::Eigen
  manif::manif
  Ceres::ceres
)

# 如果有实现文件
add_library(eststack_impl
  src/core/kabsch.cpp
  src/solution/esekf.cpp
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
- **分层解耦**: Core不依赖Problem/Solution
- **模板化**: 适应不同的状态和测量类型
- **流形支持**: 使用Manif处理李群
- **优化驱动**: 在需要时用Ceres求解非线性优化

预计开发时间：
- Phase 1-2: 2-3周
- Phase 3: 2-3周
- Phase 4: 3-4周
- Phase 5: 1-2周

**总计**: 约2-3个月可完成基础框架。
