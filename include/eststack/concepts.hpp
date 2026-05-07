// Copyright 2026 siyiovo
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//    http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

#ifndef ESTSTACK__CONCEPTS_HPP
#define ESTSTACK__CONCEPTS_HPP

// C++ standard library
#include <type_traits>
#include <concepts>
#include <optional>
#include <utility>
#include <vector>

// Eigen library
#include <Eigen/Core>

/***
 * @brief An algorithm set focus on estimation and filtering
 * @author jinhua "siyiovo" deng
 */
namespace eststack
{
    /***
     * @brief DoF of lie group concept
     * @details for DoF extraction
     */
    template <typename T>
    concept LieGroupDoF = requires {
        typename T::Scalar;
        { T::DoF } -> std::convertible_to<int>;
        requires(T::DoF > 0);
    };

    /***
     * @brief Eigen row vector concept
     * @details for Eigen fixed-size row extraction
     */
    template <typename T>
    concept EigenRow = requires {
        typename T::Scalar;
        { T::RowsAtCompileTime } -> std::convertible_to<int>;
        requires(T::RowsAtCompileTime > 0);
    };

    /***
     * @brief static dimension at compile time concept
     * @details for dimension extraction at compile time
     */
    template <typename T>
    concept DimAtCompileTime = LieGroupDoF<T> || EigenRow<T>;

    /***
     * @brief sample consensus model concept
     */
    template <typename T>
    concept SACModel = requires(T model) {
        { model.setInputSource(std::declval<const Eigen::Matrix3Xf &>()) } -> std::same_as<void>;
        { model.setInputTarget(std::declval<const Eigen::Matrix3Xf &>()) } -> std::same_as<void>;
        { model.setInlierThreshold(std::declval<float>()) } -> std::same_as<void>;
        { model.getCloudSize() } -> std::convertible_to<int>;
        { model.getSampleSize() } -> std::convertible_to<int>;
        { model.getModelSize() } -> std::convertible_to<int>;
        { model.selectInliers(std::declval<const Eigen::VectorXf &>()) } -> std::same_as<std::vector<int>>;
        { model.fit(std::declval<const std::vector<int> &>(), std::declval<Eigen::VectorXf &>()) } -> std::same_as<bool>;
    };

    /**
     * @brief loss function concept
     */
    template <typename T>
    concept LossFunction = requires(T loss, float residual, float scale) {
        { loss.weight(residual, scale) } -> std::convertible_to<float>;
        { loss.isInlier(residual, scale) } -> std::convertible_to<bool>;
    };

    /***
     * @brief estimation problem concept
     */
    template <typename T>
    concept EstProblem = requires(T prob) {
        { prob.isInitialized() } -> std::convertible_to<bool>;
        { prob.setState(std::declval<const typename T::State &>()) } -> std::same_as<void>;
        { prob.getState() } -> std::same_as<const typename T::State &>;
        { prob.setSolution(std::declval<typename T::Solution>()) } -> std::same_as<void>;
        { prob.run() } -> std::same_as<bool>;
    };

    /***
     * @brief lie group state concept
     * @details constraint from lie group base defined by [manif](https://github.com/artivis/manif/blob/devel/include/manif/impl/lie_group_base.h)
     */
    template <typename T>
    concept LieGroupState = requires(T state) {
        typename T::Scalar;
        typename T::Tangent;
        typename T::DataType;
        typename T::Jacobian;

        { T::Dim } -> std::convertible_to<int>;
        { T::DoF } -> std::convertible_to<int>;
        { T::RepSize } -> std::convertible_to<int>;

        /* plus operation constraint: SO(3) x R^3 -> SO(3) */
        { state + std::declval<typename T::Tangent>() } -> std::same_as<T>;
        /* minus operation constraint: SO(3) x SO(3) -> R^3 */
        { state - std::declval<T>() } -> std::same_as<typename T::Tangent>;
    };

    /***
     * @brief compute with autodiff concept
     */
    template <typename T>
    concept AutoComputable = std::is_same_v<typename T::ControlInput, typename T::State::Tangent>;

    /***
     * @brief transition model concept
     */
    template <typename T>
    concept TransitionModel = requires {
        typename T::State;
        typename T::Scalar;
        typename T::ControlInput;
        typename T::ProcessNoise;
        typename T::ProcessNoiseCovariance;
        typename T::StateJacobian;
        /* for dim(T_I(M)) != dim(control input), like SE_2(3) and mostly Bundle */
        typename T::NoiseJacobian;
    } && (
                                  /* autonomous: `ControlInput` = `void` */
                                  (std::is_void_v<typename T::ControlInput> && requires(T model, typename T::State x, double dt) {
                                      { model.compute(x, dt) } -> std::convertible_to<typename T::State::Tangent>;
                                      {
                                          model.computeJacobians(x, dt,
                                                                 std::declval<Eigen::Ref<typename T::StateJacobian>>(),
                                                                 std::declval<Eigen::Ref<typename T::NoiseJacobian>>())
                                      } -> std::convertible_to<bool>;
                                  }) ||
                                  /* controlled: `ControlInput` is non-`void` */
                                  (!std::is_void_v<typename T::ControlInput> && requires(T model, typename T::State x, typename T::ControlInput u, double dt) {
                                      {
                                          model.autoCompute(x, u, dt,
                                                            std::declval<Eigen::Ref<typename T::StateJacobian>>(),
                                                            std::declval<Eigen::Ref<typename T::NoiseJacobian>>())
                                      } -> std::convertible_to<typename T::State>;
                                      { model.compute(x, u, dt) } -> std::convertible_to<typename T::State::Tangent>;
                                      {
                                          model.computeJacobians(x, u, dt,
                                                                 std::declval<Eigen::Ref<typename T::StateJacobian>>(),
                                                                 std::declval<Eigen::Ref<typename T::NoiseJacobian>>())
                                      } -> std::convertible_to<bool>;
                                  }));

    /***
     * @brief measurement model concept
     */
    template <typename T>
    concept MeasurementModel = requires(T model,
                                        typename T::State x,
                                        double dt) {
        typename T::State;
        typename T::Measurement;
        typename T::MeasNoise;
        typename T::MeasNoiseCovariance;
        typename T::MeasJacobian;
        /* for dim(measurement) != dim(measurement noise) */
        typename T::NoiseJacobian;

        { model.compute(x, dt) } -> std::convertible_to<std::optional<typename T::Measurement>>;
        {
            model.computeJacobians(x, dt,
                                   std::declval<Eigen::Ref<typename T::MeasJacobian>>(),
                                   std::declval<Eigen::Ref<typename T::NoiseJacobian>>())
        } -> std::convertible_to<bool>;
    };

    /***
     * @brief kalman filter concept
     */
    template <typename T>
    concept KalmanFilter = requires(T kf) {
        typename T::State;
        typename T::StateCovariance;

        { kf.setState(std::declval<const typename T::State &>()) } -> std::same_as<void>;
        { kf.getState() } -> std::same_as<const typename T::State &>;
        { kf.setStateCovariance(std::declval<Eigen::Ref<const typename T::StateCovariance>>()) } -> std::same_as<void>;
        { kf.getStateCovariance() } -> std::same_as<const typename T::StateCovariance &>;
    };

    /***
     * @brief point cloud registration concept
     */
    template <typename T>
    concept PCR = requires(T pcr,
                           const typename T::PointCloudConstPtr &cloud,
                           const Eigen::Matrix3Xf &source_match,
                           const Eigen::Matrix3Xf &target_match,
                           typename T::Result &result) {
        typename T::PointCloud;
        typename T::PointCloudPtr;
        typename T::PointCloudConstPtr;
        typename T::Result;

        { pcr.setInputSource(cloud) } -> std::same_as<void>;
        { pcr.setInputTarget(cloud) } -> std::same_as<void>;
        { pcr.setInputCorrespondences(source_match, target_match) } -> std::same_as<void>;
        { pcr.setVoxelResolution(std::declval<float>()) } -> std::same_as<void>;
        { pcr.align() } -> std::same_as<bool>;
        { pcr.evaluate() } -> std::same_as<void>;
        { pcr.getResult(result) } noexcept -> std::same_as<void>;
    };
}

#endif //! ESTSTACK__CONCEPTS_HPP
