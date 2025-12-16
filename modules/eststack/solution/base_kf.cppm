// Copyright 2025 siyiovo
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

// export as module `BaseKF`
export module BaseKF;

// C++ standard library
#include <concepts>
#include <type_traits>

// Eigen library
#include <Eigen/Dense>

/***
 * @brief a common algorithm set for `Awakelion Robot Lab`
 */
export namespace eststack
{

    template <typename T>
    concept Scalar = std::is_floating_point_v<T>;

    template <typename T>
    concept EigenMatrix = requires {
        typename T::Scalar;
        requires std::is_arithmetic_v<typename T::Scalar>;
        { T::RowsAtCompileTime } -> std::convertible_to<int>;
        { T::ColsAtCompileTime } -> std::convertible_to<int>;
    };

    template <typename Model>
    concept StateSpaceModel = requires(const Model &m) {
        { Model::StateDim } -> std::convertible_to<int>;
        { Model::MeasDim } -> std::convertible_to<int>;

        typename Model::Scalar;
        typename Model::StateVec;
        typename Model::StateMat;
        typename Model::MeasVec;
        typename Model::MeasMat;
        typename Model::NoiseCov;
        typename Model::ObsCov;

        { m.F() } -> EigenMatrix;
        { m.H() } -> EigenMatrix;
        { m.Q() } -> EigenMatrix;
        { m.R() } -> EigenMatrix;
    };

    template <typename KF>
    concept BaseKF = requires(KF kf,
                              typename KF::Model model,
                              typename KF::Model::StateVec x0,
                              typename KF::Model::StateMat P0,
                              typename KF::Model::StateVec u,
                              typename KF::Model::MeasVec z) {
        requires StateSpaceModel<typename KF::Model>;
        { kf.reset(x0, P0) };
        { kf.predict(model, u) };
        { kf.update(model, z) };
        { kf.state() } -> std::same_as<const typename KF::Model::StateVec &>;
        { kf.cov() } -> std::same_as<const typename KF::Model::StateMat &>;
    };

    export template <StateSpaceModel M>
    class KalmanFilter
    {
    public:
        using Model = M;
        using Scalar = typename M::Scalar;
        using StateVec = typename M::StateVec;
        using StateMat = typename M::StateMat;
        using MeasVec = typename M::MeasVec;

        virtual ~KalmanFilter() = default;

        virtual void reset(const StateVec &x0, const StateMat &P0) = 0;

        virtual void predict(const Model &m, const StateVec &u) = 0;

        virtual void update(const Model &m, const MeasVec &z) = 0;

        virtual const StateVec &state() const noexcept = 0;

        virtual const StateMat &cov() const noexcept = 0;
    };

} // namespace eststack