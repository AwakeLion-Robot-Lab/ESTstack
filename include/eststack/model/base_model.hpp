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

#ifndef MODEL__BASE_MODEL_HPP
#define MODEL__BASE_MODEL_HPP

// Eigen library
#include <Eigen/Dense>

// C++ standard library
#include <concepts>
#include <type_traits>

/***
 * @brief An algorithm set focus on estimation and filtering
 * @author jinhua "siyiovo" deng
 */
namespace eststack
{
    /***
     * @brief models for problems
     */
    namespace model
    {
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
            { state + std::declval<T::Tangent>() } -> std::same_as<T>;
            /* minus operation constraint: SO(3) x SO(3) -> R^3 */
            { state - std::declval<T>() } -> std::same_as<typename T::Tangent>;
        };

        /***
         * @brief transition model concept
         */
        template <typename T>
        concept TransitionModel = requires(T model) {
            typename T::State;
            typename T::ControlInput;
            typename T::StateJacobian;
            /* for dim(T_I(M)) != dim(control input), like SE_2(3) and mostly Bundle */
            typename T::NoiseJacobian;
            typename T::ProcessNoise;

            { model.compute(std::declval<typename T::State>(), std::declval<typename T::ControlInput>()) } -> std::same_as<typename T::State>;
            { model.computeStateJacobian(std::declval<typename T::State>(), std::declval<typename T::ControlInput>()) } -> std::same_as<typename T::StateJacobian>;
            { model.computeNoiseJacobian(std::declval<typename T::State>(), std::declval<typename T::ControlInput>()) } -> std::same_as<typename T::NoiseJacobian>;
        };

        /***
         * @brief measurement model concept
         */
        template <typename T>
        concept MeasurementModel = requires(T model) {
            typename T::State;
            typename T::Measurement;
            typename T::MeasJacobian;
            /* for dim(measurement) != dim(measurement noise) */
            typename T::NoiseJacobian;
            typename T::MeasNoise;

            { model.compute(std::declval<typename T::State>()) } -> std::same_as<typename T::Measurement>;
            { model.computeMeasJacobian(std::declval<typename T::State>()) } -> std::same_as<typename T::MeasJacobian>;
            { model.computeNoiseJacobian(std::declval<typename T::State>()) } -> std::same_as<typename T::NoiseJacobian>;
        };

        /***
         * @brief base class for transition model
         * @tparam Derived derived transition model class
         * @details CRTP design for compile-time polymorphism
         */
        template <typename Derived>
        class BaseTransitionModel
        {
        public:
            using State = typename Derived::State;
            using ControlInput = typename Derived::ControlInput;
            using StateJacobian = typename Derived::StateJacobian;
            using NoiseJacobian = typename Derived::NoiseJacobian;
            using ProcessNoise = Eigen::Matrix<typename State::Scalar,
                                               NoiseJacobian::ColsAtCompileTime,
                                               NoiseJacobian::ColsAtCompileTime>;

            /***
             * @brief compute the transition model
             * @param state current state
             * @param u control input
             * @return next state
             */
            State compute(const State &state, const ControlInput &u) const
            {
                return static_cast<const Derived *>(this)->computeImpl(state, u);
            }

            /***
             * @brief compute the state jacobian df/dx about state transition w.r.t nominal state
             * @param state current state
             * @param u control input
             * @return state jacobian matrix
             */
            StateJacobian computeStateJacobian(const State &state, const ControlInput &u) const
            {
                return static_cast<const Derived *>(this)->computeStateJacobianImpl(state, u);
            }

            /***
             * @brief compute the process noise jacobian (F_w)
             * @param state current state
             * @param u control input
             * @return noise jacobian matrix
             */
            NoiseJacobian computeNoiseJacobian(const State &state, const ControlInput &u) const
            {
                return static_cast<const Derived *>(this)->computeNoiseJacobianImpl(state, u);
            }

        protected:
            /***
             * @brief default constructor
             */
            BaseTransitionModel() = default;
        };

        /***
         * @brief base class for measurement model
         * @tparam Derived derived measurement model class
         * @details CRTP design for compile-time polymorphism
         */
        template <typename Derived>
        class BaseMeasurementModel
        {
        public:
            using State = typename Derived::State;
            using Measurement = typename Derived::Measurement;
            using MeasJacobian = typename Derived::MeasJacobian;
            using NoiseJacobian = typename Derived::NoiseJacobian;
            using MeasNoise = Eigen::Matrix<typename State::Scalar,
                                            NoiseJacobian::ColsAtCompileTime,
                                            NoiseJacobian::ColsAtCompileTime>;

            /***
             * @brief compute the measurement model
             * @param state current state
             * @return expected measurement
             */
            Measurement compute(const State &state) const
            {
                return static_cast<const Derived *>(this)->computeImpl(state);
            }

            /***
             * @brief compute the measurement jacobian dh/dx about measurement w.r.t nominal state
             * @param state current state
             * @return measurement jacobian matrix
             */
            MeasJacobian computeMeasJacobian(const State &state) const
            {
                return static_cast<const Derived *>(this)->computeMeasJacobianImpl(state);
            }

            /***
             * @brief compute the measurement noise jacobian (H_w)
             * @param state current state
             * @return measurement noise jacobian matrix
             */
            NoiseJacobian computeNoiseJacobian(const State &state) const
            {
                return static_cast<const Derived *>(this)->computeNoiseJacobianImpl(state);
            }

        protected:
            /***
             * @brief default constructor
             */
            BaseMeasurementModel() = default;
        };

    }
}

#endif //! MODEL__BASE_MODEL_HPP