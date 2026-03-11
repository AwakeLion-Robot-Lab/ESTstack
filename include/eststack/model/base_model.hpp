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
            typename T::Covariance;
            /* for dim(T_I(M)) != dim(control input), like SE_2(3) */
            typename T::NoiseJacobian;

            { model.compute(std::declval<typename T::State>(), std::declval<typename T::ControlInput>()) } -> std::same_as<typename T::State>;
            { model.computeStateJacobian(std::declval<typename T::State>(), std::declval<typename T::ControlInput>()) } -> std::same_as<typename T::Covariance>;
            { model.computeNoiseJacobian(std::declval<typename T::State>(), std::declval<typename T::ControlInput>()) } -> std::same_as<typename T::NoiseJacobian>;
        };

        /***
         * @brief measurement model concept
         */
        template <typename T>
        concept MeasurementModel = requires(T model) {
            typename T::State;
            typename T::Measurement;
            typename T::Covariance;
            typename T::MeasurementJacobian;

            { model.compute(std::declval<typename T::State>()) } -> std::same_as<typename T::Measurement>;
            { model.computeMeasJacobian(std::declval<typename T::State>()) } -> std::same_as<typename T::MeasurementJacobian>;
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
            using Scalar = typename State::Scalar;
            using ControlInput = typename Derived::ControlInput;
            using Covariance = typename Derived::Covariance;
            using NoiseJacobian = typename Derived::NoiseJacobian;

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
             * @brief compute the state jacobian (F_x)
             * @param state current state
             * @param u control input
             * @return state jacobian matrix
             */
            Covariance computeStateJacobian(const State &state, const ControlInput &u) const
            {
                return static_cast<const Derived *>(this)->computeStateJacobianImpl(state, u);
            }

            /***
             * @brief compute the noise jacobian (F_w)
             * @param state current state
             * @param u control input
             * @return noise jacobian matrix
             */
            NoiseJacobian computeNoiseJacobian(const State &state, const ControlInput &u) const
            {
                return static_cast<const Derived *>(this)->computeNoiseJacobianImpl(state, u);
            }

        protected:
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
            using Covariance = typename Derived::Covariance;
            using MeasurementJacobian = typename Derived::MeasurementJacobian;

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
             * @brief compute the measurement jacobian (H)
             * @param state current state
             * @return measurement jacobian matrix
             */
            MeasurementJacobian computeMeasJacobian(const State &state) const
            {
                return static_cast<const Derived *>(this)->computeMeasJacobianImpl(state);
            }

        protected:
            BaseMeasurementModel() = default;
        };

    }
}

#endif //! MODEL__BASE_MODEL_HPP