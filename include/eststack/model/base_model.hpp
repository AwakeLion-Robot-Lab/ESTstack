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
#include <tuple>

// eststack library
#include "eststack/types.hpp"
#include "eststack/eigen_traits.hpp"

/***
 * @brief An algorithm set focus on state estimation
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
        concept TransitionModel = requires(T model) {
            typename T::State;
            typename T::Scalar;
            typename T::ControlInput;
            typename T::ProcessNoise;
            typename T::StateJacobian;
            /* for dim(T_I(M)) != dim(control input), like SE_2(3) and mostly Bundle */
            typename T::NoiseJacobian;

            {
                model.autoCompute(std::declval<typename T::State>(), std::declval<typename T::ControlInput>(), std::declval<double>(),
                                  std::declval<Eigen::Ref<typename T::StateJacobian>>(), std::declval<Eigen::Ref<typename T::NoiseJacobian>>())
            } -> std::same_as<typename T::State>;
            { model.compute(std::declval<typename T::State>(), std::declval<typename T::ControlInput>(), std::declval<double>()) } -> std::same_as<typename T::State::Tangent>;
            { model.computeJacobians(std::declval<typename T::State>(), std::declval<typename T::ControlInput>(), std::declval<double>()) } -> std::same_as<std::tuple<typename T::StateJacobian, typename T::NoiseJacobian>>;
        };

        /***
         * @brief measurement model concept
         */
        template <typename T>
        concept MeasurementModel = requires(T model) {
            typename T::State;
            typename T::Measurement;
            typename T::MeasNoise;
            typename T::MeasJacobian;
            /* for dim(measurement) != dim(measurement noise) */
            typename T::NoiseJacobian;

            { model.compute(std::declval<typename T::State>(), std::declval<double>()) } -> std::same_as<typename T::Measurement>;
            { model.computeJacobians(std::declval<typename T::State>(), std::declval<double>()) } -> std::same_as<std::tuple<typename T::MeasJacobian, typename T::NoiseJacobian>>;
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
            using ProcessNoise = typename Derived::ProcessNoise;
            using StateJacobian = typename State::Jacobian;
            using NoiseJacobian = eststack::Jacobian<State, ProcessNoise>;

            /***
             * @brief compute the transition model with autodiff
             * @param[in] x current state
             * @param[in] u control input
             * @param[out] Fx state jacobian
             * @param[out] Fw noise jacobian
             * @param[in] dt time step
             * @return priori state
             * @details explanation reference (1): https://github.com/artivis/manif/blob/devel/docs/explanation/autodiff.md
             *          explanation reference (2): docs/explanation/model.md
             */
            State autoCompute(const State &x, const ControlInput &u,
                              Eigen::Ref<StateJacobian> Fx, Eigen::Ref<NoiseJacobian> Fw, const double &dt) const
            {
                static_assert(AutoComputable<Derived>, "this is an non-auto-computable model, please modify other functions!");
                const typename State::Tangent tau = u * dt;
                return x.rplus(tau, Fx, Fw);
            }

            /***
             * @brief compute small increment via transition model
             * @param[in] x current state
             * @param[in] u control input
             * @param[in] dt time step
             * @return small increment in tangent space
             */
            State::Tangent compute(const State &x, const ControlInput &u, const double &dt) const
            {
                return static_cast<const Derived *>(this)->computeImpl(x, u, dt);
            }

            /***
             * @brief compute jacobians of the transition model
             * @param[in] x current state
             * @param[in] u control input
             * @param[in] dt time step
             * @return state jacobian and noise jacobian
             */
            auto computeJacobians(const State &x, const ControlInput &u, const double &dt) -> std::tuple<StateJacobian, NoiseJacobian>
            {
                return static_cast<const Derived *>(this)->computeJacobiansImpl(x, u, dt);
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
            using Scalar = typename State::Scalar;
            using Measurement = typename Derived::Measurement;
            using MeasNoise = typename Derived::MeasNoise;
            using MeasJacobian = eststack::Jacobian<Measurement, State>;
            using NoiseJacobian = eststack::Jacobian<Measurement, MeasNoise>;

            /***
             * @brief compute expected measurement
             * @param[in] x current state
             * @param[in] dt time step
             */
            Measurement compute(const State &x, const double &dt) const
            {
                return static_cast<const Derived *>(this)->computeImpl(x, dt);
            }

            /***
             * @brief compute jacobians of the measurement model
             * @param[in] x current state
             * @param[in] dt time step
             * @return measurement jacobian and noise jacobian
             */
            auto computeJacobians(const State &x, const double &dt) -> std::tuple<MeasJacobian, NoiseJacobian>
            {
                return static_cast<const Derived *>(this)->computeJacobiansImpl(x, dt);
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