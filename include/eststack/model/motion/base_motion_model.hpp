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

#ifndef MODEL__BASE_MOTION_MODEL_HPP
#define MODEL__BASE_MOTION_MODEL_HPP

// C++ standard library
#include <type_traits>
#include <tuple>

// Eigen library
#include <Eigen/Core>

// eststack library
#include "eststack/concepts.hpp"
#include "eststack/types.hpp"

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
         * @brief base class for transition model
         * @tparam Derived derived transition model class
         * @tparam StateT state type
         * @tparam ProcessNoiseT process noise vector type
         * @tparam ControlInputT control input type
         * @details CRTP design for compile-time polymorphism. `ControlInputT` is `void` as default
         */
        template <typename Derived, typename StateT, typename ProcessNoiseT, typename ControlInputT = void>
        class BaseTransitionModel
        {
        public:
            using State = StateT;
            using Scalar = typename State::Scalar;
            using ControlInput = ControlInputT;
            using ProcessNoise = ProcessNoiseT;
            using ProcessNoiseCovariance = eststack::Covariance<ProcessNoise>;
            using StateJacobian = typename State::Jacobian;
            using NoiseJacobian = eststack::Jacobian<State, ProcessNoise>;

            /***
             * @brief compute the transition model with autodiff
             * @param[in] x current state
             * @param[in] u control input
             * @param[in] dt time step
             * @param[out] Fx state jacobian
             * @param[out] Fw noise jacobian
             * @note ONLY with control input
             * @return priori state
             * @details explanation reference (1): https://github.com/artivis/manif/blob/devel/docs/explanation/autodiff.md
             *          explanation reference (2): docs/explanation/model.md
             */
            template <typename U = ControlInputT>
                requires(!std::is_void_v<U>)
            State autoCompute(const State &x, const U &u, double dt,
                              Eigen::Ref<StateJacobian> Fx, Eigen::Ref<NoiseJacobian> Fw) const
            {
                static_assert(eststack::AutoComputable<Derived>, "this is an non-auto-computable model, please modify other functions!");
                const typename State::Tangent tau = u * dt;
                return x.rplus(tau, Fx, Fw);
            }

            /***
             * @brief compute small increment via transition model
             * @param x current state
             * @param u control input
             * @param dt time step
             * @note with control input
             * @return small increment in tangent space
             */
            template <typename U = ControlInputT>
                requires(!std::is_void_v<U>)
            State::Tangent compute(const State &x, const U &u, double dt) const
            {
                return static_cast<const Derived *>(this)->computeImpl(x, u, dt);
            }

            /***
             * @brief compute jacobians of the transition model
             * @param x current state
             * @param u control input
             * @param dt time step
             * @note with control input
             * @return state jacobian and noise jacobian
             */
            template <typename U = ControlInputT>
                requires(!std::is_void_v<U>)
            auto computeJacobians(const State &x, const U &u, double dt) const
                -> std::tuple<StateJacobian, NoiseJacobian>
            {
                return static_cast<const Derived *>(this)->computeJacobiansImpl(x, u, dt);
            }

            /***
             * @brief compute small increment via transition model
             * @param x current state
             * @param dt time step
             * @note no control input
             * @return small increment in tangent space
             */
            template <typename U = ControlInputT>
                requires std::is_void_v<U>
            State::Tangent compute(const State &x, double dt) const
            {
                return static_cast<const Derived *>(this)->computeImpl(x, dt);
            }

            /***
             * @brief compute jacobians of the transition model
             * @param x current state
             * @param dt time step
             * @note no control input
             * @return state jacobian and noise jacobian
             */
            template <typename U = ControlInputT>
                requires std::is_void_v<U>
            auto computeJacobians(const State &x, double dt) const
                -> std::tuple<StateJacobian, NoiseJacobian>
            {
                return static_cast<const Derived *>(this)->computeJacobiansImpl(x, dt);
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
         * @tparam StateT state type (Lie group state)
         * @tparam MeasurementT measurement vector type
         * @tparam MeasNoiseT measurement noise vector type (used as dimension source)
         * @details CRTP design for compile-time polymorphism
         */
        template <typename Derived, typename StateT, typename MeasurementT, typename MeasNoiseT>
        class BaseMeasurementModel
        {
        public:
            using State = StateT;
            using Scalar = typename State::Scalar;
            using Measurement = MeasurementT;
            using MeasNoise = MeasNoiseT;
            using MeasNoiseCovariance = eststack::Covariance<MeasNoise>;
            using MeasJacobian = eststack::Jacobian<Measurement, State>;
            using NoiseJacobian = eststack::Jacobian<Measurement, MeasNoise>;

            /***
             * @brief compute expected measurement
             * @param x current state
             * @param dt time step
             */
            Measurement compute(const State &x, double dt) const
            {
                return static_cast<const Derived *>(this)->computeImpl(x, dt);
            }

            /***
             * @brief compute jacobians of the measurement model
             * @param x current state
             * @param dt time step
             * @return measurement jacobian and noise jacobian
             */
            auto computeJacobians(const State &x, double dt) const
                -> std::tuple<MeasJacobian, NoiseJacobian>
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

#endif //! MODEL__BASE_MOTION_MODEL_HPP