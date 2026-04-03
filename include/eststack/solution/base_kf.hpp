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

#ifndef SOLUTION__BASE_KF_HPP
#define SOLUTION__BASE_KF_HPP

// C++ standard library
#include <array>
#include <concepts>
#include <type_traits>
#include <utility>
#include <deque>
#include <cmath>

// Eigen library
#include <Eigen/Dense>

// ESTstack library
#include "eststack/types.hpp"
#include "eststack/eigen_traits.hpp"
#include "eststack/concepts.hpp"

/***
 * @brief An algorithm set focus on state estimation
 * @author jinhua "siyiovo" deng
 */
namespace eststack
{
    /***
     * @brief algorithms for problems
     */
    namespace solution
    {
        /***
         * @brief base class for kalman filter
         * @tparam Derived derived kalman filter class
         * @tparam StateT state type used by the filter
         * @details CRTP design for compile-time polymorphism, which is more efficient and flexible
         */
        template <typename Derived, typename StateT>
        class BaseKF
        {
        public:
            using State = StateT;
            using Scalar = typename State::Scalar;
            using StateCovariance = eststack::Covariance<State>;

            /***
             * @brief set state vector
             * @param state state vector
             */
            void setState(const State &state) noexcept
            {
                this->x_ = state;
            }

            /***
             * @brief get state vector
             */
            const State &getState() const noexcept
            {
                return this->x_;
            }

            /***
             * @brief set state covariance matrix
             * @param cov state covariance matrix
             */
            void setStateCovariance(const Eigen::Ref<const StateCovariance> &cov) noexcept
            {
                this->P_ = cov;
            }

            /***
             * @brief get state covariance matrix
             */
            const StateCovariance &getStateCovariance() const noexcept
            {
                return this->P_;
            }

            /***
             * @brief prediction step
             * @tparam TransitionModel transition model of the system
             * @tparam ControlInput control input
             * @tparam Args additional arguments, e.g. bias or white noise, etc.
             * @param F transition matrix
             * @param u control input vector
             * @param args additional arguments
             */
            template <typename TransitionModel, typename ControlInput, typename ProcessNoise, typename... Args>
            bool predict(const TransitionModel &model, const ControlInput &u, const ProcessNoise &Q, Args &&...args)
            {
                return static_cast<Derived *>(this)->predictImpl(model, u, Q, std::forward<Args>(args)...);
            }

            /***
             * @brief update step
             * @tparam MeasurementModel measurement model
             * @tparam Measurement measurement vector
             * @tparam MeasNoise measurement noise covariance
             * @tparam Args additional arguments
             * @param model measurement model
             * @param z measurement vector
             * @param R measurement noise covariance
             * @param args additional arguments
             */
            template <typename MeasurementModel, typename Measurement, typename MeasNoise, typename... Args>
            bool update(const MeasurementModel &model, const Measurement &z, const MeasNoise &R, Args &&...args)
            {
                return static_cast<Derived *>(this)->updateImpl(model, z, R, std::forward<Args>(args)...);
            }

        protected:
            /***
             * @brief default constructor
             * @details protected to prevent direct instantiation of base class, only allow derived classes to construct
             */
            BaseKF() = default;

            /***
             * @brief state vector
             * @details it represents nominal state in ESKFOM
             */
            State x_;

            /***
             * @brief state covariance matrix
             * @details it represents error covariance in ESKFOM which is defined on the tangent space of the manifold (homeomorphic to Euclidean space),
             *          since the "nominal covariance" is not well-defined on manifold
             */
            StateCovariance P_;

            /***
             * @brief priori NIS from update step
             * @details it just stores NIS while converged
             */
            std::deque<double> priori_nis_;

            /***
             * @brief maximum number of priori NIS
             */
            int max_priori_nis_num_;
        };

    }
}

#endif //! SOLUTION__BASE_KF_HPP
