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
#include <optional>
#include <utility>
#include <deque>

// Eigen library
#include <Eigen/Core>

// ESTstack library
#include "eststack/types.hpp"

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
             * @tparam ProcessNoiseCovariance process noise covariance
             * @tparam Args additional arguments, e.g. bias or white noise, etc.
             * @param F transition matrix
             * @param u control input vector
             * @param args additional arguments
             * @note with control input
             */
            template <typename TransitionModel, typename ControlInput, typename ProcessNoiseCovariance, typename... Args>
                requires(!std::is_void_v<ControlInput>)
            bool predict(const TransitionModel &model, const ControlInput &u, const ProcessNoiseCovariance &Q, Args &&...args)
            {
                return static_cast<Derived *>(this)->predictImpl(model, u, Q, std::forward<Args>(args)...);
            }

            /***
             * @brief prediction step
             * @tparam TransitionModel transition model of the system
             * @tparam ProcessNoiseCovariance process noise covariance
             * @tparam Args additional arguments, e.g. bias or white noise, etc.
             * @param F transition matrix
             * @param args additional arguments
             * @note no control input
             */
            template <typename TransitionModel, typename ProcessNoiseCovariance, typename... Args>
                requires std::is_void_v<typename TransitionModel::ControlInput>
            bool predict(const TransitionModel &model, const ProcessNoiseCovariance &Q, Args &&...args)
            {
                return static_cast<Derived *>(this)->predictImpl(model, Q, std::forward<Args>(args)...);
            }

            /***
             * @brief update step
             * @tparam MeasurementModel measurement model
             * @tparam Measurement measurement vector
             * @tparam MeasNoiseCovariance measurement noise covariance
             * @tparam Args additional arguments
             * @param model measurement model
             * @param z measurement vector
             * @param R measurement noise covariance
             * @param args additional arguments
             */
            template <typename MeasurementModel, typename Measurement, typename MeasNoiseCovariance, typename... Args>
            std::optional<double> update(const MeasurementModel &model, const Measurement &z, const MeasNoiseCovariance &R, Args &&...args)
            {
                return static_cast<Derived *>(this)->updateImpl(model, z, R, std::forward<Args>(args)...);
            }

        protected:
            /***
             * @brief default constructor
             * @param max_priori_nis_num maximum number of priori NIS to keep for convergence check
             * @details protected to prevent direct instantiation of base class, only allow derived classes to construct
             */
            BaseKF(int max_priori_nis_num = 100) : max_priori_nis_num_(max_priori_nis_num) {}

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
    } // namespace solution
} // namespace eststack

#endif //! SOLUTION__BASE_KF_HPP
