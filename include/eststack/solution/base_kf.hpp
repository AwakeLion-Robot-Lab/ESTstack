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
#include <concepts>
#include <type_traits>
#include <utility>

// Eigen library
#include <Eigen/Dense>

/***
 * @brief An algorithm set focus on estimation and filtering
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
         * @brief base class for Kalman filter
         * @tparam Derived derived kalman filter class
         * @tparam StateT state type used by the filter
         * @details CRTP design for compile-time polymorphism, which is more efficient and flexibie
         */
        template <typename Derived, typename StateT>
        class BaseKF
        {
        public:
            using State = StateT;
            using Scalar = typename State::Scalar;
            using Covariance = Eigen::Matrix<Scalar, State::DoF, State::DoF>;

            /***
             * @brief set state vector
             */
            void setState(const State &state) noexcept
            {
                x_ = state;
            }

            /***
             * @brief get state vector
             */
            const State &getState() const noexcept
            {
                return x_;
            }

            /***
             * @brief set state covariance matrix
             */
            void setCovariance(const Covariance &cov) noexcept
            {
                P_ = cov;
            }

            /***
             * @brief get state covariance matrix
             */
            const Covariance &getCovariance() const noexcept
            {
                return P_;
            }

            /***
             * @brief prediction step
             * @tparam TransitionModel transition model of the system
             * @tparam ControlInput control input
             * @tparam Args additional arguments, e.g. bias or white noise, etc.
             * @param F transition matrix
             * @param u control input vector
             * @param args additional perturbation
             * @return priori state vector
             */
            template <typename TransitionModel, typename ControlInput, typename... Args>
            const State &predict(const TransitionModel &F, const ControlInput &u, Args &&...args)
            {
                return static_cast<Derived *>(this)->predictImpl(F, u, std::forward<Args>(args)...);
            }

            /***
             * @brief update step
             * @tparam MeasurementModel measurement model
             * @tparam Measurement measurement vector
             * @tparam Args additional arguments, e.g. bias or white noise, etc.
             * @param H measurement matrix
             * @param z measurement vector
             * @param args additional arguments
             * @return posteriori state vector
             */
            template <typename MeasurementModel, typename Measurement, typename... Args>
            const State &update(const MeasurementModel &H, const Measurement &z, Args &&...args)
            {
                return static_cast<Derived *>(this)->updateImpl(H, z, std::forward<Args>(args)...);
            }

        protected:
            /***
             * @brief default constructor
             * @details protected to prevent direct instantiation of base class, only allow derived classes to construct
             */
            BaseKF() = default;

            /***
             * @brief state vector
             */
            State x_;

            /***
             * @brief covariance matrix of state
             */
            Covariance P_ = Covariance::Identity();
        };

    }
}

#endif //! SOLUTION__BASE_KF_HPP
