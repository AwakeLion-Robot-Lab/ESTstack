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

#ifndef SOLUTION__ESKFOM_HPP
#define SOLUTION__ESKFOM_HPP

// Eigen library
#include <Eigen/Dense>

// ESTstack library
#include "eststack/solution/base_kf.hpp"
#include "eststack/model/base_model.hpp"
#include "eststack/types.hpp"

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
         * @brief error state Kalman filter on manifold
         * @tparam StateT state type
         * @tparam P perturbation convention
         * @details according to (Li and Mourikis, 2012), global perturbation is more consistent with the state definition
         *          which transition matrix is more clear, here we use local perturbation as default
         * @note local perturbation on manifold as default
         */
        template <eststack::model::LieGroupState StateT,
                  eststack::Perturbation P = eststack::Perturbation::Local>
        class ESKFOM : public BaseKF<ESKFOM<StateT, P>, StateT>
        {
        public:
            EIGEN_MAKE_ALIGNED_OPERATOR_NEW;

            using Base = BaseKF<ESKFOM<StateT, P>, StateT>;
            using State = typename Base::State;
            using StateCovariance = typename Base::StateCovariance;

            /***
             * @brief default constructor
             */
            ESKFOM(const State &x0, const StateCovariance &P0)
            {
                this->setState(x0);
                this->setStateCovariance(P0);
            }

            /***
             * @brief ESKFOM prediction step implementation
             * @param model transition model providing Fx and Fw
             * @param u     control input
             * @param Q     process noise covariance
             * @param args  reserved for future extensions
             */
            template <eststack::model::TransitionModel TransitionModel, typename... Args>
            bool predictImpl(const TransitionModel &model,
                             const typename TransitionModel::ControlInput &u,
                             const typename TransitionModel::ProcessNoise &Q,
                             Args &&...args)
            {
                const auto Fx = model.computeStateJacobian(this->x_, u);
                const auto Fw = model.computeNoiseJacobian(this->x_, u);

                this->x_ = model.compute(this->x_, u);

                this->P_ = Fx * this->P_ * Fx.transpose();
                /* noalias() can fasten matrix operation while LHS do not appear in RHS */
                (this->P_).noalias() += Fw * Q * Fw.transpose();
            }

            /***
             * @brief ESKFOM update step implementation
             * @param model measurement model providing H, H_w
             * @param z     measurement vector
             * @param R     measurement noise covariance
             * @param args  reserved for future extensions
             */
            template <eststack::model::MeasurementModel MeasurementModel, typename... Args>
            bool updateImpl(const MeasurementModel &model,
                            const typename MeasurementModel::Measurement &z,
                            const typename MeasurementModel::MeasNoise &R,
                            Args &&...args)
            {
                const auto H = model.computeMeasJacobian(this->x_);
                const auto Hw = model.computeNoiseJacobian(this->x_);

                const auto y = (z - model.compute(this->x_)).eval();

                /***
                 * standard kalman gain K = PHT(HPHT + HwRHwT)^{-1}
                 * if measurement dimension >> state dimension (threshold is 6x here), (HPHT+R) is hard to get inverse matrix
                 * we can use SMW equation, a.k.a. matrix inversion lemma (please refer to chapter 2.2.15, state estimation for robotics, second edition)
                 * to get K = (P^{-1} + HT(HwRHw)^{-1} H)^{-1} HT(HwRHwT)^{-1} faster
                 */
                constexpr int MeasDim = MeasurementModel::MeasJacobian::RowsAtCompileTime;
                constexpr int StateDim = static_cast<int>(State::DoF);
                const auto K = [&]()
                {
                    if constexpr (MeasDim > 6 * StateDim)
                    {
                        /* inverse matrix of full-state measurement covariance */
                        const auto full_R_inv = (Hw * R * Hw.transpose()).inverse().eval();
                        return ((this->P_.inverse() + H.transpose() * full_R_inv * H).inverse() * H.transpose() * full_R_inv).eval();
                    }
                    else
                    {
                        const auto S_inv = (H * this->P_ * H.transpose() + Hw * R * Hw.transpose()).inverse().eval();
                        return (this->P_ * H.transpose() * S_inv).eval();
                    }
                }();

                /* compute error state and inject to nominal state */
                const typename State::Tangent dx(K * y);
                if constexpr (P == eststack::Perturbation::Global)
                    this->x_ = this->x_.lplus(dx);
                else
                    this->x_ = this->x_.rplus(dx);

                /* Joseph form covariance to keep `P_` positive semi-definite and symmetric */
                const StateCovariance IKH = StateCovariance::Identity() - K * H;
                const auto KHwRHwTKT = (K * Hw * R * Hw.transpose() * K.transpose()).eval();

                this->P_ = IKH * this->P_ * IKH.transpose();
                (this->P_).noalias() += KHwRHwTKT;

                /* last thing is reset expectation and covaraince of error state */
                /* if state includes quaternion, you need reset Jacobian = d(dx-E(dx))/d(dx) */
                const auto Ix = Eigen::Matrix<typename State::Scalar, StateDim, StateDim>::Identity();
                if constexpr (P == eststack::Perturbation::Global)
                {
                    const auto G = Ix + (0.5 * State::hat(dx));
                    this->P_ = G * this->P_ * G.transpose();
                }
                else
                {
                    const auto G = Ix - (0.5 * State::hat(dx));
                    this->P_ = G * this->P_ * G.transpose();
                }
            }
        };
    }
}

#endif //! SOLUTION__ESKFOM_HPP
