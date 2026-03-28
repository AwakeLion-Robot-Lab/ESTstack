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
#include <Eigen/Cholesky>
#include <Eigen/Dense>

// ESTstack library
#include "eststack/solution/base_kf.hpp"
#include "eststack/model/base_model.hpp"

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
         * @brief error state Kalman filter on manifold
         * @tparam StateT state type
         * @details according to (Li and Mourikis, 2012), global perturbation is more consistent with the state definition
         *          which transition matrix is more clear, here we use local perturbation as default
         * @note local perturbation on manifold as default
         */
        template <eststack::model::LieGroupState StateT>
        class ESKFOM : public BaseKF<ESKFOM<StateT>, StateT>
        {
        public:
            EIGEN_MAKE_ALIGNED_OPERATOR_NEW;

            using Base = BaseKF<ESKFOM<StateT>, StateT>;
            using State = typename Base::State;
            using Scalar = typename State::Scalar;
            using Tangent = typename State::Tangent;
            using StateCovariance = typename Base::StateCovariance;

            using Ptr = std::shared_ptr<ESKFOM>;
            using ConstPtr = std::shared_ptr<const ESKFOM>;
            using UpdateResult = typename Base::UpdateResult;

            /***
             * @brief default constructor
             */
            ESKFOM(const State &x0, const Eigen::Ref<const StateCovariance> &P0)
            {
                this->setState(x0);
                this->setStateCovariance(P0);
            }

            /***
             * @brief ESKFOM prediction step implementation
             * @param model transition model providing Fx and Fw
             * @param u     control input
             * @param Q     process noise covariance
             * @param dt    time step
             */
            template <eststack::model::TransitionModel TransitionModel, typename... Args>
            bool predictImpl(const TransitionModel &model,
                             const typename TransitionModel::ControlInput &u,
                             const typename TransitionModel::ProcessNoise &Q,
                             const double &dt)
            {
                const Tangent tau = model.compute(this->x_, u, dt);
                const auto [J_tau_dx, J_tau_i] = model.computeJacobians(this->x_, u, dt);

                typename TransitionModel::StateJacobian J_x, J_u, Fx;
                this->x_ = this->x_.rplus(tau, J_x, J_u);

                Fx.noalias() = (J_x + J_u * J_tau_dx).eval();
                const auto Fi = (J_u * J_tau_i).eval();

                this->P_ = Fx * this->P_ * Fx.transpose();
                /* noalias() can fasten matrix operation while elements in LHS do not appear in RHS */
                (this->P_).noalias() += Fi * Q * Fi.transpose();

                return true;
            }

            /***
             * @brief ESKFOM update step implementation
             * @param model measurement model providing H, H_w
             * @param z     measurement vector
             * @param R     measurement noise covariance
             * @param dt    time step
             */
            template <eststack::model::MeasurementModel MeasurementModel, typename... Args>
            UpdateResult updateImpl(const MeasurementModel &model,
                                    const typename MeasurementModel::Measurement &z,
                                    const typename MeasurementModel::MeasNoise &R, const double &dt)
            {
                const auto z_pred = model.compute(this->x_, dt);
                /* PLEASE refer to docs/explanation/model.md to figure out how to compute H */
                const auto [H, Hv] = model.computeJacobians(this->x_, dt);

                /* priori innovation */
                const auto inno = (z - z_pred).eval();
                /* priori innovation covariance */
                const auto S = (H * this->P_ * H.transpose() + Hv * R * Hv.transpose()).eval();

                /***
                 * standard kalman gain is defined as K = PHT(HPHT + HvRHvT)^{-1}
                 * if measurement dimension >> state dimension (threshold is 3 times here), (HPHT+R) is hard to get inverse matrix
                 * we can use SMW equation, a.k.a. matrix inversion lemma (please refer to chapter 2.2.15, [3])
                 * to get K = (P^{-1} + HT(HvRHv)^{-1} H)^{-1} HT(HvRHvT)^{-1} faster
                 */
                constexpr int StateDim = State::DoF;
                constexpr int MeasDim = MeasurementModel::MeasJacobian::RowsAtCompileTime;

                const auto K = [&]()
                {
                    if constexpr (MeasDim >= 3 * StateDim)
                    {
                        /* inverse matrix of full-state measurement covariance */
                        const auto full_R_inv = (Hv * R * Hv.transpose()).inverse().eval();
                        return ((this->P_.inverse() + H.transpose() * full_R_inv * H).inverse() * H.transpose() * full_R_inv).eval();
                    }
                    else
                    {
                        const auto S_inv = S.inverse().eval();
                        return (this->P_ * H.transpose() * S_inv).eval();
                    }
                }();

                /***
                 * compute error state and inject to nominal state
                 * NOTE that `tau` is a random variable (discretize to sequence),
                 * `tau` represents its expectation in injection and reset operations
                 */
                const Tangent tau(K * inno);
                this->x_ = this->x_.rplus(tau);

                /* Joseph form covariance to keep `P_` positive semi-definite and symmetric */
                const StateCovariance IKH = StateCovariance::Identity() - K * H;
                const auto J = (K * Hv * R * Hv.transpose() * K.transpose()).eval();

                this->P_ = IKH * this->P_ * IKH.transpose();
                (this->P_).noalias() += J;

                /***
                 * reset expectation and covariance of error state
                 * about why wedge operator could be replaced by small adjoint a.k.a. lie algebra adjoint, please refer to explanation/eskfom.md
                 */
                const auto Ix = Eigen::Matrix<Scalar, StateDim, StateDim>::Identity();
                const auto G = Ix - (0.5 * tau.smallAdj());
                this->P_ = G * this->P_ * G.transpose();

                /* last thing is check whether converge via NIS */
                if constexpr (MeasDim < 10)
                {
                    /***
                     * NOTE that we can use priori state N(x-, P-) at this point instead of batch estimation N(x+, P+) from all measurements
                     * details can be found in [3], chapter 5.1.2
                     */
                    const double nis = inno.transpose() * S.ldlt().solve(inno);
                    const bool converge = nis < ChiSquareTable::value<MeasDim>();
                    return {converge, nis};
                }
                else /* if measurement dimension is too high, did not calculate NIS */
                {
                    return {true, 0.0};
                }
            }
        };
    }
}

#endif //! SOLUTION__ESKFOM_HPP
