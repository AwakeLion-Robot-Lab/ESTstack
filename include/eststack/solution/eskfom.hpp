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

            using Ptr = std::shared_ptr<ESKFOM>;
            using ConstPtr = std::shared_ptr<const ESKFOM>;

            /***
             * @brief configuration options for ESKFOM
             */
            struct Options
            {
                bool use_joseph_form = true;     ///< use Joseph form for covariance update (more numerically stable)
                bool apply_reset_jacobian = true; ///< apply reset Jacobian correction after state update
            };

            /***
             * @brief default constructor
             */
            ESKFOM(const State &x0, const StateCovariance &P0, const Options &options = Options{})
                : options_(options)
            {
                this->setState(x0);
                this->setStateCovariance(P0);
            }

            /***
             * @brief get current options
             */
            const Options& getOptions() const { return options_; }

            /***
             * @brief set options
             */
            void setOptions(const Options &options) { options_ = options; }

            /***
             * @brief get the innovation (measurement residual) from last update
             * @note Only valid after calling update(). Returns the innovation y = z - h(x)
             */
            template <typename MeasurementT>
            MeasurementT getInnovation() const { return last_innovation_.template cast<typename MeasurementT::Scalar>(); }

            /***
             * @brief get the innovation covariance from last update
             * @note Only valid after calling update(). Returns S = H * P * H^T + R
             */
            template <int Dim>
            Eigen::Matrix<double, Dim, Dim> getInnovationCovariance() const
            {
                return last_innovation_cov_.template topLeftCorner<Dim, Dim>();
            }

        private:
            Options options_;
            Eigen::VectorXd last_innovation_;           ///< last innovation (for debugging/analysis)
            Eigen::MatrixXd last_innovation_cov_;       ///< last innovation covariance (for debugging/analysis)

            /***
             * @brief ESKFOM prediction step implementation
             * @param model transition model providing Fx and Fw
             * @param u     control input
             * @param Q     process noise covariance
             * @param dt    time step (default: 1.0)
             * @param args  reserved for future extensions
             */
            template <eststack::model::TransitionModel TransitionModel, typename... Args>
            bool predictImpl(const TransitionModel &model,
                             const typename TransitionModel::ControlInput &u,
                             const typename TransitionModel::ProcessNoise &Q,
                             const double &dt = 1.0,
                             Args &&...args)
            {
                /***
                 * by default, `manif` and model class provide Jacobians derived
                 * under local perturbation convention: x(k+1) = x(k) ⊕ δx_l
                 * So computeStateJacobian() returns F_raw = F_l
                 * When going to global perturbation, the error state is defined as:
                 * δx_g = Adj(x) δx_l. We can convert the local transition matrix F_l to global F_g:
                 * δx_l(k+1) = F_l * δx_l(k)
                 * Adj(x(k+1))^{-1} * δx_g(k+1) = F_l * Adj(x(k))^{-1} * δx_g(k) + F_{w,l} * w
                 * => δx_g(k+1) = [Adj(x(k+1)) * F_l * Adj(x(k))^{-1}] * δx_g(k) + [Adj(x(k+1)) * F_{w,l}] * w
                 */
                const auto Fx_raw = model.computeStateJacobian(this->x_, u, dt);
                const auto Fw_raw = model.computeNoiseJacobian(this->x_, u, dt);

                const State priori_x = model.compute(this->x_, u, dt);

                const auto Fx = [&]() -> StateCovariance
                {
                    if constexpr (P == eststack::Perturbation::Global)
                        return (priori_x.adj() * Fx_raw * this->x_.inverse().adj()).eval();
                    else
                        return Fx_raw;
                }();

                const auto Fw = [&]()
                {
                    if constexpr (P == eststack::Perturbation::Global)
                        return (priori_x.adj() * Fw_raw).eval();
                    else
                        return Fw_raw;
                }();

                this->x_ = priori_x;

                this->P_ = Fx * this->P_ * Fx.transpose();
                /* noalias() can fasten matrix operation while elements in LHS do not appear in RHS */
                (this->P_).noalias() += Fw * Q * Fw.transpose();

                return true;
            }

            /***
             * @brief ESKFOM update step implementation
             * @param model measurement model providing H, H_w
             * @param z     measurement vector
             * @param R     measurement noise covariance
             * @param dt    time step (default: 1.0)
             * @param args  reserved for future extensions
             */
            template <eststack::model::MeasurementModel MeasurementModel, typename... Args>
            bool updateImpl(const MeasurementModel &model,
                            const typename MeasurementModel::Measurement &z,
                            const typename MeasurementModel::MeasNoise &R,
                            const double &dt = 1.0,
                            Args &&...args)
            {
                /***
                 * similar to predict(), computeMeasJacobian() provides the raw measurement Jacobian
                 * H_raw = H_l derived under Local perturbation
                 * For Global perturbation, applying the relation δx_l = Adj(x)^{-1} δx_g:
                 * y = H_l * δx_l + R = H_l * Adj(x)^{-1} * δx_g + R
                 * => H_g = H_l * Adj(x)^{-1}
                 */
                const auto H_raw = model.computeMeasJacobian(this->x_, dt);
                const auto Hw = model.computeNoiseJacobian(this->x_, dt);
                const auto H = [&]()
                {
                    if constexpr (P == eststack::Perturbation::Global)
                        return (H_raw * this->x_.inverse().adj()).eval();
                    else
                        return H_raw;
                }();

                const auto y = (z - model.compute(this->x_, dt)).eval();

                // Store innovation for analysis (useful for debugging, as in manif examples)
                last_innovation_ = y;

                /***
                 * standard kalman gain is defined as K = PHT(HPHT + HwRHwT)^{-1}
                 * if measurement dimension >> state dimension (threshold is 6x here), (HPHT+R) is hard to get inverse matrix
                 * we can use SMW equation, a.k.a. matrix inversion lemma (please refer to chapter 2.2.15, [3])
                 * to get K = (P^{-1} + HT(HwRHw)^{-1} H)^{-1} HT(HwRHwT)^{-1} faster
                 */
                constexpr int MeasDim = MeasurementModel::MeasJacobian::RowsAtCompileTime;
                constexpr int StateDim = static_cast<int>(State::DoF);

                // Compute innovation covariance S = H * P * H^T + Hw * R * Hw^T
                const auto S = (H * this->P_ * H.transpose() + Hw * R * Hw.transpose()).eval();

                // Store innovation covariance for analysis (as separate computation in manif examples)
                last_innovation_cov_ = S;

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
                        return (this->P_ * H.transpose() * S.inverse()).eval();
                    }
                }();

                /***
                 * compute error state and inject to nominal state
                 * NOTE that dx is a random variable (discretize to sequence),
                 * `dx` represents its expectation in injection and reset operations
                 */
                const typename State::Tangent dx(K * y);
                if constexpr (P == eststack::Perturbation::Global)
                    this->x_ = this->x_.lplus(dx);
                else
                    this->x_ = this->x_.rplus(dx);

                /***
                 * covariance update: support both Joseph form (numerically stable) and standard form (simpler, as in manif examples)
                 * Joseph form: P = (I - KH) * P * (I - KH)^T + K * R * K^T
                 * Standard form: P = P - K * S * K^T, where S = H * P * H^T + R
                 */
                if (options_.use_joseph_form)
                {
                    /* Joseph form covariance to keep `P_` positive semi-definite and symmetric */
                    const StateCovariance IKH = StateCovariance::Identity() - K * H;
                    const auto KHwRHwTKT = (K * Hw * R * Hw.transpose() * K.transpose()).eval();

                    this->P_ = IKH * this->P_ * IKH.transpose();
                    (this->P_).noalias() += KHwRHwTKT;
                }
                else
                {
                    /* Standard form covariance update (as in manif se2_localization.cpp) */
                    this->P_ = this->P_ - K * S * K.transpose();
                }

                /***
                 * last thing is reset expectation and covariance of error state
                 * about why wedge operator could be replaced by small adjoint a.k.a. lie algebra adjoint, please refer to explanation/eskfom.md
                 * NOTE: this is optional and can be disabled for simple cases (as in manif se2_localization.cpp)
                 */
                if (options_.apply_reset_jacobian)
                {
                    const auto Ix = Eigen::Matrix<typename State::Scalar, StateDim, StateDim>::Identity();
                    if constexpr (P == eststack::Perturbation::Global)
                    {
                        const auto G = Ix + (0.5 * dx.smallAdj());
                        this->P_ = G * this->P_ * G.transpose();
                    }
                    else
                    {
                        const auto G = Ix - (0.5 * dx.smallAdj());
                        this->P_ = G * this->P_ * G.transpose();
                    }
                }

                return true;
            }
        };
    }
}

#endif //! SOLUTION__ESKFOM_HPP
