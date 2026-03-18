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

#ifndef MODEL__CT_HPP
#define MODEL__CT_HPP

// Eigen library
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <stdexcept>

// manif library
#include <manif/Bundle.h>
#include <manif/Rn.h>
#include <manif/SO2.h>

// ESTstack library
#include "eststack/model/base_model.hpp"
#include "eststack/types.hpp"

/***
 * @brief An algorithm set focus on estimation and filtering
 * @author jinhua "siyiovo" deng
 */
namespace eststack
{ /***
   * @brief models for problems
   */
    namespace model
    {
        /***
         * @brief coordinated turn motion model
         * @details state: [x, y, vx, vy, theta, omega]
         *              - x, y: target position
         *              - vx, vy: target velocity
         *              - theta: target yaw angle (SO(2) manifold)
         *              - omega: yaw angular velocity
         *          control input: [x, y, theta]
         *          process noise: [vx, vy, omega]
         */
        class CT : public BaseTransitionModel<CT>
        {
        public:
            using State = manif::Bundle<double, manif::SO2, manif::R5>;
            using ProcessNoise = Eigen::Vector3d;
            using ControlInput = Eigen::Vector3d;

            /***
             * @brief compute the transition model
             * @param x current state
             * @param u control input
             * @param dt time step
             * @return priori state
             */
            State computeImpl(const State &x, const ControlInput &u, const double &dt) const
            {
                double dx = u(0);
                double dy = u(1);
                /* well, it actually is omega, but we denote it as dtheta for consistency */
                double dtheta = u(2);

                manif::SO2<double> so2 = x.element<0>();
                auto r5_vec = x.element<1>().coeffs();
                double px = r5_vec(0);
                double py = r5_vec(1);
                double vx = r5_vec(2);
                double vy = r5_vec(3);
                double omega = r5_vec(4);

                double next_px, next_py, next_vx, next_vy, next_omega;

                // Handle singular case when omega approaches 0
                if (std::abs(omega) < 1e-5)
                {
                    next_px = px + vx * dt + dx;
                    next_py = py + vy * dt + dy;
                    next_vx = vx;
                    next_vy = vy;
                    next_omega = omega;
                }
                else
                {
                    double sin_wt = std::sin(omega * dt);
                    double cos_wt = std::cos(omega * dt);
                    double inv_w = 1.0 / omega;

                    next_px = px + (vx * sin_wt - vy * (1.0 - cos_wt)) * inv_w + dx;
                    next_py = py + (vy * sin_wt + vx * (1.0 - cos_wt)) * inv_w + dy;
                    next_vx = vx * cos_wt - vy * sin_wt;
                    next_vy = vx * sin_wt + vy * cos_wt;
                    next_omega = omega;
                }

                manif::SO2<double> next_so2(so2.angle() + omega * dt + dtheta);

                Eigen::Matrix<double, 5, 1> next_r5_vec;
                next_r5_vec << next_px, next_py, next_vx, next_vy, next_omega;
                return State(next_so2, manif::R5<double>(next_r5_vec));
            }

            /***
             * @brief compute the state jacobian (F_x)
             * @param x current state
             * @param u control input
             * @param dt time step
             * @return state jacobian matrix
             */
            StateJacobian computeStateJacobianImpl(const State &x, const ControlInput &u, const double &dt) const
            {
                auto r5_vec = x.element<1>().coeffs();
                double vx = r5_vec(2);
                double vy = r5_vec(3);
                double omega = r5_vec(4);

                StateJacobian Fx = StateJacobian::Identity();

                // tangent states order: {dtheta, dx, dy, dvx, dvy, domega}
                // Fx.row(0) represents dtheta
                Fx(0, 5) = dt; // dtheta / domega

                if (std::abs(omega) < 1e-5)
                {
                    Fx(1, 3) = dt;
                    Fx(1, 4) = 0.0;
                    Fx(1, 5) = -vy * dt * dt / 2.0;

                    Fx(2, 3) = 0.0;
                    Fx(2, 4) = dt;
                    Fx(2, 5) = vx * dt * dt / 2.0;

                    Fx(3, 3) = 1.0;
                    Fx(3, 4) = 0.0;
                    Fx(3, 5) = -vy * dt;

                    Fx(4, 3) = 0.0;
                    Fx(4, 4) = 1.0;
                    Fx(4, 5) = vx * dt;
                }
                else
                {
                    double w_dt = omega * dt;
                    double sin_wt = std::sin(w_dt);
                    double cos_wt = std::cos(w_dt);
                    double inv_w = 1.0 / omega;
                    double inv_w2 = inv_w * inv_w;

                    Fx(1, 3) = sin_wt * inv_w;
                    Fx(1, 4) = -(1.0 - cos_wt) * inv_w;
                    Fx(1, 5) = vx * (dt * cos_wt * inv_w - sin_wt * inv_w2) - vy * (dt * sin_wt * inv_w - (1.0 - cos_wt) * inv_w2);

                    Fx(2, 3) = (1.0 - cos_wt) * inv_w;
                    Fx(2, 4) = sin_wt * inv_w;
                    Fx(2, 5) = vy * (dt * cos_wt * inv_w - sin_wt * inv_w2) + vx * (dt * sin_wt * inv_w - (1.0 - cos_wt) * inv_w2);

                    Fx(3, 3) = cos_wt;
                    Fx(3, 4) = -sin_wt;
                    Fx(3, 5) = -vx * dt * sin_wt - vy * dt * cos_wt;

                    Fx(4, 3) = sin_wt;
                    Fx(4, 4) = cos_wt;
                    Fx(4, 5) = vx * dt * cos_wt - vy * dt * sin_wt;
                }

                return Fx;
            }

            /***
             * @brief compute the noise jacobian (F_w)
             * @param x current state
             * @param u control input
             * @param dt time step
             * @return noise jacobian matrix
             */
            NoiseJacobian computeNoiseJacobianImpl(const State &x, const ControlInput &u, const double &dt) const
            {
                NoiseJacobian Fw = NoiseJacobian::Zero();
                // Map the 3D process noise [noise_vx, noise_vy, noise_omega]
                // directly into the tangent state increments {dtheta, dx, dy, dvx, dvy, domega}
                Fw(3, 0) = 1.0; // dvx / dnoise_vx
                Fw(4, 1) = 1.0; // dvy / dnoise_vy
                Fw(5, 2) = 1.0; // domega / dnoise_omega
                return Fw;
            }
        };
    }
}

#endif //! MODEL__CT_HPP