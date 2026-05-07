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

#ifndef MODEL__ARMORS_HPP
#define MODEL__ARMORS_HPP

// Eigen library
#include <Eigen/Core>

// manif library
#include <manif/Bundle.h>
#include <manif/SO2.h>
#include <manif/Rn.h>

// ESTstack library
#include "eststack/model/motion/base_motion_model.hpp"

/***
 * @brief An algorithm set focus on state estimation
 * @author jinhua "siyiovo" deng
 */
namespace eststack
{ /***
   * @brief models for problems
   */
    namespace model
    {
        /***
         * @brief details of motion models, not for public use
         */
        namespace details
        {
            using ArmorsState = manif::Bundle<double, manif::R3, manif::R3, manif::SO2, manif::R1, manif::R1, manif::R1, manif::R1>;
        } // namespace details

        /***
         * @brief RM four armors motion model
         * @details state: [x, y, z, vx, vy, vz, theta, omega, r1, r2, h]
         *              - x, y, z: car center position in nozzle frame
         *              - vx, vy, vz: car velocity in nozzle frame
         *              - theta: car yaw angle (SO(2) manifold) in nozzle frame
         *              - omega: car yaw angular velocity
         *              - r1: distances from car center to first type (0,2) armor center
         *              - r2: distances from car center to second type (1,3) armor center
         *              - h: height between two continuous armors center
         *          process noise: $\[\dot{vx}, \dot{vy}, \dot{vz}, \dot{omega}\]$
         *          control input: 0 input
         *
         *          reference: https://github.com/TongjiSuperPower/sp_vision_25/blob/main/tasks/auto_aim/target.cpp#L34
         */
        class ArmorsTransistionModel final : public BaseTransitionModel<ArmorsTransistionModel, details::ArmorsState, Eigen::Vector4d>
        {
        public:
            /***
             * @brief compute small increment via transition model
             * @param x current state
             * @param dt time step
             * @return small increment in tangent space
             */
            State::Tangent computeImpl(const State &x, double dt) const
            {
                /* get velocities */
                const auto v = x.element<1>().coeffs();
                const auto omega = x.element<3>().coeffs()(0);

                /* compute small increment */
                auto tau = State::Tangent::Zero();
                /* pos = v * dt */
                tau.coeffs().segment(0, 3) = v * dt;
                /* angle = omega * dt */
                tau.coeffs()(6) = omega * dt;
                return tau;
            }

            /***
             * @brief compute jacobians of the armors transition model
             * @param[in] x current state
             * @param[in] dt time step
             * @param[out] Fx state jacobian
             * @param[out] Fw noise jacobian
             * @note no control input
             */
            bool computeJacobiansImpl(const State &x, double dt, Eigen::Ref<StateJacobian> Fx, Eigen::Ref<NoiseJacobian> Fw) const
            {
                /* compute `J_tau_dx` */
                StateJacobian J_tau_dx = StateJacobian::Zero();
                J_tau_dx.block<3, 3>(0, 3) = Eigen::Matrix3d::Identity() * dt;
                J_tau_dx(6, 7) = dt;

                /* compute `J_tau_i` — piecewise white noise model */
                auto dt2_2 = dt * dt / 2.0;
                NoiseJacobian J_tau_i = NoiseJacobian::Zero();
                J_tau_i(0, 0) = dt2_2; /* x  <- ax */
                J_tau_i(1, 1) = dt2_2; /* y  <- ay */
                J_tau_i(2, 2) = dt2_2; /* z  <- az */
                J_tau_i(3, 0) = dt;    /* vx <- ax */
                J_tau_i(4, 1) = dt;    /* vy <- ay */
                J_tau_i(5, 2) = dt;    /* vz <- az */
                J_tau_i(6, 3) = dt2_2; /* theta <- alpha */
                J_tau_i(7, 3) = dt;    /* omega <- alpha */

                Fx = J_tau_dx;
                Fw = J_tau_i;
                return true;
            }
        };

        /***
         * @brief RM four armors motion model
         * @details state and noise are same as transition model
         *          measurement: [x, y, z, theta]
         *
         *          reference: https://github.com/TongjiSuperPower/sp_vision_25/blob/main/tasks/auto_aim/target.cpp#L34
         */
        class ArmorsMeasModel final : public BaseMeasurementModel<ArmorsMeasModel, details::ArmorsState, Eigen::Vector4d, Eigen::Vector4d>
        {
        public:
            /***
             * @brief compute expected measurement
             * @param x current state
             * @param dt time step
             */
            std::optional<Measurement> computeImpl(const State &x, double dt) const
            {
                const auto pos = x.element<0>().coeffs();
                const auto yaw = x.element<2>().angle();
                return Measurement{pos(0), pos(1), pos(2), yaw};
            }

            /***
             * @brief compute jacobians of the armors measurement model
             * @param[in] x current state
             * @param[in] dt time step
             * @param[out] H measurement jacobian
             * @param[out] Hv noise jacobian
             */
            bool computeJacobiansImpl(const State &x, double dt, Eigen::Ref<MeasJacobian> H, Eigen::Ref<NoiseJacobian> Hv) const
            {
                H = MeasJacobian::Zero();
                H(0, 0) = 1.0; /* pos_x */
                H(1, 1) = 1.0; /* pos_y */
                H(2, 2) = 1.0; /* pos_z */
                H(3, 6) = 1.0; /* yaw (SO2 tangent) */

                Hv = NoiseJacobian::Identity();
                return true;
            }
        };
    } // namespace model
} // namespace eststack

#endif //! MODEL__ARMORS_HPP
