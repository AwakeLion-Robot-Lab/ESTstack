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
#include <manif/Rn.h>
#include <manif/SO2.h>

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
             * @brief compute jacobians of the transition model
             * @param[in] x current state
             * @param[in] dt time step
             * @return state jacobian and noise jacobian
             */
            auto computeJacobiansImpl(const State &x, double dt) const
                -> std::tuple<StateJacobian, NoiseJacobian>
            {
                /* compute `J_tau_dx` */
                StateJacobian J_tau_dx = StateJacobian::Zero();
                J_tau_dx.topLeftCorner(3, 3) = Eigen::Matrix3d::Identity() * dt;
                J_tau_dx(6, 7) = dt;

                /* compute `J_tau_i` */
                NoiseJacobian J_tau_i = NoiseJacobian::Zero();

                return std::make_tuple(J_tau_dx, J_tau_i);
            }
        };

        /***
         * @brief RM four armors motion model
         * @details state and process noise are same as transition model
         *          measurement: [x, y, z, theta]
         *
         *          reference: https://github.com/TongjiSuperPower/sp_vision_25/blob/main/tasks/auto_aim/target.cpp#L34
         */
        class ArmorsMeasModel final : public BaseMeasurementModel<ArmorsMeasModel, details::ArmorsState, Eigen::Vector4d, Eigen::Vector4d>
        {
            /***
             * @brief compute expected measurement
             * @param x current state
             * @param dt time step
             */
            Measurement computeImpl(const State &x, double dt) const
            {
            }

            /***
             * @brief compute jacobians of the measurement model
             * @param x current state
             * @param dt time step
             * @return measurement jacobian and noise jacobian
             */
            auto computeJacobiansImpl(const State &x, double dt) const
                -> std::tuple<MeasJacobian, NoiseJacobian>
            {
            }
        };
    } // namespace model
} // namespace eststack

#endif //! MODEL__ARMORS_HPP
