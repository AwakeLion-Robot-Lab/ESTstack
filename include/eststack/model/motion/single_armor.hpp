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

#ifndef MODEL__RM_INFANTRY_HPP
#define MODEL__RM_INFANTRY_HPP

// Eigen library
#include <Eigen/Dense>

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
         * @brief RM single armor motion model
         * @details state: [theta | x, y, z, vx, vy, vz, omega, r]
         *          - theta: target yaw angle (SO(2) manifold)
         *          - x, y, z: target position in gun frame
         *          - vx, vy, vz: target velocity in gun frame
         *          - omega: yaw angular velocity
         *          - r: radius from rotation center to armor
         *
         *          dynamics:
         *          - position: circular arc around rotation center
         *            x' = x + vx*dt - r*(sin(theta + omega*dt) - sin(theta))
         *            y' = y + vy*dt + r*(cos(theta + omega*dt) - cos(theta))
         *            z' = z + vz*dt
         *          - theta' = theta (+) omega*dt  (SO(2) group operation)
         *          - vx, vy, vz, omega, r: nearly static (process noise driven)
         */
        class SingleArmor : public BaseTransitionModel<SingleArmor>
        {
        public:
            using State = manif::Bundle<double, manif::SO2, manif::R8>;
            using ControlInput = Eigen::Matrix<double, 1, 1>;
            using NoiseInput = Eigen::Vector<double, 5>;
            using StateJacobian = eststack::Jacobian<State, State>;
            using NoiseJacobian = eststack::Jacobian<State, NoiseInput>;

            /***
             * @brief compute the transition model
             * @param state current state
             * @param u control input (dt)
             * @return next state
             */
            State computeImpl(const State &state, const ControlInput &u) const;

            /***
             * @brief compute the state jacobian
             * @param state current state
             * @param u control input (dt)
             * @return state jacobian matrix
             */
            StateJacobian computeStateJacobianImpl(const State &state, const ControlInput &u) const;

            /***
             * @brief compute the noise jacobian
             * @param state current state
             * @param u control input (dt)
             * @return noise jacobian matrix
             */
            NoiseJacobian computeNoiseJacobianImpl(const State &state, const ControlInput &u) const;
        };
    }
}

#endif //! MODEL__RM_INFANTRY_HPP
