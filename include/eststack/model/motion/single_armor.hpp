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

#ifndef MODEL__SINGLE_ARMOR_HPP
#define MODEL__SINGLE_ARMOR_HPP

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
         * @details state: [x, y, z, vx, vy, vz, theta, omega, r]
         *              - x, y, z: target position in gun frame
         *              - vx, vy, vz: target velocity in gun frame
         *              - theta: target yaw angle (SO(2) manifold)
         *              - omega: yaw angular velocity
         *              - r: radius from rotation center to armor
         *          control input: [x, y, z, theta]
         *          process noise: [vx, vy, vz, omega, r]
         */
        class SingleArmor : public BaseTransitionModel<SingleArmor>
        {
        public:
            using State = manif::Bundle<double, manif::SO2, manif::R8>;
            using ControlInput = Eigen::Vector4d;
            using ProcessNoise = Eigen::Vector<double, 5>;

            /***
             * @brief compute the transition model
             * @param x current state
             * @param u control input (dt)
             * @return next state
             */
            State computeImpl(const State &x, const ControlInput &u) const;

            /***
             * @brief compute the state jacobian
             * @param x current state
             * @param u control input (dt)
             * @return state jacobian matrix
             */
            StateJacobian computeStateJacobianImpl(const State &x, const ControlInput &u) const;

            /***
             * @brief compute the noise jacobian
             * @param x current state
             * @param u control input (dt)
             * @return noise jacobian matrix
             */
            NoiseJacobian computeNoiseJacobianImpl(const State &x, const ControlInput &u) const;
        };
    }
}

#endif //! MODEL__SINGLE_ARMOR_HPP
