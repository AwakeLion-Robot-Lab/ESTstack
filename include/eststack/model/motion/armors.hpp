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
#include <Eigen/Dense>

// manif library
#include <manif/Bundle.h>
#include <manif/Rn.h>
#include <manif/SO2.h>

// ESTstack library
#include "eststack/model/base_model.hpp"
#include "eststack/types.hpp"

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
         * @brief RM entire 4 armors motion model
         * @details state: [x, y, z, vx, vy, vz, theta, omega, r, l, h]
         *              - x, y, z: armor position in gun frame
         *              - vx, vy, vz: armor velocity in gun frame
         *              - theta: armor yaw angle (SO(2) manifold) in gun frame
         *              - omega: armor yaw angular velocity
         *              - r: distance from armor center to car center
         *              - l: minus distance of long axis and short axis of the car
         *              - h: height between two armors center
         *
         *          control input: [x, y, z, theta]
         *          process noise: [vx, vy, vz, omega]
         *
         *          reference: https://github.com/TongjiSuperPower/sp_vision_25/blob/main/tasks/auto_aim/target.cpp#L34
         */
        class Armors final : public BaseTransitionModel<Armors>
        {
        public:
            using State = manif::Bundle<double, manif::R3, manif::R3, manif::SO2, manif::R1, manif::R1, manif::R1, manif::R1>;
            /* x y z theta */
            using ControlInput = Eigen::Vector4d;
            /* sigma_vx sigma_vy sigma_vz sigma_omega */
            using ProcessNoise = Eigen::Vector4d;

            using Base = BaseTransitionModel<Armors>;
            using StateJacobian = typename Base::StateJacobian;
            using NoiseJacobian = typename Base::NoiseJacobian;

            [[deprecated]]
            State autoComputeImpl(const State &x, const ControlInput &u,
                                  Eigen::Ref<StateJacobian> Fx, Eigen::Ref<NoiseJacobian> Fw, const double &dt) const = delete;

            /***
             * @brief compute the state jacobian
             * @param x current state
             * @param u control input (dt)
             * @return state jacobian matrix
             */
            StateJacobian computeStateJacobianImpl(const State &x, const ControlInput &u, const double &dt)
            {
                StateJacobian Fx = StateJacobian::Identity();
                return Fx;
            }

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

#endif //! MODEL__ARMORS_HPP
