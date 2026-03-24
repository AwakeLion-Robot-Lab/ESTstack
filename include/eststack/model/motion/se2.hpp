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

#ifndef MODEL__MOTION_SE2_HPP
#define MODEL__MOTION_SE2_HPP

// Eigen library
#include <Eigen/Dense>

// manif library
#include <manif/SE2.h>

// ESTstack library
#include "eststack/model/base_model.hpp"
#include "eststack/types.hpp"

/***
 * @brief An algorithm set focus on estimation and filtering
 * @author jinhua "siyiovo" deng
 */
namespace eststack
{
    /***
     * @brief models for problems
     */
    namespace model
    {
        /***
         * @brief SE2 motion model
         * @details Simple motion model for SE2: X(k+1) = X(k) * exp(u)
         *          where u is a twist in se(2) tangent space
         *          This is the same as X(k+1) = X(k) ⊞ u (right-plus operation)
         */
        class MotionSE2 : public BaseTransitionModel<MotionSE2>
        {
        public:
            using State = manif::SE2d;
            using ProcessNoise = Eigen::Matrix3d;
            using ControlInput = manif::SE2Tangentd;
            using StateJacobian = Eigen::Matrix3d;
            using NoiseJacobian = Eigen::Matrix3d;

            /***
             * @brief compute the transition model
             * @param x current state
             * @param u control input (twist in se(2))
             * @param dt time step (not used in this simple model)
             * @return next state
             * @note This implements X(k+1) = X(k) * exp(u) = X(k) ⊞ u
             */
            State computeImpl(const State& x, const ControlInput& u, const double& dt) const
            {
                // Use manif's right-plus operation
                return x.rplus(u);
            }

            /***
             * @brief compute the state jacobian F_x = df/dx
             * @param x current state
             * @param u control input
             * @param dt time step (not used in this simple model)
             * @return state jacobian matrix
             * @note This computes the jacobian of the right-plus operation
             */
            StateJacobian computeStateJacobianImpl(const State& x, const ControlInput& u, const double& dt) const
            {
                // Jacobian of X ⊞ u with respect to X
                StateJacobian J_x;
                x.rplus(u, J_x);
                return J_x;
            }

            /***
             * @brief compute the noise jacobian F_w = df/dw
             * @param x current state
             * @param u control input
             * @param dt time step (not used in this simple model)
             * @return noise jacobian matrix
             * @note Process noise is applied as: X(k+1) = X(k) * exp(u + w)
             *       So we need the jacobian of exp(u) with respect to u
             */
            NoiseJacobian computeNoiseJacobianImpl(const State& x, const ControlInput& u, const double& dt) const
            {
                // Jacobian of X ⊞ u with respect to u
                StateJacobian J_u;
                manif::SE2d::Jacobian J_x_dummy;
                x.rplus(u, J_x_dummy, J_u);
                return J_u;
            }
        };
    }
}

#endif //! MODEL__MOTION_SE2_HPP
