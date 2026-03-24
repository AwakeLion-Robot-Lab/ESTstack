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

#ifndef MODEL__LANDMARK_SE2_HPP
#define MODEL__LANDMARK_SE2_HPP

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
         * @brief SE2 landmark measurement model
         * @details Measures landmark position in robot frame
         *          h(X, b) = X^{-1} * b, where X is robot pose (SE2) and b is landmark in world frame
         */
        class LandmarkSE2 : public BaseMeasurementModel<LandmarkSE2>
        {
        public:
            using State = manif::SE2d;
            using Measurement = Eigen::Vector2d;
            using MeasNoise = Eigen::Matrix2d;
            using MeasJacobian = Eigen::Matrix<double, 2, 3>;
            using NoiseJacobian = Eigen::Matrix2d;

            /***
             * @brief constructor with landmark position
             * @param landmark_position landmark position in world frame
             */
            explicit LandmarkSE2(const Eigen::Vector2d& landmark_position)
                : landmark_position_(landmark_position)
            {
            }

            /***
             * @brief compute the measurement model
             * @param x current state (robot pose)
             * @param dt time step (not used in this model)
             * @return expected measurement (landmark position in robot frame)
             */
            Measurement computeImpl(const State& x, const double& dt) const
            {
                // h(X, b) = X^{-1} * b
                // Use manif's inverse and act operations
                return x.inverse().act(landmark_position_);
            }

            /***
             * @brief compute the measurement jacobian H = dh/dx
             * @param x current state (robot pose)
             * @param dt time step (not used in this model)
             * @return measurement jacobian matrix
             * @note This computes the jacobian using the chain rule:
             *       H = J_h_xi * J_xi_x, where xi is the inverse pose
             */
            MeasJacobian computeMeasJacobianImpl(const State& x, const double& dt) const
            {
                // Compute jacobians using manif's built-in methods
                // First, get the jacobian of inverse: J_xi_x
                manif::SE2d::Jacobian J_xi_x;
                manif::SE2d x_inv = x.inverse(J_xi_x);

                // Then, get the jacobian of act: J_h_xi
                Eigen::Matrix<double, 2, 3> J_h_xi;
                x_inv.act(landmark_position_, J_h_xi);

                // Chain rule: H = J_h_xi * J_xi_x
                return J_h_xi * J_xi_x;
            }

            /***
             * @brief compute the noise jacobian H_w = dh/dw
             * @param x current state
             * @param dt time step (not used in this model)
             * @return noise jacobian matrix (identity for this model)
             */
            NoiseJacobian computeNoiseJacobianImpl(const State& x, const double& dt) const
            {
                // Measurement noise is additive: y = h(x) + w
                // So dh/dw = I
                return NoiseJacobian::Identity();
            }

            /***
             * @brief set landmark position
             * @param landmark_position landmark position in world frame
             */
            void setLandmark(const Eigen::Vector2d& landmark_position)
            {
                landmark_position_ = landmark_position;
            }

            /***
             * @brief get landmark position
             * @return landmark position in world frame
             */
            const Eigen::Vector2d& getLandmark() const
            {
                return landmark_position_;
            }

        private:
            Eigen::Vector2d landmark_position_; ///< landmark position in world frame
        };
    }
}

#endif //! MODEL__LANDMARK_SE2_HPP
