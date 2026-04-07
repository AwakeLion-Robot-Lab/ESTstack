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

#ifndef CORE__RANSAC_HPP
#define CORE__RANSAC_HPP

// C++ standard library
#include <cmath>
#include <limits>

// Eigen library
#include <Eigen/Dense>

// ESTstack library
#include "eststack/concepts.hpp"

/***
 * @brief An algorithm set focus on state estimation
 * @author jinhua "siyiovo" deng
 */
namespace eststack
{
    /***
     * @brief core algorithms for estimation
     */
    namespace core
    {
        /***
         * @brief RANSAC result
         */
        struct RANSACResult
        {
            RANSACResult(const Eigen::Isometry3d &T = Eigen::Isometry3d::Identity(),
                         int inlier_count = 0, int iterations = 0,
                         bool converged = false, double inlier_ratio = 0.0)
                : T_(T), inlier_count_(inlier_count), iterations_(iterations),
                  converged_(converged), inlier_ratio_(inlier_ratio) {}

            /***
             * @brief transformation
             */
            Eigen::Isometry3d T_;

            /***
             * @brief number of inliers
             */
            int inlier_count_;

            /***
             * @brief number of iterations performed
             */
            int iterations_;

            /***
             * @brief convergence flag
             */
            bool converged_;

            /***
             * @brief ratio of inliers to total points
             */
            double inlier_ratio_;
        };

        /***
         * @brief Random Sample Consensus for model fitting
         * @tparam model sample consensus model type
         */
        template <SACModel model>
        class RANSAC
        {
        public:
            /***
             * @brief default constructor
             * @param max_iter maximum iterations
             * @param confidence confidence level
             */
            explicit RANSAC(int max_iter = 1000, double confidence = 0.99)
                : max_iter_(max_iter), confidence_(confidence)
            {
            }

            /***
             * @brief run RANSAC optimization
             * @param source 3xN source points
             * @param target 3xN target points
             * @return RANSACResult
             */
            RANSACResult optimize(const Eigen::Matrix3Xd &source, const Eigen::Matrix3Xd &target, )
            {
            }

        private:
            /***
             * @brief maximum iterations
             */
            int max_iter_;

            /***
             * @brief confidence level
             */
            double confidence_;
        };

    } // namespace core
} // namespace eststack

#endif // CORE__RANSAC_HPP
