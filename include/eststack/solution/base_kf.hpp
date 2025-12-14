// Copyright 2025 siyiovo
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

#ifndef SOLUTION__BASE_KF_HPP
#define SOLUTION__BASE_KF_HPP

// Eigen library
#include <Eigen/Core>

/***
 * @brief a common algorithm set for `Awakelion Robot Lab`
 */
namespace eststack
{
    /***
     * @brief solutions namespace for solving specific problems
     */
    namespace solution
    {
        /***
         * @brief Base class for Kalman Filter
         */
        class BaseKF
        {
        public:
            /***
             * @brief constructor
             */
            explicit BaseKF();

            /***
             * @brief predict next state
             */
            virtual void predict() = 0;

            /***
             * @brief Update the state with measurement
             * @param z measurement vector
             */
            virtual void update(const Eigen::VectorXd &z) = 0;

        protected:
            /***
             * @brief State vector
             */
            Eigen::VectorXd state_;
        };
    } // namespace solution
} // namespace eststack

#endif // SOLUTION__BASE_KF_HPP