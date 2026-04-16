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

#ifndef CORE__ESGVI_HPP
#define CORE__ESGVI_HPP

// Eigen library
#include <Eigen/Dense>

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
         * @brief Exactly Sparse Gaussian Variational Inference for parameter estimation
         * @details ESGVI for adaptive noise covariance estimation of KF
         */
        class ESGVI
        {
        public:
            EIGEN_MAKE_ALIGNED_OPERATOR_NEW;

            /**
             * @brief default constructor
             */
            ESGVI() = default;
        };

    } // namespace core
} // namespace eststack

#endif //! CORE__ESGVI_HPP
