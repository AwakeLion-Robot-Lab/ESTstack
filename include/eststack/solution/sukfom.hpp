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

#ifndef SOLUTION__SUKFOM_HPP
#define SOLUTION__SUKFOM_HPP

// Eigen library
#include <Eigen/Dense>

// ESTstack library
#include "eststack/solution/base_kf.hpp"
#include "eststack/model/base_model.hpp"

/***
 * @brief An algorithm set focus on estimation and filtering
 * @author jinhua "siyiovo" deng
 */
namespace eststack
{
    /***
     * @brief algorithms for problems
     */
    namespace solution
    {
        /***
         * @brief scaled unscented Kalman filter on (parallelizable) manifold
         */
        template <eststack::model::LieGroupState StateT>
        class SUKFOM : public BaseKF<SUKFOM<StateT>, StateT>
        {
        public:
            using Base = BaseKF<SUKFOM<StateT>, StateT>;
            using State = typename Base::State;
            using Scalar = typename Base::Scalar;
            using Covariance = typename Base::Covariance;
        };
    }
}

#endif //! SOLUTION__SUKFOM_HPP
