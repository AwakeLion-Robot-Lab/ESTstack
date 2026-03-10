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

#ifndef SOLUTION__ESKFOM_HPP
#define SOLUTION__ESKFOM_HPP

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
         * @brief error state Kalman filter on manifold
         * @tparam StateT state type
         */
        template <eststack::model::LieGroupState StateT>
        class ESKFOM : public BaseKF<ESKFOM<StateT>, StateT>
        {
        public:
            EIGEN_MAKE_ALIGNED_OPERATOR_NEW;

            using Base = BaseKF<ESKFOM<StateT>, StateT>;
            using State = typename Base::State;
            using Scalar = typename Base::Scalar;
            using Covariance = typename Base::Covariance;

            /***
             * @brief default constructor
             */
            ESKFOM() = default;

            /***
             * @brief EKSFOM prediction step implementation
             * @details see `base_kf.hpp` for more details on the arguments
             */
            template <eststack::model::TransitionModel TransitionModel, typename... Args>
            const State &predictImpl(const TransitionModel &F, const typename TransitionModel::ControlInput &u, Args &&...args)
            {
            }

            /***
             * @brief EKSFOM update step implementation
             * @details see `base_kf.hpp` for more details on the arguments
             */
            template <eststack::model::MeasurementModel MeasurementModel, typename... Args>
            const State &updateImpl(const MeasurementModel &H, const typename MeasurementModel::Measurement &z, Args &&...args)
            {
            }
        };
    }
}

#endif //! SOLUTION__ESKFOM_HPP
