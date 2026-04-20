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

#ifndef MODEL__BASE_IMU_MODEL_HPP
#define MODEL__BASE_IMU_MODEL_HPP

// Eigen library
#include <Eigen/Core>

// manif library
#include <manif/Bundle.h>
#include <manif/SE_2_3.h>
#include <manif/Rn.h>

// ESTstack library
#include "eststack/model/motion/base_motion_model.hpp"
#include "eststack/types.hpp"

/***
 * @brief An algorithm set focus on state estimation
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
         * @brief base class for IMU model
         * @tparam Derived derived IMU model class
         * @tparam StateT state type
         * @tparam ControlInputT control input type
         * @tparam ProcessNoiseT process noise vector type
         */
        template <typename Derived, typename StateT, typename ControlInputT, typename ProcessNoiseT>
        class BaseIMUModel : public BaseTransitionModel<Derived, StateT, ControlInputT, ProcessNoiseT>
        {
        };
    } // namespace model
} // namespace eststack

#endif //! MODEL__BASE_IMU_MODEL_HPP