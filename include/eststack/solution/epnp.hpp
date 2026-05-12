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

#ifndef SOLUTION__BASE_PCR_HPP
#define SOLUTION__BASE_PCR_HPP

// C++ standard library
#include <memory>

// Eigen library
#include <Eigen/Geometry>
#include <Eigen/Core>

// manif library
#include <manif/SE3.h>

/***
 * @brief An algorithm set focus on state estimation
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
         * @brief Perspective n Points solver for 3D-2D geometry
         */
        class EPnP
        {
        };
    } // namespace solution
} // namespace eststack

#endif //! SOLUTION__BASE_PCR_HPP
