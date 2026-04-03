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

#ifndef PROBLEM__POINTCLOUD_REGISTRATION_HPP
#define PROBLEM__POINTCLOUD_REGISTRATION_HPP

// ESTstack Library
#include <eststack/problem/base_problem.hpp>

/***
 * @brief An algorithm set focus on estimation and filtering
 * @author jinhua "siyiovo" deng
 */
namespace eststack
{
    /***
     * @brief problem definitions that orchestrate filters, models and noise
     */
    namespace problem
    {
        template <typename StateT, typename SolutionT>
        class PointcloudRegistration : BaseProblem<PointcloudRegistration<StateT, SolutionT>, StateT, SolutionT>
        {
        public:
            using Base = BaseProblem<PointcloudRegistration<StateT, SolutionT>, StateT, SolutionT>;
        };
    } // namespace problem
} // namespace eststack

#endif //! PROBLEM__POINTCLOUD_REGISTRATION_HPP