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

#ifndef CORE__PROSAC_HPP
#define CORE__PROSAC_HPP

// ESTstack library
#include "eststack/core/base_sac.hpp"

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
         * @brief PROgressive SAmple Consensus for robust model fitting within good correspondences
         * @tparam Model SAmple Consensus model to fit
         */
        template <SACModel Model>
        class PROSAC final : public BaseSAC<PROSAC<Model>, Model>
        {
        };
    } // namespace core
} // namespace eststack

#endif //! CORE__PROSAC_HPP