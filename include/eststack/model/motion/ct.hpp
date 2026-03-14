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

#ifndef MODEL__CT_HPP
#define MODEL__CT_HPP

// Eigen library
#include <Eigen/Dense>
#include <Eigen/Sparse>

// manif library
#include <manif/Bundle.h>
#include <manif/Rn.h>
#include <manif/SO2.h>

// ESTstack library
#include "eststack/model/base_model.hpp"
#include "eststack/types.hpp"

/***
 * @brief An algorithm set focus on estimation and filtering
 * @author jinhua "siyiovo" deng
 */
namespace eststack
{ /***
   * @brief models for problems
   */
    namespace model
    {
        /***
         * @brief coordinated turn motion model
         */
        class CT : public BaseTransitionModel<CT>
        {
        public:
            using State = manif::Bundle<double, manif::SO2, manif::R4>;
            using ControlInput = Eigen::Vector2d;
            using NoiseInput = Eigen::Vector4d;
            using StateJacobian = eststack::Jacobian<State, State>;
            using NoiseJacobian = eststack::Jacobian<State, NoiseInput>;

            /***
             * @brief compute the transition model
             * @param state current state
             * @param u control input
             * @return next state
             */
            State computeImpl(const State &state, const ControlInput &u) const;

            /***
             * @brief compute the state jacobian (F_x)
             * @param state current state
             * @param u control input
             * @return state jacobian matrix
             */
            StateJacobian computeStateJacobianImpl(const State &state, const ControlInput &u) const;

            /***
             * @brief compute the noise jacobian (F_w)
             * @param state current state
             * @param u control input
             * @return noise jacobian matrix
             */
            NoiseJacobian computeNoiseJacobianImpl(const State &state, const ControlInput &u) const;
        };
    }
}

#endif //! MODEL__CT_HPP