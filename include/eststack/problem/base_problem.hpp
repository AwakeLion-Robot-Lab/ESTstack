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

#ifndef PROBLEM__BASE_PROBLEM_HPP
#define PROBLEM__BASE_PROBLEM_HPP

// C++ standard library
#include <concepts>
#include <type_traits>
#include <utility>
#include <memory>
#include <exception>

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
        /***
         * @brief filter-based problem concept
         * @details constrains any type that acts as a filtering problem:
         *   owns a filter, exposes State / Covariance, and can be initialized
         */
        template <typename T>
        concept EstProblem = requires(T prob) {
            { prob.isInitialized() } -> std::convertible_to<bool>;
            { prob.setSolution(std::declval<typename T::Solution>()) } -> std::same_as<void>;
            { prob.isOK() } -> std::same_as<bool>;
            { prob.run() } -> std::same_as<bool>;
        };

        /***
         * @brief base class for estimation problems
         * @tparam Derived problem class
         * @tparam StateT state type
         * @tparam SolutionT solution type
         */
        template <typename Derived, typename StateT, typename SolutionT>
        class BaseProblem
        {
        public:
            using State = StateT;
            using Solution = SolutionT;

            /***
             * @brief check if the problem has been initialized
             */
            bool isInitialized() const noexcept
            {
                return initialized_;
            }

            /***
             * @brief get state
             */
            const State &getState() const noexcept
            {
                return state_;
            }

            /***
             * @brief set state
             * @param state state
             */
            void setState(const State &state) noexcept
            {
                state_ = state;
            }

            /***
             * @brief set specific solution of the problem
             * @param sol solution
             */
            void setSolution(SolutionT &&sol)
            {
                solution_ = std::make_unique<SolutionT>(std::forward<SolutionT>(sol));
                static_cast<Derived *>(this)->setSolutionImpl(solution_.get());
                initialized_ = true;
            }

            /***
             * @brief run to solve the problem
             */
            bool run()
            {
                if (!isInitialized())
                    throw std::runtime_error("problem is not initialized!");

                return static_cast<Derived *>(this)->runImpl();
            }

        protected:
            /***
             * @brief default constructor
             */
            BaseProblem() = default;

            /***
             * @brief state vector
             */
            State state_;

            /***
             * @brief the solution of the problem
             */
            std::unique_ptr<Solution> solution_;

            /***
             * @brief initialization flag
             */
            bool initialized_{false};

            /***
             * @brief OK flag
             */
            bool ok_{false};
        };

    } // namespace problem
} // namespace eststack

#endif //! PROBLEM__BASE_PROBLEM_HPP