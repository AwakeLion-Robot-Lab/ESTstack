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

#ifndef CORE__BASE_SAC_HPP
#define CORE__BASE_SAC_HPP

// C++ standard library
#include <vector>

// ESTstack library
#include "eststack/concepts.hpp"

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
         * @brief base SAmple Consensus class
         * @tparam Derived derived SAmple Consensus class
         * @tparam Model SAmple Consensus model
         * @details refer to [PCL SAC](https://pointclouds.org/documentation/sample__consensus_2include_2pcl_2sample__consensus_2sac_8h_source.html#l00304)
         */
        template <typename Derived, SACModel Model>
        class BaseSAC
        {
        public:
            using ModelPtr = typename Model::Ptr;
            using ModelConstPtr = typename Model::ConstPtr;

            /***
             * @brief set probability
             * @param probability of choosing at least one sample free from outliers
             */
            void setProbability(float probability)
            {
                probability_ = probability;
            }

            /***
             * @brief set SAmple Consensus model
             * @param model SAmple Consensus model to fit
             */
            void setSACModel(const ModelPtr &model)
            {
                model_ = model;
            }

            /***
             * @brief set maximum iterations
             * @param max_iter maximum iterations threshold
             */
            void setMaxIterations(int max_iter) noexcept
            {
                max_iter_ = max_iter;
            }

            /***
             * @brief set thread number for parallel optimization
             */
            void setThreadNum(int thread_num)
            {
                thread_num_ = thread_num;
            }

            /***
             * @brief check if RANSAC has converged
             */
            bool hasConverged() const noexcept
            {
                return converged_;
            }

            /***
             * @brief get best fitness of the model fitting
             */
            const Eigen::VectorXf &getBestFitness() const noexcept
            {
                return best_fitness_;
            }

            /***
             * @brief get best inliers of the model fitting
             */
            const std::vector<int> &getBestInliers() const noexcept
            {
                return best_inliers_;
            }

            /***
             * @brief run RANSAC optimization
             */
            void optimize()
            {
                static_cast<Derived *>(this)->optimizeImpl();
            }

        protected:
            /***
             * @brief default constructor
             * @param probability probability of choosing at least one sample free from outliers
             * @param max_iter maximum iterations
             */
            explicit BaseSAC(float probability = 0.99f, int max_iter = 10000)
                : probability_(probability), max_iter_(max_iter),
                  thread_num_(0), converged_(false)
            {
            }

            /***
             * @brief SAmple Consensus model to fit
             */
            ModelPtr model_;

            /***
             * @brief best fitness of the model
             */
            Eigen::VectorXf best_fitness_;

            /***
             * @brief best inliers of the model fitting
             */
            std::vector<int> best_inliers_;

            /***
             * @brief probability of choosing at least one sample free from outliers
             */
            float probability_;

            /***
             * @brief maximum iterations
             */
            int max_iter_;

            /***
             * @brief thread number for parallel optimization
             */
            int thread_num_;

            /***
             * @brief flag for convergence
             */
            bool converged_;
        };

    } // namespace core
} // namespace eststack

#endif //! CORE__BASE_SAC_HPP