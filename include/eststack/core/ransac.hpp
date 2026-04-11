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

#ifndef CORE__RANSAC_HPP
#define CORE__RANSAC_HPP

// C++ standard library
#include <cmath>
#include <limits>
#include <random>
#include <algorithm>
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
         * @brief RANdom SAmple Consensus for model fitting
         * @tparam Model SAmple Consensus model
         */
        template <SACModel Model>
        class RANSAC
        {
        public:
            using ModelPtr = typename Model::Ptr;
            using ModelConstPtr = typename Model::ConstPtr;

            /***
             * @brief default constructor
             * @param probability probability of good fitting in iterations
             * @param max_iter maximum iterations
             */
            explicit RANSAC(double probability = 0.99, int max_iter = 10000)
                : probability_(probability), max_iter_(max_iter), converged_(false)
            {
                std::random_device rd;
                rng_.seed(rd());
            }

            /***
             * @brief set probability
             * @param probability of choosing at least one sample free from outliers
             */
            void setProbability(double probability)
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
                /* check model */
                if (model_ == nullptr)
                {
                    converged_ = false;
                    return;
                }

                /* check source-target dimensions */
                const int cloud_size = model_->getCloudSize();
                const int sample_size = model_->getSampleSize();
                if (cloud_size == 0 || cloud_size < sample_size)
                {
                    converged_ = false;
                    return;
                }

                /* initialization */
                Eigen::VectorXf model_coeffs(model_->getModelSize());
                int best_inliers_num = 0;
                int adaptive_max_iter = max_iter_;

                /* RANSAC loop */
                for (int iter = 0; iter < adaptive_max_iter; iter++)
                {
                    /* get random samples in specific number */
                    std::vector<int> samples = getSamples(sample_size, cloud_size);

                    /* check fitting status */
                    if (!model_->fit(samples, model_coeffs))
                        continue;

                    /* select inliers for this iteration */
                    std::vector<int> inliers = model_->selectInliers(model_coeffs);
                    int inliers_num = static_cast<int>(inliers.size());
                    /* update best inliers if better than last time */
                    if (inliers_num > best_inliers_num)
                    {
                        best_inliers_num = inliers_num;
                        best_inliers_ = std::move(inliers);
                        best_fitness_ = model_coeffs;

                        /* update adaptive maximum iterations */
                        double inlier_fraction = static_cast<double>(inliers_num) / cloud_size;
                        double p_outliers = 1.0 - std::pow(inlier_fraction, sample_size);
                        /* avoid -inf */
                        p_outliers = std::max(std::numeric_limits<double>::epsilon(), p_outliers);
                        /* avoid 0 */
                        p_outliers = std::min(1.0 - std::numeric_limits<double>::epsilon(), p_outliers);

                        /* if probability of at least one outlier < 1.0, update */
                        if (p_outliers < 1.0)
                        {
                            adaptive_max_iter = static_cast<int>(std::ceil(std::log(1.0 - probability_) / std::log(p_outliers)));
                            /***
                             * if inlier_fraction is close to 0, adaptive_max_iter become very large
                             * so constrain the upper bound is max_iter_ to avoid infinite loop
                             */
                            adaptive_max_iter = std::min(max_iter_, adaptive_max_iter);
                            /***
                             * samely, if inlier_fraction is close to 1, adaptive_max_iter become very small
                             * so constrain the lower bound is 1 to avoid zero iterations
                             * this can provide at least one more iteration
                             */
                            adaptive_max_iter = std::max(adaptive_max_iter, iter + 1);
                        }
                    }
                }

                /* after get the most inliers set, fit model finally to get best parameter */
                if (best_inliers_num > 0)
                {
                    model_->fit(best_inliers_, best_fitness_);
                    converged_ = true;
                }
                else
                {
                    converged_ = false;
                }
            }

        private:
            /***
             * @brief SAmple Consensus model to fit
             */
            ModelPtr model_;

            /***
             * @brief random generator
             */
            std::mt19937 rng_;

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
            double probability_;

            /***
             * @brief maximum iterations
             */
            int max_iter_;

            /***
             * @brief flag for convergence
             */
            bool converged_;

            /***
             * @brief get random samples for model fitting
             * @param sample_size number of samples to fit
             * @param cloud_size size of the input point cloud
             */
            std::vector<int> getSamples(int sample_size, int cloud_size)
            {
                std::vector<int> indices(cloud_size);
                std::iota(indices.begin(), indices.end(), 0);
                std::vector<int> samples;
                samples.reserve(sample_size);
                std::sample(indices.begin(), indices.end(), std::back_inserter(samples), sample_size, rng_);
                return samples;
            }
        };

    } // namespace core
} // namespace eststack

#endif //! CORE__RANSAC_HPP
