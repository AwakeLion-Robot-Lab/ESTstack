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
#include <numeric>
#include <random>
#include <algorithm>
#include <vector>

// OpenMP library
#ifdef _OPENMP
#include <omp.h>
#endif

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
         * @details refer to [PCL RANSAC](https://pointclouds.org/documentation/ransac_8hpp_source.html)
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
                int best_inliers_num = 0;
                int iter = 0;
                int adaptive_max_iter = max_iter_;
                unsigned skipped_count = 0;
                const unsigned max_skip = static_cast<unsigned>(max_iter_) * 10;

                shuffled_indices_.resize(cloud_size);
                std::iota(shuffled_indices_.begin(), shuffled_indices_.end(), 0);

                /* RANSAC loop */
#pragma omp parallel shared(best_inliers_num, adaptive_max_iter, iter, skipped_count)
                {
                    while (true)
                    {
                        /* get random samples in specific number */
                        std::vector<int> samples;
#pragma omp critical(samples)
                        {
                            samples = getSamples(sample_size);
                        }

                        /* check fitting status */
                        Eigen::VectorXf model_coeffs(model_->getModelSize());
                        if (!model_->fit(samples, model_coeffs))
                        {
                            unsigned skipped_temp;
#pragma omp atomic capture
                            skipped_temp = ++skipped_count;
                            if (skipped_temp < max_skip)
                                continue;
                            else
                                break;
                        }

                        /* select inliers for this iteration */
                        std::vector<int> inliers = model_->selectInliers(model_coeffs);
                        int inliers_num = static_cast<int>(inliers.size());

                        int best_inliers_num_local;
#pragma omp atomic read
                        best_inliers_num_local = best_inliers_num;

                        /* update best inliers if better than last time */
                        if (inliers_num > best_inliers_num_local)
                        {
#pragma omp critical(update)
                            {
                                if (inliers_num > best_inliers_num)
                                {
                                    best_inliers_num = inliers_num;
                                    best_inliers_ = inliers;
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
                                    }
                                }
                            }
                        }

                        /* check iteration at last for convergence */
                        int iter_local, adaptive_max_iter_local;

#pragma omp atomic capture
                        iter_local = ++iter;
#pragma omp atomic read
                        adaptive_max_iter_local = adaptive_max_iter;

                        if (iter_local > adaptive_max_iter_local || iter_local > max_iter_)
                            break;
                    }
                }

                /* after get the most inliers set, fit model finally to get best parameter */
                if (best_inliers_num > 0 && model_->fit(best_inliers_, best_fitness_))
                {
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
             * @brief best fitness of the model
             */
            Eigen::VectorXf best_fitness_;

            /***
             * @brief random generator
             */
            std::mt19937 rng_;

            /***
             * @brief shuffled indices for random sampling
             */
            std::vector<int> shuffled_indices_;

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
            std::vector<int> getSamples(int sample_size)
            {
                std::vector<int> samples(sample_size);

                std::size_t index_size = shuffled_indices_.size();
                for (std::size_t i = 0; i < sample_size; ++i)
                {
                    std::uniform_int_distribution<std::size_t> dist(i, index_size - 1);
                    std::swap(shuffled_indices_[i], shuffled_indices_[dist(rng_)]);
                }
                std::copy(shuffled_indices_.cbegin(), shuffled_indices_.cbegin() + sample_size, samples.begin());

                return samples;
            }
        };

    } // namespace core
} // namespace eststack

#endif //! CORE__RANSAC_HPP
