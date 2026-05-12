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
#include <algorithm>
#include <limits>
#include <random>
#include <ranges>
#include <vector>
#include <cmath>

// OpenMP library
#ifdef _OPENMP
#include <omp.h>
#endif

// ESTstack library
#include "eststack/concepts.hpp"
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
         * @brief RANdom SAmple Consensus for model fitting
         * @tparam Model SAmple Consensus model
         * @details refer to [PCL RANSAC](https://pointclouds.org/documentation/ransac_8hpp_source.html)
         */
        template <SACModel Model>
        class RANSAC : public BaseSAC<RANSAC<Model>, Model>
        {
        public:
            using Base = BaseSAC<RANSAC<Model>, Model>;
            using Base::best_fitness_;
            using Base::best_inliers_;
            using Base::converged_;
            using Base::max_iter_;
            using Base::model_;
            using Base::probability_;
            using Base::thread_num_;

            using ModelPtr = typename Base::ModelPtr;

            using Ptr = std::shared_ptr<RANSAC<Model>>;
            using ConstPtr = std::shared_ptr<const RANSAC<Model>>;

            /***
             * @brief default constructor
             * @param probability probability of choosing at least one sample free from outliers
             * @param max_iter maximum iterations
             */
            explicit RANSAC(float probability = 0.99f, int max_iter = 10000)
                : Base(probability, max_iter)
            {
                std::random_device rd;
                rng_.seed(rd());
            }

            /***
             * @brief run RANSAC optimization
             */
            void optimizeImpl()
            {
                /* check model */
                if (this->model_ == nullptr)
                {
                    this->converged_ = false;
                    return;
                }

                /* check source-target dimensions */
                const int cloud_size = this->model_->getCloudSize();
                const int sample_size = this->model_->getSampleSize();
                if (cloud_size == 0 || cloud_size < sample_size)
                {
                    this->converged_ = false;
                    return;
                }

                /* initialization */
                int best_inliers_num = 0;
                int iter = 0;
                int adaptive_max_iter = this->max_iter_;
                unsigned skipped_count = 0;
                const unsigned max_skip = static_cast<unsigned>(this->max_iter_) * 10;

                shuffled_indices_.resize(cloud_size);
                std::iota(shuffled_indices_.begin(), shuffled_indices_.end(), 0);

                int threads = this->thread_num_;
#ifdef _OPENMP
                if (threads == 0)
                    threads = omp_get_max_threads();
#endif

                /* RANSAC loop */
#pragma omp parallel if (threads > 0) num_threads(threads) shared(best_inliers_num, adaptive_max_iter, iter, skipped_count)
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
                            /* avoid infinite loop caused by skipped iterations */
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
                                    this->best_inliers_ = inliers;
                                    this->best_fitness_ = model_coeffs;

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
                                        adaptive_max_iter = static_cast<int>(std::ceil(std::log(1.0f - this->probability_) / std::log(p_outliers)));
                                        /***
                                         * if inlier_fraction is close to 0, adaptive_max_iter become very large
                                         * so constrain the upper bound is max_iter_ to avoid infinite loop
                                         */
                                        adaptive_max_iter = std::min(this->max_iter_, adaptive_max_iter);
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

                        if (iter_local > adaptive_max_iter_local || iter_local > this->max_iter_)
                            break;
                    }
                }

                /* after get the most inliers set, fit model finally to get best parameter */
                if (best_inliers_num > 0 && this->model_->fit(this->best_inliers_, this->best_fitness_))
                {
                    this->converged_ = true;
                }
                else
                {
                    this->converged_ = false;
                }
            }

        private:
            /***
             * @brief random generator
             */
            std::mt19937 rng_;

            /***
             * @brief shuffled indices for random sampling
             */
            std::vector<int> shuffled_indices_;

            /***
             * @brief get random samples for model fitting
             * @param sample_size number of samples to fit
             * @param cloud_size size of the input point cloud
             */
            inline std::vector<int> getSamples(int sample_size)
            {
                std::vector<int> samples(sample_size);

                std::size_t index_size = shuffled_indices_.size();
                for (std::size_t i = 0; i < sample_size; ++i)
                {
                    std::uniform_int_distribution<std::size_t> dist(i, index_size - 1);
                    std::ranges::swap(shuffled_indices_[i], shuffled_indices_[dist(rng_)]);
                }
                std::ranges::copy(shuffled_indices_ | std::views::take(sample_size), samples.begin());

                return samples;
            }
        };

    } // namespace core
} // namespace eststack

#endif //! CORE__RANSAC_HPP
