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
            using PointCloud = typename Model::PointCloud;
            using PointCloudPtr = typename Model::PointCloudPtr;
            using PointCloudConstPtr = typename Model::PointCloudConstPtr;

            /***
             * @brief default constructor
             * @param probability probability of good fitting in iterations
             * @param max_iter maximum iterations
             */
            explicit RANSAC(float probability = 0.99f, int max_iter = 10000)
                : probability_(probability), max_iter_(max_iter), converged_(false)
            {
                std::random_device rd;
                rng_.seed(rd());
            }

            /***
             * @brief set source cloud
             * @param source source point cloud
             */
            void setInputSource(const PointCloudConstPtr &source)
            {
                source_cloud_ = source;
            }

            /***
             * @brief set target cloud
             * @param target target point cloud
             */
            void setInputTarget(const PointCloudConstPtr &target)
            {
                target_cloud_ = target;
            }

            /***
             * @brief set probability
             */
            void setProbability(float probability) noexcept
            {
                probability_ = probability;
            }

            /***
             * @brief set max iterations
             */
            void setMaxIterations(int max_iterations) noexcept
            {
                max_iter_ = max_iterations;
            }

            const Eigen::Isometry3f &getFinalTransformation() const noexcept
            {
                return final_transformation_;
            }

            bool hasConverged() const noexcept
            {
                return converged_;
            }

            /***
             * @brief run RANSAC optimization
             */
            void optimize()
            {
                const int cloud_size = source_cloud_->size();
                if (cloud_size == 0 || cloud_size != target_cloud_->size())
                {
                    converged_ = false;
                    return;
                }

                const float eps = std::numeric_limits<float>::epsilon();

                std::vector<int> indices(cloud_size);
                std::iota(indices.begin(), indices.end(), 0);

                Eigen::VectorXf best_model;
                std::vector<int> best_inliers;

                int iteration = 0;
                int adaptive_max_iter = max_iter_;
                int best_inlier_num = 0;

                Model model;
                int sample_size = model.getSampleSize();

                while (iteration < adaptive_max_iter && iteration < max_iter_)
                {
                    /* random sample */
                    std::vector<int> sample = random_sample(indices, sample_size);

                    /* fit model */
                    if (!model.fit(sample))
                    {
                        ++iteration;
                        continue;
                    }

                    /* count inliers */
                    std::vector<int> inliers = model.countInliers();
                    int inlier_count = inliers.size();

                    /* update best */
                    if (inlier_count > best_inlier_num)
                    {
                        best_inlier_num = inlier_count;
                        best_model = model.getBestFitness();
                        best_inliers = inliers;

                        /* adaptive iteration update */
                        float frac_inliers = static_cast<float>(inlier_count) / cloud_size;
                        float p_no_outliers = 1.0f - std::pow(frac_inliers, sample_size);
                        p_no_outliers = std::max(eps, p_no_outliers);
                        p_no_outliers = std::min(1.0f - eps, p_no_outliers);

                        if (p_no_outliers < 1.0f)
                        {
                            adaptive_max_iter = static_cast<int>(std::log(1.0f - probability_) / std::log(p_no_outliers));
                            adaptive_max_iter = std::max(adaptive_max_iter, iteration + 1);
                        }
                    }

                    iteration++;
                }
            }

        private:
            /***
             * @brief source cloud
             */
            PointCloudConstPtr source_cloud_;

            /***
             * @brief target cloud
             */
            PointCloudConstPtr target_cloud_;

            /***
             * @brief random number generator
             */
            std::mt19937 rng_;

            /***
             * @brief probability of good fitting in iterations
             */
            float probability_;

            /***
             * @brief maximum iterations
             */
            int max_iter_;

            /***
             * @brief convergence flag
             */
            bool converged_;

            /***
             * @brief random sample selection
             */
            std::vector<int> random_sample(const std::vector<int> &indices, int s)
            {
                std::vector<int> sample;
                sample.reserve(s);

                std::uniform_int_distribution<int> dist(0, indices.size() - 1);

                for (int i = 0; i < s; ++i)
                {
                    int idx = dist(rng_);
                    sample.push_back(indices[idx]);
                }

                return sample;
            }
        };

    } // namespace core
} // namespace eststack

#endif //! CORE__RANSAC_HPP
