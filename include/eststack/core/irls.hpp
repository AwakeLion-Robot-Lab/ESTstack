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

#ifndef CORE__IRLS_HPP
#define CORE__IRLS_HPP

// C++ standard library
#include <cmath>
#include <limits>
#include <vector>
#include <numeric>

// Eigen library
#include <Eigen/Dense>

// ESTstack library
#include "eststack/concepts.hpp"
#include "eststack/core/kabsch.hpp"

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
         * @brief Cauchy loss function
         */
        struct CauchyLoss
        {
            /***
             * @brief assign weight for each residual
             */
            static double weight(double residual, double scale)
            {
                double scale_sq = scale * scale;
                double res_sq = residual * residual;
                return scale_sq / (res_sq + scale_sq);
            }

            /***
             * @brief check whether is inlier via residual
             */
            static bool isInlier(double residual, double scale)
            {
                return std::abs(residual) < 3.0 * scale;
            }
        };

        /***
         * @brief Huber loss function
         */
        struct HuberLoss
        {
            /***
             * @brief assign weight for each residual
             */
            static double weight(double residual, double scale)
            {
                double abs_res = std::abs(residual);
                if (abs_res <= scale)
                    return 1.0;
                else
                    return scale / abs_res;
            }

            /***
             * @brief check whether is inlier via residual
             */
            static bool isInlier(double residual, double scale)
            {
                return std::abs(residual) < 3.0 * scale;
            }
        };

        /***
         * @brief IRLS result
         */
        struct IRLSResult
        {
            IRLSResult(const Eigen::Isometry3d &T = Eigen::Isometry3d::Identity(), double cost = 0.0, double scale = 0.0,
                       int inlier_count = 0, int iterations = 0, bool converged = false)
                : T_(T), cost_(cost), scale_(scale),
                  inlier_count_(inlier_count), iterations_(iterations), converged_(converged) {}

            /***
             * @brief transformation
             */
            Eigen::Isometry3d T_;

            /***
             * @brief cost
             * @details formula: sqrt(new_weight * new_residual) - sqrt(old_weight * old_residual)
             */
            double cost_;

            /***
             * @brief scale for specific loss function
             */
            double scale_;

            /***
             * @brief number of inliers
             */
            int inlier_count_;

            /***
             * @brief number of iterations
             */
            int iterations_;

            /***
             * @brief convergence flag
             */
            bool converged_;
        };

        /***
         * @brief Iteratively Reweighted Least Squares for robust rigid transformation
         * @tparam LossFunc Loss function type
         */
        template <LossFunction LossFunc>
        class IRLS
        {
        public:
            /***
             * @brief default constructor
             * @param decay_factor scale decay factor
             * @param min_scale minimum scale threshold
             * @param cost_tolerance cost change tolerance
             * @param max_iter maximum iterations
             */
            explicit IRLS(float decay_factor = 1.3, float min_scale = 1.0,
                          float cost_tolerance = 0.01, int max_iter = 100)
                : decay_factor_(decay_factor), min_scale_(min_scale),
                  cost_tolerance_(cost_tolerance), max_iter_(max_iter)
            {
            }

            /***
             * @brief run IRLS optimization
             * @param source 3xN source points
             * @param target 3xN target points
             * @return IRLSResult
             */
            IRLSResult optimize(const Eigen::Matrix3Xd &source, const Eigen::Matrix3Xd &target)
            {
                /* check source-target dimensions */
                const int N = static_cast<int>(source.cols());
                if (N == 0 || N != target.cols())
                    return IRLSResult();

                /* initialization */
                Eigen::VectorXi inlier_indices = Eigen::VectorXi::LinSpaced(N, 0, N - 1);
                Eigen::VectorXd weights = Eigen::VectorXd::Ones(N);
                Eigen::Isometry3d trans = Eigen::Isometry3d::Identity();
                Eigen::Matrix3Xd source_current = source;
                Eigen::Matrix3Xd target_current = target;

                double scale = 0.0;
                double prev_cost = std::numeric_limits<double>::max();

                for (int iter = 1; iter <= max_iter_; ++iter)
                {
                    const int inlier_num = inlier_indices.size();

                    /***
                     * select inlier points via block operations
                     * reference: https://github.com/ShiPC-AI/TCF/blob/main/registration/irls_welsch.cpp#L16
                     */
                    Eigen::Matrix3Xd source_inlier = source_current(Eigen::all, inlier_indices);
                    Eigen::Matrix3Xd target_inlier = target_current(Eigen::all, inlier_indices);

                    /* compute transformation via weighted kabsch */
                    trans = weighted_kabsch(source_inlier, target_inlier, weights.head(inlier_num));
                    const Eigen::Matrix3d R = trans.linear();
                    const Eigen::Vector3d t = trans.translation();

                    /* compute initial residuals */
                    const Eigen::Matrix3Xd fit = (R * source_inlier).colwise() + t;
                    const Eigen::VectorXd residuals = (fit - target_inlier).colwise().norm();

                    /* initialize scale on first iteration */
                    if (iter == 1)
                    {
                        scale = residuals.maxCoeff();
                        /* early convergence, just return */
                        if (scale < min_scale_)
                            return IRLSResult(trans, 0.0, scale, inlier_num, 1, true);
                    }

                    /* compute weighted cost */
                    const double cost = weights.head(inlier_num).cwiseProduct(residuals.cwiseProduct(residuals)).sum();

                    /* check convergence */
                    if (iter > 1)
                    {
                        if (std::abs(cost - prev_cost) < cost_tolerance_ || scale < min_scale_)
                            return IRLSResult(trans, cost, scale, inlier_num, iter, true);
                    }
                    prev_cost = cost;

                    /* update inliers via scale of loss function */
                    Eigen::VectorXi new_inliers(inlier_num);
                    int new_inlier_num = 0;
                    for (int i = 0; i < inlier_num; ++i)
                    {
                        if (LossFunc::isInlier(residuals(i), scale))
                            new_inliers(new_inlier_num++) = inlier_indices(i);
                    }

                    /* no inliers remain, poor fit */
                    if (new_inlier_num == 0)
                        return IRLSResult(trans, prev_cost, scale, 0, iter, false);

                    /* resize and update inlier indices */
                    inlier_indices.conservativeResize(new_inlier_num);
                    inlier_indices = new_inliers.head(new_inlier_num);

                    /* update weights via loss function */
                    weights.conservativeResize(new_inlier_num);
                    for (int i = 0; i < new_inlier_num; ++i)
                    {
                        weights(i) = LossFunc::weight(residuals(i), scale);
                    }

                    /* update scale for next iteration */
                    scale /= decay_factor_;
                }

                /***
                 * NOTE that if not converge until max_iter_,
                 * scale must multiply by decay factor to back in,
                 * because there is no more next iteration!
                 */
                return IRLSResult(trans, prev_cost, scale * decay_factor_,
                                  inlier_indices.size(), max_iter_, false);
            }

        private:
            /***
             * @brief decay factor for scale update
             */
            float decay_factor_;

            /***
             * @brief minimum scale threshold for convergence
             */
            float min_scale_;

            /***
             * @brief cost change tolerance for convergence
             */
            float cost_tolerance_;

            /***
             * @brief maximum iterations
             */
            int max_iter_;
        };

    } // namespace core
} // namespace eststack

#endif // CORE__IRLS_HPP
