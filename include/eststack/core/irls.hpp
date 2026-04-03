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
#include <functional>

// Eigen library
#include <Eigen/Dense>

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
         * @brief Scale-Adaptive Cauchy loss function for robust estimation
         * @details weight = scale^2 / (residual^2 + scale^2)
         */
        struct SACauchyLoss
        {
            /***
             * @brief compute Cauchy weight
             * @param residual_sq squared residual
             * @param scale_sq squared scale parameter (alpha^2)
             * @return weight value
             */
            static double weight(double residual_sq, double scale_sq)
            {
                return scale_sq / (residual_sq + scale_sq);
            }

            /***
             * @brief check if residual is inlier
             * @param residual absolute residual
             * @param scale scale parameter (alpha)
             * @return true if inlier
             */
            static bool isInlier(double residual, double scale)
            {
                return std::abs(residual) < 3 * scale;
            }
        };

        /***
         * @brief Iteratively Reweighted Least Squares for robust estimation
         * @tparam LossFunc loss function type
         * @tparam StateT state type
         */
        template <typename LossFunc, DimAtCompileTime StateT>
        class IRLS
        {
        public:
            /***
             * @brief default constructor
             * @param decay_factor scale decay factor
             * @param min_scale minimum scale threshold
             * @param cost_tolerance cost convergence tolerance
             * @param max_iter maximum iterations
             */
            explicit IRLS(double decay_factor = 1.3, double min_scale = 1.0,
                          double cost_tolerance = 0.01, int max_iter = 100)
                : decay_factor_(decay_factor), min_scale_(min_scale),
                  cost_tolerance_(cost_tolerance), max_iter_(max_iter)
            {
            }

            /***
             * @brief run IRLS optimization
             * @tparam ResidualFuncType residual computation function type
             * @tparam UpdateFunc state update function type
             * @param init_state initial state estimation
             * @param measurements measurement data
             * @param compute_residuals function: (state, measurements) -> residuals vector
             * @param update_state function: (measurements, inlier_indices, weights) -> new_state
             * @return IRLSResult containing final state and convergence info
             */
            template <typename ResidualFuncType, typename UpdateFunc>
            State run(const StateT &init_state,
                     const StateT &measurements,
                     ResidualFuncType &&compute_residuals,
                     UpdateFunc &&update_state)
            {
                StateT current_state = init_state;
                double prev_cost = std::numeric_limits<double>::max();
                double scale = std::numeric_limits<double>::max();

                // Get initial residuals to determine number of points
                Eigen::VectorXd residuals = compute_residuals(current_state, measurements);
                const int n_points = residuals.size();

                // Initialize inlier indices (all points initially)
                Eigen::VectorXi inlier_indices = Eigen::VectorXi::LinSpaced(n_points, 0, n_points - 1);
                Eigen::VectorXd weights = Eigen::VectorXd::Ones(n_points);

                for (int iter = 1; iter <= max_iter_; ++iter)
                {
                    // Select inlier measurements for this iteration
                    const int n_inliers = inlier_indices.size();

                    // Step 1: Compute residuals with current state
                    residuals = compute_residuals(current_state, measurements);

                    // Step 2: Initialize scale on first iteration to max |residual|
                    if (iter == 1)
                    {
                        scale = residuals.cwiseAbs().maxCoeff();
                        if (scale < min_scale_)
                        {
                            // Perfect fit, no need to iterate
                            result.converged = true;
                            result.iterations = 0;
                            result.cost = 0.0;
                            result.scale = scale;
                            result.state = current_state;
                            result.inlier_count = n_points;
                            return result;
                        }
                    }

                    // Step 3: Compute weighted cost (energy)
                    // Only use inlier residuals for cost computation
                    double cost = 0.0;
                    for (int i = 0; i < n_inliers; ++i)
                    {
                        const int idx = inlier_indices(i);
                        cost += weights(i) * residuals(idx) * residuals(idx);
                    }

                    // Step 4: Update inlier set (truncation: residual < truncation_factor * scale)
                    int new_n_inliers = 0;
                    const double threshold = LossFunc::kTruncationFactor * scale;
                    for (int i = 0; i < n_points; ++i)
                    {
                        if (std::abs(residuals(i)) < threshold)
                        {
                            inlier_indices(new_n_inliers++) = i;
                        }
                    }

                    // Resize inlier_indices to actual count
                    inlier_indices.conservativeResize(new_n_inliers);

                    // Step 5: Update weights for inliers using loss function
                    const double scale_sq = scale * scale;
                    weights.resize(new_n_inliers);
                    for (int i = 0; i < new_n_inliers; ++i)
                    {
                        const int idx = inlier_indices(i);
                        const double residual_sq = residuals(idx) * residuals(idx);
                        weights(i) = LossFunc::weight(residual_sq, scale_sq);
                    }

                    // Step 6: Normalize weights
                    const double weight_sum = weights.sum();
                    if (weight_sum > std::numeric_limits<double>::epsilon())
                    {
                        weights /= weight_sum;
                    }

                    // Step 7: Update state using weighted inlier measurements
                    current_state = update_state(measurements, inlier_indices.head(new_n_inliers),
                                                 weights.head(new_n_inliers));

                    // Step 8: Check convergence
                    const double cost_diff = std::abs(cost - prev_cost);
                    if (cost_diff < cost_tolerance_ || scale < min_scale_)
                    {
                        result.converged = true;
                        result.iterations = iter;
                        result.cost = cost;
                        result.scale = scale;
                        result.state = current_state;
                        result.inlier_count = new_n_inliers;
                        return result;
                    }

                    // Step 9: Decay scale for next iteration (alpha = alpha / mu)
                    scale = scale / decay_factor_;
                    prev_cost = cost;
                }

                // Max iterations reached without convergence
                result.converged = false;
                result.iterations = max_iter_;
                result.cost = prev_cost;
                result.scale = scale;
                result.state = current_state;
                result.inlier_count = inlier_indices.size();
                return result;
            }

        private:
            double decay_factor_;   // mu: scale decay factor
            double min_scale_;      // gamma_min: minimum scale threshold
            double cost_tolerance_; // e_min: cost convergence tolerance
            int max_iter_;          // N_j: maximum iterations
        };

    } // namespace core
} // namespace eststack

#endif // CORE__IRLS_HPP
