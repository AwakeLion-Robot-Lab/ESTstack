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
            inline float weight(float residual, float scale)
            {
                float scale_sq = scale * scale;
                float res_sq = residual * residual;
                return scale_sq / (res_sq + scale_sq);
            }

            /***
             * @brief check whether is inlier via residual
             */
            inline bool isInlier(float residual, float scale)
            {
                return std::abs(residual) < 3.0f * scale;
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
            inline float weight(float residual, float scale)
            {
                float abs_res = std::abs(residual);
                if (abs_res <= scale)
                    return 1.0f;
                else
                    return scale / abs_res;
            }

            /***
             * @brief check whether is inlier via residual
             */
            inline bool isInlier(float residual, float scale)
            {
                return std::abs(residual) < 3.0f * scale;
            }
        };

        /***
         * @brief Iteratively Reweighted Least Squares for rigid transformation
         * @tparam LossFunc Loss function type
         */
        template <LossFunction LossFunc>
        class RigidIRLS
        {
        public:
            /***
             * @brief default constructor
             * @param decay_factor scale decay factor
             * @param min_scale minimum scale threshold
             * @param cost_tolerance cost change tolerance
             * @param max_iter maximum iterations
             */
            explicit RigidIRLS(float decay_factor = 1.3f, float min_scale = 1.0f,
                               float cost_tolerance = 0.01f, int max_iter = 100)
                : decay_factor_(decay_factor), min_scale_(min_scale),
                  cost_tolerance_(cost_tolerance), max_iter_(max_iter),
                  converged_(false)
            {
            }

            /***
             * @brief set input source points
             * @param source 3xN source point cloud
             */
            void setInputSource(const Eigen::Matrix3Xf &source)
            {
                source_ = source;
            }

            /***
             * @brief set input target points
             * @param target 3xN target point cloud
             */
            void setInputTarget(const Eigen::Matrix3Xf &target)
            {
                target_ = target;
            }

            void setMaxIterations(int max_iter) noexcept
            {
                max_iter_ = max_iter;
            }

            /***
             * @brief get final transformation matrix
             * @return const reference to final transformation
             */
            const Eigen::Isometry3f &getFinalTransformation() const noexcept
            {
                return final_transformation_;
            }

            /***
             * @brief check if optimization converged
             * @return true if converged
             */
            bool hasConverged() const noexcept
            {
                return converged_;
            }

            /***
             * @brief run IRLS optimization
             */
            void optimize()
            {
                /* check source-target dimensions */
                const int N = static_cast<int>(source_.cols());
                if (N == 0 || N != target_.cols())
                {
                    converged_ = false;
                    return;
                }

                /* initialization */
                Eigen::Isometry3f trans = Eigen::Isometry3f::Identity();
                Eigen::Matrix3Xf source_current = source_;
                Eigen::Matrix3Xf target_current = target_;
                Eigen::VectorXf weights = Eigen::VectorXf::Ones(N);
                Eigen::VectorXi inlier_indices = Eigen::VectorXi::LinSpaced(N, 0, N - 1);

                float scale = 0.0f;
                float prev_cost = std::numeric_limits<float>::max();

                LossFunc loss_func;

                for (int iter = 1; iter <= max_iter_; iter++)
                {
                    /* extract inlier points */
                    Eigen::Matrix3Xf source_temp = source_current(Eigen::all, inlier_indices);
                    Eigen::Matrix3Xf target_temp = target_current(Eigen::all, inlier_indices);
                    source_current = source_temp;
                    target_current = target_temp;

                    /* compute transform using weighted Kabsch */
                    trans = weighted_kabsch(source_current, target_current, weights);
                    Eigen::Matrix3f R = trans.linear();
                    Eigen::Vector3f t = trans.translation();

                    /* compute residuals */
                    Eigen::Matrix3Xf fit = (R * source_current).colwise() + t;
                    Eigen::VectorXf residuals = (fit - target_current).colwise().norm().transpose();

                    /* initialize scale on first iteration as max residual */
                    if (iter == 1)
                    {
                        scale = residuals.cwiseAbs().maxCoeff();
                        /* if convergence in loop 1, just return */
                        if (scale < min_scale_)
                        {
                            final_transformation_ = trans;
                            converged_ = true;
                            return;
                        }
                    }

                    /* compute current cost */
                    float cost = weights.cwiseProduct(residuals.array().square().matrix()).sum();

                    /* update next inlier indices via loss function */
                    Eigen::VectorXi flags(residuals.size());
                    for (int i = 0; i < residuals.size(); i++)
                    {
                        flags(i) = loss_func.isInlier(residuals(i), scale) ? 1 : 0;
                    }
                    inlier_indices = getNonZeroColumnIndicesFromVector(flags);

                    /* update next weights via loss function */
                    Eigen::VectorXf inlier_residuals = residuals(inlier_indices);
                    weights.resize(inlier_residuals.size());
                    for (int j = 0; j < inlier_residuals.size(); j++)
                    {
                        weights(j) = loss_func.weight(inlier_residuals(j), scale);
                    }

                    /* check convergence */
                    float cost_diff = std::abs(cost - prev_cost);
                    scale = scale / decay_factor_;
                    prev_cost = cost;

                    if (cost_diff < cost_tolerance_ || scale < min_scale_)
                    {
                        final_transformation_ = trans;
                        converged_ = true;
                        return;
                    }
                }

                final_transformation_ = trans;
                converged_ = false;
            }

        private:
            /***
             * @brief 3xN source points
             */
            Eigen::Matrix3Xf source_;

            /***
             * @brief 3xN target points
             */
            Eigen::Matrix3Xf target_;

            /***
             * @brief final transformation result
             */
            Eigen::Isometry3f final_transformation_;

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

            /***
             * @brief convergence flag
             */
            bool converged_;

            /***
             * @brief get non-zero column indices from vector
             * @param flags judgement
             * @return vector with indices of non-zero elements
             */
            inline Eigen::VectorXi getNonZeroColumnIndicesFromVector(const Eigen::VectorXi &flags)
            {
                int count = flags.count();
                Eigen::VectorXi nonzero_column(count);
                int idx = 0;
                for (int i = 0; i < flags.size(); ++i)
                {
                    if (flags(i) > 0)
                        nonzero_column(idx++) = i;
                }
                return nonzero_column;
            }
        };

    } // namespace core
} // namespace eststack

#endif // CORE__IRLS_HPP
