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
#include <memory>
#include <limits>
#include <vector>
#include <cmath>

// Eigen library
#include <Eigen/Geometry>
#include <Eigen/Core>

// ESTstack library
#include "eststack/eigen_traits.hpp"
#include "eststack/core/kabsch.hpp"
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
         * @brief Cauchy loss function
         */
        struct CauchyLoss
        {
            using Ptr = std::shared_ptr<CauchyLoss>;
            using ConstPtr = std::shared_ptr<const CauchyLoss>;

            /***
             * @brief assign weight for each residual
             */
            inline float weight(float residual, float scale) const
            {
                float scale_sq = scale * scale;
                float res_sq = residual * residual;
                return scale_sq / (res_sq + scale_sq);
            }

            /***
             * @brief check whether is inlier via residual
             */
            inline bool isInlier(float residual, float scale) const
            {
                return std::abs(residual) < 3.0f * scale;
            }
        };

        /***
         * @brief Huber loss function
         */
        struct HuberLoss
        {
            using Ptr = std::shared_ptr<HuberLoss>;
            using ConstPtr = std::shared_ptr<const HuberLoss>;

            /***
             * @brief assign weight for each residual
             */
            inline float weight(float residual, float scale) const
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
            inline bool isInlier(float residual, float scale) const
            {
                return std::abs(residual) < 3.0f * scale;
            }
        };

        /***
         * @brief Geman-McClure loss function
         */
        struct GemanMcClureLoss
        {
            using Ptr = std::shared_ptr<GemanMcClureLoss>;
            using ConstPtr = std::shared_ptr<const GemanMcClureLoss>;

            /***
             * @brief assign weight for each residual
             */
            inline float weight(float residual, float scale) const
            {
                float scale_sq = scale * scale;
                float res_sq = residual * residual;
                float denom = res_sq + scale_sq;
                return (scale_sq * scale_sq) / (denom * denom);
            }

            /***
             * @brief check whether is inlier via residual
             */
            inline bool isInlier(float residual, float scale) const
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
            EIGEN_MAKE_ALIGNED_OPERATOR_NEW

            using LossPtr = typename LossFunc::Ptr;
            using LossConstPtr = typename LossFunc::ConstPtr;

            using Ptr = std::shared_ptr<RigidIRLS<LossFunc>>;
            using ConstPtr = std::shared_ptr<const RigidIRLS<LossFunc>>;

            /***
             * @brief default constructor
             * @param decay_factor scale decay factor
             * @param min_scale minimum scale threshold
             * @param cost_tolerance cost change tolerance
             * @param max_iter maximum iterations
             */
            explicit RigidIRLS(float decay_factor = 1.3f, float min_scale = 1.0f,
                               float cost_tolerance = 0.01f, int max_iter = 100)
                : init_guess_(Eigen::Isometry3f::Identity()), decay_factor_(decay_factor),
                  min_scale_(min_scale), cost_tolerance_(cost_tolerance),
                  max_iter_(max_iter), converged_(false)
            {
            }

            /***
             * @brief set loss function
             * @param loss_func loss function
             */
            void setLossFunction(const LossConstPtr &loss_func) noexcept
            {
                loss_func_ = loss_func;
            }

            /***
             * @brief set input source point cloud
             * @param source 3xN source point cloud
             */
            void setInputSource(const Eigen::Matrix3Xf &source)
            {
                source_ = source;
            }

            /***
             * @brief set input target point cloud
             * @param target 3xN target point cloud
             */
            void setInputTarget(const Eigen::Matrix3Xf &target)
            {
                target_ = target;
            }

            /***
             * @brief set initial guess for transformation
             * @param init_guess initial guess
             */
            void setInitGuess(const Eigen::Isometry3f &init_guess)
            {
                init_guess_ = init_guess;
            }

            /***
             * @brief set maximum iterations
             * @param max_iter maximum iterations
             */
            void setMaxIterations(int max_iter) noexcept
            {
                max_iter_ = max_iter;
            }

            /***
             * @brief get final transformation matrix
             * @return final transformation
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
                /* check loss function */
                if (loss_func_ == nullptr)
                {
                    final_transformation_ = Eigen::Isometry3f::Identity();
                    converged_ = false;
                    return;
                }

                /* check source-target dimensions */
                const int N = static_cast<int>(source_.cols());
                if (N == 0 || N != target_.cols())
                {
                    final_transformation_ = Eigen::Isometry3f::Identity();
                    converged_ = false;
                    return;
                }

                /* initialization */
                Eigen::Isometry3f trans = init_guess_;
                Eigen::Matrix3Xf source_current = source_;
                Eigen::Matrix3Xf target_current = target_;
                Eigen::VectorXf weights = Eigen::VectorXf::Ones(N);
                Eigen::VectorXi inlier_indices = Eigen::VectorXi::LinSpaced(N, 0, N - 1);

                float scale = 0.0f;
                float prev_cost = std::numeric_limits<float>::max();

                for (int iter = 1; iter <= max_iter_; iter++)
                {
                    /***
                     * extract inlier points
                     * NOTE that source_current(Eigen::placeholders::all, inlier_indices) is slicing operation in Eigen,
                     * it means select all rows but specfic column that are in inlier_indices
                     * e.g. if source_current is 3x100 and inlier_indices is [0, 2, 5],
                     * then it turns out to be 3x3 matrix consisting of columns 0, 2, 5 of `source_current`
                     * slicing operation can refer to https://zhuanlan.zhihu.com/p/531476504
                     * also see `model/sac/tcf_sac_model.hpp`
                     */
                    Eigen::Matrix3Xf source_temp = source_current(Eigen::placeholders::all, inlier_indices);
                    Eigen::Matrix3Xf target_temp = target_current(Eigen::placeholders::all, inlier_indices);
                    source_current = source_temp;
                    target_current = target_temp;

                    /* compute transform via weighted kabsch */
                    if (iter > 1)
                        trans = core::weightedKabsch(source_current, target_current, weights);
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
                        flags(i) = loss_func_->isInlier(residuals(i), scale) ? 1 : 0;
                    }
                    inlier_indices = eststack::getNonZeroColumnIndicesFromVector(flags);

                    /* update next weights via loss function */
                    Eigen::VectorXf inlier_residuals = residuals(inlier_indices);
                    weights.resize(inlier_residuals.size());
                    for (int j = 0; j < inlier_residuals.size(); j++)
                    {
                        weights(j) = loss_func_->weight(inlier_residuals(j), scale);
                    }

                    /* calculate scale and previous cost for NEXT iteration, not for now! */
                    float cost_diff = std::abs(cost - prev_cost);
                    scale = scale / decay_factor_;
                    prev_cost = cost;

                    /* check convergence */
                    if (cost_diff < cost_tolerance_ || scale < min_scale_)
                    {
                        final_transformation_ = trans;
                        converged_ = true;
                        return;
                    }
                }

                /* no more iteration, had to finish */
                final_transformation_ = trans;
                converged_ = false;
            }

        private:
            /***
             * @brief loss function
             */
            LossConstPtr loss_func_;

            /***
             * @brief 3xN source point cloud
             */
            Eigen::Matrix3Xf source_;

            /***
             * @brief 3xN target point cloud
             */
            Eigen::Matrix3Xf target_;

            /***
             * @brief initial guess for transformation
             */
            Eigen::Isometry3f init_guess_;

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
        };

    } // namespace core
} // namespace eststack

#endif // CORE__IRLS_HPP
