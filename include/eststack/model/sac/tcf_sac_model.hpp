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

#ifndef MODEL__TCF_SAC_MODEL_HPP
#define MODEL__TCF_SAC_MODEL_HPP

// C++ standard library
#include <algorithm>
#include <cmath>
#include <numeric>
#include <vector>

// Eigen library
#include <Eigen/Core>
#include <Eigen/Dense>

// ESTstack library
#include "eststack/eigen_traits.hpp"
#include "eststack/model/sac/base_sac_model.hpp"

/***
 * @brief An algorithm set focus on state estimation
 * @author jinhua "siyiovo" deng
 */
namespace eststack
{
    /***
     * @brief models for problems
     */
    namespace model
    {
        /***
         * @brief one point RANSAC model based on length consistency
         * @details refer to docs/explanation/TCF.md
         */
        class OnePtSACModel final : public BaseSACModel<OnePtSACModel>
        {
        public:
            /***
             * @brief default constructor
             */
            OnePtSACModel()
            {
                /* one point can fit this model */
                sample_size_ = 1;
                /* this model parameter is index */
                model_size_ = 1;
            }

            std::vector<int> selectInliersImpl(const Eigen::VectorXf &coeffs) const
            {
                /* get sample anchor point */
                const int sample_index = static_cast<int>(std::round(coeffs(0)));
                const Eigen::Vector3f source_sample = source_.col(sample_index);
                const Eigen::Vector3f target_sample = target_.col(sample_index);

                /* compute length consistency, just watch out which one is 1xN or Nx1 */
                const Eigen::VectorXf source_length_diff = (source_.colwise() - source_sample).colwise().norm().transpose();
                const Eigen::VectorXf target_length_diff = (target_.colwise() - target_sample).colwise().norm().transpose();
                const Eigen::VectorXf length_consistency = (source_length_diff - target_length_diff).array().abs();

                /* select initial consensus set via length consistency */
                Eigen::VectorXi flags = (length_consistency.array() < threshold_).cast<int>();
                Eigen::VectorXi inlier_indices = eststack::getNonZeroColumnIndicesFromVector(flags);

                /* maintain mapping from local indices back to original point cloud */
                std::vector<int> current_indices(inlier_indices.data(), inlier_indices.data() + inlier_indices.size());

                /***
                 * now we get initial source inliers and target inliers
                 * NOTE that source_current(Eigen::placeholders::all, inlier_indices) is slicing operation in Eigen,
                 * it means select all rows but specfic column that are in inlier_indices
                 * e.g. if source_current is 3x100 and inlier_indices is [0, 2, 5],
                 * then it turns out to be 3x3 matrix consisting of columns 0, 2, 5 of `source_current`
                 * slicing operation can refer to https://zhuanlan.zhihu.com/p/531476504
                 * also see `core/irls.hpp`
                 */
                Eigen::Matrix3Xf source_inliers = source_(Eigen::placeholders::all, inlier_indices);
                Eigen::Matrix3Xf target_inliers = target_(Eigen::placeholders::all, inlier_indices);

                /***
                 * refine initial consensus set, NOTE that we will use maximal clique to solve,
                 * it just mention in [TCF issue#1](https://github.com/ShiPC-AI/TCF/issues/1)
                 * concept of maximal clique refer to https://blog.csdn.net/u010608296/article/details/119834724
                 * if you don't know graph theory, please refer to docs/explanation/TCF.md to seek an easier way to understand the idea
                 */
                int inliers_num = static_cast<int>(current_indices.size());
                /* `50` is an exp. iteration number, you can adjust */
                for (int i = 0; i <= 50; i++)
                {
                    /* we gotta calculate inliers distance matrix that distance about each inliers */
                    Eigen::MatrixXf source_dist_matrix(inliers_num, inliers_num);
                    Eigen::MatrixXf target_dist_matrix(inliers_num, inliers_num);
                    /***
                     * calculate each rows, which means distance about point j to others, e.g. row(1) = dist(1,k)
                     * it's an symmtric matrix acutually, so you can use `col(j) = source_inliers.colwise() - source_inliers.col(j)` to calculate each column too
                     */
                    for (int j = 0; j < inliers_num; j++)
                    {
                        source_dist_matrix.row(j) = (source_inliers.colwise() - source_inliers.col(j)).colwise().norm().eval();
                        target_dist_matrix.row(j) = (target_inliers.colwise() - target_inliers.col(j)).colwise().norm().eval();
                    }

                    /* calculate distance consistency of correspondences */
                    Eigen::MatrixXf inliers_distance_consistency_matrix = (source_dist_matrix - target_dist_matrix).array().abs();
                    Eigen::MatrixXi inliers_flags_matrix = (inliers_distance_consistency_matrix.array() < threshold_).cast<int>();
                    /* sqaure of the number of points in maximal clique equals to the number of real inliers approx. */
                    int new_inliers_num = static_cast<int>(std::ceil(std::sqrt(static_cast<float>(inliers_flags_matrix.sum()))));

                    /* the sum of flags about this point means how many other inliers are consistent with it */
                    Eigen::VectorXi inliers_connection = inliers_flags_matrix.rowwise().sum();
                    /* `idx_vec` just stores indices */
                    std::vector<int> idx_vec(inliers_num);
                    std::iota(idx_vec.begin(), idx_vec.end(), 0);
                    /* sort indices by their connection values in descending order */
                    std::sort(idx_vec.begin(), idx_vec.end(), [&inliers_connection](int a, int b)
                              { return inliers_connection(a) > inliers_connection(b); });
                    /* we only select top new_inliers_num indices */
                    std::vector<int> new_inliers_indices(idx_vec.begin(), idx_vec.begin() + new_inliers_num);

                    /* now we get refined correspondence */
                    Eigen::Matrix3Xf new_source_inliers(Eigen::placeholders::all, new_inliers_indices);
                    Eigen::Matrix3Xf new_target_inliers(Eigen::placeholders::all, new_inliers_indices);
                    source_inliers = std::move(new_source_inliers);
                    target_inliers = std::move(new_target_inliers);

                    /* synchronize global indices */
                    std::vector<int> next_indices;
                    next_indices.reserve(new_inliers_num);
                    for (int k = 0; k < new_inliers_num; ++k)
                        next_indices.push_back(current_indices[new_inliers_indices[k]]);
                    current_indices = std::move(next_indices);

                    /* truncate condition is also exp. */
                    if (inliers_num - new_inliers_num < 5)
                        break;
                    inliers_num = new_inliers_num;
                }

                return current_indices;
            }

            /***
             * @brief fit model with samples
             * @param samples sample indices
             * @param coeffs output model coefficients
             * @details here just pass sample indices to coeffs, it seems useless because this is NOT a real model (not like equation to solve)
             */
            bool fitImpl(const std::vector<int> &samples, Eigen::VectorXf &coeffs)
            {
                if (samples.empty())
                    return false;

                coeffs.resize(model_size_);
                coeffs(0) = static_cast<float>(samples[0]);
                return true;
            }
        };
    } // namespace model
} // namespace eststack

#endif //! MODEL__TCF_SAC_MODEL_HPP
