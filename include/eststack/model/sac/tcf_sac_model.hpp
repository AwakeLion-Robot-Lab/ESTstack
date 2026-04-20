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
#include <memory>

// Eigen library
#include <Eigen/Core>
#include <Eigen/Geometry>

// ESTstack library
#include "eststack/model/sac/base_sac_model.hpp"
#include "eststack/core/dcm_to_quat.hpp"
#include "eststack/core/kabsch.hpp"
#include "eststack/eigen_traits.hpp"

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
            using Ptr = std::shared_ptr<OnePtSACModel>;
            using ConstPtr = std::shared_ptr<const OnePtSACModel>;

            /***
             * @brief default constructor
             */
            OnePtSACModel()
            {
                /* one point can fit this model */
                this->sample_size_ = 1;
                /* this model parameter is index */
                this->model_size_ = 1;
            }

            /***
             * @brief fit model with samples
             * @param samples sample indices
             * @param[out] coeffs output model coefficients
             * @details here just pass sample index to coeff, it seems useless because this is NOT a real model (not like equation to solve)
             */
            bool fitImpl(const std::vector<int> &samples, Eigen::VectorXf &coeffs)
            {
                /* check empty */
                if (samples.empty())
                    return false;

                coeffs.resize(model_size_);
                coeffs(0) = static_cast<float>(samples[0]);
                return true;
            }

            /***
             * @brief select inliers with model coefficients
             * @param coeffs model coefficients
             * @details here just pass sample index to coeff, it seems useless because this is NOT a real model (not like equation to solve)
             */
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
                const Eigen::VectorXi flags = (length_consistency.array() < 2.0f * threshold_).cast<int>();
                const Eigen::VectorXi inlier_indices = eststack::getNonZeroColumnIndicesFromVector(flags);

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
                for (int i = 1; i <= 50; i++)
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
                    Eigen::MatrixXi inliers_flags_matrix = (inliers_distance_consistency_matrix.array() < 2.0f * threshold_).cast<int>();
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
                    Eigen::Matrix3Xf new_source_inliers = source_inliers(Eigen::placeholders::all, new_inliers_indices);
                    Eigen::Matrix3Xf new_target_inliers = target_inliers(Eigen::placeholders::all, new_inliers_indices);
                    source_inliers = std::move(new_source_inliers);
                    target_inliers = std::move(new_target_inliers);

                    /* synchronize global indices */
                    std::vector<int> next_indices;
                    next_indices.reserve(new_inliers_num);
                    for (int k = 0; k < new_inliers_num; ++k)
                        next_indices.push_back(current_indices[new_inliers_indices[k]]);
                    current_indices = std::move(next_indices);

                    /***
                     * truncate condition is also an exp., here means ratio of inliers convergence,
                     * if the number of inliers is not significantly reduced, we can stop iteration
                     */
                    if (inliers_num - new_inliers_num < 5)
                        break;
                    inliers_num = new_inliers_num;
                }

                return current_indices;
            }
        };

        /***
         * @brief two points RANSAC model based on angular discrepancy
         * @details refer to docs/explanation/TCF.md
         */
        class TwoPtSACModel final : public BaseSACModel<TwoPtSACModel>
        {
        public:
            using Ptr = std::shared_ptr<TwoPtSACModel>;
            using ConstPtr = std::shared_ptr<const TwoPtSACModel>;

            /***
             * @brief default constructor
             */
            TwoPtSACModel()
            {
                /* two points can fit this model */
                this->sample_size_ = 2;
                /* this model parameter is index */
                this->model_size_ = 2;
            }

            /***
             * @brief fit model with samples
             * @param samples sample indices
             * @param[out] coeffs output model coefficients
             * @details here just pass sample indices to coeffs, it seems useless because this is NOT a real model (not like equation to solve)
             */
            bool fitImpl(const std::vector<int> &samples, Eigen::VectorXf &coeffs)
            {
                /* check empty */
                if (samples.empty())
                    return false;

                coeffs.resize(model_size_);
                coeffs(0) = static_cast<float>(samples[0]);
                coeffs(1) = static_cast<float>(samples[1]);
                return true;
            }

            /***
             * @brief select inliers with model coefficients
             * @param coeffs model coefficients
             * @details here just pass sample indices to coeffs, it seems useless because this is NOT a real model (not like equation to solve)
             */
            std::vector<int> selectInliersImpl(const Eigen::VectorXf &coeffs) const
            {
                /* check length consistency of samples */
                const int sample_index1 = static_cast<int>(std::round(coeffs(0)));
                const int sample_index2 = static_cast<int>(std::round(coeffs(1)));

                const Eigen::Vector3f source_sample1 = source_.col(sample_index1);
                const Eigen::Vector3f target_sample1 = target_.col(sample_index1);
                const Eigen::Vector3f source_sample2 = source_.col(sample_index2);
                const Eigen::Vector3f target_sample2 = target_.col(sample_index2);

                const float source_samples_length_diff = (source_sample1 - source_sample2).norm();
                const float target_samples_length_diff = (target_sample1 - target_sample2).norm();
                const float samples_length_consistency = std::abs(source_samples_length_diff - target_samples_length_diff);
                if (samples_length_consistency > 2.0f * threshold_)
                    return {};

                /* check length consistency of two correspondences */
                const Eigen::VectorXf source_length_diff_1 = (source_.colwise() - source_sample1).colwise().norm().transpose();
                const Eigen::VectorXf target_length_diff_1 = (target_.colwise() - target_sample1).colwise().norm().transpose();
                const Eigen::VectorXf source_length_diff_2 = (source_.colwise() - source_sample2).colwise().norm().transpose();
                const Eigen::VectorXf target_length_diff_2 = (target_.colwise() - target_sample2).colwise().norm().transpose();

                const Eigen::VectorXf length_consistency_1 = (source_length_diff_1 - target_length_diff_1).array().abs();
                const Eigen::VectorXf length_consistency_2 = (source_length_diff_2 - target_length_diff_2).array().abs();
                const Eigen::VectorXi rough_flags = (length_consistency_1.array() < 2.0f * threshold_).cast<int>().cwiseMin((length_consistency_2.array() < 2.0f * threshold_).cast<int>());

                /* "3" means at least 3 points to check angular discrepancy */
                if (rough_flags.sum() < 3)
                    return {};

                /* select candidate inliers */
                const Eigen::VectorXi rough_inlier_indices = eststack::getNonZeroColumnIndicesFromVector(rough_flags);
                std::vector<int> current_indices(rough_inlier_indices.data(), rough_inlier_indices.data() + rough_inlier_indices.size());

                /* start to check angular discrepancy of each inlier */
                Eigen::VectorXi refined_flags = Eigen::VectorXi::Zero(rough_inlier_indices.size());
                for (int i = 0; i < rough_inlier_indices.size(); i++)
                {
                    /* select candidate inliers and calculate each vector of triangles */
                    int idx = rough_inlier_indices(i);
                    Eigen::Vector3f left_source_vec = source_.col(idx) - source_sample1;
                    Eigen::Vector3f left_target_vec = target_.col(idx) - target_sample1;
                    Eigen::Vector3f right_source_vec = source_.col(idx) - source_sample2;
                    Eigen::Vector3f right_target_vec = target_.col(idx) - target_sample2;

                    /***
                     * refer to TCF Fig.5, now there are two corresponding triangle
                     * very important that we coincide source-target candidate inliers, it's easier to compare augular discrepancy of two triangles
                     * that's why we can calculate the angle and check angular discrepancy
                     * cos angle = (a dot b) / (norm(a) * norm(b)) is used here
                     */
                    float left_source_vec_norm = left_source_vec.norm();
                    float left_target_vec_norm = left_target_vec.norm();
                    float right_source_vec_norm = right_source_vec.norm();
                    float right_target_vec_norm = right_target_vec.norm();

                    /* avoid sin(0) or cos(0) */
                    constexpr float denom_eps = std::numeric_limits<float>::epsilon();
                    if (left_source_vec_norm < denom_eps || right_source_vec_norm < denom_eps ||
                        left_target_vec_norm < denom_eps || right_target_vec_norm < denom_eps)
                        continue;

                    float source_cos = left_source_vec.dot(right_source_vec) / (left_source_vec_norm * right_source_vec_norm);
                    float target_cos = left_target_vec.dot(right_target_vec) / (left_target_vec_norm * right_target_vec_norm);
                    /* `std::clamp` to constrain sine/cosine value */
                    float source_angle = std::acos(std::clamp(source_cos, -1.0f, 1.0f));
                    float target_angle = std::acos(std::clamp(target_cos, -1.0f, 1.0f));
                    float angular_discrepancy = std::abs(source_angle - target_angle);

                    /* calculate angular discrepancy threshold */
                    float left_discrepancy_threshold = std::asin(std::clamp(threshold_ / left_source_vec_norm, -1.0f, 1.0f));
                    float right_discrepancy_threshold = std::asin(std::clamp(threshold_ / right_source_vec_norm, -1.0f, 1.0f));
                    float discrepancy_threshold = std::abs(left_discrepancy_threshold + right_discrepancy_threshold);

                    /* if this inlier pass the condition, mark as refined inliers */
                    if (angular_discrepancy < discrepancy_threshold)
                        refined_flags(i) = 1;
                }

                /* get local refined inliers */
                if (refined_flags.sum() == 0)
                    return {};
                const Eigen::VectorXi refined_local = eststack::getNonZeroColumnIndicesFromVector(refined_flags);

                /* synchronize global indices */
                std::vector<int> refined_inlier_indices;
                refined_inlier_indices.reserve(refined_local.size());
                for (int i = 0; i < refined_local.size(); i++)
                    refined_inlier_indices.push_back(current_indices[refined_local(i)]);

                return refined_inlier_indices;
            }
        };

        /***
         * @brief three points RANSAC model based on SVD solution of Wahba Problem
         * @details refer to docs/explanation/TCF.md
         */
        class ThreePtSACModel final : public BaseSACModel<ThreePtSACModel>
        {
        public:
            using Ptr = std::shared_ptr<ThreePtSACModel>;
            using ConstPtr = std::shared_ptr<const ThreePtSACModel>;

            /***
             * @brief default constructor
             */
            ThreePtSACModel()
            {
                /* three points can fit this model */
                this->sample_size_ = 3;
                /* this model parameter is qw qx qy qz and x y z */
                this->model_size_ = 7;
            }

            /***
             * @brief fit model with samples
             * @param samples sample indices
             * @param[out] coeffs output model coefficients
             * @details coefficients sequence: qw qx qy qz x y z
             */
            bool fitImpl(const std::vector<int> &samples, Eigen::VectorXf &coeffs)
            {
                /* check empty */
                if (samples.empty())
                    return false;

                /* check colinear or not */
                Eigen::Vector3f source_vec_1 = source_.col(samples[1]) - source_.col(samples[0]);
                Eigen::Vector3f source_vec_2 = source_.col(samples[2]) - source_.col(samples[0]);
                Eigen::Vector3f target_vec_1 = target_.col(samples[1]) - target_.col(samples[0]);
                Eigen::Vector3f target_vec_2 = target_.col(samples[2]) - target_.col(samples[0]);
                const float eps = std::numeric_limits<float>::epsilon();

                if (source_vec_1.cross(source_vec_2).norm() < eps || target_vec_1.cross(target_vec_2).norm() < eps)
                    return false;

                /* get R and t via kabsch */
                Eigen::Matrix3f source_sample = source_(Eigen::placeholders::all, samples);
                Eigen::Matrix3f target_sample = target_(Eigen::placeholders::all, samples);
                Eigen::Isometry3f T = core::kabsch(source_sample, target_sample);
                Eigen::Quaternionf q = core::fromDCM(T.linear());
                Eigen::Vector3f t = T.translation();

                /* serialize into coefficients in proper sequence */
                coeffs.resize(model_size_);
                coeffs(0) = q.w();
                coeffs(1) = q.x();
                coeffs(2) = q.y();
                coeffs(3) = q.z();
                coeffs(4) = t.x();
                coeffs(5) = t.y();
                coeffs(6) = t.z();

                return true;
            }

            /***
             * @brief select inliers with model coefficients
             * @param coeffs model coefficients
             */
            std::vector<int> selectInliersImpl(const Eigen::VectorXf &coeffs) const
            {
                /* deserialize coefficients */
                Eigen::Quaternionf q(coeffs(0), coeffs(1), coeffs(2), coeffs(3));
                Eigen::Vector3f t(coeffs(4), coeffs(5), coeffs(6));
                Eigen::Matrix3f R = q.toRotationMatrix();

                /* transform source cloud as best estimated cloud */
                Eigen::Matrix3Xf est_source = (R * source_).colwise() + t;

                /* calculate error distance as score */
                Eigen::VectorXf error_distance = (est_source - target_).colwise().norm().transpose();
                const Eigen::VectorXi inlier_flags = (error_distance.array() < threshold_).cast<int>();
                const Eigen::VectorXi inlier_indices = eststack::getNonZeroColumnIndicesFromVector(inlier_flags);

                return std::vector<int>(inlier_indices.data(), inlier_indices.data() + inlier_indices.size());
            }
        };
    } // namespace model
} // namespace eststack

#endif //! MODEL__TCF_SAC_MODEL_HPP
