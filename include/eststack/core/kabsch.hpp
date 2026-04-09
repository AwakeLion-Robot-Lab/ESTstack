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

#ifndef CORE__KABSCH_HPP
#define CORE__KABSCH_HPP

// Eigen library
#include <Eigen/Dense>
#include <Eigen/SVD>

// PCL library
#include <pcl/point_types.h>
#include <pcl/point_cloud.h>
#include <pcl/common/centroid.h>

// oneTBB library
#include <oneapi/tbb.h>

/***
 * @brief An algorithm set focus on state estimation
 * @author jinhua "siyiovo" deng
 */
namespace eststack
{
    /***
     * @brief core algorithms for state estimation
     */
    namespace core
    {
        /***
         * @brief get best rigid transformation via SVD
         * @param source source point cloud
         * @param target target point cloud
         * @return transformation from source to target
         */
        Eigen::Isometry3f kabsch(const pcl::PointCloud<pcl::PointXYZ>::Ptr &source,
                                 const pcl::PointCloud<pcl::PointXYZ>::Ptr &target)
        {
            const size_t N = source->size();
            if (N == 0)
                throw std::invalid_argument("pointcloud is empty!");
            if (N != target->size())
                throw std::invalid_argument("source/target pointcloud size mismatch!");

            /* compute centroids of each pointcloud */
            Eigen::Vector4f source_centroid4, target_centroid4;
            pcl::compute3DCentroid(*source, source_centroid4);
            pcl::compute3DCentroid(*target, target_centroid4);

            const Eigen::Vector3f source_centroid = source_centroid4.head<3>();
            const Eigen::Vector3f target_centroid = target_centroid4.head<3>();

            /* compute covariance matrix using TBB parallel reduction */
            tbb::combinable<Eigen::Matrix3f> cov([]
                                                 { return Eigen::Matrix3f::Zero(); });
            tbb::parallel_for(
                tbb::blocked_range<size_t>(0, N, 1024),
                [&](const tbb::blocked_range<size_t> &r)
                {
                    Eigen::Matrix3f &local_H = cov.local();
                    for (size_t i = r.begin(); i != r.end(); ++i)
                    {
                        /* compute translation vector between each pt to centroid both in two pointclouds */
                        const Eigen::Vector3f p = source->points[i].getVector3fMap() - source_centroid;
                        const Eigen::Vector3f q = target->points[i].getVector3fMap() - target_centroid;
                        local_H.noalias() += q * p.transpose();
                    }
                });

            Eigen::Matrix3f H = Eigen::Matrix3f::Zero();
            cov.combine_each([&](const Eigen::Matrix3f &local)
                             { H += local; });

            /* SVD for best rotation */
            Eigen::JacobiSVD<Eigen::Matrix3f> svd(H, Eigen::ComputeFullU | Eigen::ComputeFullV);
            const auto &U = svd.matrixU();
            const auto &V = svd.matrixV();

            /* constrain rotation matrix to SO(3) about det(R) = 1 */
            Eigen::Matrix3f R = V * U.transpose();
            if (U.determinant() * V.determinant() < 0)
                R.col(2) = -R.col(2); /* if reflection, flip the sign */

            /* compute translation */
            const Eigen::Vector3f t = target_centroid - R * source_centroid;

            Eigen::Isometry3f T = Eigen::Isometry3f::Identity();
            T.linear() = R;
            T.translation() = t;
            return T;
        }

        /***
         * @brief get best rigid transformation via SVD
         * @param source 3xN source points
         * @param target 3xN target points
         * @return transformation from source to target
         * @details each column of `source` and `target` is a point, and the order of points should be consistent between two sets
         */
        Eigen::Isometry3f kabsch(const Eigen::Matrix3Xf &source, const Eigen::Matrix3Xf &target)
        {
            const int N = source.cols();
            if (N == 0 || N != target.cols())
                return Eigen::Isometry3f::Identity();

            /* compute centroids */
            Eigen::Vector3f source_centroid = source.rowwise().mean();
            Eigen::Vector3f target_centroid = target.rowwise().mean();

            /* compute covariance matrix using TBB parallel reduction */
            tbb::combinable<Eigen::Matrix3f> cov([]
                                                 { return Eigen::Matrix3f::Zero(); });
            tbb::parallel_for(
                tbb::blocked_range<int>(0, N, 1024),
                [&](const tbb::blocked_range<int> &r)
                {
                    Eigen::Matrix3f &local_H = cov.local();
                    for (int i = r.begin(); i != r.end(); ++i)
                    {
                        /* compute translation vector between each pt to centroid both in two pointclouds */
                        const Eigen::Vector3f p = source.col(i) - source_centroid;
                        const Eigen::Vector3f q = target.col(i) - target_centroid;
                        local_H.noalias() += q * p.transpose();
                    }
                });

            Eigen::Matrix3f H = Eigen::Matrix3f::Zero();
            cov.combine_each([&](const Eigen::Matrix3f &local)
                             { H += local; });

            /* SVD for best rotation */
            Eigen::JacobiSVD<Eigen::Matrix3f> svd(H, Eigen::ComputeFullU | Eigen::ComputeFullV);
            const auto &U = svd.matrixU();
            const auto &V = svd.matrixV();

            /* constrain rotation matrix to SO(3) about det(R) = 1 */
            Eigen::Matrix3f R = V * U.transpose();
            if (U.determinant() * V.determinant() < 0)
                R.col(2) = -R.col(2); /* if reflection, flip the sign */

            /* compute translation */
            const Eigen::Vector3f t = target_centroid - R * source_centroid;

            Eigen::Isometry3f T = Eigen::Isometry3f::Identity();
            T.linear() = R;
            T.translation() = t;
            return T;
        }

        /***
         * @brief get best rigid transformation via SVD with weights
         * @param source source point cloud
         * @param target target point cloud
         * @param weights weights of each point
         * @return transformation from source to target
         */
        Eigen::Isometry3f weighted_kabsch(const pcl::PointCloud<pcl::PointXYZ>::Ptr &source,
                                          const pcl::PointCloud<pcl::PointXYZ>::Ptr &target,
                                          const Eigen::VectorXf &weights)
        {
            const size_t N = source->size();
            if (N == 0)
                throw std::invalid_argument("pointcloud is empty!");
            if (N != target->size())
                throw std::invalid_argument("source/target pointcloud size mismatch!");

            /* check weight sum */
            float weight_sum = weights.sum();
            /* if weight sum is zero, use uniform weights as no preference circumstance */
            if (weight_sum < std::numeric_limits<float>::epsilon())
            {
                Eigen::VectorXf uniform_weights = Eigen::VectorXf::Ones(N) / N;
                return weighted_kabsch(source, target, uniform_weights);
            }

            Eigen::VectorXf w = weights / weight_sum;

            /* compute weighted centroids */
            Eigen::Vector3f source_centroid = Eigen::Vector3f::Zero();
            Eigen::Vector3f target_centroid = Eigen::Vector3f::Zero();
            for (size_t i = 0; i < N; ++i)
            {
                source_centroid += w(i) * source->points[i].getVector3fMap();
                target_centroid += w(i) * target->points[i].getVector3fMap();
            }

            /* compute weighted covariance matrix using TBB parallel reduction */
            tbb::combinable<Eigen::Matrix3f> cov([]
                                                 { return Eigen::Matrix3f::Zero(); });
            tbb::parallel_for(
                tbb::blocked_range<size_t>(0, N, 1024),
                [&](const tbb::blocked_range<size_t> &r)
                {
                    Eigen::Matrix3f &local_H = cov.local();
                    for (size_t i = r.begin(); i != r.end(); ++i)
                    {
                        const float w_sqrt = std::sqrt(w(i));
                        const Eigen::Vector3f p = (source->points[i].getVector3fMap() - source_centroid) * w_sqrt;
                        const Eigen::Vector3f q = (target->points[i].getVector3fMap() - target_centroid) * w_sqrt;
                        local_H.noalias() += q * p.transpose();
                    }
                });

            Eigen::Matrix3f H = Eigen::Matrix3f::Zero();
            cov.combine_each([&](const Eigen::Matrix3f &local)
                             { H += local; });

            /* SVD for best rotation */
            Eigen::JacobiSVD<Eigen::Matrix3f> svd(H, Eigen::ComputeFullU | Eigen::ComputeFullV);
            const auto &U = svd.matrixU();
            const auto &V = svd.matrixV();

            /* constrain rotation matrix to SO(3) */
            Eigen::Matrix3f R = V * U.transpose();
            if (U.determinant() * V.determinant() < 0)
                R.col(2) = -R.col(2);

            /* compute translation */
            const Eigen::Vector3f t = target_centroid - R * source_centroid;

            Eigen::Isometry3f T = Eigen::Isometry3f::Identity();
            T.linear() = R;
            T.translation() = t;
            return T;
        }

        /***
         * @brief get best rigid transformation via SVD with weights
         * @param source 3xN source points
         * @param target 3xN target points
         * @param weights weights of each point
         * @details each column of `source` and `target` is a point, and the order of points should be consistent between two sets
         * @return transformation from source to target
         */
        Eigen::Isometry3f weighted_kabsch(const Eigen::Matrix3Xf &source,
                                          const Eigen::Matrix3Xf &target,
                                          const Eigen::VectorXf &weights)
        {
            const int N = source.cols();
            if (N == 0 || N != target.cols())
                return Eigen::Isometry3f::Identity();

            /* check weight sum */
            float weight_sum = weights.sum();
            /* if weight sum is zero, use uniform weights as no preference circumstance */
            if (weight_sum < std::numeric_limits<float>::epsilon())
            {
                Eigen::VectorXf uniform_weights = Eigen::VectorXf::Ones(N) / N;
                return weighted_kabsch(source, target, uniform_weights);
            }

            /* normalized weights */
            Eigen::VectorXf w = weights / weight_sum;

            /* compute weighted centroids */
            Eigen::Vector3f source_centroid = source * w;
            Eigen::Vector3f target_centroid = target * w;

            /* compute weighted covariance matrix using TBB parallel reduction */
            tbb::combinable<Eigen::Matrix3f> cov([]
                                                 { return Eigen::Matrix3f::Zero(); });
            tbb::parallel_for(
                tbb::blocked_range<int>(0, N, 1024),
                [&](const tbb::blocked_range<int> &r)
                {
                    Eigen::Matrix3f &local_H = cov.local();
                    for (int i = r.begin(); i != r.end(); ++i)
                    {
                        const float w_sqrt = std::sqrt(w(i));
                        const Eigen::Vector3f p = (source.col(i) - source_centroid) * w_sqrt;
                        const Eigen::Vector3f q = (target.col(i) - target_centroid) * w_sqrt;
                        local_H.noalias() += q * p.transpose();
                    }
                });

            Eigen::Matrix3f H = Eigen::Matrix3f::Zero();
            cov.combine_each([&](const Eigen::Matrix3f &local)
                             { H += local; });

            /* SVD for best rotation */
            Eigen::JacobiSVD<Eigen::Matrix3f> svd(H, Eigen::ComputeFullU | Eigen::ComputeFullV);
            const auto &U = svd.matrixU();
            const auto &V = svd.matrixV();

            /* constrain rotation matrix to SO(3) */
            Eigen::Matrix3f R = V * U.transpose();
            if (U.determinant() * V.determinant() < 0)
                R.col(2) = -R.col(2);

            /* compute translation */
            const Eigen::Vector3f t = target_centroid - R * source_centroid;

            Eigen::Isometry3f T = Eigen::Isometry3f::Identity();
            T.linear() = R;
            T.translation() = t;
            return T;
        }
    } // namespace core
} // namespace eststack

#endif //! CORE__KABSCH_HPP
