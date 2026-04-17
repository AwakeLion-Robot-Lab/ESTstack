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

// C++ standard library
#include <limits>
#include <exception>

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
         * @brief check weights and compute their sum
         * @tparam Index point index type
         * @param weights weights of each point
         * @param N number of points
         */
        template <typename Index>
        float checkWeights(const Eigen::VectorXf &weights, Index N)
        {
            const Eigen::Index eigen_N = static_cast<Eigen::Index>(N);
            if (weights.size() != eigen_N)
                throw std::invalid_argument("weights size mismatch!");

            /* compute sum of weights */
            float weight_sum = 0.0f;
            bool has_nonzero_weight = false;
            for (Eigen::Index i = 0; i < eigen_N; ++i)
            {
                const float weight = weights(i);
                if (weight < 0.0f)
                    throw std::invalid_argument("weights must be non-negative!");
                has_nonzero_weight = has_nonzero_weight || (weight > 0.0f);
                weight_sum += weight;
            }

            /* if have zero weights, return 0.0f */
            if (!has_nonzero_weight)
                return 0.0f;

            return weight_sum;
        }

        /***
         * @brief solve and get transformation via SVD
         * @param H covariance matrix
         * @param source_centroid centroid of source points
         * @param target_centroid centroid of target points
         */
        inline Eigen::Isometry3f solveTransform(const Eigen::Matrix3f &H,
                                                const Eigen::Vector3f &source_centroid,
                                                const Eigen::Vector3f &target_centroid)
        {
            /* SVD for best rotation */
            Eigen::JacobiSVD<Eigen::Matrix3f> svd(H, Eigen::ComputeFullU | Eigen::ComputeFullV);
            const auto &U = svd.matrixU();
            const auto &V = svd.matrixV();

            /* constrain rotation matrix to SO(3) */
            Eigen::Matrix3f D = Eigen::Matrix3f::Identity();
            if (U.determinant() * V.determinant() < 0.0f)
                D(2, 2) = -1.0f; /* if reflection, flip the sign */

            /* compute translation */
            Eigen::Isometry3f T = Eigen::Isometry3f::Identity();
            T.linear() = U * D * V.transpose();
            T.translation() = target_centroid - T.linear() * source_centroid;
            return T;
        }

        /***
         * @brief get best rigid transformation via SVD
         * @param source source point cloud
         * @param target target point cloud
         * @return transformation from source to target
         */
        inline Eigen::Isometry3f kabsch(const pcl::PointCloud<pcl::PointXYZ>::Ptr &source,
                                        const pcl::PointCloud<pcl::PointXYZ>::Ptr &target)
        {
            if (!source || !target)
                throw std::invalid_argument("pointcloud pointer is null!");

            const std::size_t N = source->size();
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
            Eigen::Matrix3f H = Eigen::Matrix3f::Zero();
            /* if points num smaller than 2048, use normal loop */
            if (N < static_cast<std::size_t>(2048))
            {
                for (std::size_t i = 0; i < N; ++i)
                {
                    const Eigen::Vector3f p = source->points[i].getVector3fMap() - source_centroid;
                    const Eigen::Vector3f q = target->points[i].getVector3fMap() - target_centroid;
                    H.noalias() += q * p.transpose();
                }
            }
            else
            {
                tbb::combinable<Eigen::Matrix3f> cov([]
                                                     { return Eigen::Matrix3f::Zero(); });
                tbb::parallel_for(
                    tbb::blocked_range<std::size_t>(0, N, static_cast<std::size_t>(1024)),
                    [&](const tbb::blocked_range<std::size_t> &r)
                    {
                        Eigen::Matrix3f &local_H = cov.local();
                        for (size_t i = r.begin(); i != r.end(); ++i)
                        {
                            const Eigen::Vector3f p = source->points[i].getVector3fMap() - source_centroid;
                            const Eigen::Vector3f q = target->points[i].getVector3fMap() - target_centroid;
                            local_H.noalias() += q * p.transpose();
                        }
                    });
                cov.combine_each([&](const Eigen::Matrix3f &local)
                                 { H += local; });
            }

            return solveTransform(H, source_centroid, target_centroid);
        }

        /***
         * @brief get best rigid transformation via SVD
         * @param source 3xN source points
         * @param target 3xN target points
         * @return transformation from source to target
         * @details each column of `source` and `target` is a point, and the order of points should be consistent between two sets
         */
        inline Eigen::Isometry3f kabsch(const Eigen::Matrix3Xf &source, const Eigen::Matrix3Xf &target)
        {
            const int N = source.cols();
            if (N == 0 || N != target.cols())
                return Eigen::Isometry3f::Identity();

            /* compute centroids */
            const Eigen::Vector3f source_centroid = source.rowwise().mean();
            const Eigen::Vector3f target_centroid = target.rowwise().mean();

            /* compute covariance matrix using TBB parallel reduction */
            Eigen::Matrix3f H = Eigen::Matrix3f::Zero();
            if (N < 2048)
            {
                for (int i = 0; i < N; ++i)
                {
                    const Eigen::Vector3f p = source.col(i) - source_centroid;
                    const Eigen::Vector3f q = target.col(i) - target_centroid;
                    H.noalias() += q * p.transpose();
                }
            }
            else
            {
                tbb::combinable<Eigen::Matrix3f> cov([]
                                                     { return Eigen::Matrix3f::Zero(); });
                tbb::parallel_for(
                    tbb::blocked_range<int>(0, N, 1024),
                    [&](const tbb::blocked_range<int> &r)
                    {
                        Eigen::Matrix3f &local_H = cov.local();
                        for (int i = r.begin(); i != r.end(); ++i)
                        {
                            const Eigen::Vector3f p = source.col(i) - source_centroid;
                            const Eigen::Vector3f q = target.col(i) - target_centroid;
                            local_H.noalias() += q * p.transpose();
                        }
                    });
                cov.combine_each([&](const Eigen::Matrix3f &local)
                                 { H += local; });
            }

            return solveTransform(H, source_centroid, target_centroid);
        }

        /***
         * @brief get best rigid transformation via SVD with weights
         * @param source source point cloud
         * @param target target point cloud
         * @param weights weights of each point
         * @return transformation from source to target
         */
        inline Eigen::Isometry3f weightedKabsch(const pcl::PointCloud<pcl::PointXYZ>::Ptr &source,
                                                const pcl::PointCloud<pcl::PointXYZ>::Ptr &target,
                                                const Eigen::VectorXf &weights)
        {
            if (!source || !target)
                throw std::invalid_argument("pointcloud pointer is null!");

            const size_t N = source->size();
            if (N == 0)
                throw std::invalid_argument("pointcloud is empty!");
            if (N != target->size())
                throw std::invalid_argument("source/target pointcloud size mismatch!");

            /* check weights */
            const float weight_sum = checkWeights(weights, N);
            if (weight_sum == 0.0f)
                return kabsch(source, target);
            const float inv_weight_sum = 1.0f / weight_sum;

            /* compute weighted centroids */
            Eigen::Vector3f source_centroid = Eigen::Vector3f::Zero();
            Eigen::Vector3f target_centroid = Eigen::Vector3f::Zero();
            for (size_t i = 0; i < N; ++i)
            {
                source_centroid += weights(i) * source->points[i].getVector3fMap();
                target_centroid += weights(i) * target->points[i].getVector3fMap();
            }
            source_centroid *= inv_weight_sum;
            target_centroid *= inv_weight_sum;

            /* compute weighted covariance matrix using TBB parallel reduction */
            Eigen::Matrix3f H = Eigen::Matrix3f::Zero();
            if (N < static_cast<size_t>(2048))
            {
                for (size_t i = 0; i < N; ++i)
                {
                    const Eigen::Vector3f p = source->points[i].getVector3fMap() - source_centroid;
                    const Eigen::Vector3f q = target->points[i].getVector3fMap() - target_centroid;
                    H.noalias() += weights(i) * q * p.transpose();
                }
            }
            else
            {
                tbb::combinable<Eigen::Matrix3f> cov([]
                                                     { return Eigen::Matrix3f::Zero(); });
                tbb::parallel_for(
                    tbb::blocked_range<size_t>(0, N, static_cast<size_t>(1024)),
                    [&](const tbb::blocked_range<size_t> &r)
                    {
                        Eigen::Matrix3f &local_H = cov.local();
                        for (size_t i = r.begin(); i != r.end(); ++i)
                        {
                            const Eigen::Vector3f p = source->points[i].getVector3fMap() - source_centroid;
                            const Eigen::Vector3f q = target->points[i].getVector3fMap() - target_centroid;
                            local_H.noalias() += weights(i) * q * p.transpose();
                        }
                    });
                cov.combine_each([&](const Eigen::Matrix3f &local)
                                 { H += local; });
            }
            H *= inv_weight_sum;

            return solveTransform(H, source_centroid, target_centroid);
        }

        /***
         * @brief get best rigid transformation via SVD with weights
         * @param source 3xN source points
         * @param target 3xN target points
         * @param weights weights of each point
         * @details each column of `source` and `target` is a point, and the order of points should be consistent between two sets
         * @return transformation from source to target
         */
        inline Eigen::Isometry3f weightedKabsch(const Eigen::Matrix3Xf &source,
                                                const Eigen::Matrix3Xf &target,
                                                const Eigen::VectorXf &weights)
        {
            const int N = source.cols();
            if (N == 0 || N != target.cols())
                return Eigen::Isometry3f::Identity();

            /* check weights */
            const float weight_sum = checkWeights(weights, N);
            if (weight_sum == 0.0f)
                return kabsch(source, target);
            const float inv_weight_sum = 1.0f / weight_sum;

            /* compute weighted centroids */
            const Eigen::Vector3f source_centroid = (source * weights) * inv_weight_sum;
            const Eigen::Vector3f target_centroid = (target * weights) * inv_weight_sum;

            /* compute weighted covariance matrix using TBB parallel reduction */
            Eigen::Matrix3f H = Eigen::Matrix3f::Zero();
            if (N < 2048)
            {
                for (int i = 0; i < N; ++i)
                {
                    const Eigen::Vector3f p = source.col(i) - source_centroid;
                    const Eigen::Vector3f q = target.col(i) - target_centroid;
                    H.noalias() += weights(i) * q * p.transpose();
                }
            }
            else
            {
                tbb::combinable<Eigen::Matrix3f> cov([]
                                                     { return Eigen::Matrix3f::Zero(); });
                tbb::parallel_for(
                    tbb::blocked_range<int>(0, N, 1024),
                    [&](const tbb::blocked_range<int> &r)
                    {
                        Eigen::Matrix3f &local_H = cov.local();
                        for (int i = r.begin(); i != r.end(); ++i)
                        {
                            const Eigen::Vector3f p = source.col(i) - source_centroid;
                            const Eigen::Vector3f q = target.col(i) - target_centroid;
                            local_H.noalias() += weights(i) * q * p.transpose();
                        }
                    });
                cov.combine_each([&](const Eigen::Matrix3f &local)
                                 { H += local; });
            }
            H *= inv_weight_sum;

            return solveTransform(H, source_centroid, target_centroid);
        }
    } // namespace core
} // namespace eststack

#endif //! CORE__KABSCH_HPP
