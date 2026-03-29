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
        Eigen::Matrix4d kabsch(const pcl::PointCloud<pcl::PointXYZ>::Ptr &source,
                               const pcl::PointCloud<pcl::PointXYZ>::Ptr &target)
        {
            const size_t N = source->size();
            if (N == 0)
                throw std::invalid_argument("pointcloud is empty!");
            if (N != target->size())
                throw std::invalid_argument("source/target pointcloud size mismatch!");

            /* compute centroids of each pointcloud */
            Eigen::Vector4d source_centroid4, target_centroid4;
            pcl::compute3DCentroid(*source, source_centroid4);
            pcl::compute3DCentroid(*target, target_centroid4);

            const Eigen::Vector3d source_centroid = source_centroid4.head<3>();
            const Eigen::Vector3d target_centroid = target_centroid4.head<3>();

            /* compute covariance matrix using TBB parallel reduction */
            tbb::combinable<Eigen::Matrix3d> cov([]
                                                 { return Eigen::Matrix3d::Zero(); });
            tbb::parallel_for(
                tbb::blocked_range<size_t>(0, N, 1024),
                [&](const tbb::blocked_range<size_t> &r)
                {
                    Eigen::Matrix3d &local_H = cov.local();
                    for (size_t i = r.begin(); i != r.end(); ++i)
                    {
                        /* compute translation vector between each pt to centroid both in two pointclouds */
                        const Eigen::Vector3d p = source->points[i].getVector3fMap().cast<double>() - source_centroid;
                        const Eigen::Vector3d q = target->points[i].getVector3fMap().cast<double>() - target_centroid;
                        local_H.noalias() += q * p.transpose();
                    }
                });

            Eigen::Matrix3d H = Eigen::Matrix3d::Zero();
            cov.combine_each([&](const Eigen::Matrix3d &local)
                             { H += local; });

            /* SVD for best rotation */
            Eigen::JacobiSVD<Eigen::Matrix3d> svd(H, Eigen::ComputeFullU | Eigen::ComputeFullV);
            const auto &U = svd.matrixU();
            const auto &V = svd.matrixV();

            /* constrain rotation matrix to SO(3) about det(R) = 1 */
            Eigen::Matrix3d R = V * U.transpose();
            if (U.determinant() * V.determinant() < 0)
                R.col(2) = -R.col(2); /* if reflection, flip the sign */

            /* compute translation */
            const Eigen::Vector3d t = target_centroid - R * source_centroid;

            Eigen::Matrix4d T = Eigen::Matrix4d::Identity();
            T.topLeftCorner<3, 3>() = R;
            T.topRightCorner<3, 1>() = t;
            return T;
        }
    }
}

#endif //! CORE__KABSCH_HPP