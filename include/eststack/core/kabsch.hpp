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

// manif library
#include <manif/SE3.h>

// PCL library
#include <pcl/point_types.h>
#include <pcl/point_cloud.h>

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
         * @brief Kabsch algorithm for point cloud registration
         * @param source source point cloud
         * @param target target point cloud
         * @return transformation from source to target
         */
        manif::SE3d kabsch(const pcl::PointCloud<pcl::PointXYZ>::Ptr &source, const pcl::PointCloud<pcl::PointXYZ>::Ptr &target);
    }
}

#endif //! CORE__KABSCH_HPP