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

#ifndef SOLUTION__BASE_PCR_HPP
#define SOLUTION__BASE_PCR_HPP

// C++ standard library
#include <concepts>

// Eigen library
#include <Eigen/Dense>

// PCL library
#include <pcl/point_cloud.h>
#include <pcl/point_types.h>

// ESTstack library
#include "eststack/types.hpp"
#include "eststack/eigen_traits.hpp"

/***
 * @brief An algorithm set focus on state estimation
 * @author jinhua "siyiovo" deng
 */
namespace eststack
{
    /***
     * @brief algorithms for problems
     */
    namespace solution
    {
        /***
         * @brief pointcloud registration result
         */
        struct PCRResult
        {
            Eigen::Isometry3d best_trans_;
            double score_;
            double inliner_fraction_;
        };

        /***
         * @brief base class for pointcloud registration algorithms
         * @tparam Derived derived PCR algorithm class
         * @tparam PointCloudT point type
         */
        template <typename Derived, typename PointT>
        class BasePCR
        {
        public:
            using Point = PointT;
            using PointCloud = pcl::PointCloud<PointT>;
            using PointCloudPtr = typename PointCloud::Ptr;
            using PointCloudConstPtr = typename PointCloud::ConstPtr;

            /***
             * @brief set source cloud
             * @param source source point cloud
             */
            void setSourceCloud(PointCloudConstPtr source)
            {
                source_cloud_ = source;
            }

            /***
             * @brief set target cloud
             * @param target target point cloud
             */
            void setTargetCloud(PointCloudConstPtr target)
            {
                target_cloud_ = target;
            }

            bool align()
            {
                return static_cast<Derived *>(this)->alignImpl();
            }

        protected:
            BasePCR() = default;

            /***
             * @brief source point cloud to be aligned
             */
            PointCloudConstPtr source_cloud_;

            /***
             * @brief target point cloud, e.g. global map
             */
            PointCloudConstPtr target_cloud_;

            /***
             * @brief registration result
             */
            PCRResult result_;
        };
    }
}

#endif //! SOLUTION__BASE_PCR_HPP