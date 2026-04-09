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
#include <cmath>

// Eigen library
#include <Eigen/Dense>

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
     * @brief models for problems
     */
    namespace model
    {
        template <typename PointT>
        class OnePtRANSACModel
        {
        public:
            using PointCloud = pcl::PointCloud<PointT>;
            using PointCloudPtr = PointCloud::Ptr;
            using PointCloudConstPtr = PointCloud::ConstPtr;

        private:
            int sample_size_ = 1;
        };

        template <typename PointT>
        class TwoPtRANSACModel
        {
        public:
            using PointCloud = pcl::PointCloud<PointT>;
            using PointCloudPtr = PointCloud::Ptr;
            using PointCloudConstPtr = PointCloud::ConstPtr;

        private:
            int sample_size_ = 2;
        };

        template <typename PointT>
        class ThreePtRANSACModel
        {
        public:
            using PointCloud = pcl::PointCloud<PointT>;
            using PointCloudPtr = PointCloud::Ptr;
            using PointCloudConstPtr = PointCloud::ConstPtr;

        private:
            int sample_size_ = 3;
        };
    } // namespace model
} // namespace eststack

#endif //! MODEL__TCF_SAC_MODEL_HPP
