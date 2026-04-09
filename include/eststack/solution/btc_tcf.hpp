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

#ifndef SOLUTION__BTC_TCF_HPP
#define SOLUTION__BTC_TCF_HPP

// Eigen library
#include <Eigen/Dense>

// PCL library
#include <pcl/point_cloud.h>
#include <pcl/point_types.h>

// ESTstack library
#include "eststack/solution/base_pcr.hpp"

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
        template <typename PointT>
        class BTC_TCF : public BasePCR<BTC_TCF<PointT>, PointT>
        {
        public:
            using Base = BasePCR<BTC_TCF<PointT>, PointT>;
            using typename Base::PointCloud;
            using typename Base::PointCloudConstPtr;
            using typename Base::PointCloudPtr;

            using Base::result_;
            using Base::source_cloud_;
            using Base::target_cloud_;

            bool alignImpl()
            {
                /* transform `pcl::PointCloud` into `Eigen::Matrix3Xf` */
                constexpr auto point_dim = sizeof(PointT);
                const auto source_matrix = source_cloud_->template getMatrixXfMap(point_dim, 8, 0).transpose();
                const auto target_matrix = target_cloud_->template getMatrixXfMap(point_dim, 8, 0).transpose();
            }
        };
    } // namespace solution
} // namespace eststack

#endif //! SOLUTION__BTC_TCF_HPP