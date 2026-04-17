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
#include <memory>

// Eigen library
#include <Eigen/Dense>

// PCL library
#include <pcl/correspondence.h>
#include <pcl/point_cloud.h>
#include <pcl/point_types.h>

// ESTstack library
#include "eststack/core/voxelset.hpp"
#include "eststack/eigen_traits.hpp"
#include "eststack/types.hpp"

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
            /***
             * @brief best rigid transformation
             */
            Eigen::Isometry3f transformation_{Eigen::Isometry3f::Identity()};

            /***
             * @brief inlier fraction
             */
            float inlier_fraction_{0.0f};

            /***
             * @brief registration score
             */
            float score_{0.0f};

            /***
             * @brief whether the registration has converged
             */
            bool converged_{false};
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
            using PointCloud = pcl::PointCloud<PointT>;
            using PointCloudPtr = typename PointCloud::Ptr;
            using PointCloudConstPtr = typename PointCloud::ConstPtr;

            /***
             * @brief set source cloud
             * @param source source point cloud
             */
            void setInputSource(const PointCloudConstPtr &source)
            {
                source_cloud_ = source;
            }

            /***
             * @brief set target cloud
             * @param target target point cloud
             */
            void setInputTarget(const PointCloudConstPtr &target)
            {
                target_cloud_ = target;
            }

            /***
             * @brief set correspondences
             * @param source_match source point coordinates of source-target correspondences
             * @param target_match target point coordinates of source-target correspondences
             */
            void setInputCorrespondences(const Eigen::Matrix3Xf &source_match, const Eigen::Matrix3Xf &target_match)
            {
                source_match_ = source_match;
                target_match_ = target_match;
            }

            /***
             * @brief set correspondences
             * @param correspondences source-target correspondences
             */
            void setInputCorrespondences(const pcl::Correspondences &correspondences)
            {
                source_match_.resize(3, correspondences.size());
                target_match_.resize(3, correspondences.size());

                for (size_t i = 0; i < correspondences.size(); ++i)
                {
                    const auto &corr = correspondences[i];
                    source_match_.col(i) = (source_cloud_->points[corr.index_query]).getVector3fMap();
                    target_match_.col(i) = (target_cloud_->points[corr.index_match]).getVector3fMap();
                }
            }

            /***
             * @brief set voxelset resolution
             * @param resolution voxelset resolution for registration evaluation
             */
            void setVoxelResolution(float resolution)
            {
                voxel_evaluator_->setResolution(resolution);
            }

            /***
             * @brief get registration result
             */
            void getResult(PCRResult &result) noexcept
            {
                result = result_;
            }

            /***
             * @brief evaluate registration result
             */
            void evaluate()
            {
                if (!result_.converged_)
                    return;

                Eigen::Isometry3f transformation = result_.transformation_;
                if (transformation.matrix().isIdentity())
                    return;

                /* get transformed source cloud */
                PointCloudPtr transformed_source_cloud(new PointCloud);
                pcl::transformPointCloud(*source_cloud_, *transformed_source_cloud, transformation);

                /* evaluate registration result via voxelset */
                this->voxel_evaluator_->convertPCLToVoxels(target_cloud_);
                const auto voxel_score = this->voxel_evaluator_->evaluateScore(*transformed_source_cloud);
                result_.inlier_fraction_ = voxel_score.inlier_fraction_;
                result_.score_ = voxel_score.rmse_;
            }

            /***
             * @brief align to get best transformation
             */
            bool align()
            {
                return static_cast<Derived *>(this)->alignImpl();
            }

        protected:
            /***
             * @brief default constructor
             */
            explicit BasePCR()
            {
                voxel_evaluator_ = std::make_unique<eststack::core::VoxelSet<PointT>>();
            }

            /***
             * @brief voxel set for registration evaluation
             */
            std::unique_ptr<eststack::core::VoxelSet<PointT>> voxel_evaluator_;

            /***
             * @brief source point cloud to be aligned
             */
            PointCloudConstPtr source_cloud_;

            /***
             * @brief target point cloud, e.g. global map
             */
            PointCloudConstPtr target_cloud_;

            /***
             * @brief source point coordinates of source-target correspondences
             */
            Eigen::Matrix3Xf source_match_;

            /***
             * @brief target point coordinates of source-target correspondences
             */
            Eigen::Matrix3Xf target_match_;

            /***
             * @brief registration result
             */
            PCRResult result_;
        };
    }
}

#endif //! SOLUTION__BASE_PCR_HPP