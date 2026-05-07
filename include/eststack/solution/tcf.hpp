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

#ifndef SOLUTION__TCF_HPP
#define SOLUTION__TCF_HPP

// C++ standard library
#include <algorithm>
#include <memory>
#include <random>
#include <vector>
#include <cmath>

// Eigen library
#include <Eigen/Geometry>
#include <Eigen/Core>

// PCL library
#include <pcl/filters/voxel_grid.h>
#include <pcl/point_cloud.h>
#include <pcl/point_types.h>

// nanoflann library
#include <nanoflann/nanoflann.hpp>

// ESTstack library
#include "eststack/model/sac/tcf_sac_model.hpp"
#include "eststack/solution/base_pcr.hpp"
#include "eststack/core/ransac.hpp"
#include "eststack/core/irls.hpp"

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
         * @brief Two-stage Consensus Filtering for real-time 3D registration
         * @tparam PointT point type
         * @details refer to [11]
         */
        template <typename PointT>
        class TCF final : public BasePCR<TCF<PointT>, PointT>
        {
        public:
            EIGEN_MAKE_ALIGNED_OPERATOR_NEW

            using Base = BasePCR<TCF<PointT>, PointT>;

            using Base::result_;
            using Base::source_cloud_;
            using Base::source_match_;
            using Base::target_cloud_;
            using Base::target_match_;
            using Base::voxel_evaluator_;

            using typename Base::PointCloud;
            using typename Base::PointCloudConstPtr;
            using typename Base::PointCloudPtr;

            using Ptr = std::shared_ptr<TCF<PointT>>;
            using ConstPtr = std::shared_ptr<const TCF<PointT>>;

            /***
             * @brief default constructor
             * @param leaf_size leaf size for voxel downsampling
             */
            explicit TCF(float leaf_size = 0.5f) : leaf_size_(leaf_size) {}

            /***
             * @brief set leaf size
             * @param leaf_size leaf size for voxel downsampling
             */
            void setLeafSize(float leaf_size) noexcept
            {
                leaf_size_ = leaf_size;
            }

            /***
             * @brief align to get best transformation
             */
            bool alignImpl()
            {
                /* reset result */
                this->result_ = eststack::solution::PCRResult();

                /* get inlier threshold via point cloud resolution */
                const float source_resolution = getPCLResolution(this->source_cloud_, this->leaf_size_);
                const float target_resolution = getPCLResolution(this->target_cloud_, this->leaf_size_);
                const float inlier_threshold = std::max(source_resolution, target_resolution);

                /* stage 1: length consistency coarse ransac */
                auto OneStageModel = std::make_shared<model::OnePtSACModel>();
                OneStageModel->setInputSource(this->source_match_);
                OneStageModel->setInputTarget(this->target_match_);
                OneStageModel->setInlierThreshold(inlier_threshold);

                auto OneStageRANSAC = std::make_shared<core::RANSAC<model::OnePtSACModel>>();
                OneStageRANSAC->setSACModel(OneStageModel);
                OneStageRANSAC->optimize();

                const auto one_stage_inliers = OneStageRANSAC->getBestInliers();
                if (one_stage_inliers.size() < 3)
                {
                    this->result_.converged_ = false;
                    return false;
                }

                Eigen::Matrix3Xf one_stage_source_inliers = this->source_match_(Eigen::placeholders::all, one_stage_inliers);
                Eigen::Matrix3Xf one_stage_target_inliers = this->target_match_(Eigen::placeholders::all, one_stage_inliers);

                /* stage 2: angular discrepancy mid ransac */
                auto TwoStageModel = std::make_shared<model::TwoPtSACModel>();
                TwoStageModel->setInputSource(one_stage_source_inliers);
                TwoStageModel->setInputTarget(one_stage_target_inliers);
                TwoStageModel->setInlierThreshold(inlier_threshold);

                auto TwoStageRANSAC = std::make_shared<core::RANSAC<model::TwoPtSACModel>>();
                TwoStageRANSAC->setSACModel(TwoStageModel);
                TwoStageRANSAC->optimize();

                const auto two_stage_inliers = TwoStageRANSAC->getBestInliers();
                if (two_stage_inliers.size() < 3)
                {
                    this->result_.converged_ = false;
                    return false;
                }

                Eigen::Matrix3Xf two_stage_source_inliers = one_stage_source_inliers(Eigen::placeholders::all, two_stage_inliers);
                Eigen::Matrix3Xf two_stage_target_inliers = one_stage_target_inliers(Eigen::placeholders::all, two_stage_inliers);

                /* stage 3: kabsch ransac for coarse rigid transformation */
                auto RigidMotionModel = std::make_shared<model::ThreePtSACModel>();
                RigidMotionModel->setInputSource(two_stage_source_inliers);
                RigidMotionModel->setInputTarget(two_stage_target_inliers);
                RigidMotionModel->setInlierThreshold(inlier_threshold);

                auto RigidMotionRANSAC = std::make_shared<core::RANSAC<model::ThreePtSACModel>>();
                RigidMotionRANSAC->setSACModel(RigidMotionModel);
                RigidMotionRANSAC->optimize();

                const auto rigid_motion_inliers = RigidMotionRANSAC->getBestInliers();
                if (rigid_motion_inliers.size() < 3)
                {
                    this->result_.converged_ = false;
                    return false;
                }

                Eigen::Matrix3Xf rigid_motion_source_inliers = two_stage_source_inliers(Eigen::placeholders::all, rigid_motion_inliers);
                Eigen::Matrix3Xf rigid_motion_target_inliers = two_stage_target_inliers(Eigen::placeholders::all, rigid_motion_inliers);

                /* get better consensus set */
                const auto coarse_trans_coeffs = RigidMotionRANSAC->getBestFitness();
                Eigen::Quaternionf q(coarse_trans_coeffs(0), coarse_trans_coeffs(1), coarse_trans_coeffs(2), coarse_trans_coeffs(3));
                Eigen::Vector3f t(coarse_trans_coeffs(4), coarse_trans_coeffs(5), coarse_trans_coeffs(6));
                Eigen::Matrix3f R = q.toRotationMatrix();

                /* projection refinement using coarse transformation */
                Eigen::Matrix3Xf coarse_transformed_source = (R * rigid_motion_source_inliers).colwise() + t;
                Eigen::VectorXf error_distance = (coarse_transformed_source - rigid_motion_target_inliers).colwise().norm().transpose();
                Eigen::VectorXi coarse_refined_flags = (error_distance.array() < inlier_threshold).cast<int>();
                Eigen::VectorXi coarse_refined_indices = eststack::getNonZeroColumnIndicesFromVector(coarse_refined_flags);

                if (coarse_refined_indices.size() < 3)
                {
                    this->result_.converged_ = false;
                    return false;
                }

                Eigen::Matrix3Xf coarse_refined_source_inliers = rigid_motion_source_inliers(Eigen::placeholders::all, coarse_refined_indices);
                Eigen::Matrix3Xf coarse_refined_target_inliers = rigid_motion_target_inliers(Eigen::placeholders::all, coarse_refined_indices);

                /* stage 4: rigid IRLS for fine rigid transformation */
                Eigen::Isometry3f init_guess = Eigen::Isometry3f::Identity();
                init_guess.linear() = R;
                init_guess.translation() = t;
                auto sa_cauchy = std::make_shared<core::CauchyLoss>();
                auto IRLS = std::make_shared<core::RigidIRLS<core::CauchyLoss>>();
                IRLS->setInputSource(coarse_refined_source_inliers);
                IRLS->setInputTarget(coarse_refined_target_inliers);
                IRLS->setInitGuess(init_guess);
                IRLS->setLossFunction(sa_cauchy);
                IRLS->optimize();

                /* return some parts of result */
                this->result_.converged_ = IRLS->hasConverged();
                this->result_.transformation_ = IRLS->getFinalTransformation();
                return this->result_.converged_;
            }

        private:
            /***
             * @brief leaf size for voxel downsampling
             */
            float leaf_size_;

            /***
             * @brief `pcl::PointCloud<PointT>` adaptor for `nanoflann`
             */
            struct PointCloudAdaptor
            {
                const PointCloud &pts;
                explicit PointCloudAdaptor(const PointCloud &p) : pts(p) {}

                inline std::size_t kdtree_get_point_count() const { return pts.size(); }

                inline float kdtree_get_pt(const std::size_t idx, const std::size_t dim) const
                {
                    if (dim == 0)
                        return pts.points[idx].x;
                    if (dim == 1)
                        return pts.points[idx].y;
                    return pts.points[idx].z;
                }

                template <class BBOX>
                bool kdtree_get_bbox(BBOX &) const { return false; }
            };

            /***
             * @brief get resolution of point cloud
             * @param cloud input point cloud
             * @param leaf_size leaf size for get resolution
             */
            float getPCLResolution(const PointCloudConstPtr &cloud, float leaf_size)
            {
                if (!cloud || cloud->size() < 2)
                    return 0.1f;

                /* if cloud size is large, down sampling first */
                PointCloudConstPtr used_cloud = cloud;
                if (cloud->size() > 10000)
                {
                    auto vg = std::make_shared<pcl::VoxelGrid<PointT>>();
                    PointCloudPtr cloud_filtered(new PointCloud);

                    vg->setInputCloud(cloud);
                    vg->setLeafSize(leaf_size, leaf_size, leaf_size);
                    vg->filter(*cloud_filtered);
                    used_cloud = cloud_filtered;
                }

                /* build `nanoflann` kd-tree */
                PointCloudAdaptor pcl_adaptor(*used_cloud);
                using KDTree = nanoflann::KDTreeSingleIndexAdaptor<nanoflann::L2_Simple_Adaptor<float, PointCloudAdaptor>, PointCloudAdaptor, 3>;
                KDTree kdtree(3, pcl_adaptor, nanoflann::KDTreeSingleIndexAdaptorParams(10));
                kdtree.buildIndex();

                /* just search for half of point cloud num */
                const int search_num = std::max(10, static_cast<int>(used_cloud->size() / 2));
                const int pt_num = static_cast<int>(used_cloud->size());
                std::vector<float> dists;
                dists.reserve(search_num);

                /* random sampling and query */
                std::random_device rd;
                std::mt19937 gen(rd());
                std::uniform_int_distribution<> dis(0, pt_num - 1);
                for (int i = 0; i < search_num; ++i)
                {
                    const int idx = dis(gen);
                    const float query_pt[3] = {used_cloud->points[idx].x,
                                               used_cloud->points[idx].y,
                                               used_cloud->points[idx].z};

                    /* query two nearest neighbors to exclude self point */
                    std::size_t indices[2];
                    float sqrt_dists[2];
                    nanoflann::KNNResultSet<float> resultSet(2);
                    resultSet.init(indices, sqrt_dists);

                    /* skip self point and take the true nearest neighbor */
                    kdtree.findNeighbors(resultSet, query_pt, nanoflann::SearchParameters());
                    float sqrt_dist = (indices[0] == static_cast<std::size_t>(idx)) ? sqrt_dists[0] : sqrt_dists[1];
                    dists.push_back(sqrt_dist);
                }

                /* get median distance */
                const std::size_t mid_pos = dists.size() / 2;
                std::nth_element(dists.begin(), dists.begin() + mid_pos, dists.end());
                return std::sqrt(dists[mid_pos]);
            }
        };
    } // namespace solution
} // namespace eststack

#endif //! SOLUTION__TCF_HPP
