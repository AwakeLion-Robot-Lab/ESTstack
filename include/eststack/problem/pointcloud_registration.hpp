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

#ifndef PROBLEM__POINTCLOUD_REGISTRATION_HPP
#define PROBLEM__POINTCLOUD_REGISTRATION_HPP

// C++ standard library
#include <type_traits>
#include <utility>

// Eigen library
#include <Eigen/Core>
#include <Eigen/Geometry>

// PCL library
#include <pcl/correspondence.h>

// ESTstack Library
#include <eststack/problem/base_problem.hpp>
#include <eststack/concepts.hpp>

/***
 * @brief An algorithm set focus on estimation and filtering
 * @author jinhua "siyiovo" deng
 */
namespace eststack
{
    /***
     * @brief problem definitions that orchestrate filters, models and noise
     */
    namespace problem
    {
        /***
         * @brief pointcloud registration problem
         * @tparam StateT state type
         * @tparam SolutionT solution type
         */
        template <typename StateT, PCR SolutionT>
        class PointCloudRegistration : public BaseProblem<PointCloudRegistration<StateT, SolutionT>, StateT, SolutionT>
        {
        public:
            using Base = BaseProblem<PointCloudRegistration<StateT, SolutionT>, StateT, SolutionT>;
            using State = typename Base::State;
            using Solution = typename Base::Solution;
            using Result = typename Solution::Result;
            using PointCloud = typename Solution::PointCloud;
            using PointCloudPtr = typename Solution::PointCloudPtr;
            using PointCloudConstPtr = typename Solution::PointCloudConstPtr;

            static_assert(std::is_assignable_v<State &, Eigen::Isometry3f>,
                          "PointCloudRegistration state must be assignable from Eigen::Isometry3f");

            /***
             * @brief default constructor
             */
            PointCloudRegistration()
            {
                this->state_ = Eigen::Isometry3f::Identity();
            }

            /***
             * @brief construct with a registration solution
             * @param sol solution
             */
            explicit PointCloudRegistration(Solution &&sol) : PointCloudRegistration()
            {
                this->setSolution(std::forward<Solution>(sol));
            }

            /***
             * @brief set source cloud
             * @param source source point cloud
             */
            void setInputSource(const PointCloudConstPtr &source)
            {
                source_cloud_ = source;
                if (this->solution_)
                    syncSolutionInputs(*this->solution_);
            }

            /***
             * @brief set target cloud
             * @param target target point cloud
             */
            void setInputTarget(const PointCloudConstPtr &target)
            {
                target_cloud_ = target;
                if (this->solution_)
                    syncSolutionInputs(*this->solution_);
            }

            /***
             * @brief set source-target correspondences
             * @param source_match source point coordinates of correspondences
             * @param target_match target point coordinates of correspondences
             */
            void setInputCorrespondences(const Eigen::Matrix3Xf &source_match, const Eigen::Matrix3Xf &target_match)
            {
                source_match_ = source_match;
                target_match_ = target_match;
                correspondence_input_ = CorrespondenceInput::Matrix;

                if (this->solution_)
                    syncSolutionInputs(*this->solution_);
            }

            /***
             * @brief set source-target correspondences by point indices
             * @param correspondences source-target correspondences
             */
            void setInputCorrespondences(const pcl::Correspondences &correspondences)
                requires requires(Solution &solution, const pcl::Correspondences &corr) {
                    { solution.setInputCorrespondences(corr) } -> std::same_as<void>;
                }
            {
                correspondences_ = correspondences;
                correspondence_input_ = CorrespondenceInput::PCL;

                if (this->solution_)
                    syncSolutionInputs(*this->solution_);
            }

            /***
             * @brief set voxelset resolution for registration evaluation
             * @param resolution voxelset resolution
             */
            void setVoxelResolution(float resolution)
            {
                voxel_resolution_ = resolution;
                has_voxel_resolution_ = true;

                if (this->solution_)
                    syncSolutionInputs(*this->solution_);
            }

            /***
             * @brief get registration result
             */
            const Result &getResult() const noexcept
            {
                return result_;
            }

            /***
             * @brief set specific solution implementation
             * @param solution solution pointer
             */
            void setSolutionImpl(Solution *solution)
            {
                if (solution)
                    syncSolutionInputs(*solution);
            }

            /***
             * @brief run pointcloud registration
             */
            bool runImpl()
            {
                this->ok_ = false;

                if (!hasValidInputs() || !this->solution_)
                    return false;

                syncSolutionInputs(*this->solution_);

                const bool aligned = this->solution_->align();
                this->solution_->evaluate();
                result_ = this->solution_->getResult();

                this->ok_ = aligned && result_.converged_;
                if (this->ok_)
                    this->state_ = result_.transformation_;

                return this->ok_;
            }

        private:
            /***
             * @brief correspondence input type
             */
            enum class CorrespondenceInput
            {
                None,
                Matrix,
                PCL
            };

            /***
             * @brief synchronize cached problem inputs into solution
             */
            void syncSolutionInputs(Solution &solution)
            {
                if (source_cloud_)
                    solution.setInputSource(source_cloud_);

                if (target_cloud_)
                    solution.setInputTarget(target_cloud_);

                if (correspondence_input_ == CorrespondenceInput::Matrix)
                    solution.setInputCorrespondences(source_match_, target_match_);
                else if constexpr (requires(Solution &sol, const pcl::Correspondences &corr) {
                                       { sol.setInputCorrespondences(corr) } -> std::same_as<void>;
                                   })
                {
                    if (correspondence_input_ == CorrespondenceInput::PCL && source_cloud_ && target_cloud_)
                        solution.setInputCorrespondences(correspondences_);
                }

                if (has_voxel_resolution_)
                    solution.setVoxelResolution(voxel_resolution_);
            }

            /***
             * @brief check whether required registration inputs are ready
             */
            bool hasValidInputs() const noexcept
            {
                if (!source_cloud_ || !target_cloud_)
                    return false;

                if (source_cloud_->empty() || target_cloud_->empty())
                    return false;

                if (correspondence_input_ == CorrespondenceInput::Matrix)
                    return source_match_.cols() == target_match_.cols();

                return true;
            }

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
             * @brief source-target correspondences by point indices
             */
            pcl::Correspondences correspondences_;

            /***
             * @brief correspondence input type
             */
            CorrespondenceInput correspondence_input_{CorrespondenceInput::None};

            /***
             * @brief voxelset resolution for registration evaluation
             */
            float voxel_resolution_{0.0f};

            /***
             * @brief whether voxelset resolution has been set
             */
            bool has_voxel_resolution_{false};

            /***
             * @brief latest registration result
             */
            Result result_;
        };
    } // namespace problem
} // namespace eststack

#endif //! PROBLEM__POINTCLOUD_REGISTRATION_HPP
