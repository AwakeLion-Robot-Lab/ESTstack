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

#ifndef CORE__VOXELSET_HPP
#define CORE__VOXELSET_HPP

// C++ standard library
#include <unordered_set>
#include <vector>
#include <limits>

// Eigen library
#include <Eigen/Dense>

// PCL library
#include <pcl/point_cloud.h>
#include <pcl/point_types.h>

/***
 * @brief An algorithm set focus on state estimation
 * @author jinhua "siyiovo" deng
 */
namespace eststack
{
    /***
     * @brief core algorithms for estimation
     */
    namespace core
    {
        /***
         * @brief score of a voxel
         */
        struct VoxelScore
        {
            /***
             * @brief inlier fraction
             */
            float inlier_fraction_;

            /***
             * @brief root mean squared error distance from inliers to voxel center
             */
            float rmse_;
        };

        /***
         * @brief hash function for voxel coordinates
         * @details refer to [22], specific url: https://github.com/koide3/hdl_global_localization/blob/master/include/hdl_global_localization/ransac/voxelset.hpp
         */
        struct HashVector3i
        {
            std::size_t operator()(const Eigen::Vector3i &x) const noexcept
            {
                size_t seed = 0;
                hash_combine(seed, x[0]);
                hash_combine(seed, x[1]);
                hash_combine(seed, x[2]);

                return seed;
            }

        private:
            /***
             * @brief combine hash values
             * @details since we just use C++ std library, we need to redefine `boost::hash_combine`,
             *          refer to https://www.boost.org/doc/libs/1_64_0/boost/functional/hash/hash.hpp
             */
            template <typename T>
            static void hash_combine(std::size_t &seed, T const &value) noexcept
            {
                std::hash<T> hasher;
                seed ^= hasher(value) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
            }
        };

        /***
         * @brief voxelize pointcloud as a hash set for fast matching and PCR score evaluation
         * @details refer to [22], specific url: https://github.com/koide3/hdl_global_localization/blob/master/include/hdl_global_localization/ransac/voxelset.hpp
         */
        template <typename PointT>
        class VoxelSet
        {
        public:
            using PointCloud = pcl::PointCloud<PointT>;
            using PointCloudPtr = typename PointCloud::Ptr;
            using PointCloudConstPtr = typename PointCloud::ConstPtr;

            /***
             * @brief default constructor
             * @param resolution voxel resolution
             */
            explicit VoxelSet(float resolution = 1.0f) : resolution_(resolution) {}

            /***
             * @brief set voxel resolution
             * @param resolution voxel resolution
             */
            void setResolution(float resolution) noexcept
            {
                resolution_ = resolution;
                voxels_.clear();
            }

            /***
             * @brief convert point cloud to voxels
             * @param cloud input point cloud
             */
            void convertPCLToVoxels(const PointCloudConstPtr &cloud)
            {
                std::vector<Eigen::Vector3i, Eigen::aligned_allocator<Eigen::Vector3i>> offsets;
                offsets.reserve(27);
                voxels_.clear();

                /* cube offset */
                for (int i = -1; i <= 1; i++)
                {
                    for (int j = -1; j <= 1; j++)
                    {
                        for (int k = -1; k <= 1; k++)
                        {
                            offsets.push_back(Eigen::Vector3i(i, j, k));
                        }
                    }
                }

                /* calculate voxel coordinate */
                for (const auto &pt : *cloud)
                {
                    Eigen::Vector4f pcl_coord = getPointCoord(pt);
                    Eigen::Vector3i voxel_coord = getVoxelCoord(pcl_coord);
                    for (const auto &offset : offsets)
                    {
                        Eigen::Vector4f center = getVoxelCenter(voxel_coord + offset);
                        if ((center - pcl_coord).squaredNorm() < resolution_ * resolution_)
                            voxels_.insert(voxel_coord + offset);
                    }
                }
            }

            /***
             * @brief evaluate matching score of a point cloud with the voxel set
             * @param cloud input point cloud
             * @return voxel score of the point cloud
             */
            VoxelScore evaluateScore(const PointCloud &cloud) const
            {
                int inliers_num = 0;
                double errors = 0.0;

                /* count inliers and accumulate errors */
                for (const auto &pt : cloud)
                {
                    Eigen::Vector4f pcl_coord = getPointCoord(pt);
                    Eigen::Vector3i voxel_coord = getVoxelCoord(pcl_coord);
                    if (voxels_.find(voxel_coord) != voxels_.end())
                    {
                        Eigen::Vector4f center = getVoxelCenter(voxel_coord);
                        ++inliers_num;
                        errors += (pcl_coord - center).squaredNorm();
                    }
                }

                VoxelScore voxel_score;
                voxel_score.inlier_fraction_ = static_cast<float>(inliers_num) / cloud.size();
                voxel_score.rmse_ = (inliers_num > 0) ? static_cast<float>(std::sqrt(errors / static_cast<double>(inliers_num))) : std::numeric_limits<float>::max();
                return voxel_score;
            }

        private:
            using Voxels = std::unordered_set<Eigen::Vector3i, HashVector3i, std::equal_to<Eigen::Vector3i>, Eigen::aligned_allocator<Eigen::Vector3i>>;

            /***
             * @brief voxel coordinates of points
             */
            Voxels voxels_;

            /***
             * @brief voxel resolution
             */
            float resolution_;

            /***
             * @brief get point coordinate in homogeneous form
             * @param pt input point
             * @details no matter what the point type is, we just get its x, y, z and set w to 1.0f for homogeneous coordinate
             */
            inline Eigen::Vector4f getPointCoord(const PointT &pt)
            {
                return Eigen::Vector4f(pt.x, pt.y, pt.z, 1.0f);
            }

            /***
             * @brief get voxel coordinate of a point
             * @param x point coordinate
             * @return voxel coordinate
             */
            inline Eigen::Vector3i getVoxelCoord(const Eigen::Vector4f &x) const
            {
                return (x.array() / resolution_ - 0.5f).floor().cast<int>().head<3>();
            }

            /***
             * @brief get voxel center coordinate of a voxel
             * @param coord voxel coordinate
             * @return voxel center coordinate
             */
            inline Eigen::Vector4f getVoxelCenter(const Eigen::Vector3i &coord) const
            {
                Eigen::Vector3f origin = (coord.template cast<float>().array() + 0.5f) * resolution_;
                Eigen::Vector3f offset = Eigen::Vector3f::Ones() * resolution_ / 2.0f;
                return (Eigen::Vector4f() << origin + offset, 1.0f).finished();
            }
        };

    } // namespace core
} // namespace eststack

#endif //! CORE__VOXELSET_HPP