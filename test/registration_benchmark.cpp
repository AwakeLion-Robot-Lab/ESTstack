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

#ifndef TEST__REGISTRATION_BENCHMARK_CPP
#define TEST__REGISTRATION_BENCHMARK_CPP

// C++ standard library
#include <filesystem>
#include <algorithm>
#include <iostream>
#include <cstdint>
#include <cstdlib>
#include <numbers>
#include <fstream>
#include <iomanip>
#include <chrono>
#include <memory>
#include <string>
#include <vector>
#include <cmath>

// Eigen library
#include <Eigen/Geometry>
#include <Eigen/Core>

// PCL library
#include <pcl/point_cloud.h>
#include <pcl/point_types.h>
#include <pcl/io/pcd_io.h>

// ESTstack library
#include "eststack/solution/tcf.hpp"
#include "utils/stack_tracer.hpp"

/***
 * @brief An algorithm set focus on estimation and filtering
 * @author jinhua "siyiovo" deng
 */
namespace eststack
{
    /***
     * @brief for test
     */
    namespace test
    {
        namespace fs = std::filesystem;

        using PointT = pcl::PointXYZ;
        using PointCloud = pcl::PointCloud<PointT>;
        using PointCloudPtr = PointCloud::Ptr;

        struct RegistrationCase
        {
            std::string name;
            fs::path source_cloud;
            fs::path target_cloud;
            fs::path matches;
            fs::path gt_pose;
        };

        struct XformError
        {
            float rotation_deg{0.0f};
            float translation_m{0.0f};
        };

        struct CaseResult
        {
            std::string case_name;
            bool converged{false};
            bool success{false};
            XformError error;
            double time_ms{0.0};
        };

        float getEnvFloat(const char *name, float default_value)
        {
            const char *value = std::getenv(name);
            if (value == nullptr)
                return default_value;

            char *end = nullptr;
            const float parsed = std::strtof(value, &end);
            return end != value ? parsed : default_value;
        }

        Eigen::Isometry3f toIsometry(const Eigen::Matrix4f &matrix)
        {
            Eigen::Isometry3f T = Eigen::Isometry3f::Identity();
            T.linear() = matrix.block<3, 3>(0, 0);
            T.translation() = matrix.block<3, 1>(0, 3);
            return T;
        }

        XformError errorBetween(const Eigen::Isometry3f &estimate, const Eigen::Isometry3f &gt)
        {
            const Eigen::Matrix3f dR = gt.linear() * estimate.linear().transpose();
            const float cos_angle = std::clamp((dR.trace() - 1.0f) * 0.5f, -1.0f, 1.0f);
            return {std::acos(cos_angle) * 180.0f / std::numbers::pi_v<float>,
                    (estimate.translation() - gt.translation()).norm()};
        }

        bool loadPose(const fs::path &path, Eigen::Matrix4f &pose)
        {
            std::ifstream ifs(path);
            if (!ifs.is_open())
                return false;

            pose.setIdentity();
            for (int r = 0; r < 4; ++r)
                for (int c = 0; c < 4; ++c)
                    if (!(ifs >> pose(r, c)))
                        return false;
            return true;
        }

        bool loadMatches(const fs::path &path, Eigen::Matrix3Xf &source_match, Eigen::Matrix3Xf &target_match)
        {
            std::ifstream ifs(path);
            if (!ifs.is_open())
                return false;

            std::vector<float> values;
            values.reserve(4096 * 6);

            float value = 0.0f;
            while (ifs >> value)
                values.push_back(value);

            if (values.empty() || values.size() % 6 != 0)
                return false;

            const int row_num = static_cast<int>(values.size() / 6);
            source_match.resize(3, row_num);
            target_match.resize(3, row_num);
            for (int i = 0; i < row_num; ++i)
            {
                const float *row = values.data() + static_cast<std::size_t>(i) * 6;
                source_match.col(i) = Eigen::Vector3f(row[0], row[1], row[2]);
                target_match.col(i) = Eigen::Vector3f(row[3], row[4], row[5]);
            }
            return true;
        }

        bool loadPointCloud(const fs::path &path, PointCloudPtr &cloud)
        {
            cloud = std::make_shared<PointCloud>();
            if (pcl::io::loadPCDFile(path.string(), *cloud) != 0)
                return false;

            cloud->points.erase(std::remove_if(cloud->points.begin(), cloud->points.end(),
                                               [](const PointT &pt)
                                               {
                                                   return !std::isfinite(pt.x) || !std::isfinite(pt.y) || !std::isfinite(pt.z);
                                               }),
                                cloud->points.end());
            cloud->width = static_cast<std::uint32_t>(cloud->points.size());
            cloud->height = 1;
            cloud->is_dense = true;
            return cloud->size() >= 3;
        }

        bool addCaseIfComplete(std::vector<RegistrationCase> &cases,
                               const fs::path &root,
                               std::string name,
                               fs::path source_cloud,
                               fs::path target_cloud,
                               fs::path matches,
                               fs::path gt_pose)
        {
            RegistrationCase c{std::move(name),
                               root / std::move(source_cloud),
                               root / std::move(target_cloud),
                               root / std::move(matches),
                               root / std::move(gt_pose)};

            if (fs::exists(c.source_cloud) && fs::exists(c.target_cloud) &&
                fs::exists(c.matches) && fs::exists(c.gt_pose))
            {
                cases.push_back(std::move(c));
                return true;
            }
            return false;
        }

        std::vector<RegistrationCase> discoverCases()
        {
            std::vector<fs::path> roots;
            if (const char *env = std::getenv("ESTSTACK_REGISTRATION_DATA_ROOT"))
                roots.emplace_back(env);
            if (const char *env = std::getenv("ESTSTACK_TCF_DATA_ROOT"))
                roots.emplace_back(env);

            roots.emplace_back("test/data");

            for (const fs::path &root : roots)
            {
                if (!fs::exists(root))
                    continue;

                std::vector<RegistrationCase> cases;
                if (addCaseIfComplete(cases, root, "KITTI/09/1391-1380",
                                      "1391_v0.3.pcd",
                                      "1380_v0.3.pcd",
                                      "1391_1380_top3.match",
                                      "1391-1380.pose"))
                    return cases;

                if (addCaseIfComplete(cases, root, "KITTI/09/1391-1380",
                                      "KITTI/09/1391_v0.3.pcd",
                                      "KITTI/09/1380_v0.3.pcd",
                                      "KITTI/09/1391_1380_top3.match",
                                      "KITTI/09/1391-1380.pose"))
                    return cases;
            }

            return {};
        }

        template <typename Runner>
        std::vector<CaseResult> runBenchmark(const std::string &algorithm_name,
                                             const std::vector<RegistrationCase> &cases,
                                             Runner &&runner)
        {
            const float success_rot_deg = getEnvFloat("ESTSTACK_REGISTRATION_SUCCESS_ROT_DEG", 0.9f);
            const float success_trans_m = getEnvFloat("ESTSTACK_REGISTRATION_SUCCESS_TRANS_M", 0.5f);
            const float eval_voxel_m = getEnvFloat("ESTSTACK_REGISTRATION_EVAL_VOXEL_M", 0.75f);

            std::vector<CaseResult> results;
            results.reserve(cases.size());

            std::cout << "\n[Registration benchmark: " << algorithm_name << "]\n"
                      << "  success threshold = " << success_rot_deg << " deg, "
                      << success_trans_m << " m\n";

            for (const RegistrationCase &c : cases)
            {
                PointCloudPtr source_cloud;
                PointCloudPtr target_cloud;
                Eigen::Matrix3Xf source_match;
                Eigen::Matrix3Xf target_match;
                Eigen::Matrix4f gt_matrix;

                if (!loadPointCloud(c.source_cloud, source_cloud))
                {
                    ADD_FAILURE() << "Failed to load source cloud: " << c.source_cloud.string();
                    continue;
                }
                if (!loadPointCloud(c.target_cloud, target_cloud))
                {
                    ADD_FAILURE() << "Failed to load target cloud: " << c.target_cloud.string();
                    continue;
                }
                if (!loadMatches(c.matches, source_match, target_match))
                {
                    ADD_FAILURE() << "Failed to load matches: " << c.matches.string();
                    continue;
                }
                if (!loadPose(c.gt_pose, gt_matrix))
                {
                    ADD_FAILURE() << "Failed to load ground-truth pose: " << c.gt_pose.string();
                    continue;
                }

                const Eigen::Isometry3f gt = toIsometry(gt_matrix);

                const auto start = std::chrono::steady_clock::now();
                const auto [converged, estimate] = runner(source_cloud, target_cloud, source_match, target_match, eval_voxel_m);
                const auto end = std::chrono::steady_clock::now();

                CaseResult result;
                result.case_name = c.name;
                result.converged = converged;
                result.error = errorBetween(estimate, gt);
                result.time_ms = std::chrono::duration<double, std::milli>(end - start).count();
                result.success = converged &&
                                 result.error.rotation_deg <= success_rot_deg &&
                                 result.error.translation_m <= success_trans_m;
                results.push_back(result);

                std::cout << "  " << std::left << std::setw(24) << c.name
                          << " success = " << (result.success ? "yes" : "no ")
                          << ", RE = " << result.error.rotation_deg << " deg"
                          << ", TE = " << result.error.translation_m << " m"
                          << ", time = " << result.time_ms << " ms\n";
            }

            return results;
        }

        void printSummary(const std::string &algorithm_name, const std::vector<CaseResult> &results)
        {
            if (results.empty())
            {
                ADD_FAILURE() << "No valid registration benchmark cases were executed.";
                return;
            }

            int success_count = 0;
            double rot_sum = 0.0;
            double trans_sum = 0.0;
            double time_sum = 0.0;

            for (const CaseResult &result : results)
            {
                success_count += result.success ? 1 : 0;
                rot_sum += result.error.rotation_deg;
                trans_sum += result.error.translation_m;
                time_sum += result.time_ms;
            }

            const double case_count = static_cast<double>(results.size());
            const double registration_recall = case_count > 0.0 ? 100.0 * success_count / case_count : 0.0;
            std::cout << "\n[Registration benchmark summary: " << algorithm_name << "]\n"
                      << "  cases = " << results.size() << "\n"
                      << "  registration recall = " << registration_recall << " %\n"
                      << "  mean rotation error (arccos) = " << rot_sum / case_count << " deg\n"
                      << "  mean translation error = " << trans_sum / case_count << " m\n"
                      << "  mean time = " << time_sum / case_count << " ms\n";
        }

        TEST(RegistrationBenchmark, TCF)
        {
            const std::vector<RegistrationCase> cases = discoverCases();
            if (cases.empty())
            {
                GTEST_SKIP() << "No registration benchmark data found. "
                             << "Set ESTSTACK_REGISTRATION_DATA_ROOT or copy TCF data under test/data/registration.";
            }

            const auto results = runBenchmark(
                "TCF", cases,
                [](const PointCloudPtr &source_cloud,
                   const PointCloudPtr &target_cloud,
                   const Eigen::Matrix3Xf &source_match,
                   const Eigen::Matrix3Xf &target_match,
                   float eval_voxel_m)
                {
                    solution::TCF<PointT> tcf;
                    tcf.setInputSource(source_cloud);
                    tcf.setInputTarget(target_cloud);
                    tcf.setInputCorrespondences(source_match, target_match);
                    tcf.setVoxelResolution(eval_voxel_m);

                    const bool converged = tcf.align();
                    return std::make_pair(converged, tcf.getResult().transformation_);
                });

            printSummary("TCF", results);
        }
    } // namespace test
} // namespace eststack

#endif //! TEST__REGISTRATION_BENCHMARK_CPP
