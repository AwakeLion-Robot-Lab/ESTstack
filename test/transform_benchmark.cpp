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

#ifndef TEST__TRANSFORM_BENCHMARK_CPP
#define TEST__TRANSFORM_BENCHMARK_CPP

// C++ standard library
#include <iomanip>
#include <utility>
#include <memory>
#include <random>
#include <vector>
#include <array>
#include <cmath>

// Eigen library
#include <Eigen/Geometry>
#include <Eigen/Core>

// PCL library
#include <pcl/point_cloud.h>
#include <pcl/point_types.h>

// ESTstack library
#include "eststack/core/dcm_to_quat.hpp"
#include "eststack/eigen_traits.hpp"
#include "eststack/core/kabsch.hpp"
#include "eststack/core/irls.hpp"
#include "utils/stack_tracer.hpp"
#include "utils/benchmarker.hpp"

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
        using eststack::core::bar_itzhack;
        using eststack::core::CauchyLoss;
        using eststack::core::fromDCM;
        using eststack::core::HuberLoss;
        using eststack::core::kabsch;
        using eststack::core::RigidIRLS;
        using eststack::core::sarabandi;
        using eststack::core::weightedKabsch;

        /* ---- tolerances ---- */
        constexpr float kQuatTolTight = 5e-5f;  /* ordinary rotations, float precision */
        constexpr float kQuatTolLoose = 1e-3f;  /* 180-deg singularity / dirty DCM */
        constexpr float kQuatNormTol = 1e-5f;   /* |q| should be 1 */
        constexpr float kXformTolRotDeg = 0.1f; /* rigid-fit rotation tolerance */
        constexpr float kXformTolTrans = 1e-2f; /* rigid-fit translation tolerance */

        /* ---- randomness ---- */
        Eigen::Quaternionf randomQuat(std::mt19937 &gen)
        {
            std::normal_distribution<float> d(0.0f, 1.0f);
            Eigen::Quaternionf q(d(gen), d(gen), d(gen), d(gen));
            return q.normalized();
        }

        Eigen::Matrix3f randomDCM(std::mt19937 &gen)
        {
            return randomQuat(gen).toRotationMatrix();
        }

        Eigen::Isometry3f randomTransform(std::mt19937 &gen,
                                          float angle_max = static_cast<float>(M_PI) / 4.0f,
                                          float trans_max = 5.0f)
        {
            std::uniform_real_distribution<float> ang(-angle_max, angle_max);
            std::normal_distribution<float> axd(0.0f, 1.0f);
            std::uniform_real_distribution<float> tr(-trans_max, trans_max);

            Eigen::Vector3f axis(axd(gen), axd(gen), axd(gen));
            axis.normalize();
            Eigen::AngleAxisf aa(ang(gen), axis);

            Eigen::Isometry3f T = Eigen::Isometry3f::Identity();
            T.linear() = aa.toRotationMatrix();
            T.translation() = Eigen::Vector3f(tr(gen), tr(gen), tr(gen));
            return T;
        }

        Eigen::Matrix3Xf randomCloud(std::mt19937 &gen, int N, float extent = 10.0f)
        {
            std::uniform_real_distribution<float> d(-extent, extent);
            Eigen::Matrix3Xf M(3, N);
            for (int c = 0; c < N; ++c)
                M.col(c) = Eigen::Vector3f(d(gen), d(gen), d(gen));
            return M;
        }

        Eigen::Matrix3Xf transformCloud(const Eigen::Matrix3Xf &X, const Eigen::Isometry3f &T)
        {
            return (T.linear() * X).colwise() + T.translation();
        }

        Eigen::Matrix3Xf addGaussianNoise(std::mt19937 &gen, const Eigen::Matrix3Xf &X, float sigma)
        {
            std::normal_distribution<float> n(0.0f, sigma);
            Eigen::Matrix3Xf Y = X;
            for (int c = 0; c < Y.cols(); ++c)
                for (int r = 0; r < 3; ++r)
                    Y(r, c) += n(gen);
            return Y;
        }

        void injectOutliers(std::mt19937 &gen, Eigen::Matrix3Xf &cloud, int N_out, float mag = 50.0f)
        {
            std::uniform_real_distribution<float> big(-mag, mag);
            for (int i = 0; i < N_out && i < cloud.cols(); ++i)
                cloud.col(i) += Eigen::Vector3f(big(gen), big(gen), big(gen));
        }

        pcl::PointCloud<pcl::PointXYZ>::Ptr toPcl(const Eigen::Matrix3Xf &m)
        {
            auto pc = std::make_shared<pcl::PointCloud<pcl::PointXYZ>>();
            pc->resize(m.cols());
            for (int i = 0; i < m.cols(); ++i)
            {
                pc->points[i].x = m(0, i);
                pc->points[i].y = m(1, i);
                pc->points[i].z = m(2, i);
            }
            return pc;
        }

        /* ---- metrics ---- */
        float angularDistance(const Eigen::Quaternionf &q1, const Eigen::Quaternionf &q2)
        {
            return q1.normalized().angularDistance(q2.normalized());
        }

        struct XformError
        {
            float rotation_deg;
            float translation;
        };

        XformError errorBetween(const Eigen::Isometry3f &a, const Eigen::Isometry3f &b)
        {
            const Eigen::Matrix3f dR = a.linear().transpose() * b.linear();
            const Eigen::AngleAxisf aa(dR);
            return {std::abs(aa.angle()) * 180.0f / static_cast<float>(M_PI),
                    (a.translation() - b.translation()).norm()};
        }

        /* ================================================================
         *                   Section 1: Kabsch (rigid fit)
         * ================================================================ */

        TEST(Kabsch, IdentityRecovery)
        {
            std::mt19937 gen(11);
            const Eigen::Matrix3Xf src = randomCloud(gen, 200);

            const Eigen::Isometry3f T = kabsch(src, src);
            const auto err = errorBetween(T, Eigen::Isometry3f::Identity());
            EXPECT_LT(err.rotation_deg, kXformTolRotDeg);
            EXPECT_LT(err.translation, kXformTolTrans);
        }

        TEST(Kabsch, CleanRigidTransformEigen)
        {
            std::mt19937 gen(22);
            const Eigen::Matrix3Xf src = randomCloud(gen, 300);
            const Eigen::Isometry3f T_true = randomTransform(gen);
            const Eigen::Matrix3Xf dst = transformCloud(src, T_true);

            const Eigen::Isometry3f T_est = kabsch(src, dst);
            const auto err = errorBetween(T_est, T_true);
            EXPECT_LT(err.rotation_deg, kXformTolRotDeg);
            EXPECT_LT(err.translation, kXformTolTrans);
        }

        TEST(Kabsch, CleanRigidTransformPcl)
        {
            std::mt19937 gen(33);
            const Eigen::Matrix3Xf src = randomCloud(gen, 300);
            const Eigen::Isometry3f T_true = randomTransform(gen);
            const Eigen::Matrix3Xf dst = transformCloud(src, T_true);

            const Eigen::Isometry3f T_est = kabsch(toPcl(src), toPcl(dst));
            const auto err = errorBetween(T_est, T_true);
            EXPECT_LT(err.rotation_deg, kXformTolRotDeg);
            EXPECT_LT(err.translation, kXformTolTrans);
        }

        TEST(Kabsch, LargeCloudParallelPath)
        {
            /* N >= 2048 triggers the TBB combinable parallel reduction. */
            std::mt19937 gen(44);
            const Eigen::Matrix3Xf src = randomCloud(gen, 4096);
            const Eigen::Isometry3f T_true = randomTransform(gen);
            const Eigen::Matrix3Xf dst = transformCloud(src, T_true);

            const Eigen::Isometry3f T_est = kabsch(src, dst);
            const auto err = errorBetween(T_est, T_true);
            EXPECT_LT(err.rotation_deg, kXformTolRotDeg);
            EXPECT_LT(err.translation, kXformTolTrans);
            /* Rotation must stay in SO(3). */
            EXPECT_NEAR(T_est.linear().determinant(), 1.0f, 1e-4f);
        }

        TEST(Kabsch, SO3ConstraintUnderReflection)
        {
            /* Construct a target related by a reflection — Kabsch must still return
             * a proper rotation (det == +1) rather than a reflection (det == -1). */
            std::mt19937 gen(55);
            const Eigen::Matrix3Xf src = randomCloud(gen, 100);
            Eigen::Matrix3Xf dst = src;
            dst.row(2) *= -1.0f; /* mirror through the xy-plane */

            const Eigen::Isometry3f T_est = kabsch(src, dst);
            EXPECT_NEAR(T_est.linear().determinant(), 1.0f, 1e-4f);
        }

        TEST(Kabsch, EigenEmptyReturnsIdentity)
        {
            const Eigen::Isometry3f T = kabsch(Eigen::Matrix3Xf(3, 0), Eigen::Matrix3Xf(3, 0));
            EXPECT_TRUE(T.matrix().isApprox(Eigen::Matrix4f::Identity()));
        }

        TEST(Kabsch, EigenSizeMismatchReturnsIdentity)
        {
            std::mt19937 gen(66);
            const Eigen::Isometry3f T = kabsch(randomCloud(gen, 10), randomCloud(gen, 20));
            EXPECT_TRUE(T.matrix().isApprox(Eigen::Matrix4f::Identity()));
        }

        TEST(Kabsch, PclNullPointerThrows)
        {
            pcl::PointCloud<pcl::PointXYZ>::Ptr null_pc;
            auto valid_pc = std::make_shared<pcl::PointCloud<pcl::PointXYZ>>();
            valid_pc->resize(1);
            EXPECT_THROW(kabsch(null_pc, valid_pc), std::invalid_argument);
            EXPECT_THROW(kabsch(valid_pc, null_pc), std::invalid_argument);
        }

        TEST(WeightedKabsch, EqualWeightsMatchesKabsch)
        {
            std::mt19937 gen(77);
            const Eigen::Matrix3Xf src = randomCloud(gen, 200);
            const Eigen::Isometry3f T_true = randomTransform(gen);
            const Eigen::Matrix3Xf dst = transformCloud(src, T_true);

            const Eigen::VectorXf w = Eigen::VectorXf::Ones(src.cols());
            const Eigen::Isometry3f T_w = weightedKabsch(src, dst, w);
            const Eigen::Isometry3f T_k = kabsch(src, dst);

            EXPECT_TRUE(T_w.matrix().isApprox(T_k.matrix(), 1e-4f));
        }

        TEST(WeightedKabsch, ZeroWeightsFallsBackToKabsch)
        {
            std::mt19937 gen(88);
            const Eigen::Matrix3Xf src = randomCloud(gen, 100);
            const Eigen::Isometry3f T_true = randomTransform(gen);
            const Eigen::Matrix3Xf dst = transformCloud(src, T_true);

            const Eigen::VectorXf w = Eigen::VectorXf::Zero(src.cols());
            const Eigen::Isometry3f T_w = weightedKabsch(src, dst, w);
            const Eigen::Isometry3f T_k = kabsch(src, dst);

            EXPECT_TRUE(T_w.matrix().isApprox(T_k.matrix(), 1e-4f));
        }

        TEST(WeightedKabsch, OutlierDownWeightingRecoversTruth)
        {
            /* Down-weight outliers to near zero; recovered transform must match truth. */
            std::mt19937 gen(99);
            constexpr int N = 300;
            constexpr int N_out = 100;

            const Eigen::Matrix3Xf src = randomCloud(gen, N);
            const Eigen::Isometry3f T_true = randomTransform(gen);
            Eigen::Matrix3Xf dst = transformCloud(src, T_true);
            injectOutliers(gen, dst, N_out);

            Eigen::VectorXf w(N);
            for (int i = 0; i < N; ++i)
                w(i) = (i < N_out) ? 0.0f : 1.0f;

            const Eigen::Isometry3f T_est = weightedKabsch(src, dst, w);
            const auto err = errorBetween(T_est, T_true);
            EXPECT_LT(err.rotation_deg, kXformTolRotDeg);
            EXPECT_LT(err.translation, kXformTolTrans);
        }

        TEST(WeightedKabsch, NegativeWeightThrows)
        {
            std::mt19937 gen(100);
            const Eigen::Matrix3Xf src = randomCloud(gen, 10);
            const Eigen::Matrix3Xf dst = src;
            Eigen::VectorXf w = Eigen::VectorXf::Ones(10);
            w(3) = -0.5f;
            EXPECT_THROW(weightedKabsch(src, dst, w), std::invalid_argument);
        }

        TEST(KabschPerformance, Benchmark)
        {
            std::mt19937 gen(2026);
            constexpr int trials = 100;
            const std::array<int, 3> sizes = {100, 1000, 4096};

            for (int N : sizes)
            {
                std::vector<std::pair<Eigen::Matrix3Xf, Eigen::Matrix3Xf>> inputs;
                inputs.reserve(trials);
                for (int t = 0; t < trials; ++t)
                {
                    Eigen::Matrix3Xf src = randomCloud(gen, N);
                    Eigen::Matrix3Xf dst = transformCloud(src, randomTransform(gen));
                    inputs.emplace_back(std::move(src), std::move(dst));
                }
                const Eigen::VectorXf w_ones = Eigen::VectorXf::Ones(N);

                std::cout << "\n[kabsch benchmark, N=" << N << ", " << trials << " trials]\n";
                benchmark("kabsch (Eigen)", trials, [&](int i)
                          {
            const auto T = kabsch(inputs[i].first, inputs[i].second);
            g_sink += T.translation().x(); });
                benchmark("weightedKabsch", trials, [&](int i)
                          {
            const auto T = weightedKabsch(inputs[i].first, inputs[i].second, w_ones);
            g_sink += T.translation().x(); });
            }
            EXPECT_TRUE(std::isfinite(static_cast<float>(g_sink)));
        }

        /* ================================================================
         *           Section 2: DCM to quaternion conversion
         * ================================================================ */

        TEST(DcmToQuat, Identity)
        {
            const Eigen::Matrix3f I = Eigen::Matrix3f::Identity();
            const Eigen::Quaternionf qI = Eigen::Quaternionf::Identity();

            EXPECT_LT(angularDistance(fromDCM(I), qI), kQuatTolTight);
            EXPECT_LT(angularDistance(sarabandi(I), qI), kQuatTolTight);
            EXPECT_LT(angularDistance(bar_itzhack(I), qI), kQuatTolTight);
        }

        TEST(DcmToQuat, CardinalAxisRotations)
        {
            /* 90-deg rotations about each axis, both signs. */
            constexpr float half_pi = static_cast<float>(M_PI) * 0.5f;
            constexpr std::array<const char *, 3> names = {"X", "Y", "Z"};

            for (int axis = 0; axis < 3; ++axis)
            {
                for (float sign : {1.0f, -1.0f})
                {
                    Eigen::Vector3f a = Eigen::Vector3f::Zero();
                    a(axis) = 1.0f;
                    const Eigen::AngleAxisf aa(sign * half_pi, a);
                    const Eigen::Matrix3f R = aa.toRotationMatrix();
                    const Eigen::Quaternionf q_ref(aa);

                    SCOPED_TRACE(std::string(sign > 0 ? "+" : "-") + "90 deg about " + names[axis]);
                    EXPECT_LT(angularDistance(fromDCM(R), q_ref), kQuatTolTight);
                    EXPECT_LT(angularDistance(sarabandi(R), q_ref), kQuatTolTight);
                    EXPECT_LT(angularDistance(bar_itzhack(R), q_ref), kQuatTolTight);
                }
            }
        }

        TEST(DcmToQuat, OneEightyDegreeRotations)
        {
            /* 180-deg rotations: dw = -1 exercises the min-branch of Shepperd-style dispatch. */
            constexpr float pi = static_cast<float>(M_PI);
            constexpr std::array<const char *, 3> names = {"X", "Y", "Z"};

            for (int axis = 0; axis < 3; ++axis)
            {
                Eigen::Vector3f a = Eigen::Vector3f::Zero();
                a(axis) = 1.0f;
                const Eigen::AngleAxisf aa(pi, a);
                const Eigen::Matrix3f R = aa.toRotationMatrix();
                const Eigen::Quaternionf q_ref(aa);

                SCOPED_TRACE(std::string("180 deg about ") + names[axis]);
                EXPECT_LT(angularDistance(fromDCM(R), q_ref), kQuatTolLoose);
                EXPECT_LT(angularDistance(sarabandi(R), q_ref), kQuatTolLoose);
                EXPECT_LT(angularDistance(bar_itzhack(R), q_ref), kQuatTolLoose);
            }
        }

        TEST(DcmToQuat, RandomRotationsAngularError)
        {
            constexpr int N = 1000;
            std::mt19937 gen(42);

            double sum[4] = {0.0, 0.0, 0.0, 0.0};
            float max_[4] = {0.0f, 0.0f, 0.0f, 0.0f};
            constexpr std::array<const char *, 4> labels = {"sarabandi", "bar_itzhack", "fromDCM", "Eigen(R)"};

            for (int i = 0; i < N; ++i)
            {
                const Eigen::Quaternionf q_ref = randomQuat(gen);
                const Eigen::Matrix3f R = q_ref.toRotationMatrix();

                const float e[4] = {
                    angularDistance(sarabandi(R), q_ref),
                    angularDistance(bar_itzhack(R), q_ref),
                    angularDistance(fromDCM(R), q_ref),
                    angularDistance(Eigen::Quaternionf(R), q_ref),
                };
                for (int k = 0; k < 4; ++k)
                {
                    sum[k] += e[k];
                    max_[k] = std::max(max_[k], e[k]);
                }
                for (int k = 0; k < 3; ++k) /* don't assert on Eigen's built-in */
                    ASSERT_LT(e[k], kQuatTolLoose) << labels[k] << " failure at sample " << i;
            }

            std::cout << "\n[angular error, N=" << N << " random rotations, radians]\n"
                      << std::scientific << std::setprecision(3);
            for (int k = 0; k < 4; ++k)
            {
                std::cout << "  " << std::left << std::setw(14) << labels[k]
                          << " mean=" << sum[k] / N << "   max=" << max_[k] << "\n";
            }
        }

        TEST(DcmToQuat, ResultIsUnitQuaternion)
        {
            std::mt19937 gen(7);
            for (int i = 0; i < 200; ++i)
            {
                const Eigen::Matrix3f R = randomDCM(gen);
                EXPECT_NEAR(sarabandi(R).norm(), 1.0f, kQuatNormTol);
                EXPECT_NEAR(bar_itzhack(R).norm(), 1.0f, kQuatNormTol);
                EXPECT_NEAR(fromDCM(R).norm(), 1.0f, kQuatNormTol);
            }
        }

        TEST(DcmToQuat, NonOrthogonalMatrix)
        {
            /* Additive noise large enough to fail isRotationMatrix but keep good conditioning. */
            std::mt19937 gen(2024);
            std::normal_distribution<float> noise(0.0f, 1e-3f);

            int robust_path_count = 0;
            for (int i = 0; i < 200; ++i)
            {
                const Eigen::Quaternionf q_ref = randomQuat(gen);
                const Eigen::Matrix3f R = q_ref.toRotationMatrix();
                Eigen::Matrix3f N;
                for (int r = 0; r < 3; ++r)
                    for (int c = 0; c < 3; ++c)
                        N(r, c) = noise(gen);
                const Eigen::Matrix3f R_dirty = R + N;

                EXPECT_LT(angularDistance(bar_itzhack(R_dirty), q_ref), 5e-3f);
                EXPECT_LT(angularDistance(fromDCM(R_dirty), q_ref), 5e-3f);
                if (!eststack::isRotationMatrix(R_dirty))
                    ++robust_path_count;
            }
            EXPECT_GT(robust_path_count, 0) << "non-orthogonal dispatch branch never exercised";
        }

        TEST(DcmToQuat, IllConditionedReturnsIdentity)
        {
            Eigen::Matrix3f R_sing;
            R_sing << 1, 0, 0,
                0, 1, 0,
                0, 0, 0;
            EXPECT_LT(angularDistance(fromDCM(R_sing), Eigen::Quaternionf::Identity()), kQuatTolTight);
        }

        TEST(DcmToQuat, DispatchMatchesSarabandiOnValidR)
        {
            std::mt19937 gen(99);
            for (int i = 0; i < 100; ++i)
            {
                const Eigen::Matrix3f R = randomDCM(gen);
                EXPECT_LT(angularDistance(fromDCM(R), sarabandi(R)), 1e-6f);
            }
        }

        TEST(DcmToQuatPerformance, Benchmark)
        {
            constexpr int N = 100000;
            std::mt19937 gen(2026);

            std::vector<Eigen::Matrix3f> rots;
            rots.reserve(N);
            for (int i = 0; i < N; ++i)
                rots.push_back(randomDCM(gen));

            std::cout << "\n[DCM->quat benchmark, N=" << N << "]\n";
            benchmark("sarabandi", N, [&](int i)
                      {
        const auto q = sarabandi(rots[i]);
        g_sink += q.w() + q.x() + q.y() + q.z(); });
            benchmark("bar_itzhack", N, [&](int i)
                      {
        const auto q = bar_itzhack(rots[i]);
        g_sink += q.w() + q.x() + q.y() + q.z(); });
            benchmark("fromDCM", N, [&](int i)
                      {
        const auto q = fromDCM(rots[i]);
        g_sink += q.w() + q.x() + q.y() + q.z(); });
            benchmark("Eigen(R)", N, [&](int i)
                      {
        const Eigen::Quaternionf q(rots[i]);
        g_sink += q.w() + q.x() + q.y() + q.z(); });
            EXPECT_TRUE(std::isfinite(static_cast<float>(g_sink)));
        }

        /* ================================================================
         *                   Section 3: IRLS
         * ================================================================ */

        TEST(IrlsLoss, CauchyWeightFormula)
        {
            CauchyLoss loss;
            /* w(r, s) = s^2 / (r^2 + s^2) */
            EXPECT_FLOAT_EQ(loss.weight(0.0f, 1.0f), 1.0f);
            EXPECT_NEAR(loss.weight(1.0f, 1.0f), 0.5f, 1e-6f);
            EXPECT_NEAR(loss.weight(3.0f, 1.0f), 0.1f, 1e-6f);
            EXPECT_FLOAT_EQ(loss.weight(2.0f, 1.0f), loss.weight(-2.0f, 1.0f));
            EXPECT_GT(loss.weight(0.5f, 1.0f), loss.weight(2.0f, 1.0f));
        }

        TEST(IrlsLoss, HuberWeightFormula)
        {
            HuberLoss loss;
            /* w = 1 within scale; scale / |r| outside */
            EXPECT_FLOAT_EQ(loss.weight(0.0f, 1.0f), 1.0f);
            EXPECT_FLOAT_EQ(loss.weight(1.0f, 1.0f), 1.0f);
            EXPECT_NEAR(loss.weight(2.0f, 1.0f), 0.5f, 1e-6f);
            EXPECT_NEAR(loss.weight(4.0f, 1.0f), 0.25f, 1e-6f);
        }

        /* isInlier shares the same threshold across loss types -> typed test */
        template <typename Loss>
        struct LossTypedTest : public ::testing::Test
        {
        };
        using LossTypes = ::testing::Types<CauchyLoss, HuberLoss>;
        TYPED_TEST_SUITE(LossTypedTest, LossTypes);

        TYPED_TEST(LossTypedTest, InlierThreshold)
        {
            TypeParam loss;
            EXPECT_TRUE(loss.isInlier(2.0f, 1.0f));
            EXPECT_FALSE(loss.isInlier(3.5f, 1.0f));
            EXPECT_TRUE(loss.isInlier(5.0f, 2.0f));
            EXPECT_FALSE(loss.isInlier(7.0f, 2.0f));
        }

        /* ---- correctness, parameterized over loss type ---- */
        template <typename Loss>
        struct IrlsTypedTest : public ::testing::Test
        {
            RigidIRLS<Loss> irls;
            void SetUp() override
            {
                irls.setLossFunction(std::make_shared<const Loss>());
            }
        };
        TYPED_TEST_SUITE(IrlsTypedTest, LossTypes);

        TYPED_TEST(IrlsTypedTest, IdentityRecovery)
        {
            std::mt19937 gen(42);
            const Eigen::Matrix3Xf src = randomCloud(gen, 200);

            this->irls.setInputSource(src);
            this->irls.setInputTarget(src);
            this->irls.optimize();

            EXPECT_TRUE(this->irls.hasConverged());
            const auto err = errorBetween(this->irls.getFinalTransformation(), Eigen::Isometry3f::Identity());
            EXPECT_LT(err.rotation_deg, kXformTolRotDeg);
            EXPECT_LT(err.translation, kXformTolTrans);
        }

        TYPED_TEST(IrlsTypedTest, CleanRigidTransform)
        {
            std::mt19937 gen(123);
            const Eigen::Matrix3Xf src = randomCloud(gen, 300);
            const Eigen::Isometry3f T_true = randomTransform(gen);
            const Eigen::Matrix3Xf dst = transformCloud(src, T_true);

            this->irls.setInputSource(src);
            this->irls.setInputTarget(dst);
            this->irls.optimize();

            EXPECT_TRUE(this->irls.hasConverged());
            const auto err = errorBetween(this->irls.getFinalTransformation(), T_true);
            EXPECT_LT(err.rotation_deg, kXformTolRotDeg);
            EXPECT_LT(err.translation, kXformTolTrans);
        }

        TYPED_TEST(IrlsTypedTest, NoisyRigidTransform)
        {
            std::mt19937 gen(456);
            const Eigen::Matrix3Xf src = randomCloud(gen, 500);
            const Eigen::Isometry3f T_true = randomTransform(gen);
            Eigen::Matrix3Xf dst = transformCloud(src, T_true);
            dst = addGaussianNoise(gen, dst, 0.05f);

            this->irls.setInputSource(src);
            this->irls.setInputTarget(dst);
            this->irls.optimize();

            const auto err = errorBetween(this->irls.getFinalTransformation(), T_true);
            EXPECT_LT(err.rotation_deg, 1.0f);
            EXPECT_LT(err.translation, 0.1f);
        }

        TYPED_TEST(IrlsTypedTest, RobustToOutliers)
        {
            std::mt19937 gen(789);
            constexpr int N = 500;
            constexpr int N_out = 100; /* 20% */

            const Eigen::Matrix3Xf src = randomCloud(gen, N);
            const Eigen::Isometry3f T_true = randomTransform(gen);
            Eigen::Matrix3Xf dst = transformCloud(src, T_true);
            dst = addGaussianNoise(gen, dst, 0.01f);
            injectOutliers(gen, dst, N_out);

            /* baseline pulled off by outliers */
            const auto err_plain = errorBetween(kabsch(src, dst), T_true);

            this->irls.setInputSource(src);
            this->irls.setInputTarget(dst);
            this->irls.optimize();
            const auto err_irls = errorBetween(this->irls.getFinalTransformation(), T_true);

            EXPECT_LT(err_irls.rotation_deg, 2.0f);
            EXPECT_LT(err_irls.translation, 0.5f);
            EXPECT_LT(err_irls.rotation_deg, err_plain.rotation_deg);
            EXPECT_LT(err_irls.translation, err_plain.translation);

            const char *tag = std::is_same_v<TypeParam, CauchyLoss> ? "Cauchy" : "Huber";
            std::cout << "\n[outlier robustness, " << tag << ", N=" << N << " w/ " << N_out << " outliers]\n"
                      << "  plain Kabsch : rot=" << err_plain.rotation_deg << " deg   trans=" << err_plain.translation << "\n"
                      << "  IRLS " << std::left << std::setw(8) << tag
                      << ": rot=" << err_irls.rotation_deg << " deg   trans=" << err_irls.translation << "\n";
        }

        TYPED_TEST(IrlsTypedTest, GoodInitialGuessStillConverges)
        {
            std::mt19937 gen(606);
            const Eigen::Matrix3Xf src = randomCloud(gen, 300);
            const Eigen::Isometry3f T_true = randomTransform(gen);
            Eigen::Matrix3Xf dst = transformCloud(src, T_true);
            dst = addGaussianNoise(gen, dst, 0.02f);

            this->irls.setInputSource(src);
            this->irls.setInputTarget(dst);
            this->irls.setInitGuess(T_true);
            this->irls.optimize();

            const auto err = errorBetween(this->irls.getFinalTransformation(), T_true);
            EXPECT_LT(err.rotation_deg, 1.0f);
            EXPECT_LT(err.translation, 0.1f);
        }

        /* ---- edge cases, loss-agnostic ---- */

        TEST(IrlsEdgeCase, EmptyInput)
        {
            RigidIRLS<CauchyLoss> irls;
            irls.setLossFunction(std::make_shared<const CauchyLoss>());
            irls.setInputSource(Eigen::Matrix3Xf(3, 0));
            irls.setInputTarget(Eigen::Matrix3Xf(3, 0));
            irls.optimize();

            EXPECT_FALSE(irls.hasConverged());
            EXPECT_TRUE(irls.getFinalTransformation().matrix().isApprox(Eigen::Matrix4f::Identity()));
        }

        TEST(IrlsEdgeCase, SizeMismatch)
        {
            std::mt19937 gen(1);
            RigidIRLS<CauchyLoss> irls;
            irls.setLossFunction(std::make_shared<const CauchyLoss>());
            irls.setInputSource(randomCloud(gen, 10));
            irls.setInputTarget(randomCloud(gen, 20));
            irls.optimize();

            EXPECT_FALSE(irls.hasConverged());
            EXPECT_TRUE(irls.getFinalTransformation().matrix().isApprox(Eigen::Matrix4f::Identity()));
        }

        TEST(IrlsEdgeCase, NoLossFunction)
        {
            std::mt19937 gen(2);
            RigidIRLS<CauchyLoss> irls;
            /* intentionally skip setLossFunction */
            irls.setInputSource(randomCloud(gen, 100));
            irls.setInputTarget(randomCloud(gen, 100));
            irls.optimize();

            EXPECT_FALSE(irls.hasConverged());
            EXPECT_TRUE(irls.getFinalTransformation().matrix().isApprox(Eigen::Matrix4f::Identity()));
        }

        TEST(IrlsPerformance, Benchmark)
        {
            std::mt19937 gen(2026);
            constexpr int N = 1000;
            constexpr int N_out = 200; /* 20% */
            constexpr int trials = 100;

            std::vector<std::pair<Eigen::Matrix3Xf, Eigen::Matrix3Xf>> inputs;
            inputs.reserve(trials);
            for (int t = 0; t < trials; ++t)
            {
                Eigen::Matrix3Xf src = randomCloud(gen, N);
                Eigen::Matrix3Xf dst = transformCloud(src, randomTransform(gen));
                dst = addGaussianNoise(gen, dst, 0.01f);
                injectOutliers(gen, dst, N_out);
                inputs.emplace_back(std::move(src), std::move(dst));
            }

            RigidIRLS<CauchyLoss> irls_c;
            irls_c.setLossFunction(std::make_shared<const CauchyLoss>());
            RigidIRLS<HuberLoss> irls_h;
            irls_h.setLossFunction(std::make_shared<const HuberLoss>());

            std::cout << "\n[IRLS benchmark, N=" << N << ", " << N_out
                      << " outliers, " << trials << " trials]\n";

            benchmark("plain Kabsch", trials, [&](int i)
                      { g_sink += kabsch(inputs[i].first, inputs[i].second).translation().x(); });
            benchmark("IRLS Cauchy", trials, [&](int i)
                      {
        irls_c.setInputSource(inputs[i].first);
        irls_c.setInputTarget(inputs[i].second);
        irls_c.optimize();
        g_sink += irls_c.getFinalTransformation().translation().x(); });
            benchmark("IRLS Huber", trials, [&](int i)
                      {
        irls_h.setInputSource(inputs[i].first);
        irls_h.setInputTarget(inputs[i].second);
        irls_h.optimize();
        g_sink += irls_h.getFinalTransformation().translation().x(); });
            EXPECT_TRUE(std::isfinite(static_cast<float>(g_sink)));
        }

    } // namespace test
} // namespace eststack

#endif //! TEST__TRANSFORM_BENCHMARK_CPP
