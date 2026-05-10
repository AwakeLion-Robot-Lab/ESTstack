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

#ifndef CORE__DCM_TO_QUAT_HPP
#define CORE__DCM_TO_QUAT_HPP

// C++ standard library
#include <cmath>

// Eigen library
#include <Eigen/Eigenvalues>
#include <Eigen/Geometry>
#include <Eigen/Core>

// ESTstack library
#include "eststack/eigen_traits.hpp"

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
         * @brief Sarabandi's method to extract precise DCM to quaternion
         * @param dcm direction cosine matrix
         * @param eta threshold to decide whether to use the fast method or the robust method
         * @return best quaternion
         * @details DCM MUST BE a valid rotation matrix, refer to [16] and https://github.com/Mayitzin/ahrs/blob/master/ahrs/common/orientation.py#L1117
         *          eta threshold selection can refer to docs/explanation/dcm_to_quat.md
         */
        Eigen::Quaternionf sarabandi(const Eigen::Matrix3f &dcm, float eta = 0.0f)
        {
            const float r11 = dcm(0, 0);
            const float r12 = dcm(0, 1);
            const float r13 = dcm(0, 2);
            const float r21 = dcm(1, 0);
            const float r22 = dcm(1, 1);
            const float r23 = dcm(1, 2);
            const float r31 = dcm(2, 0);
            const float r32 = dcm(2, 1);
            const float r33 = dcm(2, 2);

            const float dw = r11 + r22 + r33;
            const float dx = r11 - r22 - r33;
            const float dy = -r11 + r22 - r33;
            const float dz = -r11 - r22 + r33;

            float qw, qx, qy, qz;

            if (dw > eta)
            {
                qw = 0.5f * std::sqrt(1.0f + dw);
            }
            else
            {
                const float nom = (r32 - r23) * (r32 - r23) +
                                  (r13 - r31) * (r13 - r31) +
                                  (r21 - r12) * (r21 - r12);
                const float denom = 3.0f - dw;
                qw = (denom > 0.0f) ? 0.5f * std::sqrt(nom / denom) : 0.0f;
            }

            if (dx > eta)
            {
                qx = 0.5f * std::sqrt(1.0f + dx);
            }
            else
            {
                const float nom = (r32 - r23) * (r32 - r23) +
                                  (r12 + r21) * (r12 + r21) +
                                  (r31 + r13) * (r31 + r13);
                const float denom = 3.0f - dx;
                qx = (denom > 0.0f) ? 0.5f * std::sqrt(nom / denom) : 0.0f;
            }

            if (dy > eta)
            {
                qy = 0.5f * std::sqrt(1.0f + dy);
            }
            else
            {
                const float nom = (r13 - r31) * (r13 - r31) +
                                  (r12 + r21) * (r12 + r21) +
                                  (r23 + r32) * (r23 + r32);
                const float denom = 3.0f - dy;
                qy = (denom > 0.0f) ? 0.5f * std::sqrt(nom / denom) : 0.0f;
            }

            if (dz > eta)
            {
                qz = 0.5f * std::sqrt(1.0f + dz);
            }
            else
            {
                const float nom = (r21 - r12) * (r21 - r12) +
                                  (r31 + r13) * (r31 + r13) +
                                  (r23 + r32) * (r23 + r32);
                const float denom = 3.0f - dz;
                qz = (denom > 0.0f) ? 0.5f * std::sqrt(nom / denom) : 0.0f;
            }

            if (qw < 0.0f)
            {
                qw = -qw;
                qx = -qx;
                qy = -qy;
                qz = -qz;
            }
            qx *= (r32 - r23 >= 0.0f) ? 1.0f : -1.0f;
            qy *= (r13 - r31 >= 0.0f) ? 1.0f : -1.0f;
            qz *= (r21 - r12 >= 0.0f) ? 1.0f : -1.0f;

            Eigen::Quaternionf quat(qw, qx, qy, qz);
            quat.normalize();
            return quat;
        }

        /***
         * @brief Bar-Itzhack's method to extract DCM (even not orthogonal) to quaternion
         * @param dcm direction cosine matrix
         * @return best quaternion
         * @details refer to [17], [18] and docs/explanation/dcm_to_quat.md
         *          here I ONLY implement Version 3 algorithm 'cause ONLY it can handle non-orthogonal DCM
         */
        Eigen::Quaternionf bar_itzhack(const Eigen::Matrix3f &dcm)
        {
            /* compute matrix B */
            const float B00 = dcm(0, 0);
            const float B01 = dcm(0, 1);
            const float B02 = dcm(0, 2);
            const float B10 = dcm(1, 0);
            const float B11 = dcm(1, 1);
            const float B12 = dcm(1, 2);
            const float B20 = dcm(2, 0);
            const float B21 = dcm(2, 1);
            const float B22 = dcm(2, 2);

            /* compute matrix K */
            const float sigma = B00 + B11 + B22;
            const float z0 = B21 - B12;
            const float z1 = B02 - B20;
            const float z2 = B10 - B01;

            const float s00 = B00 - B11 - B22;
            const float s11 = -B00 + B11 - B22;
            const float s22 = -B00 - B11 + B22;
            const float s01 = B01 + B10;
            const float s02 = B02 + B20;
            const float s12 = B12 + B21;

            Eigen::Matrix4f K;
            K << sigma, z0, z1, z2,
                z0, s00, s01, s02,
                z1, s01, s11, s12,
                z2, s02, s12, s22;

            /* get eigenvector corresponding to maximum eigenvalue */
            Eigen::SelfAdjointEigenSolver<Eigen::Matrix4f> solver(K, Eigen::ComputeEigenvectors);
            Eigen::Vector4f q = solver.eigenvectors().col(3);

            /***
             * NOTE that Eigen::Quaternionf construct with Eigen::Vector4f is fully copy without swtiching position!
             * i.e. `Eigen::Quaternionf(q)` will construct a quaternion with x=q(0), y=q(1), z=q(2), w=q(3)
             * but normal constructor is `Eigen::Quaternionf(w, x, y, z): m_coeffs(x, y, z, w)`. see, it changes the order!
             * details refer to [Eigen 5.0.0 Quaternion documentation](https://libeigen.gitlab.io/eigen/docs-5.0/classEigen_1_1Quaternion.html),
             * just focus on constructor 2/9 and 7/9, here q(0) is real part, use constructor 2/9
             */
            Eigen::Quaternionf quat(q(0), q(1), q(2), q(3));
            quat.normalize();
            return quat;
        }

        /***
         * @brief get best quaternion from DCM, even if the DCM is not orthogonal
         * @param dcm direction cosine matrix
         * @return best quaternion
         * @details see `sarabandi` and `bar_itzhack`
         */
        Eigen::Quaternionf fromDCM(const Eigen::Matrix3f &dcm)
        {
            if (eststack::isRotationMatrix(dcm))
            {
                /* if DCM is a valid rotation matrix, use the fast method */
                return core::sarabandi(dcm);
            }
            else
            {
                /* if DCM is ill-conditioned seriously, return a default quaternion */
                if (eststack::getConditionNumber(dcm) > 1e5f)
                    return Eigen::Quaternionf::Identity();

                /* if DCM is not orthogonal but well-conditioned, use the robust method */
                return core::bar_itzhack(dcm);
            }
        }
    } // namespace core
} // namespace eststack

#endif //! CORE__DCM_TO_QUAT_HPP
