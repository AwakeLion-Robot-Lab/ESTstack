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
#include <array>
#include <cmath>

// Eigen library
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>

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
         * @brief get best quaternion from DCM, even if the DCM is not orthogonal
         * @param dcm direction cosine matrix
         * @return best quaternion
         * @details see `sarabandi` and `bar_itzhack`
         */
        Eigen::Quaterniond from_dcm(const Eigen::Matrix3d &dcm)
        {
            if (isOrthogonal(dcm))
                return sarabandi(dcm);
            else
                return bar_itzhack(dcm);
        }

        /***
         * @brief Sarabandi's method to extract precise DCM to quaternion
         * @param dcm direction cosine matrix
         * @param eta threshold to decide whether to use the fast method or the robust method
         * @return best quaternion
         * @details DCM MUST BE orthogonal, refer to [12] and https://github.com/Mayitzin/ahrs/blob/master/ahrs/common/orientation.py#L1117
         *          eta threshold selection can refer to docs/explanations/dcm_to_quat.md
         */
        Eigen::Quaterniond sarabandi(const Eigen::Matrix3d &dcm, double eta = 0.0)
        {
            const double r11 = dcm(0, 0);
            const double r12 = dcm(0, 1);
            const double r13 = dcm(0, 2);
            const double r21 = dcm(1, 0);
            const double r22 = dcm(1, 1);
            const double r23 = dcm(1, 2);
            const double r31 = dcm(2, 0);
            const double r32 = dcm(2, 1);
            const double r33 = dcm(2, 2);

            const double dw = r11 + r22 + r33;
            const double dx = r11 - r22 - r33;
            const double dy = -r11 + r22 - r33;
            const double dz = -r11 - r22 + r33;

            double qw, qx, qy, qz;

            if (dw > eta)
            {
                qw = 0.5 * std::sqrt(1.0 + dw);
            }
            else
            {
                const double nom = (r32 - r23) * (r32 - r23) +
                                   (r13 - r31) * (r13 - r31) +
                                   (r21 - r12) * (r21 - r12);
                const double denom = 3.0 - dw;
                qw = (denom > 0.0) ? 0.5 * std::sqrt(nom / denom) : 0.0;
            }

            if (dx > eta)
            {
                qx = 0.5 * std::sqrt(1.0 + dx);
            }
            else
            {
                const double nom = (r32 - r23) * (r32 - r23) +
                                   (r12 + r21) * (r12 + r21) +
                                   (r31 + r13) * (r31 + r13);
                const double denom = 3.0 - dx;
                qx = (denom > 0.0) ? 0.5 * std::sqrt(nom / denom) : 0.0;
            }

            if (dy > eta)
            {
                qy = 0.5 * std::sqrt(1.0 + dy);
            }
            else
            {
                const double nom = (r13 - r31) * (r13 - r31) +
                                   (r12 + r21) * (r12 + r21) +
                                   (r23 + r32) * (r23 + r32);
                const double denom = 3.0 - dy;
                qy = (denom > 0.0) ? 0.5 * std::sqrt(nom / denom) : 0.0;
            }

            if (dz > eta)
            {
                qz = 0.5 * std::sqrt(1.0 + dz);
            }
            else
            {
                const double nom = (r21 - r12) * (r21 - r12) +
                                   (r31 + r13) * (r31 + r13) +
                                   (r23 + r32) * (r23 + r32);
                const double denom = 3.0 - dz;
                qz = (denom > 0.0) ? 0.5 * std::sqrt(nom / denom) : 0.0;
            }

            if (qw < 0.0)
            {
                qw = -qw;
                qx = -qx;
                qy = -qy;
                qz = -qz;
            }

            qx *= (r32 - r23 >= 0.0) ? 1.0 : -1.0;
            qy *= (r13 - r31 >= 0.0) ? 1.0 : -1.0;
            qz *= (r21 - r12 >= 0.0) ? 1.0 : -1.0;

            return Eigen::Quaterniond(qw, qx, qy, qz).normalized();
        }

        /***
         * @brief Bar-Itzhack's method to extract DCM (even not orthogonal) to quaternion
         * @param dcm direction cosine matrix
         * @return best quaternion
         * @details refer to [14], here I ONLY implement Version 3 algorithm 'cause ONLY it can handle non-orthogonal DCM
         */
        Eigen::Quaterniond bar_itzhack(const Eigen::Matrix3d &dcm)
        {
            auto vec_pairs = std::array<std::pair<Eigen::Vector3d, Eigen::Vector3d>, 3>{
                std::make_pair(dcm.row(0), Eigen::Vector3d::UnitX()),
                std::make_pair(dcm.row(1), Eigen::Vector3d::UnitY()),
                std::make_pair(dcm.row(2), Eigen::Vector3d::UnitZ())};

            const double weight_val = 1.0 / 3.0;
            auto weights = std::array<double, 3>{weight_val, weight_val, weight_val};

            return steady_q_method(vec_pairs, weights);
        }

        /***
         * @brief get best quaternion via eigenvalue decomposition from 3 vector pairs
         * @param vec_pairs array of 3 vector pairs `(u_i, v_i)` where `u_i` is in body frame, `v_i` is in reference frame
         * @param weights weights for each vector pair
         * @return best quaternion
         * @details refer to [13] and docs/explanations/dcm_to_quat.md
         */
        Eigen::Quaterniond steady_q_method(const std::array<std::pair<Eigen::Vector3d, Eigen::Vector3d>, 3> &vec_pairs,
                                           const std::array<double, 3> &weights)
        {
            /* compute matrix B */
            Eigen::Matrix3d B = Eigen::Matrix3d::Zero();
            for (int i = 0; i < 3; ++i)
            {
                B += weights[i] * vec_pairs[i].first * vec_pairs[i].second.transpose(); /* 3-dim tensor, not scalar */
            }

            /* compute matrix K */
            const double sigma = B.trace();
            const Eigen::Vector3d z(B(1, 2) - B(2, 1),
                                    B(2, 0) - B(0, 2),
                                    B(0, 1) - B(1, 0));

            const double z0 = z(0);
            const double z1 = z(1);
            const double z2 = z(2);
            const double s00 = 2.0 * B(0, 0) - sigma;
            const double s11 = 2.0 * B(1, 1) - sigma;
            const double s22 = 2.0 * B(2, 2) - sigma;
            const double s01 = B(0, 1) + B(1, 0);
            const double s02 = B(0, 2) + B(2, 0);
            const double s12 = B(1, 2) + B(2, 1);

            Eigen::Matrix4d K = Eigen::Matrix4d::Zero();
            K << sigma, z0, z1, z2,
                z0, s00, s01, s02,
                z1, s01, s11, s12,
                z2, s02, s12, s22;

            /* get eigenvector corresponding to maximum eigenvalue */
            Eigen::SelfAdjointEigenSolver<Eigen::Matrix4d> solver(K, Eigen::ComputeEigenvectors);
            Eigen::Vector4d q = solver.eigenvectors().col(3);

            /* don't forget to normalize q */
            return Eigen::Quaterniond(q(0), q(1), q(2), q(3)).normalized();
        }
    } // namespace core
} // namespace eststack

#endif // CORE__DCM_TO_QUAT_HPP
