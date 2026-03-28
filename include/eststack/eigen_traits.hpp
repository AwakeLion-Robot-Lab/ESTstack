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

#ifndef ESTSTACK__EIGEN_TRAITS_HPP
#define ESTSTACK__EIGEN_TRAITS_HPP

// C++ standard library
#include <limits>
#include <exception>

// Eigen library
#include <Eigen/Dense>
#include <Eigen/Cholesky>
#include <Eigen/SVD>

/***
 * @brief An algorithm set focus on state estimation
 * @author jinhua "siyiovo" deng
 */
namespace eststack
{
    /***
     * @brief check if a matrix is symmetric
     */
    template <typename Derived>
    bool isSymmetric(const Eigen::MatrixBase<Derived> &m,
                     typename Derived::Scalar threshold = 1e-6)
    {
        if (m.rows() != m.cols())
            throw std::invalid_argument("Matrix is not square to be checked!");

        return m.isApprox(m.transpose(), threshold);
    }

    /***
     * @brief check if a symmetric matrix is positive semi-definite via Cholesky
     */
    template <typename Derived>
    bool isPositiveDefinite(const Eigen::MatrixBase<Derived> &m)
    {
        if (m.rows() != m.cols())
            throw std::invalid_argument("Matrix must be square to be checked!");

        Eigen::LLT<Eigen::Matrix<typename Derived::Scalar, Derived::RowsAtCompileTime, Derived::ColsAtCompileTime>> llt(m);
        return llt.info() == Eigen::Success;
    }

    /***
     * @brief check if a matrix is a valid covariance
     */
    template <typename Derived>
    bool isCovariance(const Eigen::MatrixBase<Derived> &m,
                      typename Derived::Scalar threshold = 1e-6)
    {
        if (isSymmetric(m, threshold))
            return isPositiveDefinite(m);
        throw std::invalid_argument("Matrix is not symmetric!");
    }

    template <typename Derived>
    bool isOrthogonal(const Eigen::MatrixBase<Derived> &m,
                      typename Derived::Scalar threshold = 1e-6)
    {
        if (m.rows() != m.cols())
            throw std::invalid_argument("Matrix must be square to be checked!");

        return (std::abs(m.determinant() - 1.0) < threshold) && (m.transpose() * m).isIdentity(threshold);
    }

    /***
     * @brief get the condition number of a matrix via SVD
     * @note L-2 norm condition number
     */
    template <typename Derived>
    typename Derived::Scalar getConditionNumber(const Eigen::MatrixBase<Derived> &m)
    {
        Eigen::JacobiSVD<Eigen::Matrix<typename Derived::Scalar, Derived::RowsAtCompileTime, Derived::ColsAtCompileTime>> svd(m, Eigen::ComputeThinU | Eigen::ComputeThinV);
        const auto singulars = svd.singularValues();
        const auto max_sigma = singulars(0);
        const auto min_sigma = singulars(singulars.size() - 1);

        if (min_sigma <= 0)
            return std::numeric_limits<typename Derived::Scalar>::infinity();
        return max_sigma / min_sigma;
    }

    //@TODO: add Posteriori CRLB to check whether kf converge

} // namespace eststack

#endif //! ESTSTACK__EIGEN_TRAITS_HPP
