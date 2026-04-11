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
     * @brief compile-time 95% confidence chi-square critical values table
     */
    struct ChiSquareTable
    {
        /***
         * @brief get 95% confidence chi-square critical value at compile-time
         */
        template <int DoF>
        static constexpr double value()
        {
            static_assert(DoF >= 1, "DoF must be greater than 0!");

            /* get value from table if low-DoF */
            if constexpr (DoF <= 10)
            {
                constexpr std::array<double, 10> table = {
                    3.841, 5.991, 7.815, 9.488, 11.070,
                    12.592, 14.067, 15.507, 16.919, 18.307};
                return table[DoF - 1];
            }
            else /* if DoF > 10, use Wilson-Hilferty approximation */
            {
                constexpr double z = 1.645;
                constexpr double a = 2.0 / (9.0 * DoF);
                /* std::sqrt may not be constexpr */
                const double term = 1.0 - a + z * std::sqrt(a);
                return DoF * term * term * term;
            }
        }
    };

    /***
     * @brief check if a matrix is symmetric
     * @param m the matrix to be checked
     * @param threshold the threshold for checking symmetry
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
     * @param m the matrix to be checked
     */
    template <typename Derived>
    bool isPositiveDefinite(const Eigen::MatrixBase<Derived> &m)
    {
        if (m.rows() != m.cols())
            throw std::invalid_argument("Matrix is not square to be checked!");

        Eigen::LLT<Eigen::Matrix<typename Derived::Scalar, Derived::RowsAtCompileTime, Derived::ColsAtCompileTime>> llt(m);
        return llt.info() == Eigen::Success;
    }

    /***
     * @brief check if a matrix is a valid covariance
     * @param m the matrix to be checked
     * @param threshold the threshold for checking symmetry and positive-definiteness
     */
    template <typename Derived>
    bool isCovariance(const Eigen::MatrixBase<Derived> &m,
                      typename Derived::Scalar threshold = 1e-6)
    {
        if (isSymmetric(m, threshold))
            return isPositiveDefinite(m);

        return false;
    }

    /***
     * @brief check if a matrix is orthogonal
     * @param m the matrix to be checked
     * @param threshold the threshold for checking orthogonality
     */
    template <typename Derived>
    bool isOrthogonal(const Eigen::MatrixBase<Derived> &m,
                      typename Derived::Scalar threshold = 1e-6)
    {
        if (m.rows() != m.cols())
            throw std::invalid_argument("Matrix must be square to be checked!");

        return (m.transpose() * m).isIdentity(threshold);
    }

    /***
     * @brief check if a matrix is a valid rotation matrix
     * @param m the matrix to be checked
     * @param ortho_threshold the threshold for checking orthogonality
     * @param det_threshold the threshold for checking determinant being 1
     */
    template <typename Derived>
    bool isRotationMatrix(const Eigen::MatrixBase<Derived> &m,
                          typename Derived::Scalar ortho_threshold = 1e-6,
                          typename Derived::Scalar det_threshold = 1e-6)
    {
        if (isOrthogonal(m, ortho_threshold))
            return (std::abs(m.determinant() - 1.0) < det_threshold);

        return false;
    }

    /***
     * @brief get the condition number of a matrix via SVD
     * @param m the matrix to be checked
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

    /***
     * @brief get non-zero column indices from vector
     * @param flags judgement
     * @return vector with indices of non-zero elements
     */
    inline Eigen::VectorXi getNonZeroColumnIndicesFromVector(const Eigen::VectorXi &flags)
    {
        int count = flags.count();
        Eigen::VectorXi nonzero_column(count);
        int idx = 0;
        for (int i = 0; i < flags.size(); ++i)
        {
            if (flags(i) > 0)
                nonzero_column(idx++) = i;
        }
        return nonzero_column;
    }

} // namespace eststack

#endif //! ESTSTACK__EIGEN_TRAITS_HPP
