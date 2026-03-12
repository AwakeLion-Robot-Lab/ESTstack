
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

#ifndef ESTSTACK__TYPES_HPP
#define ESTSTACK__TYPES_HPP

// C++ standard library
#include <concepts>

// Eigen library
#include <Eigen/Dense>

/***
 * @brief An algorithm set focus on estimation and filtering
 * @author jinhua "siyiovo" deng
 */
namespace eststack
{
    /***
     * @brief DoF of lie group concept
     * @details for DoF extraction
     */
    template <typename T>
    concept LieGroupDoF = requires {
        typename T::Scalar;
        { T::DoF } -> std::convertible_to<int>;
        requires(T::DoF > 0);
    };

    /***
     * @brief Eigen row vector concept
     * @details for Eigen fixed-size row extraction
     */
    template <typename T>
    concept EigenRow = requires {
        typename T::Scalar;
        { T::RowsAtCompileTime } -> std::convertible_to<int>;
        requires(T::RowsAtCompileTime > 0);
    };

    /***
     * @brief static dimension at compile time concept
     * @details for dimension extraction at compile time
     */
    template <typename T>
    concept DimAtCompileTime = LieGroupDoF<T> || EigenRow<T>;

    /***
     * @brief get static dimension of a type at compile time
     * @tparam T type with static dimension at compile time, either Lie group state or fixed-size Eigen vector
     */
    template <DimAtCompileTime T>
    consteval int getDimAtCompileTime()
    {
        if constexpr (LieGroupDoF<T>)
            return static_cast<int>(T::DoF);
        else
            return static_cast<int>(T::RowsAtCompileTime);
    }

    /***
     * @brief compile-time fixed-size Jacobian matrix
     * @tparam A dimension source type for rows
     * @tparam B dimension source type for columns
     * @note Jacobian<A,B> means A w.r.t. B
     */
    template <DimAtCompileTime A, DimAtCompileTime B>
    using Jacobian = Eigen::Matrix<typename A::Scalar, getDimAtCompileTime<A>(), getDimAtCompileTime<B>()>;

    /***
     * @brief compile-time fixed-size square covariance matrix
     * @tparam T dimension source type
     */
    template <DimAtCompileTime T>
    using Covariance = Jacobian<T, T>;

    /***
     * @brief perturbation convention on manifold
     * @details Global (left): y = Exp({\epsilon}^δx) ∘ x
     *          Local (right): y = x ∘ Exp({\tau}^δx)
     *          it's well known that y = Exp({\epsilon}^δx) ∘ x = x ∘ Exp({\tau}^δx)
     *          which get adjoint of manifold: Exp({\epsilon}^δx) = x ∘ Exp({\tau}^δx) ∘ x^{-1}
     */
    enum class Perturbation
    {
        Global,
        Local
    };

}

#endif //! ESTSTACK__TYPES_HPP