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

#ifndef CORE__IRLS_HPP
#define CORE__IRLS_HPP

// C++ standard library
#include <cmath>
#include <array>

// Eigen library
#include <Eigen/Dense>

// oneTBB library
#include <oneapi/tbb.h>

namespace eststack
{
    namespace core
    {
        /**
         * @brief Iteratively Reweighted Least Squares for robust regression
         * @tparam N number of residuals
         */
        template <int N>
        class IRLS
        {
        public:
            using VectorNi = Eigen::Vector<int, N>;
            using VectorNd = Eigen::Vector<double, N>;
            using MatrixNd = Eigen::Matrix<double, N, N>;

            /**
             * @brief default constructor
             * @param loss loss function instance
             * @param decay_factor factor to decay scale at each iteration
             * @param min_scale minimum scale to prevent over-shrinking
             * @param e_min minimum energy threshold for convergence
             * @param max_iter maximum number of iterations
             */
            explicit IRLS(double decay_factor = 1.3, double min_scale = 1.0, double e_min = 1e-2, int max_iter = 100)
                : decay_factor_(decay_factor), min_scale_(min_scale),
                  e_min_(e_min), max_iter_(max_iter)
            {
            }

            /**
             * @brief run IRLS optimization
             * @tparam LossFunction loss function type, e.g. Cauchy Loss, Huber Loss
             * @tparam ResidualFunction residual function type
             * @tparam Args argument types for residual function
             * @param update_weights function to update weights via robust loss function
             * @param compute_residual function to compute residuals
             * @param args arguments passed to residual function
             * @return true if converged
             */
            template <typename LossFunction, typename ResidualFunction, typename... Args>
            bool run(LossFunction &&update_weights, ResidualFunction &&compute_residual, Args &&...args)
            {
                /* initial weight vector is all ones */
                VectorNd weight_vec = VectorNd::Ones();
                /* record each inlier index in consensus set */
                VectorNi inlier_indices = VectorNi::LinSpaced(N, 0, N - 1);
                /* scale of loss function */
                double scale = 0.0;
                VectorNd residuals, weighted_residuals = VectorNd::Zero(N);

                /* iteration */
                for (int i = 1; i <= max_iter_; i++)
                {
                    /* compute residuals */
                    residuals = compute_residual(std::forward<Args>(args)...);

                    /* compute energy */
                    double e = weight_vec.cwiseProduct(residuals.cwiseAbs2()).sum();

                    /* compute weights via loss function */
                    if (i > 1)
                        weight_vec = update_weights(residuals, scale);

                    /* compute weighted residuals */
                    weighted_residuals = weight_vec.asDiagonal() * residuals;
                }
            }

        private:
            /***
             * @brief decay factor for scale update
             */
            double decay_factor_;

            /***
             * @brief minimum scale to prevent over-shrinking
             */
            double min_scale_;

            /***
             * @brief minimum energy threshold which represents change rate of residuals
             */
            double e_min_;

            /***
             * @brief max iteration of IRLS
             */
            int max_iter_;
        };
    } // namespace core
} // namespace eststack

#endif // !CORE__IRLS_HPP
