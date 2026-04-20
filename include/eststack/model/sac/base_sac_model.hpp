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

#ifndef MODEL__BASE_SAC_MODEL_HPP
#define MODEL__BASE_SAC_MODEL_HPP

// C++ standard library
#include <memory>

// Eigen library
#include <Eigen/Core>

/***
 * @brief An algorithm set focus on state estimation
 * @author jinhua "siyiovo" deng
 */
namespace eststack
{
    /***
     * @brief models for problems
     */
    namespace model
    {
        /***
         * @brief base class for SAmple Consensus model
         * @tparam Derived derived class of the base SAC model
         */
        template <typename Derived>
        class BaseSACModel
        {
        public:
            /***
             * @brief set input source points
             * @param source 3xN source points
             * @details the source points can be either source point cloud or source points coordinates of correspondences
             */
            void setInputSource(const Eigen::Matrix3Xf &source)
            {
                source_ = source;
            }

            /***
             * @brief set input target points
             * @param target 3xN target points
             * @details the target points can be either target point cloud or target points coordinates of correspondences
             */
            void setInputTarget(const Eigen::Matrix3Xf &target)
            {
                target_ = target;
            }

            /***
             * @brief set inlier threshold
             * @param threshold inlier distance threshold
             */
            void setInlierThreshold(float threshold) noexcept
            {
                threshold_ = threshold;
            }

            /***
             * @brief get cloud size
             */
            int getCloudSize() const
            {
                const int source_size = source_.cols();
                if (source_size != target_.cols())
                    return 0;
                return source_size;
            }

            /***
             * @brief get sample size
             */
            int getSampleSize() const noexcept
            {
                return sample_size_;
            }

            /***
             * @brief get model size
             */
            int getModelSize() const noexcept
            {
                return model_size_;
            }

            /***
             * @brief run model fitting
             * @param[in] samples sample indices
             * @param[out] coeffs model coefficients
             */
            bool fit(const std::vector<int> &samples, Eigen::VectorXf &coeffs)
            {
                return static_cast<Derived *>(this)->fitImpl(samples, coeffs);
            }

            /***
             * @brief select inliers
             * @param[in] coeffs model coefficients
             */
            std::vector<int> selectInliers(const Eigen::VectorXf &coeffs) const
            {
                return static_cast<const Derived *>(this)->selectInliersImpl(coeffs);
            }

        protected:
            /***
             * @brief default constructor
             */
            BaseSACModel() = default;

            /***
             * @brief 3xN source points
             * @details the source points can be either source point cloud or source points coordinates of correspondences
             */
            Eigen::Matrix3Xf source_;

            /***
             * @brief 3xN target points
             * @details the target points can be either target point cloud or target points coordinates of correspondences
             */
            Eigen::Matrix3Xf target_;

            /***
             * @brief inlier distance threshold
             */
            float threshold_{0.0f};

            /***
             * @brief minimum sample size to fit a model
             */
            int sample_size_{0};

            /***
             * @brief number of parameters in the model
             */
            int model_size_{0};
        };
    } // namespace model
} // namespace eststack

#endif //! MODEL__BASE_SAC_MODEL_HPP
