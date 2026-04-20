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

#ifndef TEST__UTILS__BENCHMARKER_HPP
#define TEST__UTILS__BENCHMARKER_HPP

// C++ standard library
#include <algorithm>
#include <chrono>
#include <iomanip>
#include <iostream>
#include <string_view>

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
        /***
         * @brief global sink to prevent dead-code elimination in benchmark loops
         * @details `inline volatile` gives a single, non-optimizable definition
         *          across test TUs.
         */
        inline volatile float g_sink = 0.0f;

        /***
         * @brief micro-benchmark helper with automatic unit selection
         * @param name    human-readable label
         * @param trials  number of iterations of `fn`
         * @param fn      callable invoked as `fn(int i)` per iteration
         * @details The callable is expected to push a value into `g_sink` so
         *          the optimizer can't remove the timed work.
         */
        template <typename F>
        inline void benchmark(std::string_view name, int trials, F &&fn)
        {
            /* warm-up: at most 1000 iterations (or all, for small `trials`) */
            const int warm = std::min(trials, 1000);
            for (int i = 0; i < warm; ++i)
                fn(i);

            const auto t0 = std::chrono::steady_clock::now();
            for (int i = 0; i < trials; ++i)
                fn(i);
            const auto t1 = std::chrono::steady_clock::now();

            const double ns = std::chrono::duration<double, std::nano>(t1 - t0).count();
            const double per = ns / trials;

            std::cout << "  " << std::left << std::setw(20) << name
                      << " total=" << std::fixed << std::setprecision(3) << ns / 1e6 << " ms"
                      << "   per-call="
                      << (per < 1e3 ? per : (per < 1e6 ? per / 1e3 : per / 1e6))
                      << (per < 1e3 ? " ns" : (per < 1e6 ? " us" : " ms")) << "\n";
        }
    } // namespace test
} // namespace eststack

#endif //! TEST__UTILS__BENCHMARKER_HPP
