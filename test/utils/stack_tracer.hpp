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

#ifndef TEST__UTILS__STACK_TRACER_HPP
#define TEST__UTILS__STACK_TRACER_HPP

// C++ standard library
#include <iostream>

// GTest library
#include <gtest/gtest.h>

// Backward-Cpp library
#include <backward.hpp>

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
         * @brief GTest event listener that appends a backward-cpp stack trace
         *        after every failed assertion.
         * @details Goes through Printer + StackTrace — the core backward-cpp API —
         *          so this file doubles as the "how to use it on demand" example.
         */
        class BackwardStackDumpListener : public ::testing::EmptyTestEventListener
        {
        public:
            void OnTestPartResult(const ::testing::TestPartResult &r) override
            {
                if (r.passed())
                    return;

                std::cerr << "\n[backward-cpp] stack trace for "
                          << (r.fatally_failed() ? "fatal" : "non-fatal")
                          << " failure at "
                          << (r.file_name() ? r.file_name() : "<unknown>")
                          << ":" << r.line_number() << "\n";

                backward::StackTrace st;
                st.load_here(32);
                /* skip OnTestPartResult itself + gtest dispatch frames */
                st.skip_n_firsts(4);

                backward::Printer p;
                p.address = true;
                p.object = false;
                p.print(st, std::cerr);
            }
        };

        /***
         * @brief register the listener at static-init time; gtest owns the pointer.
         * @note  `inline` guarantees a single registration even when the header
         *        is pulled into multiple test TUs.
         */
        [[gnu::used]]
        inline const int _register_backward_listener = []()
        {
            ::testing::UnitTest::GetInstance()
                ->listeners()
                .Append(new BackwardStackDumpListener);
            return 0;
        }();
    } // namespace test
} // namespace eststack

#endif //! TEST__UTILS__STACK_TRACER_HPP
