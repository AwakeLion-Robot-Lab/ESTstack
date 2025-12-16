set_project("ESTstack")
set_description("An algorithm set focus on estimation and filtering.")
set_version("0.0.1")
set_xmakever("2.9.8")
set_license("Apache-2.0")

set_defaultplat("linux")
set_languages("c11", "c++20")

add_rules("mode.debug", "mode.release")
if is_mode("debug") then
    set_symbols("debug")
    set_optimize("none")
    add_cxflags("-g")
elseif is_mode("release") then
    set_optimize("fastest")
    add_cxflags("-march=native", "-mtune=native", "-flto=auto")
    add_cxflags("-w")
    add_ldflags("-flto=auto")
end

if has_config("test") then
    add_requires("gtest 1.17.0", {configs = {main = true}})
end
add_requires("eigen 5.0.0", "backward-cpp v1.6", "manif 0.0.5")
add_requireconfs("manif.eigen", {override = true})

namespace("fosu-awakelion")

    target("ESTstack")
        set_kind("headeronly")
        add_headerfiles("modules/eststack/core/*.cppm")
        add_headerfiles("modules/eststack/model/motion/*.cppm")
        add_headerfiles("modules/eststack/model/imu/*.cppm")
        add_headerfiles("modules/eststack/problem/*.cppm")
        add_headerfiles("modules/eststack/solution/*.cppm")
        add_headerfiles("modules/eststack/impl/*.cppm")

    if has_config("test") then
        for _, file in ipairs(os.files("test/*.cpp")) do
            local name = path.basename(file)
            target("ESTstack-test-" .. name)
                set_kind("binary")
                set_default(false)
                set_policy("build.c++.modules", true)
                add_files(file)
                add_packages("ESTstack")
                add_packages("gtest")
                set_rundir("$(projectdir)")
                add_tests("ESTstack-test", {runargs = {"--gtest_color=yes"}})
        end
    end

namespace_end() -- namespace fosu-awakelion