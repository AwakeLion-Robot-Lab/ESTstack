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
    add_cxflags("-march=native", "-flto=auto", "-w")
    add_ldflags("-flto=auto")
end

if has_config("test") then
    add_requires("gtest 1.17.0", {configs = {main = true}})
end
add_requires("eigen 5.0.0", "pcl 1.15.1", "manif 0.0.5")
add_requires("backward-cpp v1.6")
add_requires("tbb", {system = true})                                          
add_requireconfs("manif.eigen", {override = true})  -- use eigen from xmake package

namespace("fosu-awakelion")
    target("ESTstack")
        set_kind("headeronly")
        add_includedirs("include", {public = true})
        add_headerfiles("include/eststack/**/*.hpp", {public = true})

        -- dependencies
        add_packages("eigen", "pcl", "manif", "tbb", {public = true})

    if has_config("test") then
        for _, file in ipairs(os.files("test/*.cpp")) do
            local name = path.basename(file)
            target("ESTstack-test-" .. name)
                set_kind("binary")
                set_default(false)
                add_files(file)
                add_deps("ESTstack")
                add_packages("gtest", "backward-cpp")
                set_rundir("$(projectdir)")
                add_tests("ESTstack-test", {runargs = {"--gtest_color=yes"}})
        end
    end
namespace_end() -- namespace fosu-awakelion