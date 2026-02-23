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

-- option for build ceres solver with cuda
option("with_cuda")
    set_default(false)
    set_showmenu(true)
    set_description("enable CUDA for ceres solver.")
option_end()

if has_config("test") then
    add_requires("gtest 1.17.0", {configs = {main = true}})
end
add_repositories("awakelion-xmake-repo https://github.com/AwakeLion-Robot-Lab/awakelion-xmake-repo.git")
add_requires("eigen 5.0.0", "backward-cpp v1.6", "manif 0.0.5", "awakelion-logger 1.0.2")
add_requires("openblas 0.3.30")  -- openblas for suitesparse inside ceres solver
add_requires("ceres-solver 2.2.0", {configs = {cuda = has_config("with_cuda")}})
add_requireconfs("manif.eigen", {override = true})  -- use eigen from xmake package
add_requireconfs("ceres-solver.openblas", {override = true})   -- use openblas from xmake package
add_requireconfs("ceres-solver.suitesparse.openblas", {override = true})  -- same as above

namespace("fosu-awakelion")
    target("ESTstack")
        set_kind("headeronly")
        add_includedirs("include", {public = true})
        add_headerfiles("include/eststack/**/*.hpp", {public = true})

        -- dependencies
        add_packages("eigen", "backward-cpp", "manif", "awakelion-logger", "openblas", "ceres-solver", {public = true})

    if has_config("test") then
        for _, file in ipairs(os.files("test/*.cpp")) do
            local name = path.basename(file)
            target("ESTstack-test-" .. name)
                set_kind("binary")
                set_default(false)
                add_files(file)
                add_deps("ESTstack")
                add_packages("gtest")
                set_rundir("$(projectdir)")
                add_tests("ESTstack-test", {runargs = {"--gtest_color=yes"}})
        end
    end
namespace_end() -- namespace fosu-awakelion