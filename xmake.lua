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
    add_cxflags("-march=native")
    add_cxflags("-w")
end

add_repositories("awakelion-xmake-repo https://github.com/AwakeLion-Robot-Lab/awakelion-xmake-repo.git")
add_requires("eigen 5.0.0", "awakelion-logger 1.0.0", "backward-cpp v1.6", "manif 0.0.5")
add_requireconfs("manif.eigen", {override = true})

target("test")
    set_kind("binary")
    add_files("src/*.cpp")

    add_packages("awakelion-logger")