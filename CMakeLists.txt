cmake_minimum_required(VERSION 3.0 FATAL_ERROR)
project(fPotencia VERSION 1.0.0 LANGUAGES CXX)

include(FindPkgConfig)
include(GNUInstallDirs)

find_package(GTest)
find_package(Boost REQUIRED COMPONENTS graph)
pkg_check_modules(eigen3 REQUIRED eigen3>=3.2.8)


option(BUILD_UNSUPPORTED_SOLVERS
    "Build currently unsupported (i.e., test-failing) Load Flow solvers"
    OFF)


configure_file("config.h.in" "config.h")

add_subdirectory(src)

if(${GTEST_FOUND})
    enable_testing()
    add_subdirectory(test)
endif()
