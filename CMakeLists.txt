cmake_minimum_required(VERSION 3.0 FATAL_ERROR)
project(fPotencia VERSION 1.0.0 LANGUAGES CXX)

include(FindPkgConfig)
include(GNUInstallDirs)

find_package(GTest)
find_package(Boost REQUIRED COMPONENTS graph)
pkg_check_modules(REQUIRED eigen3)

add_subdirectory(src)

if(${GTEST_FOUND})
    enable_testing()
    add_subdirectory(test)
endif()
