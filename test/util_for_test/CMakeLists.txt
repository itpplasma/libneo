cmake_minimum_required(VERSION 3.10)
add_library(util_for_test OBJECT util_for_test.f90)
target_include_directories(util_for_test PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})