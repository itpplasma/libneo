cmake_minimum_required(VERSION 3.10)

add_library(test_util SHARED test_util.f90)

add_library(test_libneo_field_util SHARED test_libneo_field_util.f90)
target_link_libraries(test_libneo_field_util PRIVATE libneo_field_base)

add_executable(test_libneo_field.x test_libneo_field.f90)
target_link_libraries(test_libneo_field.x PRIVATE libneo_field
                                                  libneo_field_base
                                                  test_util)

add_executable(test_libneo_example_field.x test_libneo_example_field.f90)
target_link_libraries(test_libneo_example_field.x PRIVATE libneo_example_field
                                                          test_libneo_field_util
                                                          test_util)

add_executable(test_libneo_biotsavart_field.x test_libneo_biotsavart_field.f90)
if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
  set(CMAKE_INSTALL_PREFIX $ENV{CODE}/local CACHE PATH "Install path prefix." FORCE)
endif()
include_directories(${CMAKE_INSTALL_PREFIX}/include)
set(LOCAL_LIB "${CMAKE_INSTALL_PREFIX}/lib")
message(STATUS "LOCAL_LIB = ${LOCAL_LIB}")
find_library(LIBBIOTSAVART biotsavart PATHS ${LOCAL_LIB})
target_link_libraries(test_libneo_biotsavart_field.x PRIVATE libneo_biotsavart_field
                                                             ${LIBBIOTSAVART}
                                                             test_libneo_field_util
                                                             test_util)
