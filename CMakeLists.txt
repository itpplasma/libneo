# This file is intended for configuration of libneo.
# List of source files are in CMakeSources.in.

# Basic cmake settings/requirements
cmake_minimum_required(VERSION 3.0)
set(CMAKE_DISABLE_IN_SOURCE_BUILD ON)

# Basic project settings.
project(libneo)

# The version number.
set (libneo_VERSION_MAJOR 0)
set (libneo_VERSION_MINOR 1)

enable_language(Fortran)
enable_testing()

# Include file for the sources.
include(CMakeSources.in)

set(CMAKE_Fortran_FLAGS " -std=f2008")

# Main goal (of the created build file): creating the library file.
add_library(libneo SHARED ${SOURCE_FILES})

#set_target_properties(libneo PROPERTIES VERSION ${libneo_VERSION})
