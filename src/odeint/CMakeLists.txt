add_library(odeint STATIC
    odeint_allroutines.f90)

if(ENABLE_OPENMP)
    find_package(OpenMP COMPONENTS Fortran)
    if(OpenMP_Fortran_FOUND)
        target_link_libraries(odeint PUBLIC OpenMP::OpenMP_Fortran)
    endif()
endif()

set_target_properties(odeint PROPERTIES
      Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
target_include_directories(odeint PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}>)
add_library(LIBNEO::odeint ALIAS odeint)
