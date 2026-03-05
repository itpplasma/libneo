add_library(odeint STATIC
    odeint_allroutines.f90)

if(ENABLE_OPENMP)
    target_link_libraries(odeint PUBLIC OpenMP::OpenMP_Fortran)
endif()

set_property(TARGET odeint PROPERTY
    INTERFACE_INCLUDE_DIRECTORIES ${CMAKE_Fortran_MODULE_DIRECTORY})
add_library(LIBNEO::odeint ALIAS odeint)
