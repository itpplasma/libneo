cmake_minimum_required(VERSION 3.10)

add_library(biotsavart SHARED biotsavart.f90)
add_library(equalize_coils_segments SHARED equalize_coils_segments.f90)
target_link_libraries(equalize_coils_segments PRIVATE biotsavart)

install(TARGETS biotsavart DESTINATION ${CMAKE_INSTALL_PREFIX}/lib)

install(FILES ${CMAKE_BINARY_DIR}/include/biotsavart.mod 
    DESTINATION ${CMAKE_INSTALL_PREFIX}/include)
