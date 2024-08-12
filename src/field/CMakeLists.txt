cmake_minimum_required(VERSION 3.10)

add_library(libneo_field SHARED
    field_base.f90
    example_field.f90
    biotsavart.f90
    biotsavart_field.f90
    field.f90)