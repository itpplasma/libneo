add_library(interpolate OBJECT
    interpolate.f90)

add_library(LIBNEO::interpolate ALIAS interpolate)
