add_library(neo_polylag OBJECT
    polylag_5.f90)

add_library(LIBNEO::polylag ALIAS neo_polylag)
