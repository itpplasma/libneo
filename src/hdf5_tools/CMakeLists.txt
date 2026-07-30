add_library(hdf5_tools INTERFACE)
target_link_libraries(hdf5_tools INTERFACE fortio::fortio)

add_library(LIBNEO::hdf5_tools ALIAS hdf5_tools)
