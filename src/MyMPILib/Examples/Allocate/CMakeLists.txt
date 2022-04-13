cmake_minimum_required(VERSION 2.6)
project (AllocateTest)

set(CMAKE_Fortran_COMPILER mpif90) 
enable_language(Fortran)

set(PROJLIBS /proj/plasma/Libs/)
### Directories of MyMPILib and MPE
set (MYMPILIB_PATH ${PROJLIBS}/MyMPILib/Build-Release)
set (MPE_PATH /afs/itp.tugraz.at/opt/mpe/1.3.0/)

include_directories("${MYMPILIB_PATH}/OBJS")
set(MPI_lib -lmpi_f90 -lmpi_f77 -lmpi -lopen-rte -lopen-pal -ldl -Wl,--export-dynamic -lnsl -lutil -lm -ldl -pthread)
find_library(MyMPILib_lib MyMPILib ${MYMPILIB_PATH})
 
add_executable(allocatetest
main.f90
)

target_link_libraries(allocatetest
        ${MyMPILib_lib}
        ${MPI_lib}
	/afs/itp.tugraz.at/opt/mpe/1.3.0/lib/libmpe.a
)

