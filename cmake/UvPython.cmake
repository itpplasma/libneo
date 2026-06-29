# Resolve Python from the uv-managed project virtualenv so the test scripts and
# f2py see the project's Python dependencies (map2disc, netCDF4, scipy, ...).
# Without this, find_package(Python) picks the system interpreter, which lacks
# them and makes every chartmap/map2disc test fail with ModuleNotFoundError.
#
# Only for a standalone libneo build: when libneo is a subproject (e.g. fetched
# by SIMPLE) the parent controls Python and runs no libneo Python tests.

find_program(UV_EXECUTABLE uv)
set(_libneo_venv "${PROJECT_SOURCE_DIR}/.venv")

if(UV_EXECUTABLE)
    message(STATUS "uv: ${UV_EXECUTABLE} -- syncing project venv (chartmap extra)")
    execute_process(
        COMMAND "${UV_EXECUTABLE}" sync --extra chartmap
        WORKING_DIRECTORY "${PROJECT_SOURCE_DIR}"
        RESULT_VARIABLE _uv_rc
        OUTPUT_VARIABLE _uv_out
        ERROR_VARIABLE _uv_out)
    if(NOT _uv_rc EQUAL 0)
        message(WARNING "uv sync failed; chartmap/map2disc tests may not run:\n${_uv_out}")
    endif()
else()
    message(STATUS "uv not found; install uv or run 'uv sync --extra chartmap' "
                   "for chartmap/map2disc tests")
endif()

# Point find_package(Python) at the venv interpreter. Set Python_EXECUTABLE
# (FORCE) so an existing build dir picks it up instead of a cached system
# interpreter; Development/NumPy still resolve from the venv's base interpreter.
if(EXISTS "${_libneo_venv}/bin/python")
    set(ENV{VIRTUAL_ENV} "${_libneo_venv}")
    set(Python_FIND_VIRTUALENV FIRST)
    set(Python_ROOT_DIR "${_libneo_venv}" CACHE PATH "libneo uv venv" FORCE)
    set(Python_EXECUTABLE "${_libneo_venv}/bin/python" CACHE FILEPATH
        "libneo uv venv interpreter" FORCE)
    message(STATUS "libneo Python venv: ${_libneo_venv}")
endif()
