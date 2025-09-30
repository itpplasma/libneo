cmake_minimum_required(VERSION 3.18)

if(NOT DEFINED PYTHON_EXECUTABLE)
  message(FATAL_ERROR "PYTHON_EXECUTABLE not set for run_test_arnoldi.cmake")
endif()
if(NOT DEFINED SETUP_SCRIPT)
  message(FATAL_ERROR "SETUP_SCRIPT not set for run_test_arnoldi.cmake")
endif()
if(NOT DEFINED ARNOLDI_EXE)
  message(FATAL_ERROR "ARNOLDI_EXE not set for run_test_arnoldi.cmake")
endif()
if(NOT DEFINED WORKDIR)
  message(FATAL_ERROR "WORKDIR not set for run_test_arnoldi.cmake")
endif()

execute_process(
  COMMAND "${PYTHON_EXECUTABLE}" "${SETUP_SCRIPT}"
  WORKING_DIRECTORY "${WORKDIR}"
  RESULT_VARIABLE setup_status
)
if(NOT setup_status EQUAL 0)
  message(FATAL_ERROR "setup_test_arnoldi.py failed with status ${setup_status}")
endif()

execute_process(
  COMMAND "${ARNOLDI_EXE}"
  WORKING_DIRECTORY "${WORKDIR}"
  RESULT_VARIABLE arnoldi_status
)
if(NOT arnoldi_status EQUAL 0)
  message(FATAL_ERROR "test_arnoldi.x failed with status ${arnoldi_status}")
endif()
