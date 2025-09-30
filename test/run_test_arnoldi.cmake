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
  OUTPUT_VARIABLE setup_stdout
  ERROR_VARIABLE setup_stderr
)
if(NOT setup_status EQUAL 0)
  message(STATUS "setup_test_arnoldi.py stdout: ${setup_stdout}")
  message(STATUS "setup_test_arnoldi.py stderr: ${setup_stderr}")
  message(FATAL_ERROR "setup_test_arnoldi.py failed with status ${setup_status}")
endif()

execute_process(
  COMMAND "${ARNOLDI_EXE}"
  WORKING_DIRECTORY "${WORKDIR}"
  RESULT_VARIABLE arnoldi_status
  OUTPUT_VARIABLE arnoldi_stdout
  ERROR_VARIABLE arnoldi_stderr
)
if(NOT arnoldi_status EQUAL 0)
  message(STATUS "test_arnoldi.x stdout: ${arnoldi_stdout}")
  message(STATUS "test_arnoldi.x stderr: ${arnoldi_stderr}")
  message(FATAL_ERROR "test_arnoldi.x failed with status ${arnoldi_status}")
endif()
