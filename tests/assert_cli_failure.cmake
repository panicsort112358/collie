if(NOT DEFINED COLLIE_EXECUTABLE OR NOT DEFINED CLI_ARGUMENTS OR
   NOT DEFINED EXPECTED_PATTERN)
  message(FATAL_ERROR "missing CLI failure-test configuration")
endif()

string(REPLACE "|" ";" arguments "${CLI_ARGUMENTS}")
execute_process(
  COMMAND "${COLLIE_EXECUTABLE}" ${arguments}
  RESULT_VARIABLE status
  OUTPUT_VARIABLE standard_output
  ERROR_VARIABLE standard_error
)

if("${status}" STREQUAL "0")
  message(
    FATAL_ERROR
    "command unexpectedly succeeded:\n${standard_output}${standard_error}"
  )
endif()

set(output "${standard_output}\n${standard_error}")
if(NOT output MATCHES "${EXPECTED_PATTERN}")
  message(
    FATAL_ERROR
    "command failed for the wrong reason (status ${status}):\n${output}"
  )
endif()
