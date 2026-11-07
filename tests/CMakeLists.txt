add_executable(
  collie_io_smoke
  collie_io_smoke.cc
)
target_link_libraries(collie_io_smoke PRIVATE CollieIO ROOT::Core)
collie_apply_compiler_options(collie_io_smoke)
add_test(NAME io.dictionary-smoke COMMAND collie_io_smoke)
set_tests_properties(
  io.dictionary-smoke
  PROPERTIES
    LABELS "smoke;io"
    TIMEOUT 60
)

add_executable(
  collie_statistical_core_smoke
  statistical_core_smoke.cc
)
target_link_libraries(collie_statistical_core_smoke PRIVATE CollieCore)
collie_apply_compiler_options(collie_statistical_core_smoke)
add_test(NAME statistics.core-smoke COMMAND collie_statistical_core_smoke)
set_tests_properties(
  statistics.core-smoke
  PROPERTIES
    LABELS "smoke;statistics"
    TIMEOUT 60
)

add_executable(
  collie_profile_likelihood_smoke
  profile_likelihood_smoke.cc
)
target_link_libraries(collie_profile_likelihood_smoke PRIVATE CollieLimit)
collie_apply_compiler_options(collie_profile_likelihood_smoke)
add_test(
  NAME statistics.profile-likelihood-smoke
  COMMAND collie_profile_likelihood_smoke
)
set_tests_properties(statistics.profile-likelihood-smoke
  PROPERTIES
    LABELS "smoke;statistics;root"
    TIMEOUT 120
)

add_executable(
  collie_clfast_parallel_reproducibility
  clfast_parallel_reproducibility.cc
)
target_link_libraries(
  collie_clfast_parallel_reproducibility
  PRIVATE CollieLimit
)
collie_apply_compiler_options(collie_clfast_parallel_reproducibility)
add_test(
  NAME statistics.clfast-parallel-reproducibility
  COMMAND collie_clfast_parallel_reproducibility
)
set_tests_properties(statistics.clfast-parallel-reproducibility
  PROPERTIES
    LABELS "statistics;parallel"
    TIMEOUT 120
)

add_test(NAME cli.help COMMAND collie --help)
set_tests_properties(cli.help PROPERTIES LABELS "smoke;cli")

add_test(
  NAME cli.rejects-unsafe-iterations
  COMMAND "${CMAKE_COMMAND}"
          "-DCOLLIE_EXECUTABLE=$<TARGET_FILE:collie>"
          "-DCLI_ARGUMENTS=limit|--output|ignored.root|--iterations|100000001|nonexistent.root"
          "-DEXPECTED_PATTERN=--iterations exceeds the safe trial-count limit"
          -P "${CMAKE_CURRENT_SOURCE_DIR}/assert_cli_failure.cmake"
)
set_tests_properties(cli.rejects-unsafe-iterations
  PROPERTIES
    LABELS "smoke;cli"
)

add_test(
  NAME cli.rejects-unsafe-workers
  COMMAND "${CMAKE_COMMAND}"
          "-DCOLLIE_EXECUTABLE=$<TARGET_FILE:collie>"
          "-DCLI_ARGUMENTS=limit|--output|ignored.root|--workers|257|nonexistent.root"
          "-DEXPECTED_PATTERN=--workers must not exceed 256"
          -P "${CMAKE_CURRENT_SOURCE_DIR}/assert_cli_failure.cmake"
)
set_tests_properties(cli.rejects-unsafe-workers
  PROPERTIES
    LABELS "smoke;cli"
)

add_executable(collie_test_fixture collie_fixture.cc)
target_link_libraries(collie_test_fixture PRIVATE CollieIO)
collie_apply_compiler_options(collie_test_fixture)

add_executable(collie_result_layout_smoke result_layout_smoke.cc)
target_link_libraries(collie_result_layout_smoke PRIVATE ROOT::Core ROOT::RIO ROOT::Tree)
collie_apply_compiler_options(collie_result_layout_smoke)

set(COLLIE_VALID_FIXTURE "${CMAKE_CURRENT_BINARY_DIR}/valid-input.root")
set(COLLIE_INVALID_FIXTURE "${CMAKE_CURRENT_BINARY_DIR}/invalid-input.root")
set(COLLIE_LIMIT_FIRST "${CMAKE_CURRENT_BINARY_DIR}/limit-first.root")
set(COLLIE_LIMIT_SECOND "${CMAKE_CURRENT_BINARY_DIR}/limit-second.root")
set(COLLIE_FIT_RESULT "${CMAKE_CURRENT_BINARY_DIR}/fit.root")
set(COLLIE_COMBINED_RESULT "${CMAKE_CURRENT_BINARY_DIR}/combined.root")

add_test(
  NAME fixture.valid-input
  COMMAND collie_test_fixture "${COLLIE_VALID_FIXTURE}"
)
set_tests_properties(fixture.valid-input
  PROPERTIES
    FIXTURES_SETUP valid_input
    LABELS "fixture;cli"
)

add_test(
  NAME fixture.invalid-input
  COMMAND collie_test_fixture "${COLLIE_INVALID_FIXTURE}" --nonfinite
)
set_tests_properties(fixture.invalid-input
  PROPERTIES
    FIXTURES_SETUP invalid_input
    LABELS "fixture;cli"
)

add_test(
  NAME cli.validate-valid-input
  COMMAND collie validate --strict "${COLLIE_VALID_FIXTURE}"
)
set_tests_properties(cli.validate-valid-input
  PROPERTIES
    FIXTURES_REQUIRED valid_input
    LABELS "integration;cli"
)

add_test(
  NAME cli.rejects-nonfinite-systematic
  COMMAND "${CMAKE_COMMAND}"
          "-DCOLLIE_EXECUTABLE=$<TARGET_FILE:collie>"
          "-DCLI_ARGUMENTS=validate|--strict|${COLLIE_INVALID_FIXTURE}"
          "-DEXPECTED_PATTERN=non-finite fractional effect"
          -P "${CMAKE_CURRENT_SOURCE_DIR}/assert_cli_failure.cmake"
)
set_tests_properties(cli.rejects-nonfinite-systematic
  PROPERTIES
    FIXTURES_REQUIRED invalid_input
    LABELS "integration;cli"
)

add_test(
  NAME cli.limit-first
  COMMAND collie limit --method fast --point 125 --precision 0
          --iterations 200 --accuracy 0.02 --seed 1234 --workers 2
          --force --output "${COLLIE_LIMIT_FIRST}" "${COLLIE_VALID_FIXTURE}"
)
set_tests_properties(cli.limit-first
  PROPERTIES
    FIXTURES_REQUIRED valid_input
    FIXTURES_SETUP limit_first
    LABELS "integration;cli"
    TIMEOUT 60
)

add_test(
  NAME cli.limit-second
  COMMAND collie limit --method fast --point 130 --precision 0
          --iterations 200 --accuracy 0.02 --seed 1234 --workers 2
          --force --output "${COLLIE_LIMIT_SECOND}" "${COLLIE_VALID_FIXTURE}"
)
set_tests_properties(cli.limit-second
  PROPERTIES
    FIXTURES_REQUIRED valid_input
    FIXTURES_SETUP limit_second
    LABELS "integration;cli"
    TIMEOUT 60
)

add_test(
  NAME cli.fit
  COMMAND collie fit --point 125 --seed 1234 --force
          --output "${COLLIE_FIT_RESULT}" "${COLLIE_VALID_FIXTURE}"
)
set_tests_properties(cli.fit
  PROPERTIES
    FIXTURES_REQUIRED valid_input
    FIXTURES_SETUP fit_result
    LABELS "integration;cli"
    TIMEOUT 120
)

add_test(
  NAME cli.combine-compatible
  COMMAND collie combine --force --output "${COLLIE_COMBINED_RESULT}"
          "${COLLIE_LIMIT_FIRST}" "${COLLIE_LIMIT_SECOND}"
)
set_tests_properties(cli.combine-compatible
  PROPERTIES
    FIXTURES_REQUIRED "limit_first;limit_second"
    FIXTURES_SETUP combined_result
    LABELS "integration;cli"
)

add_test(
  NAME cli.rejects-mixed-result-kinds
  COMMAND "${CMAKE_COMMAND}"
          "-DCOLLIE_EXECUTABLE=$<TARGET_FILE:collie>"
          "-DCLI_ARGUMENTS=combine|--output|${CMAKE_CURRENT_BINARY_DIR}/mixed.root|${COLLIE_LIMIT_FIRST}|${COLLIE_FIT_RESULT}"
          "-DEXPECTED_PATTERN=result command.*does not match"
          -P "${CMAKE_CURRENT_SOURCE_DIR}/assert_cli_failure.cmake"
)
set_tests_properties(cli.rejects-mixed-result-kinds
  PROPERTIES
    FIXTURES_REQUIRED "limit_first;fit_result"
    LABELS "integration;cli"
)

add_test(
  NAME cli.fit-layout
  COMMAND collie_result_layout_smoke "${COLLIE_FIT_RESULT}" fit
)
set_tests_properties(cli.fit-layout
  PROPERTIES
    FIXTURES_REQUIRED fit_result
    LABELS "integration;cli"
)

add_test(
  NAME cli.combined-layout
  COMMAND collie_result_layout_smoke "${COLLIE_COMBINED_RESULT}" combine
)
set_tests_properties(cli.combined-layout
  PROPERTIES
    FIXTURES_REQUIRED combined_result
    LABELS "integration;cli"
)
