cmake_minimum_required(VERSION 3.20)

project(
  Collie
  VERSION 0.4.13
  DESCRIPTION "Confidence Level Limit Evaluator"
  LANGUAGES C CXX
)

include(GNUInstallDirs)
include(CTest)
include(CMakePackageConfigHelpers)

option(COLLIE_BUILD_EXAMPLES "Build the legacy Collie command-line examples" ON)
option(
  COLLIE_ENABLE_SANITIZERS
  "Enable AddressSanitizer and UndefinedBehaviorSanitizer"
  OFF
)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")
include(CollieCompilerOptions)

find_package(
  ROOT 6.20 CONFIG REQUIRED
  COMPONENTS Core RIO Hist Graf Gpad Tree MathCore Minuit2
)
find_package(Threads REQUIRED)

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/${CMAKE_INSTALL_LIBDIR}")
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/${CMAKE_INSTALL_LIBDIR}")
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/${CMAKE_INSTALL_BINDIR}")
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

if(APPLE)
  set(COLLIE_LIBRARY_INSTALL_RPATH "@loader_path")
  set(COLLIE_EXECUTABLE_INSTALL_RPATH
      "@loader_path/../${CMAKE_INSTALL_LIBDIR}")
else()
  set(COLLIE_LIBRARY_INSTALL_RPATH "$ORIGIN")
  set(COLLIE_EXECUTABLE_INSTALL_RPATH
      "$ORIGIN/../${CMAKE_INSTALL_LIBDIR}")
endif()

add_library(CollieCore SHARED
  limit/src/BinnedModel.cc
  limit/src/RunContext.cc
)
add_library(Collie::Core ALIAS CollieCore)
set_target_properties(CollieCore PROPERTIES
  EXPORT_NAME Core
  INSTALL_RPATH "${COLLIE_LIBRARY_INSTALL_RPATH}"
)
target_include_directories(CollieCore
  PUBLIC
    "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/limit/include>"
    "$<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>"
)
target_link_libraries(CollieCore PUBLIC Threads::Threads)
collie_apply_compiler_options(CollieCore)

set(COLLIE_IO_SOURCES
  io/src/CDF_DZero_Distribution.cc
  io/src/CDF_DZero_IOfile.cc
  io/src/CDF_DZero_IOpoint.cc
  io/src/CollieChannel.cc
  io/src/CollieDistribution.cc
  io/src/CollieEventList.cc
  io/src/CollieHistogram.cc
  io/src/CollieIOFile.cc
  io/src/CollieIterator.cc
  io/src/CollieMasspoint.cc
  io/src/CollieQuickKEYS.cc
  io/src/CollieXsec.cc
)
set(COLLIE_IO_HEADERS
  io/include/CDF_DZero_Distribution.hh
  io/include/CDF_DZero_IOfile.hh
  io/include/CDF_DZero_IOpoint.hh
  io/include/CollieChannel.hh
  io/include/CollieDistribution.hh
  io/include/CollieEventList.hh
  io/include/CollieHistogram.hh
  io/include/CollieIOFile.hh
  io/include/CollieIterator.hh
  io/include/CollieMasspoint.hh
  io/include/CollieQuickKEYS.hh
  io/include/CollieXsec.hh
  io/include/InfoForEfficiencyCalculation.hh
)

add_library(CollieIO SHARED ${COLLIE_IO_SOURCES})
add_library(Collie::IO ALIAS CollieIO)
set_target_properties(CollieIO PROPERTIES
  EXPORT_NAME IO
  INSTALL_RPATH "${COLLIE_LIBRARY_INSTALL_RPATH}"
)
target_include_directories(CollieIO
  PUBLIC
    "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/io/include>"
    "$<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>"
)
target_link_libraries(CollieIO
  PUBLIC
    ROOT::Core
    ROOT::RIO
    ROOT::Hist
    ROOT::Graf
    ROOT::Gpad
    ROOT::Tree
    ROOT::MathCore
)
collie_apply_compiler_options(CollieIO)

# ROOT_GENERATE_DICTIONARY invokes rootcling and keeps all generated sources,
# rootmaps, and PCM files in the build tree.
ROOT_GENERATE_DICTIONARY(
  G__CollieIO
  ${COLLIE_IO_HEADERS}
  MODULE CollieIO
)

set(COLLIE_LIMIT_SOURCES
  limit/src/CLcompute.cc
  limit/src/CLfast.cc
  limit/src/CLfit.cc
  limit/src/CLfit2.cc
  limit/src/CLfitSS.cc
  limit/src/CLpoint.cc
  limit/src/CLsyst.cc
  limit/src/CollieLoader.cc
  limit/src/CrossSectionCalc.cc
  limit/src/CrossSectionCalc2D.cc
  limit/src/CrossSectionLimit.cc
  limit/src/CrossSectionLimit_FastApprox.cc
  limit/src/ExclusionLimit.cc
  limit/src/FileSetLoader.cc
  limit/src/FitTest.cc
  limit/src/Loader.cc
  limit/src/ProfileLH.cc
  limit/src/ProfileLH_2D.cc
  limit/src/ProfileLikelihoodEngine.cc
  limit/src/SigBkgdDist.cc
  limit/src/SigBkgdDistAdapter.cc
  limit/src/ThreeSigmaEvidence.cc
  limit/src/timeBasedSeed.cc
)

add_library(CollieLimit SHARED ${COLLIE_LIMIT_SOURCES})
add_library(Collie::Limit ALIAS CollieLimit)
set_target_properties(CollieLimit PROPERTIES
  EXPORT_NAME Limit
  INSTALL_RPATH "${COLLIE_LIBRARY_INSTALL_RPATH}"
)
target_include_directories(CollieLimit
  PUBLIC
    "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/limit/include>"
    "$<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>"
)
target_link_libraries(CollieLimit
  PUBLIC
    CollieCore
    CollieIO
    ROOT::Core
    ROOT::RIO
    ROOT::Hist
    ROOT::Minuit2
    ROOT::Tree
)
collie_apply_compiler_options(CollieLimit)

ROOT_GENERATE_DICTIONARY(
  G__CollieLimit
  limit/include/CLpoint.hh
  MODULE CollieLimit
)

add_executable(collie
  tools/collie.cpp
  tools/CollieCli.cc
)
target_include_directories(collie PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/tools")
target_link_libraries(collie PRIVATE CollieLimit CollieIO)
set_target_properties(collie PROPERTIES
  INSTALL_RPATH "${COLLIE_EXECUTABLE_INSTALL_RPATH}"
)
collie_apply_compiler_options(collie)
install(TARGETS collie RUNTIME DESTINATION "${CMAKE_INSTALL_BINDIR}")

if(COLLIE_BUILD_EXAMPLES)
  function(collie_add_example target source output_name)
    add_executable("${target}" "${source}")
    set_target_properties("${target}" PROPERTIES
      OUTPUT_NAME "${output_name}"
      SUFFIX ".exe"
      INSTALL_RPATH "${COLLIE_EXECUTABLE_INSTALL_RPATH}"
    )
    target_link_libraries("${target}" PRIVATE
      CollieLimit
      CollieIO
    )
    collie_apply_compiler_options("${target}")
    install(TARGETS "${target}" RUNTIME DESTINATION "${CMAKE_INSTALL_BINDIR}")
  endfunction()

  collie_add_example(collieIOmaker
    examples/collieIOexample.cc collieIOmaker)
  collie_add_example(collieLimitCalc
    examples/exampleLimitCalculation.cc collieLimitCalc)
  collie_add_example(collieXsecCalc
    examples/exampleXsecMeasurement.cc collieXsecCalc)
  collie_add_example(resultsFileCombiner
    examples/resultsFileCombiner.cc resultsFileCombiner)
  collie_add_example(collieIOcondenser
    examples/collieIOcondenser.cc collieIOcondenser)
  collie_add_example(collieFastApprox
    examples/exampleFastApproxLimitCalculation.cc collieFastApprox)
  collie_add_example(collieXsecFitLoader
    examples/xsecFitLoader.cc collieXsecFitLoader)
  collie_add_example(collieIOExample examples/IO.cc IO)
  collie_add_example(collieSigmas examples/Sigmas.cc Sigmas)
  collie_add_example(collieCLs examples/CLs.cc CLs)
  collie_add_example(collieIOContour examples/IO_contour.cc IO_contour)
  collie_add_example(collieCLsb examples/CLsb.cc CLsb)
  collie_add_example(collieProfileLikelihood
    examples/profile_likelihood/profile_likelihood_collie_driver.cc
    collieProfileLikelihood)
endif()

if(BUILD_TESTING)
  add_subdirectory(tests)
endif()

install(
  TARGETS CollieCore CollieIO CollieLimit
  EXPORT CollieTargets
  ARCHIVE DESTINATION "${CMAKE_INSTALL_LIBDIR}"
  LIBRARY DESTINATION "${CMAKE_INSTALL_LIBDIR}"
  RUNTIME DESTINATION "${CMAKE_INSTALL_BINDIR}"
)
install(DIRECTORY io/include/ limit/include/
  DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}"
  FILES_MATCHING PATTERN "*.hh"
)
install(
  EXPORT CollieTargets
  FILE CollieTargets.cmake
  NAMESPACE Collie::
  DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/Collie"
)
configure_package_config_file(
  cmake/CollieConfig.cmake.in
  "${CMAKE_CURRENT_BINARY_DIR}/CollieConfig.cmake"
  INSTALL_DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/Collie"
)
write_basic_package_version_file(
  "${CMAKE_CURRENT_BINARY_DIR}/CollieConfigVersion.cmake"
  VERSION "${PROJECT_VERSION}"
  COMPATIBILITY SameMajorVersion
)
install(
  FILES
    "${CMAKE_CURRENT_BINARY_DIR}/CollieConfig.cmake"
    "${CMAKE_CURRENT_BINARY_DIR}/CollieConfigVersion.cmake"
  DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/Collie"
)
