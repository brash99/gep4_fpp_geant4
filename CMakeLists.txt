cmake_minimum_required(VERSION 3.16)
project(gep4_fpp_geant4 VERSION 0.1.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

find_package(Geant4 REQUIRED ui_all vis_all)
include(${Geant4_USE_FILE})

file(GLOB FPP_SOURCES CONFIGURE_DEPENDS src/*.cc)
add_executable(gep4_fpp_sim ${FPP_SOURCES})
target_include_directories(gep4_fpp_sim PRIVATE include)
target_link_libraries(gep4_fpp_sim PRIVATE ${Geant4_LIBRARIES})

configure_file(macros/run_example.mac macros/run_example.mac COPYONLY)
configure_file(data/example_events.csv data/example_events.csv COPYONLY)

include(CTest)
if(BUILD_TESTING)
  find_package(Python3 COMPONENTS Interpreter REQUIRED)

  add_test(NAME smoke
    COMMAND gep4_fpp_sim macros/run_example.mac data/example_events.csv smoke.root)
  set_tests_properties(smoke PROPERTIES TIMEOUT 120)

  add_test(NAME smoke_single
    COMMAND gep4_fpp_sim macros/run_example.mac data/example_events.csv
      smoke_single.root --single-analyzer)
  set_tests_properties(smoke_single PROPERTIES TIMEOUT 120)
  add_test(NAME smoke_single_compact_90cm
    COMMAND gep4_fpp_sim macros/run_example.mac data/example_events.csv
      smoke_single_compact_90cm.root --single-analyzer
      --analyzer-thickness-cm 90 --compact-downstream-fpp)
  set_tests_properties(smoke_single_compact_90cm PROPERTIES TIMEOUT 120)
  add_test(NAME thickness_mechanism_integration
    COMMAND ${Python3_EXECUTABLE}
      ${CMAKE_SOURCE_DIR}/analysis/analyze_thickness_mechanisms.py
      --sample smoke 90 ${CMAKE_BINARY_DIR}/smoke_single_compact_90cm.root
      --output-stem ${CMAKE_BINARY_DIR}/smoke_mechanisms)
  set_tests_properties(thickness_mechanism_integration PROPERTIES
    DEPENDS smoke_single_compact_90cm)
  add_test(NAME track_selection_integration
    COMMAND ${Python3_EXECUTABLE} ${CMAKE_SOURCE_DIR}/analysis/compare_track_selection.py
      ${CMAKE_BINARY_DIR}/smoke_single.root ${CMAKE_BINARY_DIR}/selection_smoke)
  set_tests_properties(track_selection_integration PROPERTIES DEPENDS smoke_single)

  add_executable(fpp_digitizer_test
    tests/FPPDigitizerTest.cc
    src/FPPGeometry.cc
    src/FPPTrackFit.cc)
  target_include_directories(fpp_digitizer_test PRIVATE include)
  add_test(NAME fpp_digitizer COMMAND fpp_digitizer_test)

  add_test(NAME track_selection
    COMMAND ${Python3_EXECUTABLE} ${CMAKE_SOURCE_DIR}/tests/test_track_selection.py)

  add_test(NAME analyzing_power
    COMMAND ${Python3_EXECUTABLE} ${CMAKE_SOURCE_DIR}/tests/test_calculate_fom.py)
endif()
