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
  add_test(NAME smoke
    COMMAND gep4_fpp_sim macros/run_example.mac data/example_events.csv smoke.root)
  set_tests_properties(smoke PROPERTIES TIMEOUT 120)
endif()
