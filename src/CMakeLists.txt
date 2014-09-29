 cmake_minimum_required (VERSION 2.6)
 project (Test)

# version NUmber
set (Test_VERSION_MAJOR 1)
set (Test_VERSION_MINOR 0)
add_definitions(-DHAVE_VIENNA20=1)

set (GTEST_ROOT ~/local/gtest-1.7.0/)

configure_file (
  "${PROJECT_SOURCE_DIR}/TestConfig.h.in"
  "${PROJECT_BINARY_DIR}/TestConfig.h"
 )

include_directories("${PROJECT_BINARY_DIR}")
INCLUDE_DIRECTORIES(.)

add_library(DUPLEXHOM duplexhom.cpp)

include_directories("${PROJECT_SOURCE_DIR}/centroidalifold")
INCLUDE_DIRECTORIES("${PROJECT_SOURCE_DIR}/centroidalifold/probconsRNA")
INCLUDE_DIRECTORIES(/usr/local/Cellar/viennarna/2.1.8/include)
add_subdirectory(centroidalifold)

add_subdirectory(centroidalifold/probconsRNA)

LINK_DIRECTORIES(/usr/local/Cellar/viennarna/2.1.8/lib)

add_subdirectory(gtest)
INCLUDE_DIRECTORIES(gtest/include gtest)

add_executable(Test test.cpp pf_duplex.c)
target_link_libraries(Test centroidalifold Probcons gtest gtest_main DUPLEXHOM RNA)

option(USE_CENTROIDALIFOLD
	"Use test provided centroidalifold implementation" ON)

