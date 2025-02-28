cmake_minimum_required(VERSION 3.10)
project(Spath)

set(CMAKE_CXX_STANDARD 17)

add_library(Spath "src/spath.cpp" "src/converter.cpp")

include(GNUInstallDirs)

target_include_directories(Spath 
        PUBLIC
        "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>"
        "$<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>"
)

install(TARGETS Spath
        EXPORT Spath 
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
        INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
        CONFIGURATIONS Release
)

install(FILES include/spath.hpp include/edge.hpp include/converter.hpp DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/spath)

install(EXPORT Spath 
        FILE Spath.cmake
        NAMESPACE Spath::
        DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/Spath
)

include(CMakePackageConfigHelpers)

configure_package_config_file(${CMAKE_CURRENT_SOURCE_DIR}/Config.cmake.in
        "${CMAKE_CURRENT_BINARY_DIR}/SpathConfig.cmake"
        INSTALL_DESTINATION {CMAKE_INSTALL_LIBDIR}/cmake/Spath
)

install(FILES
        "${CMAKE_CURRENT_BINARY_DIR}/SpathConfig.cmake"     
        DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/Spath"
)


# Fetch Catch2 instead of including it as a directory
include(FetchContent)
FetchContent_Declare(
        Catch2
        GIT_REPOSITORY https://github.com/catchorg/Catch2.git
        GIT_TAG v3.8.0 
)
FetchContent_MakeAvailable(Catch2)

include(FetchContent)
FetchContent_Declare(
        jngen
        GIT_REPOSITORY https://github.com/ifsmirnov/jngen.git
        GIT_TAG master
)
FetchContent_MakeAvailable(jngen)

# Fetch Google Benchmark 
include(FetchContent)
FetchContent_Declare(
        googlebenchmark
        GIT_REPOSITORY https://github.com/google/benchmark.git
        GIT_TAG v1.9.1
)
set(BENCHMARK_ENABLE_GTEST_TESTS OFF)
set(BENCHMARK_DOWNLOAD_DEPENDENCIES ON)
set(CMAKE_BUILD_TYPE Release)

FetchContent_MakeAvailable(googlebenchmark)

add_executable(Bench "${CMAKE_CURRENT_SOURCE_DIR}/src/bench.cpp")
target_include_directories(Bench PRIVATE ${jngen_SOURCE_DIR}) # For #include "graph.h"
target_link_libraries(Bench PRIVATE benchmark::benchmark)
target_link_libraries(Bench PRIVATE Spath)



add_executable(Tests "${CMAKE_CURRENT_SOURCE_DIR}/src/tests.cpp")
target_link_libraries(Tests PRIVATE Catch2::Catch2WithMain)

target_include_directories(Tests 
        PRIVATE
        "${CMAKE_CURRENT_SOURCE_DIR}/include"
)

target_link_libraries(Tests PUBLIC Spath)