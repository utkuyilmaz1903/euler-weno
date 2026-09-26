cmake_minimum_required(VERSION 3.25)

project(euler_weno
    VERSION 0.1.0
    DESCRIPTION "1D compressible Euler solver: finite volume, WENO5, SSP-RK3"
    LANGUAGES CXX)

if(PROJECT_SOURCE_DIR STREQUAL PROJECT_BINARY_DIR)
    message(FATAL_ERROR "In-source builds are not allowed. Use a preset, e.g. `cmake --preset debug`.")
endif()

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)          # -std=c++20, not -std=gnu++20
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)  # for clangd

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    set(CMAKE_BUILD_TYPE Debug CACHE STRING "Build type" FORCE)
endif()

option(EULER_ENABLE_SANITIZERS "Enable AddressSanitizer + UBSan in Debug builds" ON)
option(EULER_BUILD_TESTS "Build the unit tests" ON)

include(cmake/ProjectOptions.cmake)

add_subdirectory(src)
add_subdirectory(app)

if(EULER_BUILD_TESTS)
    enable_testing()
    add_subdirectory(tests)
endif()
