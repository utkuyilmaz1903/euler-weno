include(FetchContent)

# Pinned by version and hash so every build (local and CI) gets identical sources.
# SYSTEM marks Catch2's headers as system headers: their warnings don't hit our -Werror.
FetchContent_Declare(Catch2
    URL      https://github.com/catchorg/Catch2/archive/refs/tags/v3.16.0.tar.gz
    URL_HASH SHA256=0957cae5821b17ce07f0833aaa52b5137643a8382203221f363a8303c109af34
    SYSTEM)
FetchContent_MakeAvailable(Catch2)

list(APPEND CMAKE_MODULE_PATH ${catch2_SOURCE_DIR}/extras)
include(Catch)

add_executable(euler_tests
    test_smoke.cpp
    test_eos.cpp
    test_numerical_flux.cpp
    test_finite_volume.cpp)

target_link_libraries(euler_tests PRIVATE
    euler::core
    euler_project_options
    Catch2::Catch2WithMain)

# Registers every TEST_CASE as its own CTest test.
catch_discover_tests(euler_tests)
