# Warning and sanitizer flags for our own targets only.
#
# Everything lives on an INTERFACE target instead of global CMAKE_CXX_FLAGS, so
# third-party code (Catch2) is not compiled with -Werror. Link it PRIVATE into
# every target we own.

add_library(euler_project_options INTERFACE)

target_compile_options(euler_project_options INTERFACE
    -Wall
    -Wextra
    -Wpedantic
    -Werror)

if(EULER_ENABLE_SANITIZERS)
    set(_euler_sanitizer_flags
        -fsanitize=address
        -fsanitize=undefined
        -fno-sanitize-recover=undefined  # abort on UB instead of printing and continuing
        -fno-omit-frame-pointer)         # readable stack traces

    target_compile_options(euler_project_options INTERFACE
        "$<$<CONFIG:Debug>:${_euler_sanitizer_flags}>")
    target_link_options(euler_project_options INTERFACE
        "$<$<CONFIG:Debug>:${_euler_sanitizer_flags}>")
endif()
