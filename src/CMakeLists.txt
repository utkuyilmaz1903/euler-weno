add_library(euler_core
    version.cpp
    eos.cpp
    numerical_flux.cpp)
add_library(euler::core ALIAS euler_core)

target_include_directories(euler_core PUBLIC ${PROJECT_SOURCE_DIR}/include)
target_link_libraries(euler_core PRIVATE euler_project_options)

# Only version.cpp sees this macro; everyone else calls euler::version().
target_compile_definitions(euler_core PRIVATE
    EULER_VERSION_STRING="${PROJECT_VERSION}")
