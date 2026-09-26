add_executable(euler_weno
    main.cpp)

target_link_libraries(euler_weno PRIVATE
    euler::core
    euler_project_options)
