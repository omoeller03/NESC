add_executable(test_physics
    test_physics.cpp
)

target_link_libraries(test_physics
    PRIVATE
        nesc_physics
)

add_test(
    NAME physics_tests
    COMMAND test_physics
)
