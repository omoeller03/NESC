cmake_minimum_required(VERSION 3.16)

project(NESC
    VERSION 0.1.0
    DESCRIPTION "Nuclear Electron Scattering Calculator"
    LANGUAGES CXX
)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

add_library(nesc_physics
    src/Nucleus.cpp
    src/Kinematics.cpp
    src/FormFactor.cpp
    src/CrossSection.cpp
)

target_include_directories(nesc_physics
    PUBLIC
        ${PROJECT_SOURCE_DIR}/include
)

target_compile_options(nesc_physics
    PRIVATE
        -Wall
        -Wextra
        -Wpedantic
)

add_executable(nesc
    src/main.cpp
)

target_link_libraries(nesc
    PRIVATE
        nesc_physics
)

include(CTest)

if(BUILD_TESTING)
    add_subdirectory(tests)
endif()
