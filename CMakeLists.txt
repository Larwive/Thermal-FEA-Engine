cmake_minimum_required(VERSION 3.16)

project(thermal_fea
    VERSION 0.1
    DESCRIPTION "Thermal finite element solver in C++"
    LANGUAGES CXX
)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

add_executable(thermal_fea
    src/main.cpp
    src/linalg/Matrix.cpp
    src/linalg/Vector.cpp
    src/mesh/Mesh.cpp
    src/solver/conjugate_gradient.cpp
    src/physics/heat_equation.cpp
)

target_compile_options(thermal_fea PRIVATE
    -Wall -Wextra -Wpedantic
)

target_include_directories(thermal_fea PUBLIC
    ${PROJECT_SOURCE_DIR}/include
)
