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

find_package(OpenMP REQUIRED)

add_library(thermal_fea
    src/linalg/Matrix.cpp
    src/linalg/Vector.cpp
    src/mesh/Mesh.cpp
    src/physics/heat_equation.cpp
    src/solver/assembly.cpp
    src/solver/conjugate_gradient.cpp
)

set_target_properties(thermal_fea PROPERTIES
  POSITION_INDEPENDENT_CODE ON
)

target_compile_options(thermal_fea PRIVATE
    -Wall -Wextra -Wpedantic
)

target_link_libraries(thermal_fea PUBLIC OpenMP::OpenMP_CXX)

target_include_directories(thermal_fea
    PUBLIC
        $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
        $<INSTALL_INTERFACE:include>
)

add_executable(thermal_fea_demo
    src/main.cpp
)

target_link_libraries(thermal_fea_demo
    PRIVATE thermal_fea
)


include(GNUInstallDirs)

install(TARGETS thermal_fea
    EXPORT thermal_feaTargets
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
)

install(DIRECTORY include/
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
)

install(TARGETS thermal_fea_demo
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
)

include(CTest)
enable_testing()

add_subdirectory(tests)