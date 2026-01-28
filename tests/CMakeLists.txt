option(THERMAL_FEA_USE_SYSTEM_CATCH2
       "Use system-provided Catch2 (Guix, distro packages)"
       OFF)

if (THERMAL_FEA_USE_SYSTEM_CATCH2)
    find_package(Catch2 3 REQUIRED)
else()
    include(FetchContent)
    FetchContent_Declare(
        Catch2
        GIT_REPOSITORY https://github.com/catchorg/Catch2.git
        GIT_TAG v3.5.4
    )
    FetchContent_MakeAvailable(Catch2)
endif()

add_executable(tests_thermal_fea
    test_linalg_vector.cpp
    test_linalg_matrix.cpp
    test_mesh.cpp
    test_solver.cpp
)

target_link_libraries(tests_thermal_fea
    PRIVATE
    Catch2::Catch2WithMain
    thermal_fea
)

add_test(
  NAME thermal_fea_tests
  COMMAND $<TARGET_FILE:tests_thermal_fea>
)
