#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>

#include "thermal_fea/linalg/Vector.hpp"

TEST_CASE("Vector norm", "[linalg]") {
    thermal_fea::linalg::Vector v(3);
    v[0] = 3.0;
    v[1] = 4.0;
    v[2] = 0.0;

    REQUIRE(v.norm() == Catch::Approx(5.0));
}