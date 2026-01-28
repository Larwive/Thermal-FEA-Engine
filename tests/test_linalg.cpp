#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>

#include "thermal_fea/linalg/Vector.hpp"

TEST_CASE("Vector norm", "[linalg]") {
    thermal_fea::linalg::Vector v(3);
    v[0] = 3.0;
    v[1] = 4.0;
    v[2] = 0.0;

    REQUIRE(v.norm() == Catch::Approx(5.0));
    
    thermal_fea::linalg::Vector w(3);
    w[0] = 0.0;
    w[1] = 0.0;
    w[2] = 0.0;

    REQUIRE(w.norm() == Catch::Approx(0.0));
    
    thermal_fea::linalg::Vector x(3);
    x[0] = 1.0;
    x[1] = -1.0;
    x[2] = 5.0;

    REQUIRE(x.norm() == Catch::Approx(std::sqrt(27.0)));
}

TEST_CASE("Dot product", "[linalg]") {
    thermal_fea::linalg::Vector v(3);
    v[0] = 3.0;
    v[1] = 4.0;
    v[2] = 0.0;
    
    thermal_fea::linalg::Vector w(3);
    w[0] = -3.0;
    w[1] = 456784782.0;
    w[2] = 0.0;

    REQUIRE(v.dot(w) == Catch::Approx(1827139119.0));
    
    thermal_fea::linalg::Vector x(3);
    x[0] = 1.0;
    x[1] = -1.0;
    x[2] = 5.0;

    REQUIRE(x.dot(w) == Catch::Approx(-456784785.0));
    
    thermal_fea::linalg::Vector y(3);
    y[0] = 1.0;
    y[1] = 1.0;
    y[2] = 1.0;

    REQUIRE(y.dot(w) == Catch::Approx(456784783.0));
    
    thermal_fea::linalg::Vector z(3);
    z[0] = 0.0;
    z[1] = 0.0;
    z[2] = 0.0;

    REQUIRE(z.dot(w) == Catch::Approx(0.0));
}
