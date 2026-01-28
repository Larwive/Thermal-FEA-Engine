#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <sstream>

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

TEST_CASE("Vector size and element access", "[linalg]") {
  thermal_fea::linalg::Vector v(5);

  REQUIRE(v.size() == 5);

  for (std::size_t i = 0; i < v.size(); ++i) {
    REQUIRE(v[i] == Catch::Approx(0.0));
  }

  v(0) = 1.5;
  v(4) = -2.5;

  REQUIRE(v[0] == Catch::Approx(1.5));
  REQUIRE(v[4] == Catch::Approx(-2.5));

  const auto &cv = v;
  REQUIRE(cv(0) == Catch::Approx(1.5));
  REQUIRE(cv(4) == Catch::Approx(-2.5));
}

TEST_CASE("Vector addition", "[linalg]") {
  thermal_fea::linalg::Vector a(3);
  thermal_fea::linalg::Vector b(3);

  a[0] = 1.0;
  a[1] = 2.0;
  a[2] = 3.0;
  b[0] = 4.0;
  b[1] = -1.0;
  b[2] = 0.5;

  auto c = a + b;

  REQUIRE(c.size() == 3);
  REQUIRE(c[0] == Catch::Approx(5.0));
  REQUIRE(c[1] == Catch::Approx(1.0));
  REQUIRE(c[2] == Catch::Approx(3.5));
}

TEST_CASE("Vector subtraction", "[linalg]") {
  thermal_fea::linalg::Vector a(3);
  thermal_fea::linalg::Vector b(3);

  a[0] = 5.0;
  a[1] = 2.0;
  a[2] = -1.0;
  b[0] = 1.0;
  b[1] = 3.0;
  b[2] = -4.0;

  auto c = a - b;

  REQUIRE(c[0] == Catch::Approx(4.0));
  REQUIRE(c[1] == Catch::Approx(-1.0));
  REQUIRE(c[2] == Catch::Approx(3.0));
}

TEST_CASE("Scalar multiplication", "[linalg]") {
  thermal_fea::linalg::Vector v(4);
  v[0] = 1.0;
  v[1] = -2.0;
  v[2] = 0.5;
  v[3] = 4.0;

  auto w = 2.0 * v;

  REQUIRE(w[0] == Catch::Approx(2.0));
  REQUIRE(w[1] == Catch::Approx(-4.0));
  REQUIRE(w[2] == Catch::Approx(1.0));
  REQUIRE(w[3] == Catch::Approx(8.0));
}

TEST_CASE("Vector addition assignment", "[linalg]") {
  thermal_fea::linalg::Vector a(3);
  thermal_fea::linalg::Vector b(3);

  a[0] = 1.0;
  a[1] = 1.0;
  a[2] = 1.0;
  b[0] = 2.0;
  b[1] = 3.0;
  b[2] = 4.0;

  auto &ref = (a += b);

  REQUIRE(&ref == &a);
  REQUIRE(a[0] == Catch::Approx(3.0));
  REQUIRE(a[1] == Catch::Approx(4.0));
  REQUIRE(a[2] == Catch::Approx(5.0));
}

TEST_CASE("Dot product properties", "[linalg]") {
  thermal_fea::linalg::Vector a(3);
  thermal_fea::linalg::Vector b(3);

  a[0] = 1.0;
  a[1] = 2.0;
  a[2] = 3.0;
  b[0] = 4.0;
  b[1] = -5.0;
  b[2] = 6.0;

  REQUIRE(a.dot(b) == Catch::Approx(b.dot(a)));
  REQUIRE(a.dot(a) == Catch::Approx(a.norm() * a.norm()));
}

TEST_CASE("Norm homogeneity", "[linalg]") {
  thermal_fea::linalg::Vector v(3);
  v[0] = 1.0;
  v[1] = -2.0;
  v[2] = 2.0;

  double alpha = -3.0;
  auto w = alpha * v;

  REQUIRE(w.norm() == Catch::Approx(std::abs(alpha) * v.norm()));
}

TEST_CASE("Output stream operator", "[linalg]") {
  thermal_fea::linalg::Vector v(3);
  v[0] = 1.0;
  v[1] = 2.5;
  v[2] = -4.0;

  std::ostringstream oss;
  oss << v;

  REQUIRE_FALSE(oss.str().empty());
}
