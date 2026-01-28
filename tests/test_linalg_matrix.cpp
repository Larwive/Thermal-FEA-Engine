#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include "thermal_fea/linalg/Matrix.hpp"
#include "thermal_fea/linalg/Vector.hpp"

using thermal_fea::linalg::Matrix;
using thermal_fea::linalg::Vector;

TEST_CASE("Matrix rows and cols", "[linalg]") {
  Matrix A(3, 4);

  REQUIRE(A.rows() == 3);
  REQUIRE(A.cols() == 4);
}

TEST_CASE("Matrix element access", "[linalg]") {
  Matrix A(2, 3);

  A(0, 0) = 1.0;
  A(0, 1) = 2.0;
  A(1, 2) = -4.5;

  REQUIRE(A(0, 0) == Catch::Approx(1.0));
  REQUIRE(A(0, 1) == Catch::Approx(2.0));
  REQUIRE(A(1, 2) == Catch::Approx(-4.5));

  const Matrix &B = A;
  REQUIRE(B(1, 2) == Catch::Approx(-4.5));
}

TEST_CASE("Matrix-vector multiplication basic", "[linalg]") {
  Matrix A(2, 3);
  Vector x(3);

  A(0, 0) = 1.0;
  A(0, 1) = 2.0;
  A(0, 2) = 3.0;
  A(1, 0) = -1.0;
  A(1, 1) = 0.0;
  A(1, 2) = 4.0;

  x[0] = 1.0;
  x[1] = -1.0;
  x[2] = 2.0;

  Vector y = A * x;

  REQUIRE(y.size() == 2);
  REQUIRE(y[0] == Catch::Approx(1.0 * 1.0 + 2.0 * (-1.0) + 3.0 * 2.0));
  REQUIRE(y[1] == Catch::Approx(-1.0 * 1.0 + 0.0 * (-1.0) + 4.0 * 2.0));
}

TEST_CASE("Matrix-vector multiplication with zero vector", "[linalg]") {
  Matrix A(3, 3);
  Vector x(3);

  for (std::size_t i = 0; i < 3; ++i)
    for (std::size_t j = 0; j < 3; ++j)
      A(i, j) = static_cast<double>(i + j + 1);

  Vector y = A * x;

  for (std::size_t i = 0; i < y.size(); ++i) {
    REQUIRE(y[i] == Catch::Approx(0.0));
  }
}

TEST_CASE("Scalar-matrix multiplication", "[linalg]") {
  Matrix A(2, 2);

  A(0, 0) = 1.0;
  A(0, 1) = -2.0;
  A(1, 0) = 3.5;
  A(1, 1) = 4.0;

  Matrix B = 2.0 * A;

  REQUIRE(B(0, 0) == Catch::Approx(2.0));
  REQUIRE(B(0, 1) == Catch::Approx(-4.0));
  REQUIRE(B(1, 0) == Catch::Approx(7.0));
  REQUIRE(B(1, 1) == Catch::Approx(8.0));
}

TEST_CASE("Linearity of matrix-vector product", "[linalg]") {
  Matrix A(2, 2);
  Vector x(2);

  A(0, 0) = 1.0;
  A(0, 1) = 2.0;
  A(1, 0) = 3.0;
  A(1, 1) = 4.0;

  x[0] = -1.0;
  x[1] = 2.0;

  Vector y1 = A * x;
  Vector y2 = (2.0 * A) * x;

  REQUIRE(y2[0] == Catch::Approx(2.0 * y1[0]));
  REQUIRE(y2[1] == Catch::Approx(2.0 * y1[1]));
}
