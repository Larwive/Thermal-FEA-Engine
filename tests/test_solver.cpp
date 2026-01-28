#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include "thermal_fea/linalg/Matrix.hpp"
#include "thermal_fea/linalg/Vector.hpp"
#include "thermal_fea/solver/conjugate_gradient.hpp"

using thermal_fea::linalg::Matrix;
using thermal_fea::linalg::Vector;
using thermal_fea::solver::conjugate_gradient;

TEST_CASE("Conjugate gradient solves 1x1 system", "[solver][cg]") {
  Matrix A(1, 1);
  Vector b(1);

  A(0, 0) = 2.0;
  b[0] = 4.0;

  Vector x = conjugate_gradient(A, b);

  REQUIRE(x.size() == 1);
  REQUIRE(x[0] == Catch::Approx(2.0).epsilon(1e-8));
}

TEST_CASE("Conjugate gradient solves diagonal system", "[solver][cg]") {
  Matrix A(3, 3);
  Vector b(3);

  A(0, 0) = 2.0;
  A(0, 1) = 0.0;
  A(0, 2) = 0.0;
  A(1, 0) = 0.0;
  A(1, 1) = 3.0;
  A(1, 2) = 0.0;
  A(2, 0) = 0.0;
  A(2, 1) = 0.0;
  A(2, 2) = 4.0;

  b[0] = 2.0;
  b[1] = 6.0;
  b[2] = 8.0;

  Vector x = conjugate_gradient(A, b);

  REQUIRE(x[0] == Catch::Approx(1.0).epsilon(1e-6));
  REQUIRE(x[1] == Catch::Approx(2.0).epsilon(1e-6));
  REQUIRE(x[2] == Catch::Approx(2.0).epsilon(1e-6));
}

TEST_CASE("Conjugate gradient solves symmetric positive definite system",
          "[solver][cg]") {
  Matrix A(2, 2);
  Vector b(2);

  // SPD matrix
  A(0, 0) = 4.0;
  A(0, 1) = 1.0;
  A(1, 0) = 1.0;
  A(1, 1) = 3.0;

  b[0] = 1.0;
  b[1] = 2.0;

  Vector x = conjugate_gradient(A, b);

  // Exact solution: x = [0.090909..., 0.636363...]
  REQUIRE(x[0] == Catch::Approx(0.09090909).epsilon(1e-6));
  REQUIRE(x[1] == Catch::Approx(0.63636363).epsilon(1e-6));
}

TEST_CASE("Conjugate gradient residual is small", "[solver][cg]") {
  Matrix A(3, 3);
  Vector b(3);

  A(0, 0) = 4.0;
  A(0, 1) = 1.0;
  A(0, 2) = 0.0;
  A(1, 0) = 1.0;
  A(1, 1) = 3.0;
  A(1, 2) = 1.0;
  A(2, 0) = 0.0;
  A(2, 1) = 1.0;
  A(2, 2) = 2.0;

  b[0] = 1.0;
  b[1] = 2.0;
  b[2] = 3.0;

  Vector x = conjugate_gradient(A, b);

  Vector r = b - A * x;

  REQUIRE(r.norm() < 1e-6);
}

TEST_CASE("Conjugate gradient returns zero for zero rhs", "[solver][cg]") {
  Matrix A(2, 2);
  Vector b(2);

  A(0, 0) = 2.0;
  A(0, 1) = 0.0;
  A(1, 0) = 0.0;
  A(1, 1) = 5.0;

  b[0] = 0.0;
  b[1] = 0.0;

  Vector x = conjugate_gradient(A, b);

  REQUIRE(x[0] == Catch::Approx(0.0));
  REQUIRE(x[1] == Catch::Approx(0.0));
}

TEST_CASE("Conjugate gradient linearity check", "[solver][cg]") {
  Matrix A(2, 2);
  Vector b1(2), b2(2);

  A(0, 0) = 3.0;
  A(0, 1) = 1.0;
  A(1, 0) = 1.0;
  A(1, 1) = 2.0;

  b1[0] = 1.0;
  b1[1] = 0.0;
  b2[0] = 0.0;
  b2[1] = 2.0;

  Vector x1 = conjugate_gradient(A, b1);
  Vector x2 = conjugate_gradient(A, b2);
  Vector x12 = conjugate_gradient(A, b1 + b2);

  Vector r = b1 + b2 - A * x12;
  REQUIRE(r.norm() < 1e-6);
}
