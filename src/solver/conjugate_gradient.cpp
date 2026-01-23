#include <cmath>

#include "thermal_fea/linalg/Matrix.hpp"
#include "thermal_fea/linalg/Vector.hpp"

namespace thermal_fea::solver {

thermal_fea::linalg::Vector
conjugate_gradient(const thermal_fea::linalg::Matrix &A,
                   const thermal_fea::linalg::Vector &b) {
  thermal_fea::linalg::Vector x(b.size());
  thermal_fea::linalg::Vector r = b - A * x;
  thermal_fea::linalg::Vector r_new(r.size());
  thermal_fea::linalg::Vector p = r;
  double alpha = 0.0;
  double beta = 0.0;

  for (int i = 0; i < 100; ++i) {
    alpha = r.dot(r) / p.dot(A * p);
    x += alpha * p;
    r_new = r - alpha * A * p;
    if (r_new.norm() < 1e-6) {
      break;
    }
    beta = r_new.dot(r_new) / r.dot(r);
    p = r_new + beta * p;
    r = r_new;
  }

  return x;
}

} // namespace thermal_fea::solver
