#pragma once

#include "thermal_fea/linalg/Matrix.hpp"
#include "thermal_fea/linalg/Vector.hpp"

namespace thermal_fea::solver {
thermal_fea::linalg::Vector
conjugate_gradient(const thermal_fea::linalg::Matrix &A,
                   const thermal_fea::linalg::Vector &b);

thermal_fea::linalg::Vector
conjugate_gradient_optimized(const thermal_fea::linalg::Matrix &A,
                             const thermal_fea::linalg::Vector &b);
} // namespace thermal_fea::solver
