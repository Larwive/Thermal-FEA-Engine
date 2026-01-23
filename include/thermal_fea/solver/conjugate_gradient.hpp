#pragma once

#include <vector>
#include "thermal_fea/linalg/Matrix.hpp"


namespace thermal_fea::solver {
std::vector<double> conjugate_gradient(const thermal_fea::linalg::Matrix &A, const std::vector<double> &b);
} // namespace thermal_fea::solver