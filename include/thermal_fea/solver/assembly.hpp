#pragma once

#include "thermal_fea/linalg/Matrix.hpp"
#include "thermal_fea/mesh/Mesh.hpp"

namespace thermal_fea::solver {

thermal_fea::linalg::Matrix K(const thermal_fea::mesh::Mesh &mesh,
                              const double k);

} // namespace thermal_fea::solver
