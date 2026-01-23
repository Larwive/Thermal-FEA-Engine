#pragma once

#include <array>

#include "thermal_fea/mesh/Element.hpp"
#include "thermal_fea/mesh/Mesh.hpp"

namespace thermal_fea::physics {

typedef std::array<std::array<double, 3>, 3> Matrix3x3;

Matrix3x3 get_Ke(const thermal_fea::mesh::Mesh &mesh,
                 const thermal_fea::mesh::Element &element, const double k);

} // namespace thermal_fea::physics
