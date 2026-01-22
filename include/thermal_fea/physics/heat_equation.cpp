#include "heat_equation.hpp"
#include <cmath>

std::array<double, 2> node_coords(const thermal_fea::mesh::Node &node) {
  return {node.x(), node.y()};
}

namespace thermal_fea::physics {

Matrix3x3 get_Ke(const thermal_fea::mesh::Mesh &mesh,
                 const thermal_fea::mesh::Element &element, const double k) {
  std::array<double, 2> coords[3];
  for (int i = 0; i < 3; ++i) {
    coords[i] = node_coords(mesh.node(element.node(i)));
  }

  double x10 = coords[1][0] - coords[0][0];
  double x20 = coords[2][0] - coords[0][0];
  double y20 = coords[2][1] - coords[0][1];
  double y10 = coords[1][1] - coords[0][1];

  double y12 = coords[1][1] - coords[2][1];
  double x21 = coords[2][0] - coords[1][0];
  double x02 = -x20;
  double y01 = -y10;

  double area = .5 * (x10 * y20 - x20 * y10);
  std::array<double, 2> gradN1 = {y12 / (2 * area), x21 / (2 * area)};
  std::array<double, 2> gradN2 = {y20 / (2 * area), x02 / (2 * area)};
  std::array<double, 2> gradN3 = {y01 / (2 * area), x10 / (2 * area)};

  std::array<double, 2> gradN[3] = {gradN1, gradN2, gradN3};

  Matrix3x3 Ke = {{0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}};
  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 3; ++j) {
      Ke[i][j] +=
          k * area * (gradN[i][0] * gradN[j][0] + gradN[i][1] * gradN[j][1]);
    }
  }

  return Ke;
}
