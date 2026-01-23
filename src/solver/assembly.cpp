#include "thermal_fea/solver/assembly.hpp"
#include "thermal_fea/linalg/Matrix.hpp"
#include "thermal_fea/mesh/Mesh.hpp"
#include "thermal_fea/physics/heat_equation.hpp"

namespace thermal_fea::solver {

thermal_fea::linalg::Matrix K(const thermal_fea::mesh::Mesh &mesh,
                              const double k) {
  thermal_fea::linalg::Matrix K(mesh.nb_nodes(), mesh.nb_nodes());

  for (const auto &element : mesh.elements()) {
    thermal_fea::physics::Matrix3x3 Ke =
        thermal_fea::physics::get_Ke(mesh, element, k);
    for (std::size_t i = 0; i < 3; ++i) {
      for (std::size_t j = 0; j < 3; ++j) {
        K(element.node_ids[i], element.node_ids[j]) += Ke[i][j];
      }
    }
  }

  return K;
}

} // namespace thermal_fea::solver
