#include <iostream>

#include <vector>

#include "thermal_fea/linalg/Matrix.hpp"
#include "thermal_fea/linalg/Vector.hpp"
#include "thermal_fea/mesh/Mesh.hpp"
#include "thermal_fea/mesh/Node.hpp"
#include "thermal_fea/solver/assembly.hpp"
#include "thermal_fea/solver/conjugate_gradient.hpp"

void set_temperature(thermal_fea::linalg::Matrix &K,
                     thermal_fea::linalg::Vector &F, std::size_t node_index,
                     double temperature) {

  // Moving contribution to the force vector.
  for (std::size_t j = 0; j < K.cols(); ++j) {
    if (j != node_index) {
      F(j) -= K(j, node_index) * temperature;
    }
  }

  // Dirichlet condition, elimination of couplings.
  for (std::size_t j = 0; j < K.cols(); ++j) {
    K(node_index, j) = 0.0;
    K(j, node_index) = 0.0;
  }
  K(node_index, node_index) = 1.0;
  F(node_index) = temperature;
}


void print_solution(const thermal_fea::mesh::Mesh &mesh,
                    const thermal_fea::linalg::Vector &u) {
  std::cout << "\nFinal temperatures:\n";
  for (std::size_t i = 0; i < mesh.nb_nodes(); ++i) {
    const auto &node = mesh.get_node(i);
    std::cout << "Node " << i
              << " : (" << node.x << ", " << node.y << ")"
              << " -> T = " << u(i) << '\n';
  }
}

int main() {

  // Four nodes in a square.
  std::vector<thermal_fea::mesh::Node> nodes(
      {thermal_fea::mesh::Node{0, 0.0, 0.0},
       thermal_fea::mesh::Node{1, 1.0, 0.0},
       thermal_fea::mesh::Node{2, 0.0, 1.0},
       thermal_fea::mesh::Node{3, 1.0, 1.0}});

  // Two triangles for the mesh. Make sure orientation is consistent !
  std::vector<thermal_fea::mesh::Element> elements(
      {thermal_fea::mesh::Element{{0, 1, 2}},
       thermal_fea::mesh::Element{{1, 3, 2}}});

  thermal_fea::mesh::Mesh mesh(nodes, elements);

  thermal_fea::linalg::Matrix K = thermal_fea::solver::K(mesh, 1.0);

  thermal_fea::linalg::Vector F(mesh.nb_nodes());

  // Imposed temperatures at nodes 0 and 3.
  set_temperature(K, F, 0, 100.0);
  set_temperature(K, F, 3, 15.0);

  thermal_fea::linalg::Vector u = thermal_fea::solver::conjugate_gradient(K, F);

  print_solution(mesh, u);
  return 0;
}
