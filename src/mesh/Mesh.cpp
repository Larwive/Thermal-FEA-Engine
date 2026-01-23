#include <cstddef>
#include <vector>

#include "thermal_fea/mesh/Element.hpp"
#include "thermal_fea/mesh/Node.hpp"
#include "thermal_fea/mesh/Mesh.hpp"

namespace thermal_fea::mesh {

std::size_t Mesh::nb_nodes() const { return nodes_.size(); }

std::size_t Mesh::nb_elements() const { return elements_.size(); }

const Node &Mesh::get_node(std::size_t index) const { return nodes_[index]; }

const Element &Mesh::get_element(std::size_t index) const {
  return elements_[index];
}

const std::vector<Node> &Mesh::nodes() const { return nodes_; }

const std::vector<Element> &Mesh::elements() const { return elements_; }

} // namespace thermal_fea::mesh
