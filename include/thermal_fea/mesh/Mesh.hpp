#pragma once

#include <cstddef>
#include <vector>

#include <thermal_fea/mesh/Element.hpp>
#include <thermal_fea/mesh/Node.hpp>

namespace thermal_fea::mesh {

class Mesh {
public:
  Mesh(const std::vector<Node> &nodes, const std::vector<Element> &elements);

private:
  std::vector<Node> nodes_;
  std::vector<Element> elements_;

public:
  std::size_t nb_nodes() const { return nodes_.size(); }

  std::size_t nb_elements() const { return elements_.size(); }

  const Node &get_node(std::size_t index) const { return nodes_[index]; }

  const Element &get_element(std::size_t index) const {
    return elements_[index];
  }

  const std::vector<Node> &nodes() const { return nodes_; }

  const std::vector<Element> &elements() const { return elements_; }
};

} // namespace thermal_fea::mesh
