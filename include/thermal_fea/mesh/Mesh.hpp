#pragma once

#include <cstddef>
#include <vector>

#include "thermal_fea/mesh/Element.hpp"
#include "thermal_fea/mesh/Node.hpp"

namespace thermal_fea::mesh {

class Mesh {
public:
  Mesh(const std::vector<Node> &nodes, const std::vector<Element> &elements):nodes_(nodes), elements_(elements){};

private:
  std::vector<Node> nodes_;
  std::vector<Element> elements_;

public:
  std::size_t nb_nodes() const;

  std::size_t nb_elements() const;

  const Node &get_node(std::size_t index) const;

  const Element &get_element(std::size_t index) const;

  const std::vector<Node> &nodes() const;

  const std::vector<Element> &elements() const;
};

} // namespace thermal_fea::mesh
