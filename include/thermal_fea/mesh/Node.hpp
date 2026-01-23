#pragma once

#include <cstddef>

namespace thermal_fea::mesh {

struct Node {
  std::size_t id;
  double x;
  double y;
};

} // namespace thermal_fea::mesh
