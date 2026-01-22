#pragma once

#include <cstddef>

namespace thermal_fea::mesh {

/*
Each element is defined by three nodes (triangle) in this P1 2D mesh. We only need the ids to reference the nodes, and it avoids implementation details (ownership and memory management).
*/

struct Element {
  std::size_t node_ids[3];
};

}