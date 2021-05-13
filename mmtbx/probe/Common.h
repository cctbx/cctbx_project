// Copyright(c) 2021, Richardson Lab at Duke
// Licensed under the Apache 2 license

#pragma once

#include <vector>
#include <scitbx/vec3.h>
#include <scitbx/array_family/shared.h>

namespace molprobity {
  namespace probe {

    /// The type to be used for a single coordinate in space
    typedef double Coord;

    /// @brief A location in space
    typedef scitbx::vec3<Coord> Point;
  }
}
