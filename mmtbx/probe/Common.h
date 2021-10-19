// Copyright(c) 2021, Richardson Lab at Duke
// Licensed under the Apache 2 license
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//     http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissionsand
// limitations under the License.

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
