// Copyright(c) 2023, Richardson Lab at Duke
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

#include "InteractionGraph.h"
#include <vector>

typedef scitbx::vec3<double> vec3;

namespace molprobity {
  namespace reduce {

bool PairsOverlap(boost::python::object const &mover1,
  scitbx::af::shared<iotbx::pdb::hierarchy::atom>  const& atoms1,
  scitbx::af::shared< scitbx::af::shared<molprobity::probe::Point> > const& positions1,
  boost::python::object const& mover2,
  scitbx::af::shared<iotbx::pdb::hierarchy::atom>  const& atoms2,
  scitbx::af::shared< scitbx::af::shared<molprobity::probe::Point> > const& positions2,
  molprobity::probe::ExtraAtomInfoMap & extraAtomInfoMap,
  double probeRad,
  boost::python::dict &atomMoverSets)
{
  // Make a vector the stores the radius for each atom in the two lists, matched based
  // on index in the list. This keeps us from having to redo these lookups multiple times
  // in the loop below.
  std::vector<double> atomRadii1(atoms1.size());
  for (int i = 0; i < atoms1.size(); i++) {
    atomRadii1[i] = extraAtomInfoMap.getMappingFor(atoms1[i]).getVdwRadius();
  }
  std::vector<double> atomRadii2(atoms2.size());
  for (int i = 0; i < atoms2.size(); i++) {
    atomRadii2[i] = extraAtomInfoMap.getMappingFor(atoms2[i]).getVdwRadius();
  }

  bool ret = false;
  for (size_t p1i = 0; p1i < positions1.size(); ++p1i) {
    scitbx::af::shared<molprobity::probe::Point> const &p1 = positions1[p1i];
    for (size_t ai1 = 0; ai1 < p1.size(); ai1++) {
      double r1 = atomRadii1[ai1];
      double x1 = p1[ai1][0];
      double y1 = p1[ai1][1];
      double z1 = p1[ai1][2];
      double limit1 = 2 * probeRad + r1;
      for (size_t p2i = 0; p2i < positions2.size(); ++p2i) {
        scitbx::af::shared<molprobity::probe::Point> const& p2 = positions2[p2i];
        for (size_t ai2 = 0; ai2 < p2.size(); ai2++) {
          double r2 = atomRadii2[ai2];
          double dx = x1 - p2[ai2][0];
          double dy = y1 - p2[ai2][1];
          double dz = z1 - p2[ai2][2];
          double dSquared = dx * dx + dy * dy + dz * dz;
          double limit = limit1 + r2;
          double limitSquared = limit * limit;
          if (dSquared <= limitSquared) {
            // Add the two movers to each other's set of movers that they overlap with.
            boost::python::object set_obj = atomMoverSets[atoms1[ai1].data->i_seq];
            boost::python::object set_type = set_obj.attr("__class__");
            set_type.attr("add")(set_obj, mover2);

            set_obj = atomMoverSets[atoms2[ai2].data->i_seq];
            set_type = set_obj.attr("__class__");
            set_type.attr("add")(set_obj, mover1);

            // We found an overlap.
            ret = true;
          }
        }
      }
    }
  }

  return ret;
}

std::string InteractionGraph_test()
{
  // The PairsOverlap() function is tested in the Python test script.

  // All tests passed.
  return "";
}

  } // end namespace reduce
} // end namespace molprobity
