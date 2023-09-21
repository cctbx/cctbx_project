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
#include "PositionReturn.h"
#include <vector>

typedef scitbx::vec3<double> vec3;

static scitbx::af::shared<iotbx::pdb::hierarchy::atom> getAtomsForMover(boost::python::object const& mover)
{
  molprobity::reduce::PositionReturn coarse =
    boost::python::extract<molprobity::reduce::PositionReturn>(mover.attr("CoarsePositions")());
  return coarse.atoms;
}

scitbx::af::shared< scitbx::af::shared<molprobity::probe::Point> > getAtomLocationsForMover(
  boost::python::object const& mover)
{
  // Find the coarse positions for this mover.
  molprobity::reduce::PositionReturn coarse =
    boost::python::extract<molprobity::reduce::PositionReturn>(mover.attr("CoarsePositions")());
  scitbx::af::shared< scitbx::af::shared<molprobity::probe::Point> > positions = coarse.positions;

  // Make a copy so that we don't change the original.
  scitbx::af::shared< scitbx::af::shared<molprobity::probe::Point> > total;
  for (size_t p = 0; p < positions.size(); ++p) {
    total.push_back(positions[p]);
  }

  // Add the fine positions for this mover.
  for (size_t c = 0; c < coarse.positions.size(); ++c) {
    molprobity::reduce::PositionReturn fine =
      boost::python::extract<molprobity::reduce::PositionReturn>(mover.attr("FinePositions")(c));
    for (size_t f = 0; f < fine.positions.size(); f++) {
      total.push_back(fine.positions[f]);
    }
  }

  return total;
}

namespace molprobity {
  namespace reduce {

AtomMoverLists::AtomMoverLists()
{
}

void AtomMoverLists::AddAtomMoverEntry(unsigned i_seq, boost::python::object mover)
{
  // Make sure that we have an entry to add it to.
  while (i_seq >= m_atomMoverLists.size()) {
    m_atomMoverLists.push_back(std::vector< boost::python::object >());
  }
  std::vector< boost::python::object > &atomMoverList = m_atomMoverLists[i_seq];

  // If we already have this, then we don't need to add it again.
  bool found = false;
  for (size_t i = 0; i < atomMoverList.size(); ++i) {
    if (atomMoverList[i].ptr() == mover.ptr()) {
      found = true;
      break;
    }
  }
  if (!found) {
    atomMoverList.push_back(mover);
  }
}

void AtomMoverLists::Clear()
{
  m_atomMoverLists.clear();
}

std::vector< boost::python::object > const& AtomMoverLists::GetAtomMoverList(unsigned i_seq) const
{
  if (i_seq >= m_atomMoverLists.size()) {
    PyErr_SetString(PyExc_RuntimeError, "Out of bounds reference in AtomMoverLists::GetAtomMoverList()");
    boost::python::throw_error_already_set();
  }
  return m_atomMoverLists[i_seq];
}

bool PairsOverlap(boost::python::object const &mover1,
  boost::python::object const& mover2,
  molprobity::probe::ExtraAtomInfoMap const &extraAtomInfoMap,
  double probeRad,
  AtomMoverLists &atomMoverLists)
{
  // Read the atoms and positions for the two movers.
  scitbx::af::shared<iotbx::pdb::hierarchy::atom> atoms1 = getAtomsForMover(mover1);
  scitbx::af::shared< scitbx::af::shared<molprobity::probe::Point> > positions1 =
    getAtomLocationsForMover(mover1);
  scitbx::af::shared<iotbx::pdb::hierarchy::atom> atoms2 = getAtomsForMover(mover2);
  scitbx::af::shared< scitbx::af::shared<molprobity::probe::Point> > positions2 =
    getAtomLocationsForMover(mover2);

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
      // Compute all portions of the distance other than the radius of the second atom.
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
            // Make sure that it is not already in the list before adding it.
            /// @todo This is the current bottleneck in interaction graph generation.
            atomMoverLists.AddAtomMoverEntry(atoms1[ai1].data->i_seq, mover2);
            atomMoverLists.AddAtomMoverEntry(atoms2[ai2].data->i_seq, mover1);

            // We found an overlap.
            ret = true;
          }
        }
      }
    }
  }

  return ret;
}

scitbx::af::shared<scitbx::af::shared<int> > FindOverlappingMoversAABB(
  scitbx::af::shared<boost::python::object> const& movers,
  molprobity::probe::ExtraAtomInfoMap const& extraAtomInfoMap,
  double probeRad)
{
  // Make a vector of bounding boxes for each mover.
  std::vector<scitbx::vec3<double> > mins(movers.size());
  std::vector<scitbx::vec3<double> > maxs(movers.size());
  for (size_t i = 0; i < movers.size(); ++i) {
    boost::python::object mover = movers[i];

    // Find the atoms and their locations for this mover.
    molprobity::reduce::PositionReturn coarse =
      boost::python::extract<molprobity::reduce::PositionReturn>(mover.attr("CoarsePositions")());
    scitbx::af::shared<iotbx::pdb::hierarchy::atom> atoms = getAtomsForMover(mover);
    scitbx::af::shared< scitbx::af::shared<molprobity::probe::Point> > total = getAtomLocationsForMover(mover);

    scitbx::vec3<double> thisMin(1e10, 1e10, 1e10);
    scitbx::vec3<double> thisMax(-1e10, -1e10, -1e10);
    for (size_t p = 0; p < total.size(); ++p) {
      // Find the bounding box that includes all dilated positions for all atoms.
      scitbx::af::shared<molprobity::probe::Point> const& pos = total[p];
      for (size_t a = 0; a < pos.size(); ++a) {
        // Find the radius of the atom, which is used to displace the bounding box.
        // We dilate by the probe radius; this will be dilated in each box, so together
        // we have 2 * probeRad.
        double r = probeRad + extraAtomInfoMap.getMappingFor(atoms[a]).getVdwRadius();

        scitbx::vec3<double> const& pt = pos[a];
        thisMin[0] = std::min(thisMin[0], pt[0] - r);
        thisMin[1] = std::min(thisMin[1], pt[1] - r);
        thisMin[2] = std::min(thisMin[2], pt[2] - r);
        thisMax[0] = std::max(thisMax[0], pt[0] + r);
        thisMax[1] = std::max(thisMax[1], pt[1] + r);
        thisMax[2] = std::max(thisMax[2], pt[2] + r);
      }
    }

    mins[i] = thisMin;
    maxs[i] = thisMax;
  }

  // For each pair of Movers whose bounding boxes overlap, add an entry to the
  // return value.
  scitbx::af::shared<scitbx::af::shared<int> > ret;
  if (movers.size() > 0) for (size_t i = 0; i < movers.size() - 1; ++i) {
    for (size_t j = i + 1; j < movers.size(); ++j) {
      if (
        mins[i][0] <= maxs[j][0] && maxs[i][0] >= mins[j][0] &&
        mins[i][1] <= maxs[j][1] && maxs[i][1] >= mins[j][1] &&
        mins[i][2] <= maxs[j][2] && maxs[i][2] >= mins[j][2]
      ) {
        scitbx::af::shared<int> overlap;
        overlap.push_back(i);
        overlap.push_back(j);
        ret.push_back(overlap);
      }
    }
  }

  return ret;
}

std::string InteractionGraph_test()
{
  // The AtomMoverLists class is tested in the Python test script.
  // The PairsOverlap() function is tested in the Python test script.
  // The FindOverlappingMoversAABB() function is tested in the Python test script.

  // All tests passed.
  return "";
}

  } // end namespace reduce
} // end namespace molprobity
