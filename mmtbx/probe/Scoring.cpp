// Copyright(c) 2021-2023, Richardson Lab at Duke
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

#include <cmath>
#include <scitbx/constants.h>
#include <algorithm>
#include "Scoring.h"
#include "DotSpheres.h"

typedef scitbx::vec3<double> vec3;

namespace molprobity {
  namespace probe {


ContactResult closest_contact(Point dot, Point atom, double atom_radius)
{
  ContactResult ret;
  Point diff = dot - atom;

  // if the dot is at the center of the atom, then pick an arbitrary point
  // on the surface and return it.
  if (diff.length_sq() == 0) {
    ret.distAboveSurface = -atom_radius;
    ret.closestContact = atom + Point(atom_radius, 0, 0);

  // Otherwise, find the point of closest approach and its distance and return
  // them.
  } else {
    double len = diff.length();
    ret.distAboveSurface = len - atom_radius;
    ret.closestContact = atom + diff * atom_radius / len;
  }
  return ret;
}

double dot2srcCenter(const Point& dot, const Point& srcLoc, double srcVDWRad, const Point& targLoc) {
  // The vector from the source pointing towards the target that is the radius of the source atom
  Point src2targVec = (targLoc - srcLoc).normalize() * srcVDWRad;
  // The point on the surface of the source atom that is closest to the target
  Point srcSurfacePoint = src2targVec + srcLoc;
  // The distance from the dot to the point on the source surface that is closest to the target
  return (srcSurfacePoint - dot).length();
}

double kissEdge2bullsEye(double ra, double rb, double rp) {
  /// @todo Describe what this is computing and how
  return 2 * ra * sqrt(rb * rp / ((ra + rb) * (ra + rp)));
}
bool annularDots(const Point& dot, const Point& srcLoc, double srcVDWRad,
  const Point& targLoc, double targVDWRad, double probeRadius) {
  return dot2srcCenter(dot, srcLoc, srcVDWRad, targLoc) > kissEdge2bullsEye(srcVDWRad, targVDWRad, probeRadius);
}


int atom_charge(iotbx::pdb::hierarchy::atom const& atom)
{
  // Get the tidy version of the charge string stripped of any extra characters.
  // Detangle the fact that the string claims to have the potential to be optional.
  // This will make the string either blank or consisting of two characters, [12][+-],
  // with the first indicating the magnitude and the second indicating the sign.
  int ret = 0;
  std::string chStr;
  boost::optional<std::string> opt = atom.charge_tidy(true);
  if (opt) { chStr = opt.get(); }

  if (chStr.size() > 0) {
    ret = chStr[0] - '0';
    if (chStr[1] == '-') { ret *= -1; }
  }

  return ret;
}

bool DotScorer::point_inside_atoms(Point const& location,
  scitbx::af::shared<iotbx::pdb::hierarchy::atom> const& atoms)
{
  for (scitbx::af::shared<iotbx::pdb::hierarchy::atom>::const_iterator e = atoms.begin();
    e != atoms.end(); e++) {
    double vdwe = m_extraInfoMap.getMappingFor(*e).getVdwRadius();
    if ((location - e->data->xyz).length_sq() < vdwe * vdwe) {
      return true;
    }
  }
  return false;
}

scitbx::af::shared<Point> DotScorer::trim_dots(iotbx::pdb::hierarchy::atom const& atom,
  scitbx::af::shared<Point> const& dots,
  scitbx::af::shared<iotbx::pdb::hierarchy::atom> const& exclude)
{
  scitbx::af::shared<Point> ret;
  for (scitbx::af::shared<Point>::const_iterator d = dots.begin(); d != dots.end(); d++) {
    Point dotLoc = atom.data->xyz + *d;
    if (!point_inside_atoms(dotLoc, exclude)) {
      ret.push_back(*d);
    }
  }
  return ret;
}

DotScorer::CheckDotResult DotScorer::check_dot(
  iotbx::pdb::hierarchy::atom const &sourceAtom,
  Point const& dotOffset, double probeRadius,
  scitbx::af::shared<iotbx::pdb::hierarchy::atom> const &interacting,
  scitbx::af::shared<iotbx::pdb::hierarchy::atom> const& exclude,
  double overlapScale)
{
  // Defaults to "no overlap" type
  CheckDotResult ret;

  // Find the extra-atom information for the source atom.
  ExtraAtomInfo const& sourceExtra = m_extraInfoMap.getMappingFor(sourceAtom);

  // If the source atom is an ion and we are ignoring ions, we're done.
  if (sourceExtra.getIsIon() && m_ignoreIonInteractions) {
    return ret;
  }

  // Find the world-space location of the dot by adding it to the location of the source atom.
  Point absoluteDotLocation = sourceAtom.data->xyz + dotOffset;

  // Check to see if the dot should be removed from consideration because it is inside an excluded atom.
  // Doing this test ahead of the neighbor-interaction test makes things faster for the long-running
  // 4fen test case.
  if (point_inside_atoms(absoluteDotLocation, exclude)) {
    return ret;
  }

  // The probe location is in the same direction as d from the source but is further away by the
  // probe radius.
  Point probLoc = sourceAtom.data->xyz + dotOffset.normalize() * (dotOffset.length() + probeRadius);

  bool isHydrogenBond = false;          ///< Are we looking at a hydrogen bond to our neighbor?
  bool tooCloseHydrogenBond = false;    ///< Are we too close to be a hydrogen bond?
  double hydrogenBondMinDist = 0;       ///< Hydrogen bond minimum distance based on the atom types (will be set below).
  bool keepDot = false;                 ///< Did we find a neighbor and we're not in an excluded atom?
  bool causeIsDummy = false;            ///< Is the cause a dummy Hydrogen?

  // Look through each atom to find the one with the smallest gap, which is the one that
  // the dot would interact with.
  for (scitbx::af::shared<iotbx::pdb::hierarchy::atom>::const_iterator b = interacting.begin();
       b != interacting.end(); b++) {

    ExtraAtomInfo const& bExtra = m_extraInfoMap.getMappingFor(*b);

    // If the potential target atom is an ion and we are ignoring ions, skip it.
    if (bExtra.getIsIon() && m_ignoreIonInteractions) {
      continue;
    }

    // See if we are too far away to interact, bail if so.
    Point const& locb = b->data->xyz;
    double vdwb = bExtra.getVdwRadius();
    double squareProbeDist = (probLoc - locb).length_sq();
    double pRadPlusVdwb = vdwb + probeRadius;
    if (squareProbeDist > pRadPlusVdwb * pRadPlusVdwb) {
      continue;
    }

    // If we're in incompatible conformations, then we can't interact.
    if (!compatible_conformations(sourceExtra.getAltLoc(), bExtra.getAltLoc())) {
      continue;
    }

    // At this point, we are within the probe radius past the edge, so we're in contention
    // to be the nearest atom.  Find the distance from the dot rather than from the probe
    // to check our actual interaction behavior.
    double dist = (absoluteDotLocation - locb).length();
    double gap = dist - vdwb;

    // See if we replace the currently-closest atom.
    if (gap < ret.gap) {

      // Figure out what kind of overlap this is based on the atom types and
      // charge status of the two atoms.
      int chargeSource = sourceExtra.getCharge();
      int chargeB = bExtra.getCharge();

      bool bothCharged = (chargeSource != 0) && (chargeB != 0);
      bool chargeComplement = bothCharged && (chargeSource * chargeB < 0);

      // See if one of the atoms is a hydrogen donor and the other can accept hydrogen bonds.
      // Then check to see whether this is a hydrogen bond and behave accordingly.
      bool couldHBond = (sourceExtra.getIsDonor() && bExtra.getIsAcceptor())
        || (sourceExtra.getIsAcceptor() && bExtra.getIsDonor());
      if (couldHBond && ((!bothCharged) || chargeComplement)) {
        isHydrogenBond = true;
        hydrogenBondMinDist = bothCharged ? m_maxChargedHydrogenOverlap : m_maxRegularHydrogenOverlap;
        tooCloseHydrogenBond = (gap < -hydrogenBondMinDist);
      } else {
        // If one of the atoms is a dummy Hydrogen, then we pretend that it does not exist
        // for non-Hydrogen-bond interactions.
        if (sourceExtra.getIsDummyHydrogen() || bExtra.getIsDummyHydrogen()) { continue; }

        // This is not a hydrogen bond.
        isHydrogenBond = tooCloseHydrogenBond = false;
      }

      // Record which atom is the closest and mark the dot to be kept because we found an
      // atom that is close enough.  Keep track of the characteristics of the collision.
      keepDot = true;
      causeIsDummy = bExtra.getIsDummyHydrogen();
      ret.gap = gap;
      ret.cause = *b;
    }
  }

  // If this dot is a keeper, fill in non-default return values.
  if (keepDot) {
    // Determine the overlap type and amount of overlap.
    if (ret.gap > 0) {
      ret.overlap = 0;
      if (m_weakHBonds && isHydrogenBond) {
        ret.overlapType = DotScorer::HydrogenBond;
      } else {
        ret.overlapType = DotScorer::NoOverlap;
      }
    } else if (isHydrogenBond) {
      ret.overlap = -overlapScale * ret.gap;
      if (tooCloseHydrogenBond) {
        // Reduce the gap magnitude by the expected hydrogen bond distance (gap is
        // negative, so we add here), compute the overlap, and report it as a clash.
        ret.gap += hydrogenBondMinDist;
        ret.overlap = -overlapScale * ret.gap;
        ret.overlapType = DotScorer::Clash;
      } else {
        ret.overlapType = DotScorer::HydrogenBond;
      }
    } else {  // ret.gap <= 0 and not a hydrogen bond
      ret.overlap = -overlapScale * ret.gap;
      ret.overlapType = DotScorer::Clash;
    }

    ret.annular = annularDots(absoluteDotLocation, sourceAtom.data->xyz, sourceExtra.getVdwRadius(),
      ret.cause.data->xyz, m_extraInfoMap.getMappingFor(ret.cause).getVdwRadius(), probeRadius);
  }

  // If the source or target is a dummy hydrogen, then we ignore dots that are too-close hydrogen bonds
  // and ones that turned out in the end not to be hydrogen bonds.
  if ((sourceExtra.getIsDummyHydrogen() || causeIsDummy)
      && (tooCloseHydrogenBond || (ret.overlapType != DotScorer::HydrogenBond))) {
    ret.overlapType = DotScorer::Ignore;
  }

  return ret;
}

DotScorer::InteractionType DotScorer::interaction_type(
  OverlapType overlapType, double gap, bool separateBadBumps) const
{
  switch (overlapType) {
    case DotScorer::NoOverlap:
      if (gap > m_contactCutoff) {
        return WideContact;
      } else {
        return CloseContact;
      }
      break;
    case DotScorer::Clash:
      if (gap > -m_bumpOverlap) {
        return SmallOverlap;
      } else if (separateBadBumps) {
        // We're checking separately for bad bump overlaps
        if (gap > -m_badBumpOverlap) {
          return Bump;
        } else {
          return BadBump;
        }
      } else {
        // We're not separately checking for bad bump overlaps
        return Bump;
      }
      break;
    case DotScorer::HydrogenBond:
      if (m_weakHBonds) {
        if (gap > 0) {
          return WeakHydrogenBond;
        } else {
          return StandardHydrogenBond;
        }
      } else {
        return StandardHydrogenBond;
      }
      break;
    case DotScorer::Ignore:
    default:
      return Invalid;
  }
}

std::string DotScorer::interaction_type_name(InteractionType t)
{
  switch (t) {
  case WideContact:
    return "wide_contact";
  case CloseContact:
    return "close_contact";
  case WeakHydrogenBond:
    return "weak_H-bond";
  case SmallOverlap:
    return "small_overlap";
  case Bump:
    return "bad_overlap";
  case BadBump:
    return "worse_overlap";
  case StandardHydrogenBond:
    return "H-bond";
  case Invalid:
    return "invalid (internal error)";
  default:
    return "unrecognized (internal error)";
  }
}

std::string DotScorer::interaction_type_short_name(InteractionType t)
{
  switch (t) {
  case WideContact:
    return "wc";
  case CloseContact:
    return "cc";
  case WeakHydrogenBond:
    return "wh";
  case SmallOverlap:
    return "so";
  case Bump:
    return "bo";
  case BadBump:
    return "wo";
  case StandardHydrogenBond:
    return "hb";
  case Invalid:
    return "invalid (internal error)";
  default:
    return "unrecognized (internal error)";
  }
}

DotScorer::ScoreDotsResult DotScorer::score_dots(
  iotbx::pdb::hierarchy::atom const &sourceAtom, double minOccupancy,
  SpatialQuery &spatialQuery, double nearbyRadius, double probeRadius,
  scitbx::af::shared<iotbx::pdb::hierarchy::atom> const &exclude,
  scitbx::af::shared<Point> const &dots, double density, bool onlyBumps,
  bool preTrimmedDots)
{
  // This method is based on AtomPositions::atomScore() from Reduce.
  // It is passed only the dots that it should score rather than excluding them
  // internally like that function does.

  // Default return has 0 value for all subscores.
  ScoreDotsResult ret;

  // Check for invalid parameters
  if ((density <= 0) || (probeRadius < 0)) {
    return ret;
  }

  // Make sure that the occupancy of the atom is high enough to score it.
  // If not, we send a valid but otherwise empty return value.
  if (sourceAtom.data->occ < minOccupancy) {
    ret.valid = true;
    return ret;
  }

  // If we're looking for other contacts besides bumps, we need to add twice the
  // probe radius to the neighbor list to ensure that we include atoms where the
  // probe spans from the source to the target.
  double twiceProbeRadius = 2 * probeRadius;
  if (!onlyBumps) {
    nearbyRadius += twiceProbeRadius;
  }

  // Find the neighboring atoms that are potentially interacting.
  // The nonzero minimum distance prevents us from selecting the source atom.
  /// @todo The nonzero distance will also prevent noticing co-located atoms,
  // which should probably be recorded as bad clashes.
  scitbx::af::shared<iotbx::pdb::hierarchy::atom> neighbors =
    spatialQuery.neighbors(sourceAtom.data->xyz, 1e-5, nearbyRadius);

  // Select only those atoms actually interacting: that have sufficient occupancy, for whom the
  // gap between the Van Der Waals surfaces is less than the probe radius, and which
  // are not in the excluded list.
  scitbx::af::shared<iotbx::pdb::hierarchy::atom> interacting;
  ExtraAtomInfo const& sourceExtra = m_extraInfoMap.getMappingFor(sourceAtom);
  for (scitbx::af::shared<iotbx::pdb::hierarchy::atom>::const_iterator a = neighbors.begin();
       a != neighbors.end(); a++) {
    ExtraAtomInfo const& aExtra = m_extraInfoMap.getMappingFor(*a);
    bool excluded = false;
    for (scitbx::af::shared<iotbx::pdb::hierarchy::atom>::const_iterator e = exclude.begin();
         e != exclude.end(); e++) {
      if (e->data.get() == a->data.get()) {
        excluded = true;
        break;
      }
    }
    if (!excluded) {
      double interactionDistance = sourceExtra.getVdwRadius() + aExtra.getVdwRadius() + twiceProbeRadius;
      if ((std::abs(a->data->occ) >= minOccupancy)
        && ((a->data->xyz - sourceAtom.data->xyz).length() <= interactionDistance)
      ) {
        interacting.push_back(*a);
      }
    }
  }

  // Construct this outside of the loop to avoid repeated construction.
  scitbx::af::shared<iotbx::pdb::hierarchy::atom> const emptyExclude;

  // Run through all of the dots and determine whether and how to score each.
  for (scitbx::af::shared<Point>::const_iterator d = dots.begin(); d != dots.end(); d++) {

    // If the dots have been pre-trimmed, we don't need to pass the excluded list.
    scitbx::af::shared<iotbx::pdb::hierarchy::atom> const* excludeToUse = &exclude;
    if (preTrimmedDots) {
      excludeToUse = &emptyExclude;
    }
    CheckDotResult score = check_dot(sourceAtom, *d, probeRadius, interacting, *excludeToUse);

    // Compute the score for the dot based on the overlap type and amount of overlap.
    // Assign it to the appropriate subscore.
    InteractionType it = interaction_type(score.overlapType, score.gap, true);
    switch (it) {

    case DotScorer::Invalid: // The dot should be ignored, so is not scored.
      break;

    case DotScorer::WideContact:  // Types where the dot is outside the target
    case DotScorer::CloseContact:
    case DotScorer::WeakHydrogenBond:
      if ((!onlyBumps) && !score.annular) {
        double scaledGap = score.gap / m_gapScale;
        ret.attractSubScore += exp(-scaledGap * scaledGap);
      }
      break;

    case DotScorer::SmallOverlap: // Types where the dots overlap and are considered clashes
    case DotScorer::Bump:
    case DotScorer::BadBump:
      {
        ret.bumpSubScore += -m_bumpWeight * score.overlap;
        // See if we should flag this atom as having a bad bump
        if (it == DotScorer::BadBump) {
          ret.hasBadBump = true;
        }
      }
      break;

    case DotScorer::StandardHydrogenBond: // Weak Hydrogen bonds are counted above
      if (!onlyBumps) {
        ret.hBondSubScore += m_hBondWeight * score.overlap;
      } else {  // In this case, we treat it as a bump
        ret.bumpSubScore += -m_bumpWeight * score.overlap;
      }
      break;

    default:
      // This should never happen.  Returns with ret invalid to indicate an error.
      std::cerr << "DotScorer::score_dots(): Internal error: Unrecognized overlap type: " <<
        static_cast<int>(score.overlapType) << std::endl;
      return ret;
    }
  }

  // Normalize the score by the density so that it does not depend on the number of dots
  // that were constructed for the atom.
  ret.bumpSubScore /= density;
  ret.hBondSubScore /= density;
  ret.attractSubScore /= density;
  ret.valid = true;
  return ret;
}

unsigned DotScorer::count_surface_dots(iotbx::pdb::hierarchy::atom const &sourceAtom,
  scitbx::af::shared<Point> const& dots,
  scitbx::af::shared<iotbx::pdb::hierarchy::atom> const& exclude)
{
  unsigned ret = 0;

  // Run through all of the dots and determine whether each is valid.
  for (scitbx::af::shared<Point>::const_iterator d = dots.begin(); d != dots.end(); d++) {

    // Find the world-space location of the dot by adding it to the location of the source atom.
    // The probe location is in the same direction as d from the source but is further away by the
    // probe radius.
    Point absoluteDotLocation = sourceAtom.data->xyz + (*d);

    // Check to see if the dot should be removed from consideration because it is also inside an excluded atom.
    bool keepDot = true;
    for (scitbx::af::shared<iotbx::pdb::hierarchy::atom>::const_iterator e = exclude.begin();
      e != exclude.end(); e++) {
      double vdwe = m_extraInfoMap.getMappingFor(*e).getVdwRadius();
      if ((absoluteDotLocation - e->data->xyz).length_sq() <= vdwe * vdwe) {
        keepDot = false;
        break;
      }
    }
    // Increment if the dot was not excluded.
    if (keepDot) {
      ret++;
    }
  }

  return ret;
}

bool DotScorer::compatible_conformations(std::string const& a1, std::string const& a2)
{
  if (a1.size() == 0 || a1[0] == ' ' || a2.size() == 0 || a2[0] == ' ') {
    return true;
  }
  return a1 == a2;
}

//===========================================================================================================
// Testing code below here

/// @brief Returns true of the two floating-point numbers are nearly equal
static bool closeTo(double a, double b) {
  return fabs(a - b) < 1e-10;
}

std::string DotScorer::test()
{
  /// @todo Check the annular-dots behavior.

  // Check compatible_conformations() for all combinations of altlocs.
  {
    if (!compatible_conformations("", "")) {
      return "DotScorer::test(): Compatible conformations failed for both nul";
    }
    if (!compatible_conformations("A", "")) {
      return "DotScorer::test(): Compatible conformations failed for A and nul";
    }
    if (!compatible_conformations("", "A")) {
      return "DotScorer::test(): Compatible conformations failed for nul and A";
    }
    if (!compatible_conformations(" ", " ")) {
      return "DotScorer::test(): Compatible conformations failed for both blank";
    }
    if (!compatible_conformations("A", " ")) {
      return "DotScorer::test(): Compatible conformations failed for A and blank";
    }
    if (!compatible_conformations(" ", "A")) {
      return "DotScorer::test(): Compatible conformations failed for blank and A";
    }
    if (!compatible_conformations("A", "A")) {
      return "DotScorer::test(): Compatible conformations failed for both A";
    }
    if (compatible_conformations("A", "B")) {
      return "DotScorer::test(): Compatible conformations falsely succeeded for A and B";
    }
  }

  // Check trim_dots(), which along with check_dot() will also check point_inside_atoms().
  {
    double targetRad = 1.5, sourceRad = 1.0;
    DotSphere ds(sourceRad, 200);
    unsigned int atomSeq = 0;
    scitbx::af::shared<Point> kept;

    // Construct a single target atom, including its extra info
    iotbx::pdb::hierarchy::atom a;
    a.set_xyz(vec3(0, 0, 0));
    a.set_occ(1);
    a.data->i_seq = atomSeq++;
    scitbx::af::shared<iotbx::pdb::hierarchy::atom> atoms;
    atoms.push_back(a);
    SpatialQuery sq(atoms);
    ExtraAtomInfo e(targetRad, true);
    scitbx::af::shared<ExtraAtomInfo> infos;
    infos.push_back(e);

    // Construct a source atom, including its extra info.
    // This will be a hydrogen but not a donor.
    iotbx::pdb::hierarchy::atom source;
    source.set_occ(1);
    source.data->i_seq = atomSeq++;
    ExtraAtomInfo se(sourceRad);
    atoms.push_back(source);
    infos.push_back(se);

    // Construct the scorer to be used with the specified bond gaps.
    DotScorer as(ExtraAtomInfoMap(atoms, infos));

    // Construct an exclusion list and add the target atom.
    scitbx::af::shared<iotbx::pdb::hierarchy::atom> exclude;
    exclude.push_back(a);

    // Check the source atom not overlapping with the target.
    source.set_xyz(vec3(targetRad + sourceRad + 0.1, 0, 0));
    kept = as.trim_dots(source, ds.dots(), exclude);
    if (kept.size() != ds.dots().size()) {
      return "DotScorer::test(): Unexpected non-excluded dot count for non-overlapping case";
    }

    // Check the source atom completely overlapping with the target.
    source.set_xyz(vec3(0, 0, 0));
    kept = as.trim_dots(source, ds.dots(), exclude);
    if (kept.size() != 0) {
      return "DotScorer::test(): Unexpected nonzero non-excluded dot count for fully-overlapping case";
    }

    // Check the source atom partially overlapping with the target.
    source.set_xyz(vec3(targetRad, 0, 0));
    kept = as.trim_dots(source, ds.dots(), exclude);
    if ((kept.size() == 0) || (kept.size() >= ds.dots().size())) {
      return "DotScorer::test(): Unexpected non-excluded dot count for partially-overlapping case";
    }
  }

  // Test the check_dot() function for atoms in different alternate conformations
  {
    double targetRad = 1.5, sourceRad = 1.0, probeRad = 0.25;
    DotSphere ds(sourceRad, 200);
    unsigned int atomSeq = 0;

    // Construct and fill the SpatialQuery information
    // with a vector of a single target atom, including its extra info.
    iotbx::pdb::hierarchy::atom a;
    a.set_xyz(vec3(0, 0, 0));
    a.set_occ(1.0);
    a.data->i_seq = atomSeq++;
    scitbx::af::shared<iotbx::pdb::hierarchy::atom> atoms;
    atoms.push_back(a);
    SpatialQuery sq(atoms);
    ExtraAtomInfo e(targetRad, true);

    // Construct an empty exclusion list.
    scitbx::af::shared<iotbx::pdb::hierarchy::atom> exclude;

    // Construct a source atom, including its extra info.
    // This will be a hydrogen but not a donor.
    iotbx::pdb::hierarchy::atom source;
    source.set_occ(1.0);
    source.data->i_seq = atomSeq++;
    ExtraAtomInfo se(sourceRad);
    atoms.push_back(source);
    ScoreDotsResult res;

    // Check the source atom just overlapping with the target.
    source.set_xyz(vec3(targetRad + sourceRad - 0.1, 0, 0));

    std::vector<std::string> altlocs;
    altlocs.push_back("");
    altlocs.push_back(" ");
    altlocs.push_back("A");
    altlocs.push_back("B");
    for (size_t i = 0; i < altlocs.size(); i++) {
      e.setAltLoc(altlocs[i]);

      for (size_t j = 0; j < altlocs.size(); j++) {
        se.setAltLoc(altlocs[j]);

        // Make the extra atom info map match the new state.
        scitbx::af::shared<ExtraAtomInfo> infos;
        infos.push_back(e);
        infos.push_back(se);

        // Construct the scorer with updated information.
        DotScorer as(ExtraAtomInfoMap(atoms, infos));

        res = as.score_dots(source, 1, sq, sourceRad + targetRad,
          probeRad, exclude, ds.dots(), ds.density(), false);
        if (!res.valid) {
          return "DotScorer::test(): Could not score dots for alternate conformation test case";
        }
        bool zero = res.totalScore() == 0;
        if (zero != !DotScorer::compatible_conformations(e.getAltLoc(), se.getAltLoc())) {
          return "DotScorer::test(): Unexpected nonzero value for alternate conformation test case '" + e.getAltLoc() +
            "' vs. '" + se.getAltLoc() + "'";
        }
      }
    }
  }

  // Test the check_dot() function to make sure that it gets correct interaction types for
  // all ranges of interaction.  Also check the cause.
  {
    double targetRad = 1.5, sourceRad = 1.0, probeRad = 0.25;
    unsigned int atomSeq = 0;

    // Construct and fill the interaction list with a vector of a single target atom
    iotbx::pdb::hierarchy::atom a;
    a.set_xyz(vec3(0,0,0));
    a.set_occ(1);
    a.data->i_seq = atomSeq++;
    scitbx::af::shared<iotbx::pdb::hierarchy::atom> atoms;
    atoms.push_back(a);
    ExtraAtomInfo e(targetRad);
    scitbx::af::shared<ExtraAtomInfo> infos;
    infos.push_back(e);
    scitbx::af::shared<iotbx::pdb::hierarchy::atom> interacting;
    interacting.push_back(a);

    // Construct the source atom, including its extra info.
    iotbx::pdb::hierarchy::atom source;
    source.set_occ(1);
    source.data->i_seq = atomSeq++;
    ExtraAtomInfo se(sourceRad);
    atoms.push_back(source);
    infos.push_back(se);

    // Construct an empty exclusion list.
    scitbx::af::shared<iotbx::pdb::hierarchy::atom> exclude;

    // Construct the scorer to be used.
    DotScorer as(ExtraAtomInfoMap(atoms, infos));

    // Check a dot that is to -X of the source atom by its radius against the
    // target for various positions of the source atom.
    CheckDotResult res;
    source.set_xyz(vec3( sourceRad + targetRad - 0.6,0,0 ));
    res = as.check_dot(source, Point(-sourceRad, 0, 0), probeRad, interacting, exclude);
    if (res.overlapType != DotScorer::Clash) {
      return "DotScorer::test(): Did not find clash when expected for dot_score()";
    }
    if (res.cause.data != a.data) {
      return "DotScorer::test(): Did not find expected cause for dot_score()";
    }
    if (as.interaction_type(res.overlapType, res.gap, true) != DotScorer::BadBump) {
      return "DotScorer::test(): Did not find WorseOverlap when expected for dot_score()";
    }
    if (as.interaction_type(res.overlapType, res.gap, false) != DotScorer::Bump) {
      return "DotScorer::test(): Found BadOverlap when not expected for dot_score()";
    }

    source.set_xyz(vec3( sourceRad + targetRad - 0.45,0,0 ));
    res = as.check_dot(source, Point(-sourceRad, 0, 0), probeRad, interacting, exclude);
    if (res.overlapType != DotScorer::Clash) {
      return "DotScorer::test(): Did not find clash when expected for dot_score()";
    }
    if (as.interaction_type(res.overlapType, res.gap) != DotScorer::Bump) {
      return "DotScorer::test(): Did not find BadOverlap when expected for dot_score()";
    }

    source.set_xyz(vec3( sourceRad + targetRad - 0.01,0,0 ));
    res = as.check_dot(source, Point(-sourceRad, 0, 0), probeRad, interacting, exclude);
    if (res.overlapType != DotScorer::Clash) {
      return "DotScorer::test(): Did not find small clash when expected for dot_score()";
    }
    if (as.interaction_type(res.overlapType, res.gap) != DotScorer::SmallOverlap) {
      return "DotScorer::test(): Did not find SmallOverlap when expected for dot_score()";
    }

    source.set_xyz(vec3( sourceRad + targetRad + 0.01,0,0 ));
    res = as.check_dot(source, Point(-sourceRad, 0, 0), probeRad, interacting, exclude);
    if (res.overlapType != DotScorer::NoOverlap) {
      return "DotScorer::test(): Did not find no overlap when expected for dot_score()";
    }
    if (as.interaction_type(res.overlapType, res.gap) != DotScorer::CloseContact) {
      return "DotScorer::test(): Did not find CloseContact when expected for dot_score()";
    }

    source.set_xyz(vec3( sourceRad + targetRad + 0.26,0,0 ));
    res = as.check_dot(source, Point(-sourceRad, 0, 0), probeRad, interacting, exclude);
    if (res.overlapType != DotScorer::NoOverlap) {
      return "DotScorer::test(): Did not find no overlap when expected for dot_score()";
    }
    if (as.interaction_type(res.overlapType, res.gap) != DotScorer::WideContact) {
      return "DotScorer::test(): Did not find WideContact when expected for dot_score()";
    }

    // Check so far away that there won't be any nearby atoms.
    source.set_xyz(vec3( sourceRad + targetRad + 10.0,0,0 ));
    res = as.check_dot(source, Point(-sourceRad, 0, 0), probeRad, interacting, exclude);
    if (res.overlapType != DotScorer::Ignore) {
      return "DotScorer::test(): Did not find ignore when expected for dot_score()";
    }
    if (as.interaction_type(res.overlapType, res.gap) != DotScorer::Invalid) {
      return "DotScorer::test(): Did not find Invalid when expected for dot_score()";
    }
  }

  // Construct test cases with all combinations of charges and extra information, holding the
  // radii of the neighbor atom and probe atom constant.  Do this in combination with adding or
  // not adding an excluded atom that completely covers the neighbor atom.
  // Run tests against all of these cases to ensure that the behavior is as expected in each case.
  // This tests the case logic within the code.
  // The target radius has to be large enough to get a bad bump even for hydrogen bond cases.
  double targetRad = 1.5, sourceRad = 1.0, probeRad = 0.25;

  // Construct the dot sphere to be used.
  DotSphere ds(sourceRad, 200);

  static const char* chargesArray[] = {"--", "-", "", "+", "++"};
  std::vector<std::string> charges;
  for (size_t i = 0; i < sizeof(chargesArray)/sizeof(chargesArray[0]); i++) {
    charges.push_back(chargesArray[i]);
  }
  // Make a vector of Booleans that covers false and true so that we can use
  // a for-loop iterator to cover both cases and set all combinations of Boolean variables.
  std::vector<bool> bools; bools.push_back(false); bools.push_back(true);
  for (std::vector<std::string>::const_iterator targetCharge = charges.begin();
       targetCharge != charges.end(); targetCharge++) {
  for (std::vector<std::string>::const_iterator sourceCharge = charges.begin();
       sourceCharge != charges.end(); sourceCharge++) {
   for (std::vector<bool>::const_iterator targetAccept = bools.begin();
        targetAccept != bools.end(); targetAccept++) {
    for (std::vector<bool>::const_iterator sourceAccept = bools.begin();
        sourceAccept != bools.end(); sourceAccept++) {
     for (std::vector<bool>::const_iterator targetDonor = bools.begin();
          targetDonor != bools.end(); targetDonor++) {
      for (std::vector<bool>::const_iterator sourceDonor = bools.begin();
           sourceDonor != bools.end(); sourceDonor++) {
       for (std::vector<bool>::const_iterator sourceDummy = bools.begin();
            sourceDummy != bools.end(); sourceDummy++) {
        for (std::vector<bool>::const_iterator targetDummy = bools.begin();
             targetDummy != bools.end(); targetDummy++) {
         for (std::vector<bool>::const_iterator sourceIon = bools.begin();
              sourceIon != bools.end(); sourceIon++) {
          for (std::vector<bool>::const_iterator targetIon = bools.begin();
               targetIon != bools.end(); targetIon++) {
           for (std::vector<bool>::const_iterator ignoreIons = bools.begin();
                ignoreIons != bools.end(); ignoreIons++) {
            for (std::vector<bool>::const_iterator onlyBumps = bools.begin();
                 onlyBumps != bools.end(); onlyBumps++) {
              for (std::vector<bool>::const_iterator excludeAtom = bools.begin();
                   excludeAtom != bools.end(); excludeAtom++) {

                //================================================================
                // Test the scoring for various cases to ensure that they all behave as expected
                unsigned int atomSeq = 0;

                // Construct and fill the SpatialQuery information
                // with a vector of a single target atom, including its extra info looked up by
                // its i_seq value.
                iotbx::pdb::hierarchy::atom a;
                a.set_charge(targetCharge->c_str());
                a.set_xyz(vec3( 0,0,0 ));
                a.set_occ(1);
                if (*targetIon) {
                  a.set_element("CU");
                } else {
                  a.set_element("O");
                }
                a.data->i_seq = atomSeq++;
                scitbx::af::shared<iotbx::pdb::hierarchy::atom> atoms;
                atoms.push_back(a);
                SpatialQuery sq(atoms);
                ExtraAtomInfo e(targetRad, *targetAccept, *targetDonor, *targetDummy, *targetIon,
                  atom_charge(a));
                scitbx::af::shared<ExtraAtomInfo> infos;
                infos.push_back(e);

                // Construct the source atom, including its extra info.
                iotbx::pdb::hierarchy::atom source;
                source.set_charge(sourceCharge->c_str());
                source.set_occ(1);
                if (*sourceIon) {
                  source.set_element("CU");
                } else {
                  source.set_element("O");
                }
                source.data->i_seq = atomSeq++;
                ExtraAtomInfo se(sourceRad, *sourceAccept, *sourceDonor, *sourceDummy, *sourceIon,
                  atom_charge(source));
                atoms.push_back(source);
                infos.push_back(se);

                // Construct the scorer to be used.
                DotScorer as(ExtraAtomInfoMap(atoms, infos), 0.25, 10.0, 4.0, 0.6, 0.8, 0.4, 0.5, 0.25, false,
                  *ignoreIons);

                // Determine our hydrogen-bond state
                bool compatibleCharge = atom_charge(source) * atom_charge(a) <= 0;
                bool compatible = (*sourceDonor && *targetAccept) || (*sourceAccept && *targetDonor);
                bool hBond = compatibleCharge && compatible;

                // If we have an excluded atom, we should always get no values or bumping.
                // Skip the remainder of the tests in this case
                scitbx::af::shared<iotbx::pdb::hierarchy::atom> exclude;
                if (*excludeAtom) {
                  // Describe the extra atom to the system, including its extra info.
                  iotbx::pdb::hierarchy::atom ea;
                  ea.set_xyz(vec3( 0,0,0 ));
                  ea.set_occ(1);
                  ea.data->i_seq = atomSeq++;
                  ExtraAtomInfo ex(targetRad + 0.2, *targetAccept, *targetDonor, *targetDummy);
                  atoms.push_back(ea);
                  infos.push_back(ex);
                  exclude.push_back(ea);

                  // We added an atom, so we need a new DotScorer
                  DotScorer as2(ExtraAtomInfoMap(atoms, infos));

                  // Even when we have a close clash, we should get no response.
                  source.set_xyz(vec3( sourceRad,0,0 ));
                  ScoreDotsResult res = as2.score_dots(source, 1, sq, sourceRad + targetRad,
                    probeRad, exclude, ds.dots(), ds.density(), *onlyBumps);
                  if (!res.valid) {
                    return "DotScorer::test(): Could not score dots for excluded-atom case";
                  }
                  if ((res.totalScore() != 0) || res.hasBadBump) {
                    // Acceptor target atoms are excluded from calculations when the source atom
                    // is a dummy Hydrogen.
                    if (!(*sourceDummy && *targetAccept)) {
                      return "DotScorer::test(): Got unexpected result for excluded-atom case";
                    }
                  }

                  // Skip the rest of the tests for this case.
                  continue;
                }

                // If we have an ion and we are ignoring ions, we should always get no bumping.
                if ((*ignoreIons) && (*sourceIon || *targetIon)) {
                  if (!hBond) {
                    // Even when we have a close clash, we should get no response.
                    source.set_xyz(vec3(sourceRad, 0, 0));
                    ScoreDotsResult res = as.score_dots(source, 1, sq, sourceRad + targetRad,
                      probeRad, exclude, ds.dots(), ds.density(), *onlyBumps);
                    if (!res.valid) {
                      return "DotScorer::test(): Could not score dots for ion case";
                    }
                    if ((res.bumpSubScore != 0) || res.hasBadBump) {
                      return "DotScorer::test(): Got unexpected result for ion case";
                    }
                  }
                  // Skip the rest of the tests for this case.
                  continue;
                }

                // If we have a dummy hydrogen and we cannot be a hydrogen-bond pair,
                // we should always get no bumping.
                if (*sourceDummy || *targetDummy) {
                  if (!hBond) {
                    // Even when we have a close clash, we should get no response.
                    source.set_xyz(vec3( sourceRad,0,0 ));
                    ScoreDotsResult res = as.score_dots(source, 1, sq, sourceRad + targetRad,
                      probeRad, exclude, ds.dots(), ds.density(), *onlyBumps);
                    if (!res.valid) {
                      return "DotScorer::test(): Could not score dots for dummy hydrogen case";
                    }
                    if ((res.bumpSubScore != 0) || res.hasBadBump) {
                      return "DotScorer::test(): Got unexpected result for dummy hydrogen case";
                    }
                  }
                  // Skip the rest of the tests for this case.
                  continue;
                }

                // When we get so close that the source atom radius touches the center of the target,
                // we should get bad bumps in all cases.
                {
                  source.set_xyz(vec3( sourceRad,0,0 ));
                  ScoreDotsResult res = as.score_dots(source, 1, sq, sourceRad + targetRad,
                    probeRad, exclude, ds.dots(), ds.density(), *onlyBumps);
                  if (!res.valid) {
                    return "DotScorer::test(): Could not score dots for bad-bump case";
                  }
                  if (!res.hasBadBump) {
                    return "DotScorer::test(): Got no bad bump for bad-bump case";
                  }
                }

                // When we are only checking for bumps, we should get no interaction when the
                // atoms are not touching.  Otherwise, slight interaction.
                {
                  source.set_xyz(vec3( sourceRad + targetRad + 0.001,0,0 ));
                  ScoreDotsResult res = as.score_dots(source, 1, sq, sourceRad + targetRad,
                    probeRad, exclude, ds.dots(), ds.density(), *onlyBumps);
                  if (!res.valid) {
                    return "DotScorer::test(): Could not score dots for bump-only test case";
                  }
                  if (*onlyBumps) {
                    if (res.totalScore() != 0) {
                      return "DotScorer::test(): Got value when not expected for bump-only test case";
                    }
                  }
                  else {
                    if (res.totalScore() == 0) {
                      return "DotScorer::test(): Got no value when one expected for non-bump-only test case";
                    }
                  }
                }

                // When we are only checking for bumps, even hydrogen bonds should be counted as bumps.
                if (*onlyBumps) {
                  source.set_xyz(vec3( sourceRad + targetRad - 0.1,0,0 ));
                  ScoreDotsResult res = as.score_dots(source, 1, sq, sourceRad + targetRad,
                    probeRad, exclude, ds.dots(), ds.density(), *onlyBumps);
                  if (!res.valid) {
                    return "DotScorer::test(): Could not score dots for bump-only test case";
                  }
                  if (res.hBondSubScore != 0) {
                    return "DotScorer::test(): Got unexpected hydrogen bond score for bump-only test case";
                  }
                  if (res.bumpSubScore >= 0) {
                    return "DotScorer::test(): Got unexpected bump score for bump-only test case";
                  }
                }

              }
            }
           }
          }
         }
        }
       }
      }
     }
    }
   }
  }
  }

  // Test behavior of weak hydrogen bonds and their interaction with dummy Hydrogens.
  for (std::vector<bool>::const_iterator weakHBonds = bools.begin();
       weakHBonds != bools.end(); weakHBonds++) {
    for (std::vector<bool>::const_iterator sourceDummy = bools.begin();
         sourceDummy != bools.end(); sourceDummy++) {
      for (std::vector<bool>::const_iterator targetDummy = bools.begin();
           targetDummy != bools.end(); targetDummy++) {
        // Test the scoring for various cases to ensure that they all behave as expected
        unsigned int atomSeq = 0;

        // Construct and fill the SpatialQuery information
        // with a vector of a single target atom, including its extra info looked up by
        // its i_seq value.
        // Make the charges always compatible (one + one -) and specify the acceptor and
        // donor states based on whether the atoms are dummy hydrogens -- if so they are
        // donors and if not they are acceptors.
        iotbx::pdb::hierarchy::atom a;
        a.set_charge("-");
        a.set_xyz(vec3( 0,0,0 ));
        a.set_occ(1);
        a.data->i_seq = atomSeq++;
        scitbx::af::shared<iotbx::pdb::hierarchy::atom> atoms;
        atoms.push_back(a);
        SpatialQuery sq(atoms);
        ExtraAtomInfo e(targetRad, !*targetDummy, *targetDummy, *targetDummy, false, atom_charge(a));
        scitbx::af::shared<ExtraAtomInfo> infos;
        infos.push_back(e);

        // Construct the source atom, including its extra info.
        iotbx::pdb::hierarchy::atom source;
        source.set_charge("+");
        source.set_occ(1);
        source.data->i_seq = atomSeq++;
        ExtraAtomInfo se(sourceRad, !*sourceDummy, *sourceDummy, *sourceDummy, false, atom_charge(source));
        atoms.push_back(source);
        infos.push_back(se);

        // Empty excluded-atom list.
        scitbx::af::shared<iotbx::pdb::hierarchy::atom> exclude;

        // Construct the scorer to be used.
        DotScorer as(ExtraAtomInfoMap(atoms, infos), 0.25, 10.0, 4.0, 0.6, 0.8, 0.4, 0.5, 0.25, *weakHBonds);

        // Set the atoms to be slightly separated and then score the dots.
        source.set_xyz(vec3( sourceRad + targetRad + 0.1,0,0 ));
        ScoreDotsResult res = as.score_dots(source, 1, sq, sourceRad + targetRad,
          probeRad, exclude, ds.dots(), ds.density(), false);
        if (!res.valid) {
          return "DotScorer::test(): Could not score dots for weak hydrogen-bond test case";
        }

        // We should get nonzero scores when we're using weak hydrogen bonds and we
        // have exactly one of the two atoms being dummies.  We should also get them
        // when we're not using any dummies.
        bool expectNonzero = *weakHBonds && (*sourceDummy != *targetDummy);
        expectNonzero |= (!*sourceDummy && !*targetDummy);
        if (expectNonzero != (res.totalScore() != 0)) {
          return "DotScorer::test(): Got unexpected score for weak hydrogen-bond test case";
        }
      }
    }
  }

  // Test a polar hydrogen that is a donor overlapping with a target that is not an
  // acceptor (another polar hydrogen) to ensure we get a net negative result.
  {
      double targetRad = 1.05, sourceRad = 1.005, probeRad = 0.25;
      DotSphere ds(sourceRad, 200);
      unsigned int atomSeq = 0;

      // Construct and fill the SpatialQuery information
      // with a vector of a single target atom, including its extra info looked up by
      // its i_seq value.
      iotbx::pdb::hierarchy::atom a;
      a.set_xyz(vec3(0, 0, 0));
      a.set_occ(1);
      a.data->i_seq = atomSeq++;
      scitbx::af::shared<iotbx::pdb::hierarchy::atom> atoms;
      atoms.push_back(a);
      SpatialQuery sq(atoms);
      ExtraAtomInfo e(targetRad, false, true);
      scitbx::af::shared<ExtraAtomInfo> infos;
      infos.push_back(e);

      // Construct the source atom, including its extra info.
      iotbx::pdb::hierarchy::atom source;
      source.set_occ(1);
      source.data->i_seq = atomSeq++;
      ExtraAtomInfo se(sourceRad, false, true);
      atoms.push_back(source);
      infos.push_back(se);

      // Construct an empty exclusion list.
      scitbx::af::shared<iotbx::pdb::hierarchy::atom> exclude;

      // Construct the scorer to be used.
      DotScorer as(ExtraAtomInfoMap(atoms, infos));

      // Test the source atom with an overlap of magnitude 0.3.
      double gap = -0.3;
      source.set_xyz(vec3(targetRad + sourceRad + gap, 0, 0));
      ScoreDotsResult res = as.score_dots(source, 1, sq, sourceRad + targetRad,
          probeRad, exclude, ds.dots(), ds.density(), false);
      if (!res.valid) {
          return "DotScorer::test(): Could not score dots for hydrogen bumping case";
      }
      if ((res.bumpSubScore >= 0) || (res.totalScore() >= 0) || (res.attractSubScore <= 0)) {
          return "DotScorer::test(): Unexpected scores for hydrogen bumping case";
      }
  }

  // Sweep an atom from just touching to far away and make sure the attraction
  // curve is monotonically decreasing to 0 and does not reach 0 until the atom
  // is two probe radii away.
  {
    double targetRad = 1.5, sourceRad = 1.0, probeRad = 0.25;
    DotSphere ds(sourceRad, 200);
    unsigned int atomSeq = 0;

    // Construct and fill the SpatialQuery information
    // with a vector of a single target atom, including its extra info looked up by
    // its i_seq value.
    iotbx::pdb::hierarchy::atom a;
    a.set_xyz(vec3( 0,0,0 ));
    a.set_occ(1);
    a.data->i_seq = atomSeq++;
    scitbx::af::shared<iotbx::pdb::hierarchy::atom> atoms;
    atoms.push_back(a);
    SpatialQuery sq(atoms);
    ExtraAtomInfo e(targetRad);
    scitbx::af::shared<ExtraAtomInfo> infos;
    infos.push_back(e);

    // Construct the source atom, including its extra info.
    iotbx::pdb::hierarchy::atom source;
    source.set_occ(1);
    source.data->i_seq = atomSeq++;
    ExtraAtomInfo se(sourceRad);
    atoms.push_back(source);
    infos.push_back(se);

    // Construct an empty exclusion list.
    scitbx::af::shared<iotbx::pdb::hierarchy::atom> exclude;

    // Construct the scorer to be used.
    DotScorer as(ExtraAtomInfoMap(atoms, infos));

    // Sweep the source atom
    double lastAttract = 1e10;
    bool foundNonzero = false;
    for (double gap = 0; gap < 10; gap += 0.1) {
      source.set_xyz(vec3( targetRad + sourceRad + gap, 0, 0 ));
      ScoreDotsResult res = as.score_dots(source, 1, sq, sourceRad + targetRad,
        probeRad, exclude, ds.dots(), ds.density(), false);
      if (!res.valid) {
        return "DotScorer::test(): Could not score dots for swept-distance case";
      }
      if ((gap < 2 * probeRad) && (res.totalScore() == 0)) {
        return "DotScorer::test(): Got zero score for swept-distance case when within range";
      }
      if ((res.attractSubScore != res.totalScore()) || (res.attractSubScore > lastAttract)) {
        return "DotScorer::test(): Non-monotonic scores for swept-distance case";
      }
      lastAttract = res.attractSubScore;
      if (lastAttract != 0) { foundNonzero = true; }
    }
    if (!foundNonzero) {
      return "DotScorer::test(): No nonzero scores for swept-distance case";
    }
    if (lastAttract != 0) {
      return "DotScorer::test(): Non-empty last score for swept-distance case";
    }
  }

  // Test the setting of weights for the various subscores.
  {
    double targetRad = 1.5, sourceRad = 1.0, probeRad = 0.25;
    DotSphere ds(sourceRad, 200);
    unsigned int atomSeq = 0;

    // Construct and fill the SpatialQuery information
    // with a vector of a two target atoms, including extra info looked up by
    // their i_seq values.  One is an acceptor and the other is not.
    scitbx::af::shared<iotbx::pdb::hierarchy::atom> atoms;
    scitbx::af::shared<ExtraAtomInfo> infos;

    { // The first is a bond acceptor at the origin
      iotbx::pdb::hierarchy::atom a;
      a.set_xyz(vec3( 0,0,0 ));
      a.set_occ(1);
      a.data->i_seq = atomSeq++;
      atoms.push_back(a);
      ExtraAtomInfo e(targetRad, true);
      infos.push_back(e);
    }

    { // The second is not a bond acceptor and is a bit less than the source atom away from the edge of the first.
      iotbx::pdb::hierarchy::atom a;
      a.set_xyz(vec3( 2 * targetRad + 2*(sourceRad*0.8),0,0 ));
      a.set_occ(1);
      a.data->i_seq = atomSeq++;
      atoms.push_back(a);
      ExtraAtomInfo e(targetRad, false);
      infos.push_back(e);
    }
    SpatialQuery sq(atoms);

    // Construct the source atom, including its extra info.
    // It is a hydrogen donor and is located halfway
    // between the two atoms.
    iotbx::pdb::hierarchy::atom source;
    source.set_xyz(vec3( targetRad + sourceRad * 0.8 ,0,0 ));
    source.set_occ(1);
    source.data->i_seq = atomSeq++;
    ExtraAtomInfo se(sourceRad, false, true);
    atoms.push_back(source);
    infos.push_back(se);

    // Construct an empty exclusion list.
    scitbx::af::shared<iotbx::pdb::hierarchy::atom> exclude;

    // Results holds a vector of six values, three for the gap, bump, and bond weights used in
    // each run and three for the attract, bump, and hBond scores measured with these weights.
    std::vector< std::vector<double> > results;

    // Run tests with a hydrogen source molecule that can H-bond with one of two target atoms so that
    // we get results in all three of the score results when we compare it.
    for (double wGap = 0.25; wGap < 10; wGap += 3) {
      for (double wBump = 0.25; wBump < 10; wBump += 3) {
        for (double wBond = 0.25; wBond < 10; wBond += 3) {

          // Construct the scorer to be used.
          DotScorer as(ExtraAtomInfoMap(atoms, infos),wGap, wBump, wBond);

          // Get the dot results and store all of the params and results in a vector
          ScoreDotsResult res = as.score_dots(source, 1, sq, sourceRad + targetRad,
            probeRad, exclude, ds.dots(), ds.density(), false);
          if (!res.valid) {
            return "DotScorer::test(): Could not score dots for weight scaling cases";
          }
          std::vector<double> rowValueArray;
          rowValueArray.push_back(wGap);
          rowValueArray.push_back(wBump);
          rowValueArray.push_back(wBond);
          rowValueArray.push_back(res.attractSubScore);
          rowValueArray.push_back(res.bumpSubScore);
          rowValueArray.push_back(res.hBondSubScore);
          std::vector<double> row;
          for (size_t i = 0; i < rowValueArray.size(); i++) {
            row.push_back(rowValueArray[i]);
          }
          results.push_back(row);
        }
      }
    }

    // Ensure that the ratio of the scores matches the ratio of the weights for each entry between
    // the first row and all other rows.
    for (size_t i = 1; i < results.size(); i++) {
      // The attraction score is not a linear weighting of the
      // dot scores, so we can only check those for equality or inequality.
      bool paramClose = closeTo(1, results[0][0] / results[i][0]);
      bool resultClose = closeTo(1, results[0][3 + 0] / results[i][3 + 0]);
      if (paramClose != resultClose) {
        return "DotScorer::test(): Inconsistent gap ratio for weight scaling case";
      }

      // Check the other two scores for linearity.
      for (size_t j = 1; j < 3; j++) {
        double paramRatio = results[0][j] / results[i][j];
        double resultRatio = results[0][3 + j] / results[i][3 + j];
        if (!closeTo(paramRatio, resultRatio)) {
          return "DotScorer::test(): Unequal ratio for weight scaling case";
        }
      }
    }

  }

  // Test the setting of bond-gap distances.
  {
    double bumpOverlap = 0.2, badOverlap = bumpOverlap + 0.2;
    double maxRegularHydrogenOverlap = badOverlap + 0.2;
    double maxChargedHydrogenOverlap = maxRegularHydrogenOverlap + 0.2;

    double targetRad = 1.5, sourceRad = 1.0, probeRad = 0.25;
    DotSphere ds(sourceRad, 200);
    unsigned int atomSeq = 0;

    // Construct and fill the SpatialQuery information
    // with a vector of a single target atom, including its extra info looked up by
    // its i_seq value.  It will be an acceptor and it will have a negative charge
    iotbx::pdb::hierarchy::atom a;
    a.set_xyz(vec3( 0,0,0 ));
    a.set_occ(1);
    a.set_charge("-");
    a.data->i_seq = atomSeq++;
    scitbx::af::shared<iotbx::pdb::hierarchy::atom> atoms;
    atoms.push_back(a);
    SpatialQuery sq(atoms);
    ExtraAtomInfo e(targetRad, true, false, false, false, atom_charge(a));
    scitbx::af::shared<ExtraAtomInfo> infos;
    infos.push_back(e);

    // Construct an empty exclusion list.
    scitbx::af::shared<iotbx::pdb::hierarchy::atom> exclude;

    // Construct a source atom, including its extra info.
    // This will be a hydrogen but not a donor to check
    // for the standard bad-bump result.
    // Test it against both sides of the bad-bump line to see if it responds correctly.
    {
      iotbx::pdb::hierarchy::atom source;
      source.set_occ(1);
      source.data->i_seq = atomSeq++;
      ExtraAtomInfo se(sourceRad);
      atoms.push_back(source);
      infos.push_back(se);
      ScoreDotsResult res;

      // Construct the scorer to be used with the specified bond gaps.
      DotScorer as(ExtraAtomInfoMap(atoms, infos), 0.25, 10.0, 4.0, maxRegularHydrogenOverlap,
        maxChargedHydrogenOverlap, bumpOverlap, badOverlap);

      // Check the source atom against outside and inside the gap
      source.set_xyz(vec3( targetRad + sourceRad - badOverlap + 0.1, 0, 0 ));
      res = as.score_dots(source, 1, sq, sourceRad + targetRad,
        probeRad, exclude, ds.dots(), ds.density(), false);
      if (!res.valid) {
        return "DotScorer::test(): Could not score dots for badOverlap setting case";
      }
      if (res.hasBadBump) {
        return "DotScorer::test(): Bad bump found when not expected for badOverlap setting case";
      }

      source.set_xyz(vec3( targetRad + sourceRad - badOverlap - 0.1, 0, 0 ));
      res = as.score_dots(source, 1, sq, sourceRad + targetRad,
        probeRad, exclude, ds.dots(), ds.density(), false);
      if (!res.valid) {
        return "DotScorer::test(): Could not score dots for badOverlap setting case";
      }
      if (!res.hasBadBump) {
        return "DotScorer::test(): Bad bump not found when expected for badOverlap setting case";
      }
    }

    // Construct a source atom, including its extra info.
    // This will be an uncharged hydrogen donor to check
    // for the non-charged hydrogen bad-bump result.
    // Test it against both sides of the bad-bump line to see if it responds correctly.
    {
      iotbx::pdb::hierarchy::atom source;
      source.set_occ(1);
      source.set_charge("0");
      source.data->i_seq = atomSeq++;
      ExtraAtomInfo se(sourceRad, false, true, false, false, atom_charge(source));
      atoms.push_back(source);
      infos.push_back(se);
      ScoreDotsResult res;

      // Construct the scorer to be used with the specified bond gaps.
      DotScorer as(ExtraAtomInfoMap(atoms, infos), 0.25, 10.0, 4.0, maxRegularHydrogenOverlap,
        maxChargedHydrogenOverlap, bumpOverlap, badOverlap);

      // Check the source atom against outside and inside the gap
      source.set_xyz(vec3(targetRad + sourceRad - maxRegularHydrogenOverlap + 0.1, 0, 0));
      res = as.score_dots(source, 1, sq, sourceRad + targetRad,
        probeRad, exclude, ds.dots(), ds.density(), false);
      if (!res.valid) {
        return "DotScorer::test(): Could not score dots for maxRegularHydrogenOverlap bump setting case";
      }
      if (res.bumpSubScore != 0) {
        return "DotScorer::test(): Bump found when not expected for maxRegularHydrogenOverlap setting case";
      }

      source.set_xyz(vec3(targetRad + sourceRad - maxRegularHydrogenOverlap - 0.1, 0, 0));
      res = as.score_dots(source, 1, sq, sourceRad + targetRad,
        probeRad, exclude, ds.dots(), ds.density(), false);
      if (!res.valid) {
        return "DotScorer::test(): Could not score dots for maxRegularHydrogenOverlap bump setting case";
      }
      if (res.bumpSubScore >= 0) {
        return "DotScorer::test(): Bump not found when expected for maxRegularHydrogenOverlap setting case";
      }

      // Check the source atom against badly outside and inside the gap
      source.set_xyz(vec3( targetRad + sourceRad - maxRegularHydrogenOverlap - badOverlap + 0.1, 0, 0 ));
      res = as.score_dots(source, 1, sq, sourceRad + targetRad,
        probeRad, exclude, ds.dots(), ds.density(), false);
      if (!res.valid) {
        return "DotScorer::test(): Could not score dots for maxRegularHydrogenOverlap setting case";
      }
      if (res.hasBadBump) {
        return "DotScorer::test(): Bad bump found when not expected for maxRegularHydrogenOverlap setting case";
      }

      source.set_xyz(vec3( targetRad + sourceRad - maxRegularHydrogenOverlap - badOverlap - 0.1, 0, 0 ));
      res = as.score_dots(source, 1, sq, sourceRad + targetRad,
        probeRad, exclude, ds.dots(), ds.density(), false);
      if (!res.valid) {
        return "DotScorer::test(): Could not score dots for maxRegularHydrogenOverlap setting case";
      }
      if (!res.hasBadBump) {
        return "DotScorer::test(): Bad bump not found when expected for maxRegularHydrogenOverlap setting case";
      }
    }

    // Construct a source atom, including its extra info.
    // This will be a charged hydrogen donor to check
    // for the charged hydrogen bad-bump result.
    // Test it against both sides of the bad-bump line to see if it responds correctly.
    {
      iotbx::pdb::hierarchy::atom source;
      source.set_occ(1);
      source.set_charge("+");
      source.data->i_seq = atomSeq++;
      ExtraAtomInfo se(sourceRad, false, true, false, false, atom_charge(source));
      atoms.push_back(source);
      infos.push_back(se);
      ScoreDotsResult res;

      // Construct the scorer to be used with the specified bond gaps.
      DotScorer as(ExtraAtomInfoMap(atoms, infos), 0.25, 10.0, 4.0, maxRegularHydrogenOverlap,
        maxChargedHydrogenOverlap, bumpOverlap, badOverlap);

      // Check the source atom against outside and inside the gap
      source.set_xyz(vec3( targetRad + sourceRad - maxChargedHydrogenOverlap - badOverlap + 0.1, 0, 0 ));
      res = as.score_dots(source, 1, sq, sourceRad + targetRad,
        probeRad, exclude, ds.dots(), ds.density(), false);
      if (!res.valid) {
        return "DotScorer::test(): Could not score dots for maxChargedHydrogenOverlap setting case";
      }
      if (res.hasBadBump) {
        return "DotScorer::test(): Bad bump found when not expected for maxChargedHydrogenOverlap setting case";
      }

      source.set_xyz(vec3( targetRad + sourceRad - maxChargedHydrogenOverlap - badOverlap - 0.1, 0, 0 ));
      res = as.score_dots(source, 1, sq, sourceRad + targetRad,
        probeRad, exclude, ds.dots(), ds.density(), false);
      if (!res.valid) {
        return "DotScorer::test(): Could not score dots for maxChargedHydrogenOverlap setting case";
      }
      if (!res.hasBadBump) {
        return "DotScorer::test(): Bad bump not found when expected for maxChargedHydrogenOverlap setting case";
      }
    }
  }

  // Test the control of occupancy level
  {
    double targetRad = 1.5, sourceRad = 1.0, probeRad = 0.25;
    DotSphere ds(sourceRad, 200);
    unsigned int atomSeq = 0;

    // Construct and fill the SpatialQuery information
    // with a vector of a single target atom, including its extra info looked up by
    // its i_seq value.  It will have an occupancy of 0.5.
    iotbx::pdb::hierarchy::atom a;
    a.set_xyz(vec3( 0,0,0 ));
    a.set_occ(0.5);
    a.data->i_seq = atomSeq++;
    scitbx::af::shared<iotbx::pdb::hierarchy::atom> atoms;
    atoms.push_back(a);
    SpatialQuery sq(atoms);
    ExtraAtomInfo e(targetRad, true);
    scitbx::af::shared<ExtraAtomInfo> infos;
    infos.push_back(e);

    // Construct an empty exclusion list.
    scitbx::af::shared<iotbx::pdb::hierarchy::atom> exclude;

    // Construct a source atom, including its extra info.
    // This will be a hydrogen but not a donor.
    // It will have an occupancy of 0.25.
    iotbx::pdb::hierarchy::atom source;
    source.set_occ(0.25);
    source.data->i_seq = atomSeq++;
    ExtraAtomInfo se(sourceRad);
    atoms.push_back(source);
    infos.push_back(se);
    ScoreDotsResult res;

    // Construct the scorer to be used with the specified bond gaps.
    DotScorer as(ExtraAtomInfoMap(atoms, infos));

    // Check the source atom just overlapping with the target.
    source.set_xyz(vec3( targetRad + sourceRad - 0.1, 0, 0 ));

    // Occupancy values to check and whether the result should be zero:
    static const double occupancyArray[] = { 1.0, 0.75, 0.5, 0.25, 0.1 };
    static const double expectZeroArray[] = { true, true, true, false, false };
    std::vector<double> occupancies;
    std::vector<bool> expectZero;
    for (size_t i = 0; i < sizeof(occupancyArray)/sizeof(occupancyArray[0]); i++) {
      occupancies.push_back(occupancyArray[i]);
      expectZero.push_back(expectZeroArray[i]);
    }
    for (size_t i = 0; i < occupancies.size(); i++) {
      res = as.score_dots(source, occupancies[i], sq, sourceRad + targetRad,
        probeRad, exclude, ds.dots(), ds.density(), false);
      if (!res.valid) {
        return "DotScorer::test(): Could not score dots for occupancy test case";
      }
      bool zero = res.totalScore() == 0;
      if (zero != expectZero[i]) {
        return "DotScorer::test(): Unexpected zero state for occupancy test case";
      }
    }
  }

  // Check with invalid parameters.
  {
    double targetRad = 1.5, sourceRad = 1.0, probeRad = 0.25;
    DotSphere ds(sourceRad, 200);
    unsigned int atomSeq = 0;

    // Construct and fill the SpatialQuery information
    // with a vector of a single target atom, including its extra info.
    iotbx::pdb::hierarchy::atom a;
    a.set_xyz(vec3( 0,0,0 ));
    a.set_occ(1);
    a.data->i_seq = atomSeq++;
    scitbx::af::shared<iotbx::pdb::hierarchy::atom> atoms;
    atoms.push_back(a);
    SpatialQuery sq(atoms);
    ExtraAtomInfo e(targetRad, true);
    scitbx::af::shared<ExtraAtomInfo> infos;
    infos.push_back(e);

    // Construct an empty exclusion list.
    scitbx::af::shared<iotbx::pdb::hierarchy::atom> exclude;

    // Construct a source atom, including its extra info.
    // This will be a hydrogen but not a donor.
    iotbx::pdb::hierarchy::atom source;
    source.set_occ(1);
    source.data->i_seq = atomSeq++;
    ExtraAtomInfo se(sourceRad);
    atoms.push_back(source);
    infos.push_back(se);
    ScoreDotsResult res;

    // Construct the scorer to be used with the specified bond gaps.
    DotScorer as(ExtraAtomInfoMap(atoms, infos));

    // Check the source atom just overlapping with the target.
    source.set_xyz(vec3( targetRad + sourceRad - 0.1, 0, 0 ));

    res = as.score_dots(source, 1, sq, sourceRad + targetRad,
      -0.1, exclude, ds.dots(), ds.density(), false);
    if (res.valid) {
      return "DotScorer::test(): Unexpected valid result for probeRadius < 0 case";
    }
    res = as.score_dots(source, 1, sq, sourceRad + targetRad,
      probeRad, exclude, ds.dots(), 0, false);
    if (res.valid) {
      return "DotScorer::test(): Unexpected valid result for density <= 0 case";
    }
  }

  // Check surface-dot counting.
  {
    double targetRad = 1.5, sourceRad = 1.0;
    DotSphere ds(sourceRad, 200);
    unsigned int atomSeq = 0;
    unsigned count;

    // Construct a single target atom, including its extra info
    iotbx::pdb::hierarchy::atom a;
    a.set_xyz(vec3( 0,0,0 ));
    a.set_occ(1);
    a.data->i_seq = atomSeq++;
    scitbx::af::shared<iotbx::pdb::hierarchy::atom> atoms;
    atoms.push_back(a);
    SpatialQuery sq(atoms);
    ExtraAtomInfo e(targetRad, true);
    scitbx::af::shared<ExtraAtomInfo> infos;
    infos.push_back(e);

    // Construct a source atom, including its extra info.
    // This will be a hydrogen but not a donor.
    iotbx::pdb::hierarchy::atom source;
    source.set_occ(1);
    source.data->i_seq = atomSeq++;
    ExtraAtomInfo se(sourceRad);
    atoms.push_back(source);
    infos.push_back(se);

    // Construct the scorer to be used with the specified bond gaps.
    DotScorer as(ExtraAtomInfoMap(atoms, infos));

    // Construct an exclusion list and add the target atom.
    scitbx::af::shared<iotbx::pdb::hierarchy::atom> exclude;
    exclude.push_back(a);

    // Check the source atom not overlapping with the target.
    source.set_xyz(vec3( targetRad + sourceRad + 0.1, 0, 0 ));
    count = as.count_surface_dots(source, ds.dots(), exclude);
    if (count != ds.dots().size()) {
      return "DotScorer::test(): Unexpected surface dot count for non-overlapping case";
    }

    // Check the source atom completely overlapping with the target.
    source.set_xyz(vec3( 0, 0, 0 ));
    count = as.count_surface_dots(source, ds.dots(), exclude);
    if (count != 0) {
      return "DotScorer::test(): Unexpected nonzero surface dot count for fully-overlapping case";
    }

    // Check the source atom partially overlapping with the target.
    source.set_xyz(vec3( targetRad, 0, 0 ));
    count = as.count_surface_dots(source, ds.dots(), exclude);
    if ((count == 0) || (count >= ds.dots().size())) {
      return "DotScorer::test(): Unexpected surface dot count for partially-overlapping case";
    }
  }

  return "";
}

std::string Scoring_test()
{
  std::string ret;
  ContactResult res;

  // Test the atom-charge code.
  {
    static const char* chargesArray[] = { "--", "-", "", "+", "++", "+2", "-1", "0" };
    static const int expectedChargeArray[] = { -2, -1, 0, 1, 2, 2, -1, 0 };
    std::vector<std::string> charges;
    std::vector<int> expectedCharge;
    for (size_t i = 0; i < sizeof(chargesArray) / sizeof(chargesArray[0]); i++) {
      charges.push_back(chargesArray[i]);
      expectedCharge.push_back(expectedChargeArray[i]);
    }

    for (size_t i = 0; i < charges.size(); i++) {
      iotbx::pdb::hierarchy::atom a;
      a.set_charge(charges[i].c_str());
      if (atom_charge(a) != expectedCharge[i]) {
        return "Scoring_test: atom_charge() failed";
      }
    }
  }

  // Check the ExtraAtomInfoMap
  {
    double targetRad = 1.5, sourceRad = 1.0;
    unsigned int atomSeq = 0;

    // Construct a single target atom, including its extra info
    iotbx::pdb::hierarchy::atom a;
    a.set_xyz(vec3(0, 0, 0));
    a.set_occ(1);
    a.data->i_seq = atomSeq++;
    scitbx::af::shared<iotbx::pdb::hierarchy::atom> atoms;
    atoms.push_back(a);
    ExtraAtomInfo e(targetRad, true);
    scitbx::af::shared<ExtraAtomInfo> infos;
    infos.push_back(e);

    // Construct a source atom, including its extra info.
    iotbx::pdb::hierarchy::atom source;
    source.set_occ(1);
    source.data->i_seq = atomSeq++;
    ExtraAtomInfo se(sourceRad);
    atoms.push_back(source);
    infos.push_back(se);

    // Construct an ExtraAtomInfoMap and verify that its lookups work.
    ExtraAtomInfoMap eam(atoms, infos);
    if (eam.getMappingFor(a).getVdwRadius() != targetRad) {
      return "Scoring_test(): Unexpected original target radius from ExtraAtomInfoMap";
    }
    if (eam.getMappingFor(a).getIsAcceptor() != true) {
      return "DotScorer::test(): Unexpected target acceptor status from ExtraAtomInfoMap";
    }
    if (eam.getMappingFor(source).getVdwRadius() != sourceRad) {
      return "Scoring_test(): Unexpected source radius from ExtraAtomInfoMap";
    }

    // Modify the value of one of the atoms and verify that the lookup still works.
    ExtraAtomInfo e2(targetRad);
    eam.setMappingFor(a, e2);
    if (eam.getMappingFor(a).getIsAcceptor() != false) {
      return "Scoring_test(): Unexpected changed target acceptor status from ExtraAtomInfoMap";
    }
  }

  // Test the DotScorer class
  ret = DotScorer::test();
  if (ret.size() > 0) {
    return std::string("Scoring_test: DotScorer failed: ") + ret;
  }

  // Test that the distance from a location to itself is negative radius, and that
  // the projection is the radius away from the origin.
  Point ones(1, 1, 1);
  res = closest_contact(ones, ones, 1);
  if (!closeTo(res.distAboveSurface, -1)) {
    return "Scoring_test: closest_contact() returned bad distance for same point";
  }
  if (!closeTo((ones - res.closestContact).length(), 1)) {
    return "Scoring_test: closest_contact() returned bad closest contact for same point";
  }

  // Test that the projection distance is the same for all points on the cardinal directions
  // and diagonal directions around a sphere.
  for (int x = -1; x <= 1; x++) {
    for (int y = -1; y <= 1; y++) {
      for (int z = -1; z <= 1; z++) {
        if (x != 0 || y != 0 || z != 0) {
          Point test(x, y, z);
          test += ones;
          double rad = 0.5;
          res = closest_contact(test, ones, rad);
          if (!closeTo(rad, (res.closestContact - ones).length())) {
            return "Scoring_test: closest_contact() returned bad closest contact for surrounding point";
          }
        }
      }
    }
  }

  // Test that the distance scales as expected as we move away from the center of an atom
  for (int x = 0; x < 100; x++) {
    Point test(x, 0, 0);
    test += ones;
    double rad = 25;
    res = closest_contact(test, ones, rad);
    if (!closeTo(x - rad, res.distAboveSurface)) {
      return "Scoring_test: closest_contact() returned bad distance for point";
    }
  }

  /// @todo Test the overlap scale factor

  // All tests passed.
  return "";
}


} // end namespace probe
} // end namespace molprobity
