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

#pragma once

#include "SpatialQuery.h"
#include <iotbx/pdb/hierarchy.h>
#include <vector>
#include <map>
#include <boost/python.hpp>

namespace molprobity {
  namespace probe {

    //=====================================================================================================
    // Helper functions outside the class.

    /// @brief Structure to hold the results from a call to closest_contact()
    class ContactResult {
    public:
      Point   closestContact;     ///< The point on the radius of the tested sphere closest to the dot
      double  distAboveSurface;   ///< Distance that the dot is above the tested sphere (negative for inside)
    };

    /// @brief Find the point of closest contact and distance above/below atomic surface
    /// @param [in] dot The dot whose location is to be projected onto the atom
    /// @param [in] atom The center of the atom
    /// @param [in] atom_radius The radius of the atom
    /// @return The point of closest contact on the atom's surface to the dot and the
    ///       signed distance that the dot is above that surface, with negative indicating
    ///       that it is inside the atom's surface.
    ContactResult closest_contact(Point dot, Point atom, double atom_radius);

    /// @brief Return the signed-integer charge of the atom.
    /// @param atom Atom whose charge is to be determined.
    /// return Integer representing the charge: -2, -1, 0, 1, 2
    int atom_charge(iotbx::pdb::hierarchy::atom const& atom);

    // functions used to restrict annular rings of good dots around clashes

    /// @brief The distance from a dot to the point on the source surface that is closest to the target
    /// @param [in] dot The dot being considered
    /// @param [in] srcLoc The center of the source atom
    /// @param [in] srcVDWRad The van Der Waals radius of the source atom in Angstroms
    /// @param [in] targLoc The center of the target atom
    /// @return The distance from a dot to the point on the source surface that is closest to the target
    double dot2srcCenter(const Point& dot, const Point& srcLoc, double srcVDWRad, const Point& targLoc);

    /// @brief This is the distance from the point on the source atom closest to the target to the edge
    ///        of contact when the source and target are in optimal contact?
    ///        @todo Figure out for sure what this is computing.
    /// @param [in] ra Radius of the source atom
    /// @param [in] rb Radius of the target atom
    /// @param [in] rp Probe radius
    double kissEdge2bullsEye(double ra, double rb, double rp);

    /// @brief A dot is annular if it is further from the center of contact than edge of the overlap
    /// region is at optimum contact.
    ///
    /// This checks to make sure that dots that would not have contributed to a good score at optimium
    /// contact are not considered to contribute to a good score when the atoms are overlapping.
    /// This should be run on a dot that is not in contact to see if it is within the acceptable region
    /// or not (if not, it is annular).  When the source and target atom just touch ("kiss") at their
    /// edges, the probe radius will describe a circle on the source atom.  A dot is annular if it is
    /// further from the kiss location than this ring.  This can only happen for the condition where
    /// the atoms are overlapping because no dot will be constructed beyond this boundary for atoms
    /// that are just touching or not touching.
    /// @return true is returned to indicate that this dot is annular and should not be counted.
    bool annularDots(const Point& dot, const Point& srcLoc, double srcVDWRad,
      const Point& targLoc, double targVDWRad, double probeRadius);

    //=====================================================================================================
    /// @brief Class to hold data values for an atom beyond those present in the hierarchy::atom class itself
    /// that are needed by the Probe calculations.  These must be filled in by the client, perhaps by calling
    /// the mmtbx/probe/Helpers.py::getExtraAtomInfo() function.

    class ExtraAtomInfo {
    public:
      /// @brief Constructor with default parameters
      ExtraAtomInfo(double vdwRadius = 0, bool isAcceptor = false, bool isDonor = false,
        bool isDummyHydrogen = false, bool isIon = false, int charge = 0, std::string altLoc = " ")
        : m_vdwRadius(vdwRadius), m_isAcceptor(isAcceptor), m_isDonor(isDonor)
        , m_isDummyHydrogen(isDummyHydrogen), m_isIon(isIon), m_charge(charge)
        , m_altLoc(altLoc) {}
      /// @brief Constructor from another ExtraAtomInfo
      ExtraAtomInfo(const ExtraAtomInfo &e)
        : m_vdwRadius(e.m_vdwRadius), m_isAcceptor(e.m_isAcceptor), m_isDonor(e.m_isDonor)
        , m_isDummyHydrogen(e.m_isDummyHydrogen), m_isIon(e.m_isIon), m_charge(e.m_charge)
        , m_altLoc(e.m_altLoc) {}

      /// @brief Get and set methods
      double getVdwRadius() const { return m_vdwRadius; }
      void setVdwRadius(double val) { m_vdwRadius = val; }

      bool getIsAcceptor() const { return m_isAcceptor; }
      void setIsAcceptor(bool val) { m_isAcceptor = val; }

      bool getIsDonor() const { return m_isDonor; }
      void setIsDonor(bool val) { m_isDonor = val; }

      bool getIsDummyHydrogen() const { return m_isDummyHydrogen; }
      void setIsDummyHydrogen(bool val) { m_isDummyHydrogen = val; }

      bool getIsIon() const { return m_isIon; }
      void setIsIon(bool val) { m_isIon = val; }

      int getCharge() const { return m_charge; }
      void setCharge(int val) { m_charge = val; }

      std::string getAltLoc() const { return m_altLoc; }
      void setAltLoc(std::string val) { m_altLoc = val; }

      /// @brief == operator is required so that we can wrap the standard vector operators in Boost::Python
      bool operator ==(ExtraAtomInfo const& o) {
        return ((getVdwRadius() == o.getVdwRadius())
          && (getIsAcceptor() == o.getIsAcceptor())
          && (getIsDonor() == o.getIsDonor())
          && (getIsDummyHydrogen() == o.getIsDummyHydrogen())
          && (getIsIon() == o.getIsIon())
          && (getCharge() == o.getCharge())
          && (getAltLoc() == o.getAltLoc())
          );
      }
      bool operator !=(ExtraAtomInfo const& o) { return !(*this == o); }

    protected:
      double m_vdwRadius;      ///< van Der Waals radius of the atom
      bool m_isAcceptor;       ///< Does this accept hydrogen bonds?
      bool m_isDonor;          ///< Is this a donor hydrogen?
      bool m_isDummyHydrogen;  ///< These are inserted on Oxygens that are waters to provide potential
                               ///  hydrogen bonds to nearby acceptors.
      bool m_isIon;            ///< Is this an ion?
      int m_charge;            ///< The integer charge of the atom: -2, -1, 0, 1, 2
      std::string m_altLoc;    ///< The alternate location indicator for the atom
    };

    //=====================================================================================================
    /// @brief Class to map from an atom to the ExtraAtomInfo associated with it.

    class ExtraAtomInfoMap {
    public:
      /// @brief Constructor creates the map from vectors of atoms and ExtraAtomInfo.
      /// @param [in] atoms Vector of atoms that will be mapped from.  Must include all of the ones added
      ///             to the SpatialQuery structure so that anyone using its outputs will find an entry.
      ///             NOTE: This uses the i_seq values of atoms to index the vector, so there must be a unique
      ///             i_seq value for each atom (including Phantom hydrogens) and the i_seq must not
      ///             change between setting and getting.
      /// @param [in] extraInfo Vector of extra information pertaining to each atom.  Must be the same
      ///             length as the atoms vector.
      ExtraAtomInfoMap(scitbx::af::shared<iotbx::pdb::hierarchy::atom> const &atoms
        , scitbx::af::shared<ExtraAtomInfo> const &extraInfo)
      {
        if (atoms.size() == extraInfo.size()) {
          // Pre-allocated the vector to a reasonable size -- it will not be correct if we
          // have skipped i_seq values, but should be close.
          m_extraInfo.reserve(atoms.size());
          // Insert the values
          for (size_t i = 0; i < atoms.size(); i++) {
            setMappingFor(atoms[i], extraInfo[i]);
          }
        }
      }

      /// @brief Get and set methods
      ExtraAtomInfo  const &getMappingFor(iotbx::pdb::hierarchy::atom const &atom) const
      {
        if (atom.data->i_seq >= m_extraInfo.size()) {
          PyErr_SetString(PyExc_RuntimeError, "Out of bounds reference in ExtraAtomInfoMap::getMappingFor()");
          boost::python::throw_error_already_set();
        }
        return m_extraInfo[atom.data->i_seq];
      }
      void setMappingFor(iotbx::pdb::hierarchy::atom const &atom, ExtraAtomInfo const &info)
      {
        unsigned i_seq = atom.data->i_seq;
        if (m_extraInfo.size() < i_seq + 1) {
          m_extraInfo.resize(i_seq + 1);
        }
        m_extraInfo[i_seq] = info;
      }

    protected:
      // Vector indexed by i_seq
      std::vector<ExtraAtomInfo> m_extraInfo;
    };

    //=====================================================================================================
    /// @brief Class to handle scoring dots given sets of atoms.
    class DotScorer {
    public:
      /// @brief Constructor stores the non-bonded parameters to be used to determine atom features.
      /// @param [in] extraInfoMap maps from atoms to ExtraAtomInfo.
      /// @param [in] gapScale Scale factor to apply to gap between atoms (gap is divided by this)
      /// @param [in] bumpWeight Factor to apply when atoms are in bumping overlap
      /// @param [in] hBondWeight Factor to apply to hydrogen-bond overlaps
      /// @param [in] maxRegularHydrogenOverlap How much overlap can there be between a hydrogen
      ///             and the atom it is hydrogen-bonded to before we call it a clash when the atoms
      ///             are not both charged.  It must go badBumpOverlap beyond this before we call
      ///             it a bad clash.
      /// @param [in] maxChargedHydrogenOverlap How much overlap can there be between a hydrogen
      ///             and the atom it is hydrogen-bonded to before we call it a clash when the
      ///             atoms are both charged.  It must go badBumpOverlap beyond this before we call
      ///             it a bad clash.
      /// @param [in] bumpOverlap Atoms that overlap more than this will cause a clash.  This is a
      ///             positive number indicating how much overlap.
      /// @param [in] badBumpOverlap Atoms that overlap more than this will cause bad bump to be flagged.
      ///             This number can be set very large to cause no bumps to be flagged as bad bumps.  This is a
      ///             positive number indicating how much overlap.
      /// @param [in] contactCutoff Atoms that are nearer than this will be considered to be in near
      ///             contact, with atoms further in far contact.  This should be at least as large
      ///             as the probe radius.
      /// @param [in] weakHBonds Include weak hydrogen bonds (dots that are in hydrogen bonds but which
      ///             are outside contact).  These are categorized as a different InteractionType but
      ///             are counted as hydrogen bonds when scoring.
      /// @param [in] ignoreIonInteractions Ignore ions when computing interactions.  An ion source
      ///             atom will have no interactions and an ion will not be considered in the
      ///             target atom calculations (unless it is listed as an excluder).  This is because
      ///             as of 4/20/2022 Probe does not properly handle ionic interactions.
      /// @todo Consider moving the probe radius into the constructor parameters
      DotScorer(ExtraAtomInfoMap extraInfoMap
        , double gapScale = 0.25
        , double bumpWeight = 10.0
        , double hBondWeight = 4.0
        , double maxRegularHydrogenOverlap = 0.6
        , double maxChargedHydrogenOverlap = 0.8
        , double bumpOverlap = 0.4
        , double badBumpOverlap = 0.5
        , double contactCutoff = 0.25
        , bool weakHBonds = false
        , bool ignoreIonInteractions = false
      ) : m_extraInfoMap(extraInfoMap)
        , m_gapScale(gapScale), m_bumpWeight(bumpWeight), m_hBondWeight(hBondWeight)
        , m_maxRegularHydrogenOverlap(maxRegularHydrogenOverlap)
        , m_maxChargedHydrogenOverlap(maxChargedHydrogenOverlap)
        , m_bumpOverlap(bumpOverlap)
        , m_badBumpOverlap(badBumpOverlap)
        , m_contactCutoff(contactCutoff)
        , m_weakHBonds(weakHBonds)
        , m_ignoreIonInteractions(ignoreIonInteractions)
      {}

      /// @brief Tells whether the specified location is inside any of a list of atoms.
      ///
      /// The original Probe code does internal checks for Phantom Hydrogens, but we handle that
      /// in the calling routine by properly adjusting the list of atoms to be excluded.
      /// @param [in] location The location to check
      /// @param [in] atoms The list of atoms to check
      /// @return True if the location is inside any of the atoms, false otherwise
      bool point_inside_atoms(Point const &location, scitbx::af::shared<iotbx::pdb::hierarchy::atom> const &atoms);

      /// @brief Trim down a list of dots to only those that are not inside any of a list of atoms.
      /// @param [in] atom The center of the source atom to check
      /// @param [in] dots Vector of dot offsets from the center
      /// @param [in] exclude The list of atoms to check
      /// @return A vector of dots that are not inside any of the atoms
      scitbx::af::shared<Point> trim_dots(iotbx::pdb::hierarchy::atom const &atom,
        scitbx::af::shared<Point> const &dots,
        scitbx::af::shared<iotbx::pdb::hierarchy::atom> const &exclude);

      /// @brief Enumeration listing the basic types of overlap a dot can have with an atom.
      /// The values mean: NoOverlap => dot outside atom, Clash => dot inside atom and not hydrogen bonding
      /// (including too-close hydrogen), HydrogenBond => Hydrogen bond, Ignore = this dot was inside
      /// an excluded atom or had no neighboring atoms so should be ignored.
      enum OverlapType { Ignore = -2, Clash = -1, NoOverlap = 0, HydrogenBond = 1 };

      /// @brief Structure to hold the results from a call to check_dot()
      class CheckDotResult {
      public:
        CheckDotResult() : overlapType(DotScorer::Ignore), overlap(0), gap(1e100), annular(false) {};
        OverlapType overlapType;            ///< What kind of overlap, if any, was found
        iotbx::pdb::hierarchy::atom cause;  ///< Cause of the overlap, if overlapType != Ignore
        double  overlap;                    ///< Amount of overlap assigned to source atom if there is a clash
        double  gap;                        ///< Gap distance (overlap may only be a fraction of this).
        /// Was this an annular dot, which is further from the point on the source atom that is closest
        /// to the target atom than a point on the tangent edge of the target atom where the ray just
        /// grazes the surface of the atom.  @todo The math being done here is a bit opaque, so the
        /// origin of the second measurement is not clear.  The intent is to prevent over-spreading of
        /// dots beyond what would have been present when the atoms just touched, which prevents
        /// additional good contact outside the clash; see page 1717 of J. Mol. Biol. 285 (1999).
        bool    annular;
      };

      /// @brief Score an individual dot against a specific set of interacting atoms unless within an excluded atom
      /// @param [in] sourceAtom Atom that the dot is offset with respect to.
      /// @param [in] dotOffset Offset from the center of the source atom to the dot
      /// @param [in] probeRadius Radius of the probe rolled between the two potentially-contacting atoms
      ///             If this is < 0, an invalid result will be returned.
      /// @param [in] interacting The atoms that are to be checked because they are close enough to sourceAtom
      ///             to interact with it.  This list should not include any excluded atoms because spurious
      ///             interactions may occur along the edge of the excluded atom when a dot is close to it but
      ///             not inside.
      /// @param [in] excluded Atoms that are to be excluded from contact, for example this could be a list
      ///             of atoms bonded to sourceAtom.  If the dot is inside an excluded atom, it will not be
      ///             considered even if it is overlapping with an interacting atom.
      /// @param [in] overlapScale: The fraction of overlap to assign to each of the two atoms, scaling the
      ///             spike drawn for each.  The default value of 0.5 draws half of the spike for one atom
      ///             and the other half for the other.
      CheckDotResult check_dot(iotbx::pdb::hierarchy::atom const &sourceAtom,
        Point const& dotOffset, double probeRadius,
        scitbx::af::shared<iotbx::pdb::hierarchy::atom> const& interacting,
        scitbx::af::shared<iotbx::pdb::hierarchy::atom> const& exclude,
        double overlapScale = 0.5);

      /// @brief Enumeration listing the detailed class of interaction a dot can have with an atom.
      /// This is a weak enumeration rather than a class so that it can by type-cast back into integers
      /// and used as an index.  These index values are chosen to match those used by the orignal Probe.
      enum InteractionType {
        WideContact = 0, CloseContact = 1, WeakHydrogenBond = 2, SmallOverlap = 3,
        Bump = 4, BadBump = 5, StandardHydrogenBond = 6, Invalid = -1
      };

      /// @brief Determine the type of interaction that is represented by a CheckDotResult.
      /// @param [in] overlapType Value returned in the CheckDotResult.
      /// @param [in] gap Value returned in the CheckDotResult.
      /// @param [in] separateBadBumps Whether we should return bad bumps (true) or group
      ///             them in with bumps (false).
      /// @return InteractionType of the contact, Invalid if the OverlapType was Ignore or
      ///         a non-negative number indicating the index of the interaction type in the
      ///         original Probe code.
      InteractionType interaction_type(OverlapType overlapType, double gap, bool separateBadBumps = false) const;

      /// @brief String name for the interaction type, meaningful to the user.
      /// @param [in] t Type to report on.
      /// @return String name suitable for printing.
      static std::string interaction_type_name(InteractionType t);

      /// @brief Short string name for the interaction type, meaningful to the user.
      /// @param [in] t Type to report on.
      /// @return String name suitable for printing.
      static std::string interaction_type_short_name(InteractionType t);

      /// @brief Structure to hold the results from a call to score_dots()
      class ScoreDotsResult {
      public:
        ScoreDotsResult() : valid(false), bumpSubScore(0), hBondSubScore(0), attractSubScore(0), hasBadBump(false) {};
        bool    valid;              ///< False if this information has not yet been computed.
        double  bumpSubScore;       ///< Portion of the score due to bumps
        double  hBondSubScore;      ///< Portion of the score due to hydrogen bonds
        double  attractSubScore;    ///< Portion of the score due to non-bumping attraction
        bool    hasBadBump;         ///< Did this atom have a bad bump for any of its dots?
        /// @brief Sum of all of the sub-scores, the total score
        double  totalScore() const { return bumpSubScore + hBondSubScore + attractSubScore; }
      };

      /// @brief Determine the bump and hydrogen-bond subscores for a vector of dots on an atom
      /// @param [in] sourceAtom Atom that the dots are surrounding.
      /// @param [in] minOccupancy The minimum occupancy of the atom to be scored (applies to both
      ///             source and target; if either is below this, it will not be scored.
      ///             If the source is below this, the return value will be marked valid but will
      ///             have 0 in all of its interactions.
      /// @param [in] spatialQuery Structure to ask for neighbors of the atom.  This must contain
      ///             only atoms that are to be considered; those that are in the same conformation
      ///             or in all conformations.
      /// @param [in] nearbyRadius Maximum distance that an atom can be away and still be a neighbor.
      ///             This should NOT include consideration of the probe radius, which will be added
      ///             inside this function when it is required (depending on whether onlyBumps is
      ///             set).
      /// @param [in] probeRadius Radius of the probe rolled between the two potentially-contacting atoms
      ///             If this is < 0, an invalid result will be returned.
      /// @param [in] excluded Atoms that are to be excluded from contact, for example this could be a list
      ///             of atoms bonded to sourceAtom.
      /// @param [in] dots Vector of dots to compare.  Each is added to the sourceAtom origin.
      /// @param [in] density Density of the dots on the probe sphere, used to normalize results.
      ///             If this is <= 0, an invalid result will be returned.
      /// @param [in] onlyBumps If true, ignore near touches and count even hydrogen bonds as bumps.
      /// @param [in] preTrimmedDots If true, the dots have already been trimmed to only those that
      ///             are not inside excluded atoms. The excluded atoms are still used to determine
      ///             which neighboring atoms should be excluded, but they are not passed on to
      ///             the check_dot() function.
      /// @return Normalized sum of scores, also broken down by hydrogen bond vs. bump scores.
      ScoreDotsResult score_dots(iotbx::pdb::hierarchy::atom const &sourceAtom, double minOccupancy,
        SpatialQuery &spatialQuery, double nearbyRadius, double probeRadius,
        scitbx::af::shared<iotbx::pdb::hierarchy::atom> const &exclude,
        scitbx::af::shared<Point> const &dots, double density, bool onlyBumps,
        bool preTrimmedDots = false);

      /// @brief Count how many surface dots on an atom are not within an excluded atom.
      /// @param [in] sourceAtom Atom that the dot is offset with respect to.
      /// @param [in] dots Vector of dots to compare.  Each is added to the sourceAtom origin.
      /// @param [in] excluded Atoms that are to be excluded from contact, for example this could be a list
      ///             of atoms bonded to sourceAtom.  If the dot is inside an excluded atom, it will not be
      ///             counted.
      /// @return Number of surface dots that are not inside an excluded atom.
      unsigned count_surface_dots(iotbx::pdb::hierarchy::atom const &sourceAtom,
        scitbx::af::shared<Point> const& dots,
        scitbx::af::shared<iotbx::pdb::hierarchy::atom> const& exclude);

      //===========================================================================
      // Seldom-used methods below here.

      /// @brief See if the two atoms are in compatible conformations.
      /// @param [in] a1 First atom altloc to compare
      /// @param [in] a2 Second atom altloc to compare
      /// @return True if the atoms are in compatible conformations, false otherwise.
      static bool compatible_conformations(std::string const& a1, std::string const& a2);

      /// @brief Test method to verify that the class is behaving as intended.
      /// @return Empty string on success, string telling what went wrong on failure.
      static std::string test();

    protected:
      /// Parameters stored from constructor.
      ExtraAtomInfoMap m_extraInfoMap;
      double m_gapScale;
      double m_bumpWeight;
      double m_hBondWeight;
      double m_maxRegularHydrogenOverlap;
      double m_maxChargedHydrogenOverlap;
      double m_bumpOverlap;
      double m_badBumpOverlap;
      double m_contactCutoff;
      bool m_weakHBonds;
      bool m_ignoreIonInteractions;
    };

    /// @todo Figure out what all of the things needed by Probe (as opposed to Reduce) are.

    //=====================================================================================================

    /// @brief Test function to verify that all classes are behaving as intended.
    /// @return Empty string on success, string telling what went wrong on failure.
    std::string Scoring_test();

  }
}
