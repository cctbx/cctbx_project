// $Id$
/* Copyright (c) 2001 The Regents of the University of California through
   E.O. Lawrence Berkeley National Laboratory, subject to approval by the
   U.S. Department of Energy. See files COPYRIGHT.txt and
   cctbx/LICENSE.txt for further details.

   Revision history:
     2001 Sep 13: SpaceGroupType -> SpaceGroupInfo (R.W. Grosse-Kunstleve)
     2001 Jul 02: Merged from CVS branch sgtbx_special_pos (rwgk)
     Apr 2001: SourceForge release (R.W. Grosse-Kunstleve)
 */

#ifndef CCTBX_SGTBX_COORDINATES_H
#define CCTBX_SGTBX_COORDINATES_H

#include <vector>
#include <complex>
#include <cctbx/uctbx.h>
#include <cctbx/coordinates.h>
#include <cctbx/miller.h>
#include <cctbx/sgtbx/groups.h>

namespace sgtbx {

  namespace detail {

    template <class T>
    TrVec getUnitShifts(const boost::array<T, 3>& Delta)
    {
      TrVec result(1);
      for(std::size_t i=0;i<3;i++) {
        if (Delta[i] >= 0.) result[i] = static_cast<int>(Delta[i] + 0.5);
        else                result[i] = static_cast<int>(Delta[i] - 0.5);
      }
      return result;
    }

    template <class T>
    double modShortLength2(const uctbx::UnitCell& uc,
                           const boost::array<T, 3>& Diff) {
      return uc.Length2(fractional<T>(Diff).modShort());
    }

    void SetUniqueOps(const SpaceGroup& SgOps,
                      const RTMx& SpecialOp,
                      std::vector<RTMx>& UniqueOps);

  } // namespace detail

  //! Container for special position snap parameters.
  class SpecialPositionSnapParameters {
    public:
      //! Define special position snap parameters.
      /*! The definitions are used in a constructor of
          class SpecialPosition. The parameters are used in a
          numerically robust algorithm that first determines the site
          symmetry (point group) of a special position, and then moves
          ("snaps") the input coordinates to the exact location of the
          special position.
          <p>
          If the distance between symmetry mates is less than or equal
          to MinMateDistance, the site is moved to a special position.
          <p>
          If, after moving the site to a special position, the distance
          between symmetry mates is still less than or equal to
          MinMateDistance, an exception is raised if
          MustBeWellBehaved == true. This condition is usually
          the consequence of an input error: the unit cell is too
          small relative to MinMateDistance.
          <p>
          If MustBeWellBehaved == false, no exception is raised
          because of this condition. SpecialPosition::isWellBehaved()
          can be used to query the status.
          <p>
          For efficiency, the UnitCell object and the SpaceGroup object are
          only copied by reference.  This is, these objects must exist
          as long as the SpecialPositionSnapParameters are in use.
          <p>
          See also: SpecialPositionTolerances
       */
      SpecialPositionSnapParameters(const uctbx::UnitCell& uc,
                                    const SpaceGroup& SgOps,
                                    bool MustBeWellBehaved = true,
                                    double MinMateDistance = 0.5)
        : m_UnitCell(uc),
          m_SgOps(SgOps),
          m_MustBeWellBehaved(MustBeWellBehaved),
          m_MinMateDistance2(MinMateDistance * MinMateDistance) {}
    private:
      friend class SpecialPosition;
      friend class WyckoffTable;
      const uctbx::UnitCell& m_UnitCell;
      const SpaceGroup& m_SgOps;
      bool m_MustBeWellBehaved;
      double m_MinMateDistance2;
  };

  template <class T> class SymEquivCoordinates;

  //! Container for special position tolerance parameters.
  class SpecialPositionTolerances {
    public:
      //! Define special position tolerance parameters.
      /*! The definitions are used in a constructor of
          class SymEquivCoordinates. The parameters are used in an
          algorithm that uses simple distance calculations to determine
          if a site is on a special position. This algorithm is
          slightly faster then the alternative algorithm
          (see SpecialPositionSnapParameters). However, the
          input coordinates are not moved to the exact location
          of the special position. Therefore the value for
          Tolerance should in general be very small.
          <p>
          If the distance between symmetry mates is less than or equal
          to Tolerance, the site is considered to be on a special
          position.
          <p>
          If the distance between symmetry mates is less than
          MinimumDistance, but not less than or equal to
          Tolerance, an exception is raised.
          <p>
          The simple distance calculation is not guaranteed to be
          numerically stable.  As a safeguard it is asserted that the
          number of symmetry equivalent coordinates is a factor of the
          space group multiplicity.
          <p>
          To avoid numerical instabilities and potentially inaccurate
          results, MinimumDistance should be strictly greater than
          Tolerance. Only under this condition are the results
          guaranteed to be correct.
          <p>
          For efficiency, the UnitCell object and the SpaceGroup object are
          only copied by reference.  This is, these objects must exist
          as long as the SpecialPositionTolerances are in use.
       */
      SpecialPositionTolerances(const uctbx::UnitCell& uc,
                                const SpaceGroup& SgOps,
                                double MinimumDistance = 0.5,
                                double Tolerance = 0.01)
        : m_UnitCell(uc),
          m_SgOps(SgOps),
          m_MinimumDistance2(MinimumDistance * MinimumDistance),
          m_Tolerance2(Tolerance * Tolerance) {
        cctbx_assert(m_MinimumDistance2 >= m_Tolerance2);
      }
    private:
# if !(defined(BOOST_MSVC) && BOOST_MSVC <= 1200) // 1200 == VC++ 6.0
      template <class T> friend class SymEquivCoordinates;
# else
      friend class SymEquivCoordinates<float>;
      friend class SymEquivCoordinates<double>;
# endif
      const uctbx::UnitCell& m_UnitCell;
      const SpaceGroup& m_SgOps;
      double m_MinimumDistance2;
      double m_Tolerance2;
  };

  namespace detail { class RT_PointGroup; }

  //! Implementation of a robust special position algorithm.
  class SpecialPosition {
    public:
      //! Apply the robust special position algorithm to X.
      /*! The SpecialPositionSnapParameters are used in a numerically
          robust algorithm that first determines the point group of the
          site symmetry of X. If the site point group symmetry is
          higher than 1, the exact location of the special position
          (SnapPosition()) is determined.
          See SpecialPositionSnapParameters and
          SpecialPositionTolerances for details.
          <p>
          If auto_expand == true, expand() is called at the
          end of the constructor.
          <p>
          If determinePointGroupType == true, the point group type of
          the site symmetry is determined and can later be accessed via
          getPointGroupType().
          <p>
          See also: WyckoffTable, SymEquivCoordinates,
                    examples/python/generate_hklf.py
       */
      SpecialPosition(const SpecialPositionSnapParameters& params,
                      const fractional<double>& X,
                      bool auto_expand = false,
                      bool determinePointGroupType = false);
      //! Retrieve the original coordinates (X in the constructor).
      inline const fractional<double>& OriginalPosition() const {
        return m_OriginalPosition;
      }
      //! Exact location of the special position.
      inline const fractional<double>& SnapPosition() const {
        return m_SnapPosition;
      }
      //! Distance squared between OriginalPosition() and SnapPosition().
      inline double DistanceMoved2() const {
        return m_Parameters.m_UnitCell.Distance2(m_SnapPosition,
                                                 m_OriginalPosition);
      }
      //! Distance between OriginalPosition() and SnapPosition().
      inline double DistanceMoved() const {
        return std::sqrt(DistanceMoved2());
      }
      //! Shortest distance squared between symmetry mates of SnapPosition().
      inline double ShortestDistance2() const {
        return m_ShortestDistance2;
      }
      //! Shortest distance between symmetry mates of SnapPosition().
      inline double ShortestDistance() const {
        return std::sqrt(m_ShortestDistance2);
      }
      //! Test if ShortestDistance() > MinMateDistance.
      /*! See SpecialPositionSnapParameters for details.
       */
      inline bool isWellBehaved() const {
        return m_ShortestDistance2 > m_Parameters.m_MinMateDistance2;
      }
      //! Multiplicity (number of distinct symmetry mates) of SnapPosition().
      inline int M() const { return m_M; }
      //! Special position operation.
      /*! This operation is used to compute SnapPosition() from the
          input coordinates X and satisfies the following two
          conditions:
          <ul>
          <li>SpecialOp() * X = SnapPosition()
          <li>SpecialOp() * SnapPosition() = SnapPosition()
          </ul>
       */
      inline const RTMx& SpecialOp() const { return m_SpecialOp; }
      /*! \brief Access the site symmetry point group type that has been
          determined in the constructor.
       */
      /*! An exception is thrown if the point group type has not been
          determined in the constructor (determinePointGroupType == false).
       */
      tables::MatrixGroup::Code getPointGroupType() const;
      //! Expand the special position symmetry operation.
      /*! The SpecialOp() is multiplied with all symmetry operations.
          The unique results are stored in an internal list which
          can be accessed with the operator()().
          <p>
          expand() must be called before an instance of SpecialPosition
          can be used as an argument in the constructor of
          SymEquivCoordinates.
       */
      inline void expand() {
        detail::SetUniqueOps(m_Parameters.m_SgOps, m_SpecialOp, m_UniqueOps);
        cctbx_assert(m_UniqueOps.size() == m_M);
      }
      //! Test if expand() has been called.
      /*! See also: CheckExpanded()
       */
      inline bool isExpanded() const { return m_UniqueOps.size() != 0; }
      //! Assert that expand() has been called.
      /*! An exception is thrown if this assertion fails.<br>
          See also: isExpanded()
       */
      inline void CheckExpanded() const {
        if (!isExpanded()) {
          throw error(
          "Unique operations not initialized. Use expand() to initialize.");
        }
      }
      //! Access the i'th element of the list generated by expand().
      /*! i must be in the range [0,M()[. No range checking is
          performed for maximal performance.
       */
      inline const RTMx& operator[](std::size_t i) const {
        return m_UniqueOps[i];
      }
      //! Access the i'th element of the list generated by expand().
      /*! i must be in the range [0,M()[. An exception (cctbx::error_index)
          is thrown if i is out of range.
       */
      inline const RTMx& operator()(std::size_t i) const {
        CheckExpanded();
        if (i >= m_M) throw error_index();
        return m_UniqueOps[i];
      }
    private:
      friend class WyckoffTable;
      void BuildSpecialOp(detail::RT_PointGroup& SiteSymmetry);
      const SpecialPositionSnapParameters& m_Parameters;
      const fractional<double> m_OriginalPosition;
      fractional<double> m_SnapPosition;
      double m_ShortestDistance2;
      int m_M;
      RTMx m_SpecialOp;
      tables::MatrixGroup::Code m_PointGroupType;
      std::vector<RTMx> m_UniqueOps;
  };

  //! Information for Wyckoff positions.
  /*! See also: WyckoffTable, WyckoffMapping, SymEquivCoordinates
   */
  class WyckoffPosition {
    public:
      //! Constructor. For internal use only.
      WyckoffPosition(int M, const char Letter, const RTMx& SpecialOp)
        : m_M(M), m_Letter(Letter), m_SpecialOp(SpecialOp) {}
      //! The multiplicity of the Wyckoff position.
      inline int M() const { return m_M; }
      //! The Wyckoff letter according to the International Tables.
      /*! The Wyckoff letter "alpha" is represented by the character "@".
       */
      inline char Letter() const { return m_Letter; }
      //! A representative special position operation.
      inline const RTMx& SpecialOp() const { return m_SpecialOp; }
      //! Expand the representative special position operation.
      /*! The SpecialOp() is multiplied with all symmetry operations.
          The unique results are stored in an internal list which
          can be accessed with the operator()().
          <p>
          expand() must be called before an instance of WyckoffPosition
          can be used as an argument in the constructor of
          SymEquivCoordinates.
          <p>
          See also: WyckoffTable::expand()
       */
      inline void expand(const SpaceGroup& SgOps) {
        detail::SetUniqueOps(SgOps, m_SpecialOp, m_UniqueOps);
        cctbx_assert(m_UniqueOps.size() == m_M);
      }
      //! Test if expand() has been called.
      /*! See also: CheckExpanded()
       */
      inline bool isExpanded() const { return m_UniqueOps.size() != 0; }
      //! Assert that expand() has been called.
      /*! An exception is thrown if this assertion fails.<br>
          See also: isExpanded()
       */
      inline void CheckExpanded() const {
        if (!isExpanded()) {
          throw error(
          "Unique operations not initialized. Use expand() to initialize.");
        }
      }
      //! Access the i'th element of the list generated by expand().
      /*! i must be in the range [0,M()[. No range checking is
          performed for maximal performance.
       */
      inline const RTMx& operator[](std::size_t i) const {
        return m_UniqueOps[i];
      }
      //! Access the i'th element of the list generated by expand().
      /*! i must be in the range [0,M()[. An exception (cctbx::error_index)
          is thrown if i is out of range.
       */
      inline const RTMx& operator()(std::size_t i) const {
        CheckExpanded();
        if (i >= m_M) throw error_index();
        return m_UniqueOps[i];
      }
    private:
      int m_M;
      char m_Letter;
      RTMx m_SpecialOp;
      std::vector<RTMx> m_UniqueOps;
  };

  //! Pair of WyckoffPosition and a symmetry operation.
  /*! The WyckoffPositions are defined by a representative special
      position operation. These are tabulated according to the
      conventions of the International Tables for Crystallography,
      Volume A, 1983. A representative special position operation
      only applies to one particular point, axis, or plane.
      To determine the Wyckoff position of an aritrary point,
      a symmetry operation (Mapping()) is constructed that maps the
      input coordinates to the point, axis, or plane for which
      the representative operation is defined.
      <p>
      Both the rotation part and the translation part of Mapping() are
      unique for a given position. In particular, the modulus operation
      that can be applied to regular symmetry operations cannot be
      applied to Mapping().
      <p>
      See also: WyckoffTable, WyckoffPosition, SymEquivCoordinates
   */
  class WyckoffMapping {
    public:
      //! Constructor. For internal use only.
      WyckoffMapping() : m_WP(0), m_Mapping(RTMx(0, 0)) {}
      //! Constructor. For internal use only.
      WyckoffMapping(const WyckoffPosition& WP, const RTMx& Mapping)
        : m_WP(&WP), m_Mapping(Mapping) {}
      //! Reference to entry in WyckoffTable.
      /*! For efficiency, the entry is not copied. This is, the
          WyckoffTable must exist as long as this reference
          is used.
       */
      inline const WyckoffPosition& WP() const { return *m_WP; }
      //! Symmetry operation that maps to the representative Wyckoff position.
      /*! See class details and WyckoffTable::getWyckoffMapping().
       */
      inline const RTMx& Mapping() const { return m_Mapping; }
      //! Exact location of the representative Wyckoff position.
      /*! Formula used: WP().SpecialOp() * Mapping() * X
          <p>
          The result can be used directly with the representative
          Wyckoff position to compute SymEquivCoordinates. This
          is useful for repeated, efficient computation of
          symmetry mates.
       */
      inline fractional<double>
      snap_to_representative(const fractional<double>& X) const {
        return (*m_WP).SpecialOp() * (m_Mapping * X);
      }
      //! Exact location of the special position.
      /*! This is equivalent to SpecialPosition::SnapPosition().
          <p>
          Formula used: Mapping().inverse() * WP().SpecialOp() * Mapping() * X
       */
      inline fractional<double>
      snap(const fractional<double>& X) const {
        return    m_Mapping.inverse_with_cancel()
               * ((*m_WP).SpecialOp() * (m_Mapping * X));
      }
    private:
      const WyckoffPosition* m_WP;
      RTMx m_Mapping;
  };

  //! Table of Wyckoff positions.
  /*! This class represents the information for all Wyckoff positions
      for a given space group.
      <p>
      See also: WyckoffPosition, WyckoffMapping,
                SpecialPosition, SymEquivCoordinates
      <p>
      An example for the use of this class is given in the file
      examples/python/generate_hklf.py.
   */
  class WyckoffTable {
    public:
      //! Constructor.
      /*! The Wyckoff positions for the 230 reference settings are
          tabulated. SgInfo::CBOp() is used to transform the
          tabulated settings to the given setting.
          <p>
          If auto_expand == true, SgInfo::SgOps() is
          used to expand the representative special position
          symmetry operations to lists of unique operations.
          <p>
          See also: expand(), WyckoffPosition::expand()
       */
      WyckoffTable(const SpaceGroupInfo& SgInfo, bool auto_expand = false);
      //! Call expand() for all Wyckoff positions.
      /*! See also: WyckoffPosition::expand()
       */
      void expand(const SpaceGroup& SgOps);
      //! The number of Wyckoff positions for the given space group.
      /*! This number varies between 1 and 27.
       */
      inline std::size_t N() const { return m_Operations.size(); }
      //! Return a reference to the i'th Wyckoff position.
      /*! i must be in the range [0,N()[. No range checking is
          performed for maximal performance.
          <p>
          The general position has the index 0. The Wyckoff position
          with the letter "a" has the index N()-1.
       */
      inline const WyckoffPosition& operator[](std::size_t i) const {
        return m_Operations[i];
      }
      //! Return a reference to the i'th Wyckoff position.
      /*! i must be in the range [0,N()[. An exception (cctbx::error_index)
          is thrown if i is out of range.
          <p>
          The general position has the index 0. The Wyckoff position
          with the letter "a" has the index N()-1.
       */
      inline const WyckoffPosition& operator()(std::size_t i) const {
        if (i >= N()) throw error_index();
        return m_Operations[i];
      }
      //! Return a reference to the Wyckoff position with the given Letter.
      /*! An exception is thrown if the Letter is not valid.
          <p>
          The Wyckoff letter "alpha" must be given as the character "@".
          <p>
          Looking up a Wyckoff position by letter is the slowest form
          of access. For repeated access, use LookupIndex() to map
          the Wyckoff letter to an index.
          For a given space group type, there is a setting independent
          one-to-one correspondence between letters and indices.
       */
      inline const WyckoffPosition& operator()(char Letter) const {
        return m_Operations[LookupIndex(Letter)];
      }
      //! Look up the Wyckoff position index for a given Wyckoff Letter.
      /*! An exception is thrown if the Letter is not valid.
          <p>
          The Wyckoff letter "alpha" must be given as the character "@".
          <p>
          For a given space group type, there is a setting independent
          one-to-one correspondence between letters and indices.
       */
      std::size_t LookupIndex(char Letter) const;
      //! Determine the Wyckoff position using SpecialPosition::SpecialOp().
      /*! Use this overload for maximum reliability and flexibility.
          In general, determining SpecialPosition::SpecialOp() followed
          by using this overload is also faster than the alternative.
          Usage:<pre>
          uctbx::UnitCell uc = ...;
          sgtbx::SpaceGroup SgOps = ...;
          fractional<double> X = ...;
          SpecialPositionSnapParameters SnapParameters(uc, SgOps);
          SpecialPosition SP = SpecialPosition(SnapParameters, X);
          WyckoffMapping WM = getWyckoffMapping(SP);</pre>
       */
      const WyckoffMapping
      getWyckoffMapping(const SpecialPosition& SP) const;
      //! Determine the Wyckoff position of X.
      /*! This overload is included mainly for debugging purposes.
          It is recommended to use the alternative algorithm.
       */
      const WyckoffMapping
      getWyckoffMapping(const uctbx::UnitCell& uc,
                        const SpaceGroup& SgOps,
                        const fractional<double>& X,
                        double SnapRadius = 0.5) const;
    private:
      void InitializeOperations(const SpaceGroupInfo& SgInfo);
      std::vector<WyckoffPosition> m_Operations;
  };

  //! Container for symmetry equivalent (atomic) coordinates.
  template <class T>
  class SymEquivCoordinates {
    public:
      //! Compute symmetry equivalent coordinates using a robust algorithm.
      /*! See class SpecialPositionSnapParameters for details.
       */
      SymEquivCoordinates(const SpecialPositionSnapParameters& params,
                          const fractional<T>& X)
      {
        SpecialPosition SP(params, X, true);
        for(std::size_t i=0;i<SP.M();i++) {
          m_Coordinates.push_back(SP[i] * SP.SnapPosition());
        }
      }
      //! Compute symmetry equivalent coordinates using a robust algorithm.
      /*! See class SpecialPosition for details.
          <p>
          SpecialPosition::expand() has to be called for SP before SP
          can be used in this constructor.
       */
      SymEquivCoordinates(const SpecialPosition& SP)
      {
        SP.CheckExpanded();
        for(std::size_t i=0;i<SP.M();i++) {
          m_Coordinates.push_back(SP[i] * SP.SnapPosition());
        }
      }
      //! Compute symmetry equivalent coordinates using a WyckoffMapping.
      /*! See class WyckoffMapping for details.
          <p>
          WyckoffPosition::expand() has to be called for SP before SP
          can be used in this constructor. This is normally achieved
          by calling WyckoffTable::expand().
       */
      SymEquivCoordinates(const WyckoffMapping& WM,
                          const fractional<T>& X)
      {
        const WyckoffPosition& WP = WM.WP();
        WP.CheckExpanded();
        fractional<T> Xr = WM.snap_to_representative(X);
        for(std::size_t i=0;i<WP.M();i++) {
          m_Coordinates.push_back(WP[i] * Xr);
        }
      }
      //! Compute symmetry equivalent coordinates using a WyckoffPosition.
      /*! See class WyckoffMapping for details. If X is known to be
          close to the representative Wyckoff position, the
          WyckoffMapping::Mapping() is not needed. The reference
          to the tabulated Wyckoff position can be used directly.
          <p>
          This constructor should only be used with great care:
          If X is not close to the representative Wyckoff position, the
          symmetry equivalent coordinates will in general be incorrect.
          The representative Wyckoff position in turn depends on the
          change-of-basis matrix of SgType used in the constructor of
          WyckoffTable. This is, it is important that TidyCBOp = true
          when calling getSpaceGroupType().
       */
      SymEquivCoordinates(const WyckoffPosition& WP,
                          const fractional<T>& X)
      {
        WP.CheckExpanded();
        for(std::size_t i=0;i<WP.M();i++) {
          m_Coordinates.push_back(WP[i] * X);
        }
      }
      /*! \brief Compute symmetry equivalent coordinates using simple
          distance calculations.
       */
      /*! See class SpecialPositionTolerances for details.
          <p>
          To ensure numerical stability, this constructor should only
          be used with coordinates that have been preprocessed
          and are known to be well defined. In general it is best
          to use the robust algorithm. Unless there is a large number
          of special positions in a structure, the robust algorithm
          is only marginally slower.
       */
      SymEquivCoordinates(const SpecialPositionTolerances& params,
                          const fractional<T>& X)
      {
        T Tolerance2 = params.m_Tolerance2;
        T MinimumDistance2 = params.m_MinimumDistance2;
        m_Coordinates.push_back(X);
        for(int i=1;i<params.m_SgOps.OrderZ();i++) {
          fractional<T> SX = params.m_SgOps(i) * X;
          T Delta2 = getShortestDistance2(params.m_UnitCell, SX);
          if (Delta2 >= Tolerance2) {
            if (Delta2 < MinimumDistance2) {
              throw error(
              "Special position not well defined."
              " Use SpecialPositionSnapParameters.");
            }
            else {
              m_Coordinates.push_back(SX);
            }
          }
        }
        if (params.m_SgOps.OrderZ() % m_Coordinates.size() != 0) {
          throw error(
          "Numerical instability. Use SpecialPositionSnapParameters.");
        }
      }
      /*! \brief Compute symmetry equivalent coordinates without
          treatment of special positions.
       */
      /*! The symmetry operations are applied to X. Duplicates on
          special positions are not removed. The multiplicty M() will
          always be equal to SpaceGroup::OrderZ(). This algorithm is
          therefore very fast and is suitable for structures with most
          atoms on a general position (e.g. protein structures).
          The true multiplicity of X can be tabulated and used
          as a weight in structure factor calculations.
       */
      SymEquivCoordinates(const SpaceGroup& SgOps,
                          const fractional<T>& X)
      {
        m_Coordinates.push_back(X);
        for(int i=1;i<SgOps.OrderZ();i++) {
          fractional<T> SX = SgOps(i) * X;
          m_Coordinates.push_back(SX);
        }
      }
      //! Number of symmetry equivalent coordinates (multiplicity).
      inline int M() const { return m_Coordinates.size(); }
      //! Return the i'th symmetry equivalent coordinate.
      /*! i must be in the range [0,M()[. No range checking is
          performed for maximal performance.
       */
      inline const fractional<T>&
      operator[](std::size_t i) const {
        return m_Coordinates[i];
      }
      //! Return the i'th symmetry equivalent coordinate.
      /*! An exception is thrown if i is out of range.
       */
      inline const fractional<T>&
      operator()(std::size_t i) const {
        if (i >= M()) throw error_index();
        return m_Coordinates[i];
      }
      //! Shortest distance squared between the symmetry mates of X and Y.
      /*! Determine the shortest distance between Y and the symmetry
          mates in the internal table.
       */
      T getShortestDistance2(const uctbx::UnitCell& uc,
                             const fractional<T>& Y) const
      {
        T result = detail::modShortLength2(uc, Y - m_Coordinates[0]);
        for(std::size_t i=1;i<m_Coordinates.size();i++) {
          T Delta2 = detail::modShortLength2(uc, Y - m_Coordinates[i]);
          if (result > Delta2)
              result = Delta2;
        }
        return result;
      }
      //! Shortest distance between the symmetry mates of X and Y.
      /*! Determine the shortest distance between Y and the symmetry
          mates in the internal table.
       */
      inline T getShortestDistance(const uctbx::UnitCell& uc,
                                   const fractional<T>& Y) const {
        return std::sqrt(getShortestDistance2(uc, Y));
      }
      //! Compute Sum(exp(2 pi i H X)) for all symmetry equivalent X.
      /*! This sum is a sub-expression in the structure factor
          calculation. See file examples/python/generate_hklf.py.
       */
      std::complex<T> StructureFactor(const Miller::Index& H) const
      {
        using cctbx::constants::pi;
        std::complex<T> F(0., 0.);
        for(std::size_t i=0;i<M();i++) {
          T phase = 2. * pi * (H * m_Coordinates[i]);
          F += std::complex<T>(std::cos(phase), std::sin(phase));
        }
        return F;
      }

    private:
      std::vector<fractional<T> > m_Coordinates;
  };

} // namespace sgtbx

#endif // CCTBX_SGTBX_COORDINATES_H
