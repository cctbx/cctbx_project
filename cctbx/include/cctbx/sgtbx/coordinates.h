// $Id$
/* Copyright (c) 2001 The Regents of the University of California through
   E.O. Lawrence Berkeley National Laboratory, subject to approval by the
   U.S. Department of Energy. See files COPYRIGHT.txt and
   cctbx/LICENSE.txt for further details.

   Revision history:
     2001 Oct 12: SpecialPosition -> SiteSymmetry (R.W. Grosse-Kunstleve)
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

namespace cctbx { namespace sgtbx {

  namespace detail {

    void SetUniqueOps(const SpaceGroup& SgOps,
                      const RTMx& SpecialOp,
                      std::vector<RTMx>& UniqueOps);

  } // namespace detail

  //! Container for special position snap parameters.
  class SpecialPositionSnapParameters {
    public:
      // NO DEFAULT CONSTRUCTOR
      //! Define special position snap parameters.
      /*! The definitions are used in a constructor of
          class SiteSymmetry. The parameters are used in a
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
          because of this condition. SiteSymmetry::isWellBehaved()
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
#ifndef DOXYGEN_SHOULD_SKIP_THIS
      friend class SiteSymmetry;
      friend class WyckoffTable;
#endif // DOXYGEN_SHOULD_SKIP_THIS
      const uctbx::UnitCell& m_UnitCell;
      const SpaceGroup& m_SgOps;
      bool m_MustBeWellBehaved;
      double m_MinMateDistance2;
  };

  template <class FloatType> class SymEquivCoordinates;

  //! Container for special position tolerance parameters.
  class SpecialPositionTolerances {
    public:
      // NO DEFAULT CONSTRUCTOR
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
          number of symmetrically equivalent coordinates is a factor of the
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
#ifndef DOXYGEN_SHOULD_SKIP_THIS
# if (defined(BOOST_MSVC) && BOOST_MSVC <= 1200) || defined(__MWERKS__)
                                        // 1200 == VC++ 6.0
      friend class SymEquivCoordinates<float>;
      friend class SymEquivCoordinates<double>;
# else
      template <class FloatType> friend class SymEquivCoordinates;
# endif
#endif // DOXYGEN_SHOULD_SKIP_THIS
      const uctbx::UnitCell& m_UnitCell;
      const SpaceGroup& m_SgOps;
      double m_MinimumDistance2;
      double m_Tolerance2;
  };

  namespace detail {

    class RT_PointGroup {
      public:
        typedef std::vector<RTMx> vec_type;
        vec_type Matrices;
        bool invalid;
        RT_PointGroup() : invalid(false) {}
        void reset(const RTMx& M);
        void add(const RTMx& M);
        void expand(const RTMx& M);
        bool try_expand(const RTMx& M);
        RTMx accumulate() const;
    };

  } // namespace detail

  //! Robust algorithm for the determination of site-symmetries.
  class SiteSymmetry {
    public:
      //! Default constructor. Calling certain methods may cause crashes!
      SiteSymmetry() : m_Parameters(0) {}
      //! Determine the site symmetry of X.
      /*! The SpecialPositionSnapParameters are used in a numerically
          robust algorithm that first determines the point group of the
          site symmetry of X. If the site symmetry is higher than 1,
          the exact location of the special position (SnapPosition())
          is determined. See SpecialPositionSnapParameters and
          SpecialPositionTolerances for details.
          <p>
          If auto_expand == true, expand() is called at the
          end of the constructor.
          <p>
          See also: WyckoffTable, SymEquivCoordinates,
                    examples/python/generate_hklf.py
       */
      SiteSymmetry(const SpecialPositionSnapParameters& params,
                   const fractional<double>& X,
                   bool auto_expand = false);
      //! Retrieve the original coordinates (X in the constructor).
      const fractional<double>& OriginalPosition() const {
        return m_OriginalPosition;
      }
      //! Exact location of the special position.
      const fractional<double>& SnapPosition() const {
        return m_SnapPosition;
      }
      //! Distance squared between OriginalPosition() and SnapPosition().
      double DistanceMoved2() const {
        return m_Parameters->m_UnitCell.Distance2(m_SnapPosition,
                                                  m_OriginalPosition);
      }
      //! Distance between OriginalPosition() and SnapPosition().
      double DistanceMoved() const {
        return std::sqrt(DistanceMoved2());
      }
      //! Shortest distance squared between symmetry mates of SnapPosition().
      double ShortestDistance2() const {
        return m_ShortestDistance2;
      }
      //! Shortest distance between symmetry mates of SnapPosition().
      double ShortestDistance() const {
        return std::sqrt(m_ShortestDistance2);
      }
      //! Test if ShortestDistance() > MinMateDistance.
      /*! See SpecialPositionSnapParameters for details.
       */
      bool isWellBehaved() const {
        return m_ShortestDistance2 > m_Parameters->m_MinMateDistance2;
      }
      //! Multiplicity (number of distinct symmetry mates) of SnapPosition().
      int M() const { return m_M; }
      //! Special position operation.
      /*! This operation is used to compute SnapPosition() from the
          input coordinates X and satisfies the following two
          conditions:
          <ul>
          <li>SpecialOp() * X = SnapPosition()
          <li>SpecialOp() * SnapPosition() = SnapPosition()
          </ul>
       */
      const RTMx& SpecialOp() const { return m_SpecialOp; }
      //! Determine the site symmetry point group type.
      tables::MatrixGroup::Code PointGroupType() const;
      //! Apply special position operator to coordinates X.
      template <class FloatType>
      const fractional<FloatType>
      ApplySpecialOp(const fractional<FloatType>& X) const {
        return m_SpecialOp * X;
      }
      /*! \brief Test if given anisotropic displacement parameters
          Ustar are compatible with site symmetry.
       */
      /*! The expression
          <p>
          R Ustar R_transposed == Ustar
          <p>
          is evaluated for all rotation parts R of the site
          symmetry.
       */
      template <class FloatType>
      bool
      isCompatibleUstar(const af::tiny<FloatType, 6>& Ustar,
                        FloatType tolerance = 1.e-6) const
      {
        FloatType scaled_tolerance = 0.;
        for(std::size_t j=0;j<6;j++) {
          FloatType x = Ustar[j];
          if (x < 0.) x = -x;
          if (scaled_tolerance < x) scaled_tolerance = x;
        }
        scaled_tolerance *= tolerance;
        af::tiny<FloatType, 9>
        U = MatrixLite::CondensedSymMx33_as_FullSymMx33(
          Ustar, type_holder<FloatType>());
        for (std::size_t i=0;i<m_PointGroup.Matrices.size();i++) {
          af::tiny<FloatType, 9>
          R = m_PointGroup.Matrices[i].Rpart().as_array(FloatType());
          af::tiny<FloatType, 9>
          RURt = MatrixLite::FullTensorTransformation(R, U);
          af::tiny<FloatType, 6>
          Up = MatrixLite::FullSymMx33_as_CondensedSymMx33(
            RURt, type_holder<FloatType>());
          if (!(af::approx_equal_scaled(Ustar, Up, scaled_tolerance) == true)){
            return false;
          }
        }
        return true;
      }
      /*! \brief Check if given anisotropic displacement parameters
          Ustar are compatible with site symmetry.
       */
      /*! Similar to isCompatibleUstar(), but an exception
          is thrown if the Ustar tensor is incompatible with
          the site symmetry.
       */
      template <class FloatType>
      void
      CheckUstar(const af::tiny<FloatType, 6>& Ustar,
                 double tolerance = 1.e-6) const
      {
        if (!isCompatibleUstar(Ustar, tolerance)) {
          throw error(
            "Ustar tensor is incompatible with site symmetry.");
        }
      }
      /*! \brief Average symmetry copies of Ustar tensor to obtain a
          tensor that satisfies the symmetry constraints.
       */
      /*! The averaged tensor is equivalent to beta_inv
          of Giacovazzo, Fundamentals of Crystallography 1992,
          p. 189.
       */
      template <class FloatType>
      af::tiny<FloatType, 6>
      AverageUstar(const af::tiny<FloatType, 6>& Ustar) const
      {
        af::tiny<FloatType, 9>
        U = MatrixLite::CondensedSymMx33_as_FullSymMx33(
          Ustar, type_holder<FloatType>());
        af::tiny<FloatType, 9> SumRURt;
        SumRURt.fill(0.);
        for (std::size_t i=0;i<m_PointGroup.Matrices.size();i++) {
          af::tiny<FloatType, 9>
          R = m_PointGroup.Matrices[i].Rpart().as_array(FloatType());
          af::tiny<FloatType, 9>
          RURt = MatrixLite::FullTensorTransformation(R, U);
          SumRURt = SumRURt + RURt;
        }
        return MatrixLite::FullSymMx33_as_CondensedSymMx33(
          SumRURt, type_holder<FloatType>())
            / FloatType(m_PointGroup.Matrices.size());
      }
      //! Expand the special position symmetry operation.
      /*! The SpecialOp() is multiplied with all symmetry operations.
          The unique results are stored in an internal list which
          can be accessed with the operator()().
          <p>
          expand() must be called before an instance of SiteSymmetry
          can be used as an argument in the constructor of
          SymEquivCoordinates.
       */
      void expand() {
        detail::SetUniqueOps(m_Parameters->m_SgOps, m_SpecialOp, m_UniqueOps);
        cctbx_assert(m_UniqueOps.size() == m_M);
      }
      //! Test if expand() has been called.
      /*! See also: CheckExpanded()
       */
      bool isExpanded() const { return m_UniqueOps.size() != 0; }
      //! Assert that expand() has been called.
      /*! An exception is thrown if this assertion fails.
          <p>
          See also: isExpanded()
       */
      void CheckExpanded() const {
        if (!isExpanded()) {
          throw error(
          "Unique operations not initialized. Use expand() to initialize.");
        }
      }
      //! Access the i'th element of the list generated by expand().
      /*! i must be in the range [0,M()[. No range checking is
          performed for maximal performance.
       */
      const RTMx& operator[](std::size_t i) const {
        return m_UniqueOps[i];
      }
      //! Access the i'th element of the list generated by expand().
      /*! i must be in the range [0,M()[. An exception (cctbx::error_index)
          is thrown if i is out of range.
       */
      const RTMx& operator()(std::size_t i) const {
        CheckExpanded();
        if (i >= m_M) throw error_index();
        return m_UniqueOps[i];
      }
    private:
#ifndef DOXYGEN_SHOULD_SKIP_THIS
      friend class WyckoffTable;
#endif // DOXYGEN_SHOULD_SKIP_THIS
      void BuildSpecialOp();
      const SpecialPositionSnapParameters* m_Parameters;
      const fractional<double> m_OriginalPosition;
      detail::RT_PointGroup m_PointGroup;
      fractional<double> m_SnapPosition;
      double m_ShortestDistance2;
      int m_M;
      RTMx m_SpecialOp;
      std::vector<RTMx> m_UniqueOps;
  };

  //! Information for Wyckoff positions.
  /*! See also: WyckoffTable, WyckoffMapping, SymEquivCoordinates
   */
  class WyckoffPosition {
    public:
      //! Default constructor. Some data members are not initialized!
      WyckoffPosition() {}
      //! Constructor. For internal use only.
      WyckoffPosition(int M, const char Letter, const RTMx& SpecialOp)
        : m_M(M), m_Letter(Letter), m_SpecialOp(SpecialOp) {}
      //! The multiplicity of the Wyckoff position.
      int M() const { return m_M; }
      //! The Wyckoff letter according to the International Tables.
      /*! The Wyckoff letter "alpha" is represented by the character "@".
       */
      char Letter() const { return m_Letter; }
      //! A representative special position operation.
      const RTMx& SpecialOp() const { return m_SpecialOp; }
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
      void expand(const SpaceGroup& SgOps) {
        detail::SetUniqueOps(SgOps, m_SpecialOp, m_UniqueOps);
        cctbx_assert(m_UniqueOps.size() == m_M);
      }
      //! Test if expand() has been called.
      /*! See also: CheckExpanded()
       */
      bool isExpanded() const { return m_UniqueOps.size() != 0; }
      //! Assert that expand() has been called.
      /*! An exception is thrown if this assertion fails.<br>
          See also: isExpanded()
       */
      void CheckExpanded() const {
        if (!isExpanded()) {
          throw error(
          "Unique operations not initialized. Use expand() to initialize.");
        }
      }
      //! Access the i'th element of the list generated by expand().
      /*! i must be in the range [0,M()[. No range checking is
          performed for maximal performance.
       */
      const RTMx& operator[](std::size_t i) const {
        return m_UniqueOps[i];
      }
      //! Access the i'th element of the list generated by expand().
      /*! i must be in the range [0,M()[. An exception (cctbx::error_index)
          is thrown if i is out of range.
       */
      const RTMx& operator()(std::size_t i) const {
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
      //! Default constructor. Calling certain methods may cause crashes!
      WyckoffMapping() : m_WP(0), m_Mapping(RTMx(0, 0)) {}
      //! Constructor. For internal use only.
      WyckoffMapping(const WyckoffPosition& WP, const RTMx& Mapping)
        : m_WP(&WP), m_Mapping(Mapping) {}
      //! Reference to entry in WyckoffTable.
      /*! For efficiency, the entry is not copied. This is, the
          WyckoffTable must exist as long as this reference
          is used.
       */
      const WyckoffPosition& WP() const { return *m_WP; }
      //! Symmetry operation that maps to the representative Wyckoff position.
      /*! See class details and WyckoffTable::getWyckoffMapping().
       */
      const RTMx& Mapping() const { return m_Mapping; }
      //! Exact location of the representative Wyckoff position.
      /*! Formula used: WP().SpecialOp() * Mapping() * X
          <p>
          The result can be used directly with the representative
          Wyckoff position to compute SymEquivCoordinates. This
          is useful for repeated, efficient computation of
          symmetry mates.
       */
      fractional<double>
      snap_to_representative(const fractional<double>& X) const {
        return m_WP->SpecialOp() * (m_Mapping * X);
      }
      //! Exact location of the special position.
      /*! This is equivalent to SiteSymmetry::SnapPosition().
          <p>
          Formula used: Mapping().inverse() * WP().SpecialOp() * Mapping() * X
       */
      fractional<double>
      snap(const fractional<double>& X) const {
        return    m_Mapping.inverse_with_cancel()
               * (m_WP->SpecialOp() * (m_Mapping * X));
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
                SiteSymmetry, SymEquivCoordinates
      <p>
      An example for the use of this class is given in the file
      examples/python/generate_hklf.py.
   */
  class WyckoffTable {
    public:
      //! Default constructor. Some data members are not initialized!
      WyckoffTable() {}
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
      explicit
      WyckoffTable(const SpaceGroupInfo& SgInfo, bool auto_expand = false);
      //! Call expand() for all Wyckoff positions.
      /*! See also: WyckoffPosition::expand()
       */
      void expand(const SpaceGroup& SgOps);
      //! The number of Wyckoff positions for the given space group.
      /*! This number varies between 1 and 27.
       */
      std::size_t N() const { return m_Operations.size(); }
      //! Return a reference to the i'th Wyckoff position.
      /*! i must be in the range [0,N()[. No range checking is
          performed for maximal performance.
          <p>
          The general position has the index 0. The Wyckoff position
          with the letter "a" has the index N()-1.
       */
      const WyckoffPosition& operator[](std::size_t i) const {
        return m_Operations[i];
      }
      //! Return a reference to the i'th Wyckoff position.
      /*! i must be in the range [0,N()[. An exception (cctbx::error_index)
          is thrown if i is out of range.
          <p>
          The general position has the index 0. The Wyckoff position
          with the letter "a" has the index N()-1.
       */
      const WyckoffPosition& operator()(std::size_t i) const {
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
      const WyckoffPosition& operator()(char Letter) const {
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
      //! Determine the Wyckoff position using SiteSymmetry::SpecialOp().
      /*! Use this overload for maximum reliability and flexibility.
          In general, determining SiteSymmetry::SpecialOp() followed
          by using this overload is also faster than the alternative.
          Usage:<pre>
          uctbx::UnitCell uc = ...;
          sgtbx::SpaceGroup SgOps = ...;
          fractional<double> X = ...;
          SpecialPositionSnapParameters SnapParameters(uc, SgOps);
          SiteSymmetry SS = SiteSymmetry(SnapParameters, X);
          WyckoffMapping WM = getWyckoffMapping(SS);</pre>
       */
      const WyckoffMapping
      getWyckoffMapping(const SiteSymmetry& SS) const;
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

  //! Container for symmetrically equivalent (atomic) coordinates.
  template <class FloatType>
  class SymEquivCoordinates {
    public:
      //! Default constructor. Some data members are not initialized!
      SymEquivCoordinates() {}
      /*! \brief Compute symmetrically equivalent coordinates
          using a robust algorithm.
       */
      /*! See class SpecialPositionSnapParameters for details.
       */
      SymEquivCoordinates(const SpecialPositionSnapParameters& params,
                          const fractional<FloatType>& X)
      {
        SiteSymmetry SS(params, X, true);
        for(std::size_t i=0;i<SS.M();i++) {
          m_Coordinates.push_back(SS[i] * SS.SnapPosition());
        }
      }
      /*! \brief Compute symmetrically equivalent coordinates
          using a robust algorithm.
       */
      /*! See class SiteSymmetry for details.
          <p>
          SiteSymmetry::expand() has to be called for SS before SS
          can be used in this constructor.
       */
      explicit
      SymEquivCoordinates(const SiteSymmetry& SS)
      {
        SS.CheckExpanded();
        for(std::size_t i=0;i<SS.M();i++) {
          m_Coordinates.push_back(SS[i] * SS.SnapPosition());
        }
      }
      //! Compute symmetrically equivalent coordinates using a WyckoffMapping.
      /*! See class WyckoffMapping for details.
          <p>
          WyckoffPosition::expand() has to be called for WM.WP() before WM
          can be used in this constructor. This is normally achieved
          by calling WyckoffTable::expand().
       */
      SymEquivCoordinates(const WyckoffMapping& WM,
                          const fractional<FloatType>& X)
      {
        const WyckoffPosition& WP = WM.WP();
        WP.CheckExpanded();
        fractional<FloatType> Xr = WM.snap_to_representative(X);
        for(std::size_t i=0;i<WP.M();i++) {
          m_Coordinates.push_back(WP[i] * Xr);
        }
      }
      //! Compute symmetrically equivalent coordinates using a WyckoffPosition.
      /*! See class WyckoffMapping for details. If X is known to be
          close to the representative Wyckoff position, the
          WyckoffMapping::Mapping() is not needed. The reference
          to the tabulated Wyckoff position can be used directly.
          <p>
          This constructor should only be used with great care:
          If X is not close to the representative Wyckoff position, the
          symmetrically equivalent coordinates will in general be incorrect.
          The representative Wyckoff position in turn depends on the
          change-of-basis matrix of SgType used in the constructor of
          WyckoffTable. This is, it is important that TidyCBOp = true
          when calling getSpaceGroupType().
       */
      SymEquivCoordinates(const WyckoffPosition& WP,
                          const fractional<FloatType>& X)
      {
        WP.CheckExpanded();
        for(std::size_t i=0;i<WP.M();i++) {
          m_Coordinates.push_back(WP[i] * X);
        }
      }
      /*! \brief Compute symmetrically equivalent coordinates using simple
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
                          const fractional<FloatType>& X)
      {
        FloatType Tolerance2 = params.m_Tolerance2;
        FloatType MinimumDistance2 = params.m_MinimumDistance2;
        m_Coordinates.push_back(X);
        for(int i=1;i<params.m_SgOps.OrderZ();i++) {
          fractional<FloatType> SX = params.m_SgOps(i) * X;
          FloatType Delta2 = getShortestDistance2(params.m_UnitCell, SX);
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
      /*! \brief Compute symmetrically equivalent coordinates without
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
                          const fractional<FloatType>& X)
      {
        m_Coordinates.push_back(X);
        for(int i=1;i<SgOps.OrderZ();i++) {
          fractional<FloatType> SX = SgOps(i) * X;
          m_Coordinates.push_back(SX);
        }
      }
      //! Number of symmetrically equivalent coordinates (multiplicity).
      int M() const { return m_Coordinates.size(); }
      //! Return the coordinates of the i'th symmetrically equivalent position.
      /*! i must be in the range [0,M()[. No range checking is
          performed for maximal performance.
       */
      const fractional<FloatType>&
      operator[](std::size_t i) const {
        return m_Coordinates[i];
      }
      //! Return the coordinates of the i'th symmetrically equivalent position.
      /*! An exception is thrown if i is out of range.
       */
      const fractional<FloatType>&
      operator()(std::size_t i) const {
        if (i >= M()) throw error_index();
        return m_Coordinates[i];
      }
      //! Shortest difference vector between the symmetry mates of X and Y.
      /*! Determine the shortest difference vector between Y and the symmetry
          mates in the internal table. The result is the shortest difference
          vector Y - symmetry mate of X under application of unit shifts.
       */
      fractional<FloatType>
      getShortestDifference(const uctbx::UnitCell& uc,
                            const fractional<FloatType>& Y) const
      {
        FloatType shortest_dist2(0);
        fractional<FloatType> shortest_diff;
        for(std::size_t i=0;i<m_Coordinates.size();i++) {
          fractional<FloatType> diff = Y - m_Coordinates[i];
          diff = diff.modShort();
          FloatType dist2 = uc.Length2(diff);
          if (shortest_dist2 > dist2 || i == 0) {
              shortest_dist2 = dist2;
              shortest_diff = diff;
          }
        }
        return shortest_diff;
      }
      //! Shortest distance squared between the symmetry mates of X and Y.
      /*! Determine the shortest distance between Y and the symmetry
          mates in the internal table.
       */
      FloatType getShortestDistance2(const uctbx::UnitCell& uc,
                                     const fractional<FloatType>& Y) const
      {
        FloatType result = uc.modShortLength2(Y - m_Coordinates[0]);
        for(std::size_t i=1;i<m_Coordinates.size();i++) {
          FloatType Delta2 = uc.modShortLength2(Y - m_Coordinates[i]);
          if (result > Delta2)
              result = Delta2;
        }
        return result;
      }
      //! Shortest distance between the symmetry mates of X and Y.
      /*! Determine the shortest distance between Y and the symmetry
          mates in the internal table.
       */
      FloatType getShortestDistance(const uctbx::UnitCell& uc,
                                    const fractional<FloatType>& Y) const {
        return std::sqrt(getShortestDistance2(uc, Y));
      }
      //! Compute Sum(exp(2 pi i H X)) for all symmetrically equivalent X.
      /*! This sum is a sub-expression in the structure factor
          calculation.
       */
      std::complex<FloatType> StructureFactor(const Miller::Index& H) const
      {
        using cctbx::constants::pi;
        std::complex<FloatType> F(0., 0.);
        for(std::size_t i=0;i<M();i++) {
          FloatType phase = 2. * pi * (H * m_Coordinates[i]);
          F += std::complex<FloatType>(std::cos(phase), std::sin(phase));
        }
        return F;
      }

    private:
      std::vector<fractional<FloatType> > m_Coordinates;
  };

}} // namespace cctbx::sgtbx

#endif // CCTBX_SGTBX_COORDINATES_H
