// $Id$
/* Copyright (c) 2001 The Regents of the University of California through
   E.O. Lawrence Berkeley National Laboratory, subject to approval by the
   U.S. Department of Energy. See files COPYRIGHT.txt and
   cctbx/LICENSE.txt for further details.

   Revision history:
     2001 Sep 13: SpaceGroupType -> SpaceGroupInfo (R.W. Grosse-Kunstleve)
     2001 Jul 02: Merged from CVS branch sgtbx_special_pos (rwgk)
     2001 May 31: merged from CVS branch sgtbx_type (R.W. Grosse-Kunstleve)
     Apr 2001: SourceForge release (R.W. Grosse-Kunstleve)
 */

#ifndef CCTBX_SGTBX_GROUPS_H
#define CCTBX_SGTBX_GROUPS_H

#include <iostream>
#include <vector>
#include <map>
#include <complex>
#include <cctbx/coordinates.h>
#include <cctbx/sgtbx/symbols.h>
#include <cctbx/sgtbx/matrix.h>
#include <cctbx/sgtbx/change_basis.h>
#include <cctbx/sgtbx/lattice_tr.h>
#include <cctbx/sgtbx/utils.h>
#include <cctbx/sgtbx/miller.h>
#include <cctbx/sgtbx/tables.h>
#include <cctbx/uctbx.h>
#include <cctbx/adptbx.h>

namespace cctbx { namespace sgtbx {

  class TrOps {
    friend class SpaceGroup;
    private:
      std::vector<TrVec> m_Vects;
      bool add(const TrVec& NewTr);
    public:
      explicit
      TrOps(int BF = STBF) : m_Vects() {
        m_Vects.push_back(TrVec(BF));
      }
      void reset(int BF = STBF) {
        m_Vects.clear();
        m_Vects.push_back(TrVec());
      }
      const std::vector<TrVec>& Vects() const { return m_Vects; }
      std::vector<TrVec>& Vects() { return m_Vects; }
      int nVects() const { return m_Vects.size(); }
      bool expand(const TrVec& NewTr);

      TrOps ChangeBasis(const ChOfBasisOp& CBOp) const;

      TrVec& operator[](int i) { return m_Vects[i]; }
      const TrVec& operator[](int i) const { return m_Vects[i]; }

      char getConventionalCentringTypeSymbol() const;
      ChOfBasisOp getConventionalZ2POp(int RBF = CRBF,
                                       int TBF = CTBF) const;
      TrVec TidyT(const TrVec& T) const;
  };

  class SpaceGroupInfo; // forward declaration

  //! Space Group Operations.
  /*! This is the central class of the sgtbx. It is composed of the
      following main components:<br>
      <ol>
      <li>A list of nLTr() centring vectors.
      <li>A variable fInv() and the translation part InvT() of an inversion
          operation, if present. fInv() == 1 if no inversion operation
          exists, fInv() == 2 otherwise.
      <li>A list of nSMx() = OrderZ()/(nLTr()*fInv()) representative
          symmetry matrices.
      </ol>
   */
  class SpaceGroup {
    private:
      bool  m_NoExpand;
      int   nLSL;
      int   nSSL;
      int   m_fInv;
      int   m_nSMx;
      TrOps m_LTr;
      TrVec m_InvT;
      bool  m_isTidy;
      boost::array<RTMx, nMaxReprRotMx> m_SMx;
      void addInv(const TrVec& NewInvT);
      void addSMx(const RTMx& NewSMx);
    public:
      //! Reset the symmetry to P1.
      void reset();
      //! Default constructor. Symmetry is set to P1.
      /*! With NoExpand == true, group multiplication will not be
          carried out. This option is for internal use only.
       */
      explicit
      SpaceGroup(bool NoExpand = false) : m_NoExpand(NoExpand) {
        reset();
      }
      //! Initialize with symmetry encoded by a Hall symbol.
      /*! The Hall symbol parser can be instructed to be more
          or less Pedantic.<br>
          With NoCType == true ("no centring type symbol"), the parser
          can also be used to interpret matrix symbols only. This
          option is for internal use only.<br>
          If an error occurs, an exception is thrown. parse_string can
          be investigated to locate the input character that triggered
          the error.
       */
      explicit
      SpaceGroup(parse_string& HSym,
            bool Pedantic = false, bool NoCType = false,
            bool NoExpand = false);
      //! Initialize with symmetry encoded by a Hall symbol.
      /*! Identical to the constructor that takes a parse_string
          as the first argument. However, if an exception is thrown
          there is no way to locate the input character that triggered
          the error.
       */
      explicit
      SpaceGroup(const std::string& HSym,
            bool Pedantic = false, bool NoCType = false,
            bool NoExpand = false);
      /*! Identical to the constructor that takes a parse_string
          as the first argument. However, if an exception is thrown
          there is no way to locate the input character that triggered
          the error.
       */
      explicit
      SpaceGroup(const char* HSym,
            bool Pedantic = false, bool NoCType = false,
            bool NoExpand = false);
      //! Initizalize with the Hall symbol in the SpaceGroupSymbols object.
      explicit
      SpaceGroup(const SpaceGroupSymbols& SgSymbols);
      //! Add a lattice translation vector to the space group.
      /*! Group multiplication is automatically performed
          (unless NoExpand == true).
       */
      void expandLTr(const TrVec& NewLTr);
      //! Add a centre of inversion to the space group.
      /*! The matrix added is (-I, NewInvT). This is, the rotation part
          -I is implied, and the translation part is NewInvT.<br>
          Group multiplication is automatically performed
          (unless NoExpand == true).
       */
      void expandInv(const TrVec& NewInvT);
      //! Add a Seitz matrix to the space group.
      /*! Group multiplication is automatically performed
          (unless NoExpand == true).
       */
      void expandSMx(const RTMx& NewSMx);
      //! Add a vector of Seitz matrices to the space group.
      /*! See also: expandSMx()
       */
      template <typename RTMxVectorType>
      void expandSMxVector(const RTMxVectorType& SMxVector) {
        for(std::size_t i=0;i<SMxVector.size();i++) {
          expandSMx(SMxVector[i]);
        }
      }

      //! Add lattice translation vectors to the space group.
      /*! The lattice translation vectors corresponding to the
          conventional centring type symbol are determined from
          a lookup table and added to the space group.<br>
          Group multiplication is automatically performed
          (unless NoExpand == true).
       */
      int expandConventionalCentringType(char Symbol);

      //! Parse a Hall Symbol and add the encoded symmetry matrices.
      /*! Similar to the constructors that take a Hall symbol as the
          first argument. However, the new symmetry matrices are added
          to the existing ones.
       */
      int ParseHallSymbol(parse_string& HSym,
                          bool Pedantic = false, bool NoCType = false);
      //! Parse a Hall Symbol and add the encoded symmetry matrices.
      /*! The change-of-basis operator is parsed but not applied.<br>
          For internal use only.
       */
      int ParseHallSymbolCBOp(parse_string& HSym, ChOfBasisOp& CBOp,
                              bool Pedantic = false, bool NoCType = false);

      //! Apply a change-of-basis operator to the symmetry operations.
      /*! The transformed space group is returned as a new object.<br>
          An exception is thrown if the change-of-basis operator is
          invalid or if the new symmetry matrices can not be
          represented as integer matrices with the base factors used.
       */
      SpaceGroup ChangeBasis(const ChOfBasisOp& CBOp) const;

      //! Rotation base factor of Seitz Matrices.
      int RBF() const { return m_SMx[0].Rpart().BF(); }
      //! Translation base factor of Seitz Matrices.
      int TBF() const { return m_SMx[0].Tpart().BF(); }

      //! Number of lattice translations.
      int nLTr() const { return m_LTr.nVects(); }
      //! Flag for centre of inversion.
      /*! fInv() == 1 if no inversion operation exists,<br>
          fInv() == 2 otherwise.
       */
      int fInv() const { return m_fInv; }
      //! Translation part for centre of inversion.
      TrVec InvT(bool tidy = false) const {
        if (tidy == false) return m_InvT;
        if (!m_InvT.isValid()) return m_InvT;
        return m_LTr.TidyT(m_InvT);
      }
      //! Number of representative Seitz matrices.
      int nSMx() const { return m_nSMx; }
      //! Order of the point-group = fInv() * nSMx().
      int OrderP() const { return m_fInv * m_nSMx; }
      //! Order of the space-group = nLTr() * fInv() * nSMx().
      int OrderZ() const { return m_LTr.nVects() * m_fInv * m_nSMx; }

      //! Access to the list of lattice translation vectors.
      const TrVec& LTr(int i) const { return m_LTr[i]; }
      //! Access to the list of representative Seitz matrices.
      RTMx& operator[](int i) { return m_SMx[i]; }
      //! Access to the list of representative Seitz matrices.
      const RTMx& operator[](int i) const { return m_SMx[i]; }
      //! Return a symmetry operation.
      /*! Usage:<pre>
          SpaceGroup s(...);
          for (int iLTr = 0; iLTr < s.nLTr(); iLTr++)
            for (int iInv = 0; iInv < s.fInv(); iInv++)
              for (int iSMx = 0; iSMx < s.nSMx(); iSMx++)
                RTMx M = s(iLTr, iInv, iSMx);</pre>
          The symmetry operations are generated from the representative
          list of Seitz matrices. The lattice translation with the index
          iLTr is added to the translation part. If iInv == 1, the
          representative Seitz matrix is pre-multiplied by (-I, InvT).<br>
          Note that the translation part of the returned Seitz matrix
          is not modified.
       */
      RTMx operator()(int iLTr, int iInv, int iSMx) const;
      //! Return a symmetry operation.
      /*! Usage:<pre>
          SpaceGroup s(...);
          for (int iLIS = 0; iLIS < s.OrderZ(); iLIS++)
            RTMx M = s(iLIS);</pre>
          The symmetry operations are generated from the representative
          list of Seitz matrices. Internally iLIS is interpreted as:<br>
          iLIS = ((iLTr * fInv()) + iInv) * nSMx() + iSMx<br>
          The comments for operator()(int iLTr, int iInv, int iSMx) apply.
       */
      RTMx operator()(int iLIS) const;

      //! Tidy the lists of representative symmetry operations in place.
      /*! The list of lattice translations is sorted in a certain order.
          If there is a centre of inversion, a certain normalized
          translation part is selected for InvT(), and a proper
          rotation part is selected for all Seitz matrices in the
          representative list. The list of representative Seitz
          matrices is sorted, and the translation parts are
          normalized.
          <p>
          After application of makeTidy(), a given space group
          representation will always result in exactly the same
          internal representation.
       */
      void makeTidy();
      //! Test for equality.
      /*! Internally, makeTidy() is used, followed by essentially a
          byte-wise comparison of the objects.<br>
          Each SpaceGroup object maintains an internal flag
          indicating whether or not makeTidy() was applied already.
          If an SpaceGroup object is repeatedly used in a test for equality,
          the test will therefore be significantly faster if makeTidy()
          is applied outside the loop.
       */
      friend bool operator==(const SpaceGroup& lhs, const SpaceGroup& rhs);
      //! Negation of test for equality.
      friend bool operator!=(const SpaceGroup& lhs,
                                    const SpaceGroup& rhs) {
        return !(lhs == rhs);
      }

      //! Determine operator for centred->primitive basis transformation.
      /*! Determine a change-of-basis operator that transforms the
          symmetry from a centered to a primitive setting.<br>
          For the conventional lattice centring types (P, A, B, C, I, R, H, F)
          a standard change-of-basis operator is determined by a lookup in a
          table. For unconventional lattice centring types, a change-of-basis
          operator is constructed with the algorithm published in
          <A HREF="http://journals.iucr.org/a/issues/1999/02/02/au0146/"
          ><I>Acta Cryst.</I> 1999, <B>A55</B>:383-395</A>.<br>
          The change-of-basis operator can be used as the argument to
          the ChangeBasis() member function:<pre>
          SpaceGroup centred(...);
          SpaceGroup primitive = centred.ChangeBasis(centred.getZ2POp());</pre>
       */
      ChOfBasisOp getZ2POp(int RBF = CRBF, int TBF = CTBF) const;
      //! Construct operator for centred->primitive basis transformation.
      /*! Similar to getZ2POp(), but does not use standard change-of-basis
          operators for conventional lattice centring types.<br>
          For internal use only.
       */
      ChOfBasisOp ConstructZ2POp(int RBF = CRBF, int TBF = CTBF) const;
      //! Determine the conventional centring type symbol (P,A,B,C,I,R,H,F).
      /*! If the lattice translation vectors do not correspond to
          any of the conventional centring types, the null character
          is returned.
       */
      char getConventionalCentringTypeSymbol() const {
        return m_LTr.getConventionalCentringTypeSymbol();
      }

      //! Test for chirality.
      /*! A space group is chiral if all its symmetry operations
          have a positive rotation-part type (1, 2, 3, 4, 6).<br>
          If there are symmetry operations with negative rotation-part
          types (-1, -2=m, -3, -4, -6) the space group is not chiral.<br>
          There are exactly 65 chiral space groups.<br>
          Note that proteins always crystallize in a chiral space
          group.
       */
      bool isChiral() const;
      //! Test for a centre of inversion.
      bool isCentric() const { return m_fInv == 2; }
      //! Test if centre of inversion is at the origin.
      bool isOriginCentric() const {
        return isCentric() && InvT(true).isNull();
      }

      //! Test if reflection with given Miller index is systematically absent.
      bool isSysAbsent(const Miller::Index& H) const;
      //! Test if reflection with given Miller index is centric.
      /*! A reflection with the Miller index H is "centric" if
          there is a symmetry operation with rotation part R such
          that H*R = -H.<br>
          The phase of a centric reflection is restricted to two phase
          angels (modulo pi).<br>
          If the phase restriction is also needed later in the
          calculation, do not use isCentric(), but use
          getPhaseRestriction() instead.
       */
      bool isCentric(const Miller::Index& H) const;
      //! Determine the phase restriction for the given Miller index.
      /*! See class PhaseRestriction.
       */
      PhaseRestriction getPhaseRestriction(const Miller::Index& H) const;
      //! See class PhaseRestriction::isValidPhase_rad().
      bool isValidPhase_rad(const Miller::Index& H, double phi,
                                   double tolerance = 1.e-5) const {
        return getPhaseRestriction(H).isValidPhase_rad(phi, tolerance);
      }
      //! See class PhaseRestriction::isValidPhase_deg().
      bool isValidPhase_deg(const Miller::Index& H, double phi,
                                   double tolerance = 1.e-5) const {
        return getPhaseRestriction(H).isValidPhase_deg(phi, tolerance);
      }
      //! Determine the reflection multiplicity for the given Miller index.
      /*! The multiplicity is defined as the number of symmetrically
          equivalent but distinct reflections.<br>
          If FriedelFlag == true, a centre of inversion is added to the
          list of symmetry matrices considered in the determination of
          the multiplicity.
       */
      int multiplicity(const Miller::Index& H, bool FriedelFlag) const;
      //! Determine "epsilon" for the given Miller index.
      /*! The factor epsilon counts the number of times a Miller
          index H is mapped onto itself by symmetry. This factor
          is used for "statistical averaging" and in direct methods
          formulae.<br>
          See also: sgtbx::SymEquivIndices::epsilon()
       */
      int epsilon(const Miller::Index& H) const;
      //! Generate list of symmetrically equivalent reflections.
      /*! See class SymEquivMillerIndices
       */
      SymEquivMillerIndices
      getEquivMillerIndices(const Miller::Index& H) const;

      //! Structure factor without Debye-Waller factor.
      /*! Sum of exp(2 pi j H S X) over all symmetry operations S.
          j is the imaginary number.
       */
      template <class FloatType>
      std::complex<FloatType>
      StructureFactor(const Miller::Index& H,
                      const cctbx::fractional<FloatType> X) const
      {
        using cctbx::constants::pi;
        std::complex<FloatType> F(0., 0.);
        for (int s = 0; s < nSMx(); s++) {
          Miller::Index HR = H * m_SMx[s].Rpart();
          FloatType HRX = HR * X;
          TrVec T = m_SMx[s].Tpart();
          for (int i = 0; i < fInv(); i++) {
            if (i) {
              HRX = -HRX;
              T = m_InvT - T;
            }
            for (int l = 0; l < nLTr(); l++) {
              FloatType HT = FloatType(H * (T + LTr(l))) / TBF();
              FloatType phase = 2. * pi * (HRX + HT);
              F += std::complex<FloatType>(std::cos(phase), std::sin(phase));
            }
          }
        }
        return F;
      }
      //! Structure factor with isotropic Debye-Waller factor given Uiso.
      /*! Sum of exp(2 pi j H S X) over all symmetry operations S,
          multiplied by cctbx::adptbx::DebyeWallerFactorUiso().
          j is the imaginary number.
          <p>
          This function is provided for symmetry with the
          structure factor calculation given anisotropic
          displacement parameters.
       */
      template <class FloatType>
      std::complex<FloatType>
      StructureFactor(const uctbx::UnitCell& uc,
                      const Miller::Index& H,
                      const cctbx::fractional<FloatType> X,
                      double Uiso) const
      {
        return
          StructureFactor(H, X) * adptbx::DebyeWallerFactorUiso(uc, H, Uiso);
      }
      //! Structure factor with isotropic Debye-Waller factor given Uiso.
      /*! Sum of exp(2 pi j H S X) over all symmetry operations S,
          multiplied by cctbx::adptbx::DebyeWallerFactorUiso().
          j is the imaginary number.
          <p>
          This function is provided for symmetry with the
          structure factor calculation given anisotropic
          displacement parameters.
       */
      template <class FloatType>
      std::complex<FloatType>
      StructureFactor(double stol2,
                      const Miller::Index& H,
                      const cctbx::fractional<FloatType> X,
                      double Uiso) const
      {
        return
          StructureFactor(H, X) * adptbx::DebyeWallerFactorUiso(stol2, Uiso);
      }
      //! Structure factor with anisotropic Debye-Waller factor given Ustar.
      /*! Sum of cctbx::adptbx::DebyeWallerFactorUstar() * exp(2 pi j H S X)
          over all symmetry operations S.
          j is the imaginary number.
       */
      template <class FloatType>
      std::complex<FloatType>
      StructureFactor(const Miller::Index& H,
                      const cctbx::fractional<FloatType> X,
                      const boost::array<FloatType, 6>& Ustar) const
      {
        using cctbx::constants::pi;
        std::complex<FloatType> F(0., 0.);
        for (int s = 0; s < nSMx(); s++) {
          Miller::Index HR = H * m_SMx[s].Rpart();
          FloatType HRX = HR * X;
          TrVec T = m_SMx[s].Tpart();
          std::complex<FloatType> Fs(0., 0.);
          for (int i = 0; i < fInv(); i++) {
            if (i) {
              HRX = -HRX;
              T = m_InvT - T;
            }
            for (int l = 0; l < nLTr(); l++) {
              FloatType HT = FloatType(H * (T + LTr(l))) / TBF();
              FloatType phase = 2. * pi * (HRX + HT);
              Fs += std::complex<FloatType>(std::cos(phase), std::sin(phase));
            }
          }
          F += Fs * adptbx::DebyeWallerFactorUstar(HR, Ustar);
        }
        return F;
      }

      //! Check if a unit cell is compatible with the symmetry operations.
      /*! This function is designed to work together with the UnitCell
          class in the uctbx:<pre>
        sgtbx::SpaceGroup sg = ... // define space group
        uctbx::UnitCell uc = ... // define unit cell
        bool OK = sg.isCompatibleMetricalMatrix(uc.getMetricalMatrix());</pre>
          A given unit cell is compatible with a given space group
          representation if the following relation holds for all
          rotation matrices R of the space group:<p>
          Transpose[R].G.R == G<p>
          G is the metrical matrix for the unit cell.<br>
          The tolerance compensates for rounding errors.<br>
          See also: isCompatibleUnitCell()
       */
      bool isCompatibleMetricalMatrix(const uctbx::Mx33& G,
                                      double tolerance = 1.e-4) const;
      //! Check if a unit cell is compatible with the symmetry operations.
      /*! Similar to isCompatibleMetricalMatrix(), but an exception
          is thrown if the unit cell is incompatible with the
          symmetry operations.<br>
          See also: CheckUnitCell()
       */
      void CheckMetricalMatrix(const uctbx::Mx33& G,
                               double tolerance = 1.e-4) const;
      //! Check if a unit cell is compatible with the symmetry operations.
      /*! Similar to isCompatibleMetricalMatrix(),
          but the unit cell is passed instead of the metrical matrix.
          The metrical matrix is computed from the given unit cell.<br>
          The tolerance compensates for rounding errors.
       */
      bool isCompatibleUnitCell(const uctbx::UnitCell& uc,
                                       double tolerance = 1.e-4) const {
        return isCompatibleMetricalMatrix(uc.getMetricalMatrix(), tolerance);
      }
      //! Check if a unit cell is compatible with the symmetry operations.
      /*! Similar to isCompatibleUnitCell(), but an exception
          is thrown if the unit cell is incompatible with the
          symmetry operations.
       */
      void CheckUnitCell(const uctbx::UnitCell& uc,
                         double tolerance = 1.e-4) const {
        CheckMetricalMatrix(uc.getMetricalMatrix(), tolerance);
      }

      //! The translation parts of the symmetry operations are set to 0.
      /*! If DiscardZ = true, the lattice translation vectors are not
          modified. If AddInv = true, a centre of inversion is added
          at the origin.
          <p>
          See also:
            BuildDerivedPattersonGroup(),
            BuildDerivedPointGroup(),
            BuildDerivedLaueGroup()
       */
      SpaceGroup BuildDerivedGroup(bool DiscardZ, bool AddInv) const;

      //! Build the corresponding Patterson space group.
      /*! The translation parts of the symmetry operations are set to 0.
          However, the lattice translation vectors are not modified.
       */
      SpaceGroup BuildDerivedPattersonGroup() const {
        return BuildDerivedGroup(false, true);
      }

      //! Build the corresponding point group.
      /*! The translation parts of the symmetry operations are set to 0,
          and the lattice translation vectors are discarded.
       */
      SpaceGroup BuildDerivedPointGroup() const {
        return BuildDerivedGroup(true, false);
      }

      //! Build the corresponding Laue group.
      /*! The translation parts of the symmetry operations are set to 0,
          the lattice translation vectors are discarded,
          and a centre of inversion is added at the origin.
       */
      SpaceGroup BuildDerivedLaueGroup() const {
        return BuildDerivedGroup(true, true);
      }

      //! Histogram of rotation part types.
      /*! The histogram is accumulated in a std::map. The key
          is the rotation part type (-6,-4,-3,-2,-1,1,2,3,4,6),
          and the value is the number of times this rotation part
          type occurs in the list of representative symmetry
          matrices.<br>
          This function is used by the algorithm for the determination
          of the space group type and its usefulness is probably
          limited to this context.
       */
      std::map<int, int> CountRtypes() const;
      //! Determine the point group type.
      /*! The code returned is a matrix group code. There are
          exactly 32 possible return values, corresponding to
          the 32 crystallographic point group types.
       */
      tables::MatrixGroup::Code getPointGroupType() const;
      //! Match given symmetry operations with the 656 tabulated settings.
      /*! This is a light-weight alternative to using
          class SpaceGroupInfo, but will only work for the
          656 tabulated settings.
          <p>
          This algorithm is not particularly optimized. Therefore
          the runtimes for MatchTabulatedSettings() and
          instantiating class SpaceGroupInfo are in general comparable.
          The main purpose of this algorithm is therefore the
          retrieval of conventional Hermann-Mauguin symbols.
       */
      SpaceGroupSymbols MatchTabulatedSettings() const;
      //! Convenience method for instantiating class SpaceGroupInfo.
      SpaceGroupInfo Info() const;
      /*! \brief Refine gridding such that each grid point is
          mapped onto another grid point by all symmetry operations.
       */
      template <typename GridTupleType>
      GridTupleType refine_gridding(const GridTupleType& grid) const {
        GridTupleType prev_grid = grid;
        GridTupleType ref_grid = grid;
        for (;;) {
          for(int i=0;i<OrderZ();i++) {
            ref_grid = operator()(i).refine_gridding(ref_grid);
          }
          if (prev_grid == ref_grid) break;
          prev_grid = ref_grid;
        }
        return ref_grid;
      }
      //! Refine gridding starting with grid 1,1,1.
      /*! See also: other overload.
       */
      boost::array<int, 3> refine_gridding() const {
        return refine_gridding(array<int, 3>(1, 1, 1));
      }
  };

  //! iostream output operator for class SpaceGroup.
  /*! Provided mainly for debugging purposes.
   */
  std::ostream& operator<<(std::ostream& os, const SpaceGroup& SgOps);

  //! Compute and use change-of-basis matrix to reference setting.
  /*! A space group type is characterized by the space group number
      according to the International Tables for Crystallography,
      Volume A, and a change-of-basis matrix that transforms
      the coordinates of a given setting to a reference setting.
   */
  class SpaceGroupInfo {
    public:
      //! Default constructor.
      /*! Equivalent to SpaceGroupInfo(SpaceGroup())
       */
      SpaceGroupInfo() : m_SgOps(), m_SgNumber(1), m_CBOp() {}
      //! Determine the space group type.
      /*! The algorithm for the determination of the space group
          type is published in
          <A HREF="http://journals.iucr.org/a/issues/1999/02/02/au0146/"
          ><I>Acta Cryst.</I> 1999, <B>A55</B>:383-395</A>.
          The result is the space group number according to the
          International Tables for Crystallography, and a
          change-of-basis matrix that transforms the given space group
          to its reference setting. RBF and TBF are the rotation base
          factor and the translation base factor, respectively, for the
          change-of-basis matrix.<p>
          In general, there are several change-of-basis matrices that
          could be used. If TidyCBOp == true, the operations of the
          affine normalizer are used to determine a "standard"
          change-of-basis matrix. This is, the change-of-basis matrix
          is independent of the order of the given symmetry operations.
          In principle, the space group number and the "standard"
          change-of-basis matrix could be used to test for the
          equivalence of two groups of symmetry operations (but the
          algorithm used for operator== is more efficient).
          On average, the additional time required for the
          determination of the "standard" change-of-basis matrix
          is about 70% of the time for the main determination
          of the space group type. This is, the TidyCBOp == true
          option is relatively expensive. However, the absolute
          runtime is only a small fraction of a second.
       */
      explicit
      SpaceGroupInfo(const SpaceGroup& SgOps,
                     bool TidyCBOp = true,
                     int RBF = CRBF, int TBF = CTBF);
      //! Access to space group passed to the constructor.
      const SpaceGroup& SgOps() const { return m_SgOps; }
      //! Space group number according to the International Tables.
      int SgNumber() const { return m_SgNumber; }
      //! Change-of-basis operator.
      const ChOfBasisOp& CBOp() const { return m_CBOp; }
      //! Get the additional generators of the Euclidean normalizer.
      /*! See International Tables for Crystallography Volume A,
          1983, Table 15.3.2. The generators are tabulated for
          reference settings and transformed to the given setting
          using CBOp().
          <p>
          getK2L = true requests the additional generator from
          column "Inversion through a centre at" of Table 15.3.2.
          <p>
          getL2N = true requests the additional generators from
          column "Further generators" of Table 15.3.2.
          <p>
          See also: class StructureSeminvariant
       */
      std::vector<RTMx>
      getAddlGeneratorsOfEuclideanNormalizer(bool getK2L, bool getL2N) const;
      /*! \brief Add the additional generators of the Euclidean normalizer
          to the space group.
       */
      /*! See also: getAddlGeneratorsOfEuclideanNormalizer(),
                    SpaceGroup::expandSMxVector()
       */
      SpaceGroup
      expandAddlGeneratorsOfEuclideanNormalizer(bool useK2L, bool useL2N)
      {
        SpaceGroup result = SgOps();
        result.expandSMxVector(
          getAddlGeneratorsOfEuclideanNormalizer(useK2L, useL2N));
        return result;
      }
      //! Test for the 22 (11 pairs) enantiomorphic space groups.
      /*! A space group G is enantiomorphic if G and -I.G.-I have
          two different space group types. I is the unit matrix.<br>
          There are 11 pairs of enantiomorphic space groups,
          for example P41 and P43.
          <p>
          The notion of enantiomorphic space groups is connected
          to the notion of <i>affine</i> space group types. There
          are 219 affine space group types, compared to the
          230 conventional space group types. A pair of
          enantiomorphic space groups corresponds to a single
          affine space group type.
          <p>
          See also: getChangeOfHandOp()
       */
      bool isEnantiomorphic() const;
      //! Determine a change-of-hand matrix.
      /*! If the space group is centro-symmetric, the change-of-hand
          matrix is the identity matrix.
          <p>
          If the space group belongs to one of the 22 enantiomorphic
          space group types, the change-of-hand matrix is determined as
          a centre of inversion that is located at the origin of the
          reference setting.
          <p>
          If the space group is not enantiomorphic and not
          centro-symmetric, the change-of-hand matrix is determined as
          a centre of inversion of the Euclidean normalizer.
          <p>
          The change-of-hand matrix can be used to transform the
          symmetry operations to obtain the enantiomorph symmetry
          operations (use SpaceGroup::ChangeBasis()), and to transform
          fractional coordinates. For example:<pre>
          SpaceGroup SgOps = ...;
          SpaceGroupInfo SgInfo(SgOps);
          ChOfBasisOp CHOp = SgInfo.getChangeOfHandOp();
          SpaceGroup EnantiomorphSgOps = SgOps.ChangeBasis(CHOp);
          fractional<double> X = ...;
          fractional<double> EnatiomorphX = CHOp(X);</pre>
          <p>
          See also: getAddlGeneratorsOfEuclideanNormalizer(),
                    isEnantiomorphic(), SpaceGroup::ChangeBasis()
       */
      ChOfBasisOp getChangeOfHandOp() const;
      //! Build a Hall symbol for the given symmetry operations.
      /*! For a given group of symmetry operations, there are in
          general several plausible choices for a Hall symbol.
          The Hall symbol returned by this algorithm is derived
          from the Hall symbol for the reference setting. A change-of-basis
          operator is attached if necessary (i.e. "P 2 (y,z,x)").
          <p>
          If TidyCBOp == true, the returned Hall symbol
          is a reproducible representation of a
          given setting that is independent of the order of
          the symmetry operations.
       */
      std::string BuildHallSymbol(bool TidyCBOp = true) const;
      //! Determine conventional Hermann-Mauguin symbol or Hall symbol.
      /*! First, SgOps().MatchTabulatedSettings() is called. If the given
          symmetry operations correspond to one of the 656 tabulated
          settings, the extended Hermann-Mauguin symbol is returned
          (e.g. "P n n n :2").  Otherwise the string "Hall: " + the
          result of BuildHallSymbol() is returned (e.g. "Hall:  P 2 2
          -1n").
          <p>
          The result of BuildLookupSymbol() can be used as the input
          for the constructor of class SpaceGroupSymbols.
       */
      std::string BuildLookupSymbol() const;
    private:
      SpaceGroup m_SgOps;
      int m_SgNumber;
      ChOfBasisOp m_CBOp;
  };

  inline SpaceGroupInfo SpaceGroup::Info() const {
    return SpaceGroupInfo(*this);
  }

}} // namespace cctbx::sgtbx

#endif // CCTBX_SGTBX_GROUPS_H
