// $Id$
/* Copyright (c) 2001 The Regents of the University of California through
   E.O. Lawrence Berkeley National Laboratory, subject to approval by the
   U.S. Department of Energy. See files COPYRIGHT.txt and
   cctbx/LICENSE.txt for further details.

   Revision history:
     2001 Jul 02: Merged from CVS branch sgtbx_special_pos (rwgk)
     2001 May 31: merged from CVS branch sgtbx_type (R.W. Grosse-Kunstleve)
     Apr 2001: SourceForge release (R.W. Grosse-Kunstleve)
 */

#ifndef CCTBX_SGTBX_GROUPS_H
#define CCTBX_SGTBX_GROUPS_H

#include <iostream>
#include <vector>
#include <map>
#include <cctbx/sgtbx/symbols.h>
#include <cctbx/sgtbx/matrix.h>
#include <cctbx/sgtbx/change_basis.h>
#include <cctbx/sgtbx/lattice_tr.h>
#include <cctbx/sgtbx/utils.h>
#include <cctbx/sgtbx/miller.h>
#include <cctbx/sgtbx/tables.h>
#include <cctbx/uctbx.h>

namespace sgtbx {

  class TrOps {
    friend class SgOps;
    private:
      std::vector<TrVec> m_Vects;
      bool add(const TrVec& NewTr);
    public:
      inline TrOps() : m_Vects() { m_Vects.push_back(TrVec()); }

      inline void reset() { m_Vects.clear(); m_Vects.push_back(TrVec()); }
      inline const std::vector<TrVec>& Vects() const { return m_Vects; }
      inline std::vector<TrVec>& Vects() { return m_Vects; }
      inline int nVects() const { return m_Vects.size(); }
      bool expand(const TrVec& NewTr);

      TrOps ChangeBasis(const ChOfBasisOp& CBOp) const;

      inline TrVec& operator[](int i) { return m_Vects[i]; }
      inline const TrVec& operator[](int i) const { return m_Vects[i]; }

      char getConventionalCentringTypeSymbol() const;
      ChOfBasisOp getConventionalZ2POp(int RBF = CRBF,
                                       int TBF = CTBF) const;
      TrVec TidyT(const TrVec& T) const;
  };

  //! Result type for SgOps::getSpaceGroupType.
  /*! A space group type is characterized by the space group number
      according to the International Tables for Crystallography,
      Volume A, and a change-of-basis matrix that transforms
      the coordinates of a given setting to a reference setting.
   */
  class SpaceGroupType {
    public:
      //! Constructor. For internal use only.
      inline SpaceGroupType(int SgNumber, const ChOfBasisOp& CBOp)
        : m_SgNumber(SgNumber), m_CBOp(CBOp) {
        cctbx_assert(0 < SgNumber && SgNumber <= 230);
        cctbx_assert(CBOp.M().Rpart().det() > 0);
        cctbx_assert(CBOp.InvM().Rpart().det() > 0);
      }
      //! Space group number according to the International Tables.
      inline int SgNumber() const { return m_SgNumber; }
      //! Change-of-basis operator.
      inline const ChOfBasisOp& CBOp() const { return m_CBOp; }
    private:
      int m_SgNumber;
      ChOfBasisOp m_CBOp;
  };

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
  class SgOps {
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
      inline SgOps(bool NoExpand = false) : m_NoExpand(NoExpand) {
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
      SgOps(parse_string& HSym,
            bool Pedantic = false, bool NoCType = false,
            bool NoExpand = false);
      //! Initialize with symmetry encoded by a Hall symbol.
      /*! Identical to the constructor that takes a parse_string
          as the first argument. However, if an exception is thrown
          there is no way to locate the input character that triggered
          the error.
       */
      SgOps(const std::string& HSym,
            bool Pedantic = false, bool NoCType = false,
            bool NoExpand = false);
      /*! Identical to the constructor that takes a parse_string
          as the first argument. However, if an exception is thrown
          there is no way to locate the input character that triggered
          the error.
       */
      SgOps(const char* HSym,
            bool Pedantic = false, bool NoCType = false,
            bool NoExpand = false);
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
      SgOps ChangeBasis(const ChOfBasisOp& CBOp) const;

      //! Rotation base factor of Seitz Matrices.
      inline int RBF() const { return m_SMx[0].Rpart().BF(); }
      //! Translation base factor of Seitz Matrices.
      inline int TBF() const { return m_SMx[0].Tpart().BF(); }

      //! Number of lattice translations.
      inline int nLTr() const { return m_LTr.nVects(); }
      //! Flag for centre of inversion.
      /*! fInv() == 1 if no inversion operation exists,<br>
          fInv() == 2 otherwise.
       */
      inline int fInv() const { return m_fInv; }
      //! Translation part for centre of inversion.
      inline TrVec InvT(bool tidy = false) const {
        if (tidy == false) return m_InvT;
        if (!m_InvT.isValid()) return m_InvT;
        return m_LTr.TidyT(m_InvT);
      }
      //! Number of representative Seitz matrices.
      inline int nSMx() const { return m_nSMx; }
      //! Order of the point-group = fInv() * nSMx().
      inline int OrderP() const { return m_fInv * m_nSMx; }
      //! Order of the space-group = nLTr() * fInv() * nSMx().
      inline int OrderZ() const { return m_LTr.nVects() * m_fInv * m_nSMx; }

      //! Access to the list of lattice translation vectors.
      inline const TrVec& LTr(int i) const { return m_LTr[i]; }
      //! Access to the list of representative Seitz matrices.
      inline RTMx& operator[](int i) { return m_SMx[i]; }
      //! Access to the list of representative Seitz matrices.
      inline const RTMx& operator[](int i) const { return m_SMx[i]; }
      //! Return a symmetry operation.
      /*! Usage:<pre>
          SgOps s(...);
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
          SgOps s(...);
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
          translation part is selected. The list of representative
          Seitz matrices is sorted, and the translation parts are
          normalized.<br>
          After application of makeTidy(), a given space group
          representation will always result in exactly the same
          internal representation.
       */
      void makeTidy();
      //! Test for equality.
      /*! Internally, makeTidy() is used, followed by essentially a
          byte-wise comparison of the objects.<br>
          Each SgOps object maintains an internal flag
          indicating whether or not makeTidy() was applied already.
          If an SgOps object is repeatedly used in a test for equality,
          the test will therefore be significantly faster if makeTidy()
          is applied outside the loop.
       */
      friend bool operator==(const SgOps& lhs, const SgOps& rhs);
      //! Negation of test for equality.
      inline friend bool operator!=(const SgOps& lhs, const SgOps& rhs) {
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
          SgOps centred(...);
          SgOps primitive = centred.ChangeBasis(centred.getZ2POp());</pre>
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
      inline char getConventionalCentringTypeSymbol() const {
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
      //! Test for the 22 (11 pairs) enantiomorphic space groups.
      /*! A space group G is enantiomorphic if G and -I.G.-I have
          two different space group types. I is the unit matrix.<br>
          There are 11 pairs of enantiomorphic space groups,
          for example P41 and P43.<br>
          The notion of enantiomorphic space groups is connected
          to the notion of <i>affine</i> space group types. There
          are only 219 affine space group types, compared to the
          230 conventional space group types. An pair of
          enantiomorphic space groups corresponds to a single
          affine space group type.
          <p>
          See also: getEnantiomorphSgOps()
       */
      bool isEnantiomorphic() const;
      //! Return the enantiomorph symmetry operations.
      /*! If the given space group belongs to one of the 22
          enantiomorphic space groups (11 pairs), the enantiomorph symmetry
          operations are determined by using -x,-y,-z as a
          change-of-basis matrix. Otherwise the symmetry operations
          are copied unmodified.
          <p>
          See also: isEnantiomorphic()
       */
      SgOps getEnantiomorphSgOps() const;
      //! Test for a centre of inversion.
      inline bool isCentric() const { return m_fInv == 2; }

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
      inline bool isValidPhase_rad(const Miller::Index& H, double phi,
                                   double tolerance = 1.e-5) const {
        return getPhaseRestriction(H).isValidPhase_rad(phi, tolerance);
      }
      //! See class PhaseRestriction::isValidPhase_deg().
      inline bool isValidPhase_deg(const Miller::Index& H, double phi,
                                   double tolerance = 1.e-5) const {
        return getPhaseRestriction(H).isValidPhase_deg(phi, tolerance);
      }
      //! Determine the reflection multiplicity for the given Miller index.
      /*! The multiplicity is defined as the number of symmetry
          equivalent but distinct reflections.<br>
          If FriedelSym == true, a centre of inversion is added to the
          list of symmetry matrices considered in the determination of
          the multiplicity.
       */
      int multiplicity(const Miller::Index& H, bool FriedelSym) const;
      //! Determine "epsilon" for the given Miller index.
      /*! The factor epsilon counts the number of times a Miller
          index H is mapped onto itself by symmetry. This factor
          is used for "statistical averaging" and in direct methods
          formulae.<br>
          See also: sgtbx::SymEquivIndices::epsilon()
       */
      int epsilon(const Miller::Index& H) const;
      //! Generate list of symmetry equivalent reflections.
      /*! See class SymEquivMillerIndices
       */
      SymEquivMillerIndices
      getEquivMillerIndices(const Miller::Index& H) const;
      //! Determine "cut parameters" for building Miller indices.
      /*! When building (or generating) a large list of Miller indices,
          it is useful to restrict the loop over all possible indices
          to 1/2, 1/4, or 1/8 of reciprocal space, if possible.<br>
          The getCutParameters() function returns parameters that can be
          used in the following way:<pre>
          Miller::Index Hmax = ... // determine max hkl
          Miller::Vec3 CutP = getCutParameters(false);
          Miller::Index Hmin;
          for(int i=0;i<3;i++) Hmin[i] = CutP[i] * Hmax[i];
          Miller::Index H;
          for (H[0] = Hmin[0]; H[0] <= Hmax[0]; H[0]++)
            for (H[1] = Hmin[1]; H[1] <= Hmax[1]; H[1]++)
              for (H[2] = Hmin[2]; H[2] <= Hmax[2]; H[2]++)
                // use H</pre>
          This is, each element of CutP is either -1 or 0. A value
          of 0 indicates that the corresponding negative half-space
          can be omitted in the loop over possible indices.<br>
          If FriedelSym == true, a centre of inversion is added to the
          list of symmetry matrices considered in the determination of
          the cut parameters.
       */
      Miller::Vec3 getCutParameters(bool FriedelSym) const;
      //! Determine a representative ("asymmetric") Miller index.
      /*! A light-weight, general alternative to using contiguous
          asymmetric units.<br>
          See sgtbx::SymEquivMillerIndices::getMasterIndex()
       */
      Miller::MasterIndex
      getMasterIndex(const Miller::Index& H,
                     bool Pretty = false) const;
      //! Determine a representative ("asymmetric") Miller index.
      /*! Similar to getMasterIndex(H, Pretty), but the master index
          is only selected from the active region defined by the
          cut parameters. This is useful for building large lists
          of indices quickly.<br>
          See also: getCutParameters(),
          sgtbx::SymEquivMillerIndices::getMasterIndex()
       */
      Miller::MasterIndex
      getMasterIndex(const Miller::Index& H, const Miller::Vec3& CutP,
                     bool Pretty = false) const;

      //! Check if a unit cell is compatible with the symmetry operations.
      /*! This function is designed to work together with the UnitCell
          class in the uctbx:<pre>
        sgtbx::SgOps sg = ... // define space group
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
      inline bool isCompatibleUnitCell(const uctbx::UnitCell& uc,
                                       double tolerance = 1.e-4) const {
        return isCompatibleMetricalMatrix(uc.getMetricalMatrix(), tolerance);
      }
      //! Check if a unit cell is compatible with the symmetry operations.
      /*! Similar to isCompatibleUnitCell(), but an exception
          is thrown if the unit cell is incompatible with the
          symmetry operations.
       */
      inline void CheckUnitCell(const uctbx::UnitCell& uc,
                         double tolerance = 1.e-4) const {
        CheckMetricalMatrix(uc.getMetricalMatrix(), tolerance);
      }

      //! The translation parts of the symmetry operations are set to 0.
      /*! The lattice translation vectors are not modified.<br>
          This function is used by the algorithm for the determination
          of the space group type and its usefulness is probably
          limited to this context.
       */
      SgOps toZPointGroup() const;

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
          change-of-basis matrix.<br>
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
          runtime is only in the order of 1/10 of a second.
       */
      SpaceGroupType getSpaceGroupType(bool TidyCBOp = true,
                                       int RBF = CRBF, int TBF = CTBF) const;
      //! Build a Hall symbol for the given symmetry operations.
      /*! For a given group of symmetry operations, there are in
          general several plausible choices for a Hall symbol.
          The Hall symbol returned by this algorithm is derived
          from the Hall symbol for the reference setting. A change-of-basis
          operator is attached, if necessary (i.e. "P 2 (y,z,x)").<br>
          See getSpaceGroupType() for a discussion of TidyCBOp.<br>
          This algorithm involves the determination of the space group
          type. If the space group type has been determined already,
          the result can be re-used to avoid run-time overhead.
       */
      std::string BuildHallSymbol(const SpaceGroupType& SgType,
                                  bool TidyCBOp = true) const;
      //! Build a Hall symbol for the given symmetry operations.
      /*! See previous overload for this function.
       */
      std::string BuildHallSymbol(bool TidyCBOp = true) const;
      //! Match given symmetry operations with the 530 tabulated settings.
      /*! This is a light-weight alternative to getSpaceGroupType(),
          but will only work for the 530 tabulated settings.<br>
          This algorithm is not particularly optimized. Therefore
          the run-times for MatchTabulatedSettings() and
          getSpaceGroupType() are in general comparable.
          The main purpose of this algorithm is therefore the
          retrieval of conventional Hermann-Mauguin symbols.
       */
      SpaceGroupSymbols MatchTabulatedSettings() const;
      //! Determine conventional Hermann-Mauguin symbol or Hall symbol.
      /*! First, MatchTabulatedSettings() is called. If the given
          symmetry operations correspond to one of the 530 tabulated
          settings, the extended Hermann-Mauguin symbol is returned
          (e.g. "P n n n :2").  Otherwise the string "Hall: " + the
          result of BuildHallSymbol() is returned (e.g. "Hall:  P 2 2
          -1n").<br>
          The result of BuildLookupSymbol() can be used as the input
          for the constructor of class SpaceGroupSymbols.<br>
          This algorithm involves the determination of the space group
          type. If the space group type has been determined already,
          the result can be re-used to avoid run-time overhead.
       */
      std::string BuildLookupSymbol(const SpaceGroupType& SgType) const;
      //! Determine conventional Hermann-Mauguin symbol or Hall symbol.
      /*! See previous overload for this function.
       */
      std::string BuildLookupSymbol() const;
  };

  //! iostream output operator for class SgOps.
  /*! Provided mainly for debugging purposes.
   */
  std::ostream& operator<<(std::ostream& os, const SgOps& sgo);

} // namespace sgtbx

#endif // CCTBX_SGTBX_GROUPS_H
