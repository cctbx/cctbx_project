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

#ifndef CCTBX_SGTBX_MILLER_H
#define CCTBX_SGTBX_MILLER_H

#include <complex>
#include <cctbx/miller.h>

namespace cctbx {
  namespace Miller {

    //! Helper class for symmetry equivalent Miller indices.
    class SymEquivIndex {
      private:
        int   m_TBF;
        Index m_HR;
        int   m_HT;
      public:
        //! Constructor. For internal use only.
        inline SymEquivIndex(const Index& HR, int HT, int TBF)
          : m_HR(HR), m_HT(HT) , m_TBF(TBF) { }
        //! Product of Miller index and rotation part of symmetry operation.
        inline const Index& HR() const { return m_HR; }
        //! Product of Miller index and translation part of symmetry operation.
        inline int HT() const { return m_HT; }
        //! Translation base factor.
        /*! This is the factor by which HT() is multiplied.
         */
        inline int TBF() const { return m_TBF; }

        //! Phase of this index in radians, given phase of input Miller index.
        /*! Formula used:<br>
            this_phi = phi - (2 * pi * HT()) / TBF();
         */
        template <class T>
        T Phase_rad(const T& phi) const {
          using cctbx::constants::pi;
          return phi - (2. * pi * HT()) / TBF();
        }
        //! Phase of this index in degrees, given phase of input Miller index.
        /*! Formula used:<br>
            this_phi = phi - (2 * 180 * HT()) / TBF();
         */
        template <class T>
        T Phase_deg(const T& phi) const {
          return phi - (2. * 180. * HT()) / TBF();
        }
        //! Complex value for this index, given complex value for input index.
        /*! Formula used:<br>
            this_F = F * exp(-2 * pi * j * HT() / TBF());<br>
            where j is the imaginary number.
         */
        template <class T>
        std::complex<T> ShiftPhase(const std::complex<T>& F) const {
          using cctbx::constants::pi;
          T theta = (-2. * pi * HT()) / TBF();
          return F * std::polar(1., theta);
        }
    };

    //! Master Miller index class.
    /*! See documentation for sgtbx::SymEquivIndices::getMasterIndex().<br>
        The MasterIndex class supports the following data layouts:<p>
        1. Assuming Friedel symmetry (no anomalous signal), only one
        value is stored for a Friedel pair:<pre>
         h  k  l  F</pre>
        The values associated with h,k,l and -h,-k,-l are assumed to be
        equal, and the phases are related by the equation phi(h,k,l) =
        -phi(-h,-k,-l).<br>
        In this case, MasterIndex.H() is intended for use as the
        representative ("asymmetric") index in a list.<p>
        2. No Friedel symmetry (i.e. in presence of an anomalous signal),
        two columns of data:<pre>
         h  k  l  F+  F-</pre>
        Both Friedel mates are associated with the same index in a list.
        The Miller index for F+ is (h, k, l) and the implied Miller index
        for F- is (-h, -k, -l).<br>
        In this case, MasterIndex.H() is intended for use as the
        representative ("asymmetric") index in a list.
        MasterIndex.iMate() is intended for selecting the data column.<p>
        3. No Friedel symmetry (i.e. in presence of an anomalous signal),
        one column of data:<pre>
         h  k  l  F
        -h -k -l  F</pre>
        There is a separate entry for each Friedel mate in a list.<br>
        In this case, MasterIndex.HR() is intended for use as the
        representative ("asymmetric") index in a list.
     */
    class MasterIndex {
      private:
        int   m_TBF;
        bool  m_haveCutP;
        Vec3  m_CutP;
        bool  m_Pretty;
        Index m_H;
        bool  m_iMate;
        Index m_HR;
        int   m_HT;

        inline bool isInActiveArea(const Miller::Index& H) {
          for(int i=0;i<3;i++) if (m_CutP[i] == 0 && H[i] < 0) return false;
          return true;
        }

      public:
        //! For internal use only.
        inline MasterIndex()
          : m_TBF(0), m_haveCutP(false) {}
        //! For internal use only.
        inline MasterIndex(bool Pretty)
          : m_TBF(0), m_haveCutP(false), m_Pretty(Pretty) {}
        //! For internal use only.
        inline MasterIndex(const Vec3& CutP, bool Pretty)
          : m_TBF(0), m_haveCutP(true), m_CutP(CutP), m_Pretty(Pretty) {}

        //! For internal use only.
        inline void Update(const Miller::Index& TrialH,
                           int iMate,
                           const Miller::SymEquivIndex& HS) {
          if (m_haveCutP && !isInActiveArea(TrialH)) return;
          if (m_TBF == 0
              || (m_Pretty ? TrialH < m_H // operator< is very slow
                           : Miller::hashCompare()(TrialH, m_H))) {
            m_H = TrialH;
            m_iMate = iMate;
            m_HR = HS.HR();
            m_HT = HS.HT();
            m_TBF = HS.TBF();
          }
        }

        //! Test if cut parameters were used in determining the master index.
        inline bool haveCutP() const { return m_haveCutP; }
        //! The cut parameters used in determining the master index.
        /*! If haveCutP() == false, the cut parameters are not defined.
         */
        inline const Vec3& CutP() const { return m_CutP; }
        //! Test if the Pretty flag was used in determining the master index.
        /*! If Pretty == true, a "nice" master index for human readable
            listings was requested.
         */
        inline bool Pretty() const { return m_Pretty; }

        //! The master index for data layouts 1 and 2. See class details.
        inline const Index& H() const { return m_H; }
        //! Selection of F+ or F- data column. See class details.
        /*! The values returned are 0 for F+, and 1 for F-.
         */
        inline int iMate() const { return m_iMate; }
        //! The master index for data layout 3. See class details.
        /*! This index is the product of the input Miller index for which
            the master index was computed, and the rotation part of a
            symmetry operation.
         */
        inline const Index& HR() const { return m_HR; }
        //! Phase shift for HR() with respect to the input Miller index.
        /*! Low-level information for computing the phase of the master
            index, given the phase of the input Miller index. Note that
            high-level functions are also available (e.g. ShiftPhase()).<br>
            HT() is multiplied by a base factor TBF() in order to obtain
            an integer value.
         */
        inline int HT() const { return m_HT; }
        //! Translation base factor.
        /*! This is the factor by which HT() is multiplied.
         */
        inline int TBF() const { return m_TBF; }

        //! Phase of master index in radians, given phase of input Miller index.
        /*! Formula used:<br>
            master_phi = phi - (2 * pi * HT()) / TBF();<br>
            if (FriedelSym && iMate()) master_phi = -master_phi;
         */
        template <class T>
        T Phase_rad(const T& phi, bool FriedelSym) const {
          using cctbx::constants::pi;
          T master_phi = phi - (2. * pi * HT()) / TBF();
          if (FriedelSym && iMate()) return -master_phi;
          return master_phi;
        }
        //! Phase of master index in degrees, given phase of input Miller index.
        /*! Formula used:<br>
            master_phi = phi - (2 * 180 * HT()) / TBF();<br>
            if (FriedelSym && iMate()) master_phi = -master_phi;
         */
        template <class T>
        T Phase_deg(const T& phi, bool FriedelSym) const {
          T master_phi = phi - (2. * 180. * HT()) / TBF();
          if (FriedelSym && iMate()) return -master_phi;
          return master_phi;
        }
        //! Complex value for master index, given complex value for input index.
        /*! Formula used:<br>
            master_F = F * exp(-2 * pi * j * HT() / TBF());<br>
            where j is the imaginary number.
         */
        template <class T>
        std::complex<T> ShiftPhase(const std::complex<T>& F,
                                   bool FriedelSym) const {
          using cctbx::constants::pi;
          T theta = (-2. * pi * HT()) / TBF();
          std::complex<T> master_F = F * std::polar(1., theta);
          if (FriedelSym && iMate()) return std::conj(master_F);
          return master_F;
        }
    };

  } // namespace Miller
} // namespace cctbx

namespace sgtbx {

  inline Miller::Index operator*(const RotMx& lhs, const Miller::Index& rhs) {
    return Miller::Index(
      lhs[0] * rhs[0] + lhs[3] * rhs[1] + lhs[6] * rhs[2],
      lhs[1] * rhs[0] + lhs[4] * rhs[1] + lhs[7] * rhs[2],
      lhs[2] * rhs[0] + lhs[5] * rhs[1] + lhs[8] * rhs[2]);
  }

  inline Miller::Index operator*(const Miller::Index& lhs, const RotMx& rhs) {
    return rhs * lhs;
  }

  inline int operator*(const Miller::Index& lhs, const TrVec& rhs) {
    int result = 0;
    for(int i=0;i<3;i++) result += lhs[i] * rhs[i];
    return result;
  }

  inline int HT_mod_1(const Miller::Index& H, const TrVec& T) {
    return sgtbx::modPositive(H * T, T.BF());
  }

  //! class for the high-level handling of centric reflections.
  /*! A reflection with the Miller index H is "centric" if
      there is a symmetry operation with rotation part R such
      that H*R = -H.<br>
      The phase of a centric reflection is restricted to two phase
      angels (modulo pi).
   */
  class PhaseRestriction {
    private:
      int m_HT;
      int m_TBF;
    public:
      //! For internal use only.
      inline PhaseRestriction(int HT, int TBF) : m_HT(HT), m_TBF(TBF) {}

      //! Test if there actually is a phase restriction.
      /*! See class details.
       */
      inline bool isCentric() const { return m_HT >= 0; }

      //! Phase shift H*T (mod 1) corresponding to H*R = -H.
      /*! Low-level information for computing the restricted phases.
          Note that high-level functions are also available (e.g.
          HT_deg()).<br>
          HT() is multiplied by a base factor TBF() in order to obtain
          an integer value.
       */
      inline int HT() const { return m_HT; }
      //! Translation base factor.
      /*! This is the factor by which HT() is multiplied.
       */
      inline int TBF() const { return m_TBF; }

      //! Compute the phase restriction for the given Period.
      /*! Formula used:<br>
          Period * HT() / TBF();<br>
          The return value is -1 if the phase is not restricted,
          and >= 0 and < Period otherwise.<br>
          See also: HT_rad(), HT_deg()
       */
      inline double HT(double Period) const {
        if (!isCentric()) return -1.;
        return Period * static_cast<double>(m_HT)
                      / static_cast<double>(m_TBF);
      }
      //! Compute the phase restriction in radians.
      /*! The return value is -1 if the phase is not restricted,
          and >= 0 and < pi otherwise.
       */
      inline double HT_rad() const { return HT(cctbx::constants::pi); }
      //! Compute the phase restriction in degrees.
      /*! The return value is -1 if the phase is not restricted,
          and >= 0 and < 180 otherwise.
       */
      inline double HT_deg() const { return HT(180.); }

      //! Test if phase phi (with given Period) is compatible with restriction.
      /*! Period is the period of the restricted phase, e.g. 180 if
          the phase is measured in degrees.<br>
          The tolerance compensates for rounding errors.<br>
          See also: isValidPhase_rad(), isValidPhase_deg()
       */
      bool isValidPhase(double Period, double phi, double tolerance) const;
      //! Test if phase phi (in radians) is compatible with restriction.
      /*! The tolerance compensates for rounding errors.
       */
      inline bool isValidPhase_rad(double phi,
                                   double tolerance = 1.e-5) const {
        return isValidPhase(cctbx::constants::pi, phi, tolerance);
      }
      //! Test if phase phi (in degrees) is compatible with restriction.
      /*! The tolerance compensates for rounding errors.
       */
      inline bool isValidPhase_deg(double phi,
                                   double tolerance = 1.e-5) const {
        return isValidPhase(180., phi, tolerance);
      }
  };

  //! class for the handling of symmetry equivalent Miller indices.
  /*! This class is exclusively used to represent the results
      of SgOps::getEquivMillerIndices().<br>
      The Miller index used in the call to getEquivMillerIndices
      is referred to as the "input Miller index."
   */
  class SymEquivMillerIndices {
    friend class SgOps;
    private:
      int m_TBF;
      int m_OrderP;
      int m_HT_Restriction;
      std::vector<Miller::SymEquivIndex> m_List;
      inline SymEquivMillerIndices(int TBF, int OrderP)
        : m_TBF(TBF), m_OrderP(OrderP), m_HT_Restriction(-1), m_List() {}
      void add(const Miller::SymEquivIndex& SEI);
      void setMasterIndex(Miller::MasterIndex& Master) const;
    public:
      //! The phase restriction (if any) for the input Miller index.
      /*! See class PhaseRestriction.
       */
      inline PhaseRestriction getPhaseRestriction() const {
        return PhaseRestriction(m_HT_Restriction, m_TBF); }
      //! Test if reflection with input Miller index is centric.
      /*! A reflection with the Miller index H is "centric" if
          there is a symmetry operation with rotation part R such
          that H*R = -H.<br>
          See also: class PhaseRestriction
       */
      inline bool isCentric() const { return m_HT_Restriction >= 0; }
      //! Number of symmetry equivalent Miller indices.
      /*! Note that this is not in general equal to the multiplicity.<br>
          See also: M()
       */
      inline int N() const { return m_List.size(); }
      //! Multiplicity of the input Miller index.
      /*! For acentric reflections and in the presence of Friedel symmetry
          (no anomalous signal), the multiplicity is twice the number
          of symmetry equivalent Miller indices N().<br>
          For centric reflections or in the absence of Friedel symmetry
          (i.e. in the presence of an anomalous signal), the multiplicity
          is equal to the number of symmetry equivalent Miller indices
          N().<br>
          See also: N()
       */
      inline int M(bool FriedelSym) const {
        if (FriedelSym && !isCentric()) return 2 * N();
        return N();
      }
      //! Flag for Friedel mates == M(FriedelSym)/N().
      /*! Useful for looping over all symmetry equivalent reflections
          including Friedel mates (if FriedelSym == true).<br>
          See also: operator()
       */
      inline int fMates(bool FriedelSym) const {
        if (FriedelSym && !isCentric()) return 2;
        return 1;
      }
      //! Determine "epsilon" for the given Miller index.
      /*! The factor epsilon counts the number of times a Miller
          index H is mapped onto itself by symmetry. This factor
          is used for "statistical averaging" and in direct methods
          formulae.<br>
          Note that epsilon is directly related to the number
          of symmetry equivalent indices N():<br>
          epsilon == sgtbx::SgOps::OrderP() / N()
       */
      inline int epsilon() const { return m_OrderP / N(); }

      //! Low-level access to the N() symmetry equivalent Miller indices.
      /*! See also: operator()
       */
      inline const Miller::SymEquivIndex& operator[](int iList) const {
        return m_List[iList];
      }
      //! Medium-level access to the symmetry equivalent Miller indices.
      /*! Intended use:<pre>
          sgtbx::SgOps sg = ... // define space group
          Miller::Index H = ... // define input Miller index.
          bool FriedelSym = ... // define Friedel symmetry.
          sgtbx::SymEquivMillerIndices SEMI = getEquivMillerIndices(H);
          for (int iList = 0; iList < SEMI.N(); iList++)
            for (int iMate = 0; iMate < SEMI.fMate(FriedelSym); iMate++)
              Miller::Index EquivH = SEMI(iMate, iList);
          </pre>
          Note that it is possible and often more convenient to have a
          one-deep loop with M() iterations.<br>
          See also: operator()(int iIL)
       */
      Miller::Index operator()(int iMate, int iList) const;
      //! High-level access to the symmetry equivalent Miller indices.
      /*! Intended use:<pre>
          sgtbx::SgOps sg = ... // define space group
          Miller::Index H = ... // define input Miller index.
          bool FriedelSym = ... // define Friedel symmetry.
          sgtbx::SymEquivMillerIndices SEMI = getEquivMillerIndices(H);
          for (int iIL = 0; iIL < SEMI.M(FriedelSym); iIL++)
            Miller::Index EquivH = SEMI(iIL);
          </pre>
       */
      Miller::Index operator()(int iIL) const;

      //! Test if phase phi (in radians) is compatible with restriction.
      /*! The tolerance compensates for rounding errors.<br>
          See also: class PhaseRestriction
       */
      inline bool isValidPhase_rad(double phi,
                                   double tolerance = 1.e-5) const {
        return getPhaseRestriction().isValidPhase_rad(phi, tolerance);
      }
      //! Test if phase phi (in degrees) is compatible with restriction.
      /*! The tolerance compensates for rounding errors.<br>
          See also: class PhaseRestriction
       */
      inline bool isValidPhase_deg(double phi,
                                   double tolerance = 1.e-5) const {
        return getPhaseRestriction().isValidPhase_deg(phi, tolerance);
      }

      //! Determine a representative ("asymmetric") Miller index.
      /*! A light-weight, general alternative to using contiguous
          asymmetric units.<br>
          The symmetry equivalent indices of the input index H are
          computed and sorted in a certain order. The first
          index in the sorted list is returned along with
          information needed to compute the phase shift between
          the original index and the master index.<br>
          The sort order is dependent on the Pretty parameter:<br>
          With Pretty == false, the comparison function used for
          sorting the indices is simple and fast. This is useful
          for building a large list of indices quickly.<br>
          With Pretty == true, the comparison function used for
          sorting the indices favors "nice" indices intended for
          human-readable listings. This comparison function is
          significantly slower than the one used if Pretty == false.
       */
      Miller::MasterIndex getMasterIndex(bool Pretty = false) const;
      //! Determine a representative ("asymmetric") Miller index.
      /*! Similar to getMasterIndex(Pretty), but the master index
          is only selected from the active region defined by the
          cut parameters. This is useful for building large lists
          of indices quickly.<br>
          See also: sgtbx::SgOps::getCutParameters()
       */
      Miller::MasterIndex getMasterIndex(const Miller::Vec3& CutP,
                                         bool Pretty = false) const;
  };

} // namespace sgtbx

#endif // CCTBX_SGTBX_MILLER_H
