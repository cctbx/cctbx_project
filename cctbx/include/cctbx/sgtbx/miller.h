// $Id$
/* Copyright (c) 2001 The Regents of the University of California through
   E.O. Lawrence Berkeley National Laboratory, subject to approval by the
   U.S. Department of Energy. See files COPYRIGHT.txt and
   cctbx/LICENSE.txt for further details.

   Revision history:
     2001 Oct 27: Redesign: AsymIndex (rwgk)
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

    //! Class for symmetrically equivalent Miller indices.
    class SymEquivIndex {
      public:
        //! Default constructor. Some data members are not initialized!
        SymEquivIndex() {}
        //! Constructor.
        /*! HR is the product of the input index H and the
            rotation part of a symmetry operation.
            <p>
            HT is the product of the input index H and the
            translation part of the symmetry operation,
            multiplied by the translation base-factor TBF
            in order to obtain an integer value for HT.
         */
        SymEquivIndex(const Index& HR, int HT, int TBF)
          : m_HR(HR), m_HT(HT) , m_TBF(TBF) {}
        //! Product of Miller index and rotation part of symmetry operation.
        const Index& HR() const { return m_HR; }
        //! Product of Miller index and translation part of symmetry operation.
        int HT() const { return m_HT; }
        //! Translation base factor.
        /*! This is the factor by which HT() is multiplied.
         */
        int TBF() const { return m_TBF; }

        /*! \brief Phase for equivalent index in radians, given phase for
            input Miller index.
         */
        /*! Formula used:<br>
            phi_eq = phi_in - (2 * pi * HT()) / TBF();
         */
        template <class FloatType>
        FloatType phase_eq_rad(const FloatType& phi_in) const {
          using cctbx::constants::pi;
          return phi_in - (2. * pi * HT()) / TBF();
        }
        /*! \brief Phase for input index in radians, given phase for
            equivalent Miller index.
         */
        /*! Formula used:<br>
            phi_in = phi_eq + (2 * pi * HT()) / TBF();
         */
        template <class FloatType>
        FloatType phase_in_rad(const FloatType& phi_eq) const {
          using cctbx::constants::pi;
          return phi_eq + (2. * pi * HT()) / TBF();
        }
        /*! \brief Phase for equivalent index in degrees, given phase for
            input Miller index.
         */
        /*! Formula used:<br>
            phi_eq = phi_in - (2 * pi * HT()) / TBF();
         */
        template <class FloatType>
        FloatType phase_eq_deg(const FloatType& phi_in) const {
          return phi_in - (2. * 180. * HT()) / TBF();
        }
        /*! \brief Complex value for equivalent index, given complex
            value for input index.
         */
        /*! Formula used:<br>
            f_eq = f_in * exp(-2 * pi * j * HT() / TBF());<br>
            where j is the imaginary number.
         */
        template <class FloatType>
        std::complex<FloatType>
        complex_eq(const std::complex<FloatType>& f_in) const {
          using cctbx::constants::pi;
          FloatType theta = (-2. * pi * HT()) / TBF();
          return f_in * std::polar(1., theta);
        }
      protected:
        Index m_HR;
        int   m_HT;
        int   m_TBF;
    };

  } // namespace Miller

  namespace sgtbx {

#ifndef DOXYGEN_SHOULD_SKIP_THIS

    inline Miller::Index
    operator*(const Miller::Index& lhs, const RotMx& rhs) {
      return Miller::Index(
        lhs[0] * rhs[0] + lhs[1] * rhs[3] + lhs[2] * rhs[6],
        lhs[0] * rhs[1] + lhs[1] * rhs[4] + lhs[2] * rhs[7],
        lhs[0] * rhs[2] + lhs[1] * rhs[5] + lhs[2] * rhs[8]);
    }

    inline int operator*(const Miller::Index& lhs, const TrVec& rhs) {
      int result = 0;
      for(int i=0;i<3;i++) result += lhs[i] * rhs[i];
      return result;
    }

    inline int HT_mod_1(const Miller::Index& H, const TrVec& T) {
      return sgtbx::modPositive(H * T, T.BF());
    }

    class SpaceGroup; // forward declaration

#endif // DOXYGEN_SHOULD_SKIP_THIS

    //! XXX
    class sys_absent_test
    {
      public:
        //! XXX
        sys_absent_test() {}
        //! XXX
        sys_absent_test(const SpaceGroup& sgops, const Miller::Index& h);
        //! XXX
        int ht_restriction() const { return ht_restriction_; }
        //! XXX
        bool is_sys_absent() const { return ht_restriction_ == -2; }
        //! XXX
        bool is_centric() const { return ht_restriction_ >= 0; }
      protected:
        int ht_restriction_;
    };

    //! class for the high-level handling of centric reflections.
    /*! A reflection with the Miller index H is "centric" if
        there is a symmetry operation with rotation part R such
        that H*R = -H.<br>
        The phase of a centric reflection is restricted to two phase
        angels (modulo pi).
     */
    class PhaseRestriction {
      public:
        //! Default constructor. Some data members are not initialized!
        PhaseRestriction() {}
        //! For internal use only.
        PhaseRestriction(int HT, int TBF) : m_HT(HT), m_TBF(TBF) {}

        //! Test if there actually is a phase restriction.
        /*! See class details.
         */
        bool isCentric() const { return m_HT >= 0; }

        //! Phase shift H*T (mod 1) corresponding to H*R = -H.
        /*! Low-level information for computing the restricted phases.
            Note that high-level functions are also available (e.g.
            HT_deg()).<br>
            HT() is multiplied by a base factor TBF() in order to obtain
            an integer value.
         */
        int HT() const { return m_HT; }
        //! Translation base factor.
        /*! This is the factor by which HT() is multiplied.
         */
        int TBF() const { return m_TBF; }

        //! Compute the phase restriction for the given Period.
        /*! Formula used:<br>
            Period * HT() / TBF();<br>
            The return value is -1 if the phase is not restricted,
            and >= 0 and < Period otherwise.<br>
            See also: HT_rad(), HT_deg()
         */
        double HT(double Period) const {
          if (!isCentric()) return -1.;
          return Period * double(m_HT) / double(m_TBF);
        }
        //! Compute the phase restriction in radians.
        /*! The return value is -1 if the phase is not restricted,
            and >= 0 and < pi otherwise.
         */
        double HT_rad() const { return HT(cctbx::constants::pi); }
        //! Compute the phase restriction in degrees.
        /*! The return value is -1 if the phase is not restricted,
            and >= 0 and < 180 otherwise.
         */
        double HT_deg() const { return HT(180.); }

        /*! \brief Test if phase phi (with given Period) is
            compatible with restriction.
         */
        /*! Period is the period of the restricted phase, e.g. 180 if
            the phase is measured in degrees.<br>
            The tolerance compensates for rounding errors.<br>
            See also: isValidPhase_rad(), isValidPhase_deg()
         */
        bool isValidPhase(double Period, double phi, double tolerance) const;
        //! Test if phase phi (in radians) is compatible with restriction.
        /*! The tolerance compensates for rounding errors.
         */
        bool isValidPhase_rad(double phi,
                                     double tolerance = 1.e-5) const {
          return isValidPhase(cctbx::constants::pi, phi, tolerance);
        }
        //! Test if phase phi (in degrees) is compatible with restriction.
        /*! The tolerance compensates for rounding errors.
         */
        bool isValidPhase_deg(double phi,
                                     double tolerance = 1.e-5) const {
          return isValidPhase(180., phi, tolerance);
        }
      private:
        int m_HT;
        int m_TBF;
    };

    //! class for the handling of symmetrically equivalent Miller indices.
    /*! This class is exclusively used to represent the results
        of SpaceGroup::getEquivMillerIndices().<br>
        The Miller index used in the call to getEquivMillerIndices
        is referred to as the "input Miller index."
     */
    class SymEquivMillerIndices
    {
      public:
        //! Default constructor. Some data members are not initialized!
        SymEquivMillerIndices() {}
        //! The phase restriction (if any) for the input Miller index.
        /*! See class PhaseRestriction.
         */
        PhaseRestriction getPhaseRestriction() const {
          return PhaseRestriction(m_HT_Restriction, m_TBF); }
        //! Test if reflection with input Miller index is centric.
        /*! A reflection with the Miller index H is "centric" if
            there is a symmetry operation with rotation part R such
            that H*R = -H.<br>
            See also: class PhaseRestriction
         */
        bool isCentric() const { return m_HT_Restriction >= 0; }
        //! Number of symmetrically equivalent Miller indices.
        /*! Note that this is not in general equal to the multiplicity.<br>
            See also: M()
         */
        int N() const { return m_List.size(); }
        //! Multiplicity of the input Miller index.
        /*! For acentric reflections and in the presence of Friedel symmetry
            (no anomalous signal), the multiplicity is twice the number
            of symmetrically equivalent Miller indices N().<br>
            For centric reflections or in the absence of Friedel symmetry
            (i.e. in the presence of an anomalous signal), the multiplicity
            is equal to the number of symmetrically equivalent Miller indices
            N().<br>
            See also: N()
         */
        int M(bool FriedelFlag) const {
          if (FriedelFlag && !isCentric()) return 2 * N();
          return N();
        }
        //! Flag for Friedel mates == M(FriedelFlag)/N().
        /*! Useful for looping over all symmetrically equivalent reflections
            including Friedel mates (if FriedelFlag == true).<br>
            See also: operator()
         */
        int fMates(bool FriedelFlag) const {
          if (FriedelFlag && !isCentric()) return 2;
          return 1;
        }
        //! Determine "epsilon" for the given Miller index.
        /*! The factor epsilon counts the number of times a Miller
            index H is mapped onto itself by symmetry. This factor
            is used for "statistical averaging" and in direct methods
            formulae.<br>
            Note that epsilon is directly related to the number
            of symmetrically equivalent indices N():<br>
            epsilon == sgtbx::SpaceGroup::OrderP() / N()
         */
        int epsilon() const { return m_OrderP / N(); }

        /*! \brief Low-level access to the N() symmetrically
            equivalent Miller indices.
         */
        /*! See also: operator()
         */
        const Miller::SymEquivIndex& operator[](int iList) const {
          return m_List[iList];
        }
        //! Medium-level access to the symmetrically equivalent Miller indices.
        /*! Intended use:<pre>
            sgtbx::SpaceGroup SgOps = ... // define space group
            Miller::Index H = ... // define input Miller index.
            bool FriedelFlag = ... // define Friedel symmetry.
            sgtbx::SymEquivMillerIndices SEMI = SgOps.getEquivMillerIndices(H);
            for (int iList = 0; iList < SEMI.N(); iList++)
              for (int iMate = 0; iMate < SEMI.fMates(FriedelFlag); iMate++)
                Miller::Index EquivH = SEMI(iMate, iList);
            </pre>
            Note that it is possible and often more convenient to have a
            one-deep loop with M() iterations.<br>
            See also: operator()(int iIL)
         */
        Miller::Index operator()(int iMate, int iList) const;
        //! High-level access to the symmetrically equivalent Miller indices.
        /*! Intended use:<pre>
            sgtbx::SpaceGroup SgOps = ... // define space group
            Miller::Index H = ... // define input Miller index.
            bool FriedelFlag = ... // define Friedel symmetry.
            sgtbx::SymEquivMillerIndices SEMI = SgOps.getEquivMillerIndices(H);
            for (int iIL = 0; iIL < SEMI.M(FriedelFlag); iIL++)
              Miller::Index EquivH = SEMI(iIL);
            </pre>
         */
        Miller::Index operator()(int iIL) const;

        /*! \brief Shift phase of structure factor f_in corresponding to
            input Miller index to structure factor corresponding to
            operator()(iMate, iList).
         */
#if !(defined(BOOST_MSVC) && BOOST_MSVC <= 1200)
        template <class FloatType>
        std::complex<FloatType>
#else
        inline
        std::complex<double>
#endif
        complex_eq(
          int iMate, int iList,
#if !(defined(BOOST_MSVC) && BOOST_MSVC <= 1200)
          const std::complex<FloatType>& f_in
#else
          const std::complex<double>& f_in
#endif
          ) const;

        /*! \brief Shift phase of structure factor F corresponding to
            input Miller index to structure factor corresponding to
            operator()(iIL).
         */
#if !(defined(BOOST_MSVC) && BOOST_MSVC <= 1200)
        template <class FloatType>
        std::complex<FloatType>
#else
        inline
        std::complex<double>
#endif
        complex_eq(
          int iIL,
#if !(defined(BOOST_MSVC) && BOOST_MSVC <= 1200)
          const std::complex<FloatType>& f_in
#else
          const std::complex<double>& f_in
#endif
          ) const;

        //! Test if phase phi (in radians) is compatible with restriction.
        /*! The tolerance compensates for rounding errors.<br>
            See also: class PhaseRestriction
         */
        bool isValidPhase_rad(double phi,
                              double tolerance = 1.e-5) const {
          return getPhaseRestriction().isValidPhase_rad(phi, tolerance);
        }
        //! Test if phase phi (in degrees) is compatible with restriction.
        /*! The tolerance compensates for rounding errors.<br>
            See also: class PhaseRestriction
         */
        bool isValidPhase_deg(double phi,
                              double tolerance = 1.e-5) const {
          return getPhaseRestriction().isValidPhase_deg(phi, tolerance);
        }
      private:
        int m_TBF;
        int m_OrderP;
        int m_HT_Restriction;
        std::vector<Miller::SymEquivIndex> m_List;
        SymEquivMillerIndices(int TBF, int OrderP)
          : m_TBF(TBF), m_OrderP(OrderP), m_HT_Restriction(-1), m_List() {}
        void add(const Miller::SymEquivIndex& SEI);
#ifndef DOXYGEN_SHOULD_SKIP_THIS
        friend class SpaceGroup;
#endif // DOXYGEN_SHOULD_SKIP_THIS
        struct iIL_decomposition
        {
          iIL_decomposition(int iMate_, int iList_)
            : iMate(iMate_), iList(iList_)
          {}
          int iMate, iList;
        };
        iIL_decomposition decompose_iIL(int iIL) const;
    };

#if !(defined(BOOST_MSVC) && BOOST_MSVC <= 1200)
    template <class FloatType>
    std::complex<FloatType>
#else
    inline
    std::complex<double>
#endif
    SymEquivMillerIndices::complex_eq(
      int iMate, int iList,
#if !(defined(BOOST_MSVC) && BOOST_MSVC <= 1200)
      const std::complex<FloatType>& f_in
#else
      const std::complex<double>& f_in
#endif
      ) const
    {
      if (iMate == 0) {
        return m_List[iList].complex_eq(f_in);
      }
      return std::conj(m_List[iList].complex_eq(f_in));
    }

#if !(defined(BOOST_MSVC) && BOOST_MSVC <= 1200)
    template <class FloatType>
    std::complex<FloatType>
#else
    inline
    std::complex<double>
#endif
    SymEquivMillerIndices::complex_eq(
      int iIL,
#if !(defined(BOOST_MSVC) && BOOST_MSVC <= 1200)
      const std::complex<FloatType>& f_in
#else
      const std::complex<double>& f_in
#endif
      ) const
    {
      iIL_decomposition d = decompose_iIL(iIL);
      return complex_eq(d.iMate, d.iList, f_in);
    }

  } // namespace sgtbx
} // namespace cctbx

#endif // CCTBX_SGTBX_MILLER_H
