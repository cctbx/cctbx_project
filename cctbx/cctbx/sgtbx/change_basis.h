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

#ifndef CCTBX_SGTBX_CHANGE_BASIS_H
#define CCTBX_SGTBX_CHANGE_BASIS_H

#include <cctbx/coordinates.h>
#include <cctbx/uctbx.h>
#include <cctbx/sgtbx/matrix.h>

namespace cctbx { namespace sgtbx {

  //! Change-of-basis (transformation) operator.
  /*! For ease of use, a change-of-basis matrix M() and its inverse InvM()
      are grouped by this class.
   */
  class ChOfBasisOp {
    public:
      //! Initialize the change-of-basis operator with M and InvM.
      /*! The input matrices are NOT checked for consistency.
       */
      ChOfBasisOp(const RTMx& M, const RTMx& InvM)
        : m_M(M), m_InvM(InvM) { }
      //! Initialize the change-of-basis operator with M.
      /*! The inverse matrix InvM is computed by inverting M.
          An exception is thrown if M is not invertible.
       */
      explicit ChOfBasisOp(const RTMx& M)
        : m_M(M), m_InvM(M.inverse()) { }
      /*! \brief Initialize the change-of-basis operator with
          matrix given as xyz string.
       */
      /*! The inverse matrix is computed by inversion.
          An exception is thrown if the given matrix is not invertible.
          <p>
          See also: constructor of class RTMx
       */
      ChOfBasisOp(const std::string& StrXYZ,
                         const char* StopChars = "",
                         int RBF = CRBF, int TBF = CTBF)
        : m_M(RTMx(StrXYZ, StopChars, RBF, TBF)), m_InvM(m_M.inverse()) {}
      //! Initialize the change-of-basis operator with unit matrices.
      /*! The unit matrices are initialized with the rotation base factor
          RBF and the translation base factor TBF.
       */
      explicit ChOfBasisOp(int RBF = CRBF, int TBF = CTBF)
        : m_M(RBF, TBF), m_InvM(RBF, TBF) { }
      //! Test if the change-of-basis operator is valid.
      /*! A ChOfBasisOp is valid only if the rotation base factor and the
          translation base factor of both M() and InvM() are not zero.
       */
      bool isValid() const {
        return m_M.isValid() && m_InvM.isValid();
      }
      //! Return a new change-of-basis operator with unit matrices.
      /*! The new matrices inherit the rotation and translation base factors.
       */
      ChOfBasisOp Unit() const {
        return ChOfBasisOp(m_M.Unit(), m_InvM.Unit());
      }
      //! Test if the change-of-basis operator contains unit matrices.
      bool isUnit() const {
        return m_M.isUnit() && m_InvM.isUnit();
      }
      //! Return a new copy with the base factors RBF and TBF.
      /*! An exception is thrown if the elements cannot be scaled to
          the new base factors.<br>
          RBF or TBF == 0 indicates that the corresponding old base
          factor is retained.
       */
      ChOfBasisOp newBaseFactors(int RBF, int TBF) {
        return ChOfBasisOp(m_M.newBaseFactors(RBF, TBF),
                           m_InvM.newBaseFactors(RBF, TBF));
      }
      //! Return a new copy of the operator, but with the base factors of CBOp.
      /*! An exception is thrown if the elements cannot be scaled to
          the new base factors.
       */
      ChOfBasisOp newBaseFactors(const ChOfBasisOp& CBOp) {
        return ChOfBasisOp(m_M.newBaseFactors(CBOp.M()),
                           m_InvM.newBaseFactors(CBOp.InvM()));
      }
      //! Return the change-of-basis matrix.
      const RTMx& M() const { return m_M; }
      //! Return the inverse of the change-of-basis matrix.
      const RTMx& InvM() const { return m_InvM; }
      //! Return M() for Inv == false, and InvM() for Inv == true.
      const RTMx& select(bool Inv) {
        if (Inv) return m_InvM;
        return m_M;
      }
      //! Return a new copy with M() and InvM() swapped.
      ChOfBasisOp swap() const { return ChOfBasisOp(m_InvM, m_M); }
      //! Apply modulus operation such that 0 <= x < TBF().
      /*! The operation is applied to the elements of the
          translation vectors of M() and InvM(). The vectors are
          modified in place.
       */
      void modPositiveInPlace() {
        m_M.modPositive();
        m_InvM.modPositive();
      }
      //! Apply modulus operation such that -TBF()/2+1 < x <= TBF()/2.
      /*! The operation is applied to the elements of the
          translation vectors of M() and InvM(). The vectors are
          modified in place.
       */
      void modShortInPlace() {
        m_M.modShort();
        m_InvM.modShort();
      }
      //! M().Rpart() * R * InvM().Rpart(), for R with rotation base factor 1.
      /*! The rotation base factor of the result is 1.
       */
      RotMx operator()(const RotMx& R) const;
      //! M() * RT * InvM(), for RT with rotation base factor 1.
      /*! Similar to apply(RT), but faster. The translation base factor
          of the result is equal to the translation base factor of RT.
       */
      RTMx operator()(const RTMx& RT) const;
      //! M() * RT * InvM(), for general RT.
      /*! Similar to opertor()(). RT may have any rotation base factor
          or translation base factor. The base factors of the result
          are made as small as possible.
          <p>
          See also: RTMx::multiply, RTMx::cancel()
       */
      RTMx apply(const RTMx& RT) const;
      //! Transform unit cell parameters.
      /*! Equivalent to uc.ChangeBasis(InvM().Rpart()).
          See cctbx::uctbx::UnitCell::ChangeBasis.
          <p>
          To transform in the other direction, use swap() followed by
          apply().
       */
      uctbx::UnitCell apply(const uctbx::UnitCell& uc) const {
        return uc.ChangeBasis(InvM().Rpart());
      }
      //! M() * (RotMx(SignI)|T) * InvM()
      TrVec operator()(const TrVec& T, int SignI) const;
      //! Transform fractional coordinates: M() * X
      template <class FloatType>
      fractional<FloatType>
      operator()(const fractional<FloatType>& X) const { return m_M * X; }
      //! M() = CBOp.M() * M(); InvM() = InvM() * CBOp.InvM();
      void update(const ChOfBasisOp& CBOp) {
        m_M = (CBOp.M() * m_M).newBaseFactors(CBOp.M());
        m_InvM = (m_InvM * CBOp.InvM()).newBaseFactors(CBOp.InvM());
      }
      //! M() = (I|Shift) * M(); InvM() = InvM() * (I|-Shift);
      void update(const TrVec& Shift) {
        // (I|S)*(R|T) = (R|T+S)
        m_M = RTMx(m_M.Rpart(), m_M.Tpart() + Shift);
        // (R|T)*(I|-S) = (R|T-R*S)
        m_InvM = RTMx(
          m_InvM.Rpart(), m_InvM.Tpart()
            - (m_InvM.Rpart() * Shift).newBaseFactor(m_InvM.Tpart().BF()));
      }
    private:
      RTMx m_M;
      RTMx m_InvM;
  };

  //! Multiplication of change-of-basis operators.
  inline ChOfBasisOp operator*(const ChOfBasisOp& lhs,
                               const ChOfBasisOp& rhs) {
    return ChOfBasisOp(
      (lhs.M() * rhs.M()).newBaseFactors(lhs.M()),
      (rhs.InvM() * lhs.InvM()).newBaseFactors(lhs.InvM()));
  }

}} // namespace cctbx::sgtbx

#endif // CCTBX_SGTBX_CHANGE_BASIS_H
