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
#include <cctbx/sgtbx/matrix.h>

namespace sgtbx {

  //! Change-of-basis (transformation) operator.
  /*! For ease of use, a change-of-basis matrix M() and its inverse InvM()
      are grouped by this class.
   */
  class ChOfBasisOp {
    private:
      RTMx Mx;
      RTMx InvMx;
    public:
      //! Initialize the change-of-basis operator with M and InvM.
      /*! The input matrices are NOT checked for consistency.
       */
      inline ChOfBasisOp(const RTMx& M, const RTMx& InvM)
        : Mx(M), InvMx(InvM) { }
      //! Initialize the change-of-basis operator with M.
      /*! The inverse matrix InvM is computed by inverting M.
          An exception is thrown if M is not invertible.
       */
      explicit inline ChOfBasisOp(const RTMx& M)
        : Mx(M), InvMx(M.inverse()) { }
      //! Initialize the change-of-basis operator with unit matrices.
      /*! The unit matrices are initialized with the rotation base factor
          RBF and the translation base factor TBF.
       */
      explicit inline ChOfBasisOp(int RBF = CRBF, int TBF = CTBF)
        : Mx(RBF, TBF), InvMx(RBF, TBF) { }
      //! Test if the change-of-basis operator is valid.
      /*! A ChOfBasisOp is valid only if the rotation base factor and the
          translation base factor of both M() and InvM() are not zero.
       */
      inline bool isValid() const {
        return Mx.isValid() && InvMx.isValid();
      }
      //! Return a new change-of-basis operator with unit matrices.
      /*! The new matrices inherit the rotation and translation base factors.
       */
      inline ChOfBasisOp Unit() const {
        return ChOfBasisOp(Mx.Unit(), InvMx.Unit());
      }
      //! Test if the change-of-basis operator contains unit matrices.
      inline bool isUnit() const {
        return Mx.isUnit() && InvMx.isUnit();
      }
      //! Return a new copy with the base factors RBF and TBF.
      /*! An exception is thrown if the elements cannot be scaled to
          the new base factors.<br>
          RBF or TBF == 0 indicates that the corresponding old base
          factor is retained.
       */
      inline ChOfBasisOp newBaseFactors(int RBF, int TBF) {
        return ChOfBasisOp(Mx.newBaseFactors(RBF, TBF),
                           InvMx.newBaseFactors(RBF, TBF));
      }
      //! Return a new copy of the operator, but with the base factors of CBOp.
      /*! An exception is thrown if the elements cannot be scaled to
          the new base factors.
       */
      inline ChOfBasisOp newBaseFactors(const ChOfBasisOp& CBOp) {
        return ChOfBasisOp(Mx.newBaseFactors(CBOp.M()),
                           InvMx.newBaseFactors(CBOp.InvM()));
      }
      //! Return the change-of-basis matrix.
      inline const RTMx& M() const { return Mx; }
      //! Return the inverse of the change-of-basis matrix.
      inline const RTMx& InvM() const { return InvMx; }
      //! Return M() for Inv == false, and InvM() for Inv == true.
      inline const RTMx& select(bool Inv) {
        if (Inv) return InvMx;
        return Mx;
      }
      //! Return a new copy with M() and InvM() swapped.
      inline ChOfBasisOp swap() const { return ChOfBasisOp(InvMx, Mx); }
      //! Apply modulus operation such that 0 <= x < TBF().
      /*! The operation is applied to the elements of the
          translation vectors of M() and InvM(). The vectors are
          modified in place.
       */
      inline void modPositiveInPlace() {
        Mx.modPositive();
        InvMx.modPositive();
      }
      //! Apply modulus operation such that -TBF()/2+1 < x <= TBF()/2.
      /*! The operation is applied to the elements of the
          translation vectors of M() and InvM(). The vectors are
          modified in place.
       */
      inline void modShortInPlace() {
        Mx.modShort();
        InvMx.modShort();
      }
      //! M().Rpart() * R * InvM().Rpart(), for R with rotation base factor 1.
      /*! The rotation base factor of the result is 1.
       */
      RotMx operator()(const RotMx& R) const;
      //! M() * RT * InvM(), for RT with rotation base factor 1.
      /*! Similar to apply(), but faster. The translation base factor
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
      RTMx apply(const RTMx& RT);
      //! M() * (RotMx(SignI)|T) * InvM()
      TrVec operator()(const TrVec& T, int SignI) const;
      //! Transform fractional coordinates: M() * X
      template <class T>
      inline fractional<T>
      operator()(const fractional<T>& X) const { return Mx * X; }
      //! M() = CBOp.M() * M(); InvM() = InvM() * CBOp.InvM();
      inline void update(const ChOfBasisOp& CBOp) {
        Mx = (CBOp.M() * Mx).newBaseFactors(CBOp.M());
        InvMx = (InvMx * CBOp.InvM()).newBaseFactors(CBOp.InvM());
      }
      //! M() = (I|Shift) * M(); InvM() = InvM() * (I|-Shift);
      inline void update(const TrVec& Shift) {
        // (I|S)*(R|T) = (R|T+S)
        Mx = RTMx(Mx.Rpart(), Mx.Tpart() + Shift);
        // (R|T)*(I|-S) = (R|T-R*S)
        InvMx = RTMx(
          InvMx.Rpart(), InvMx.Tpart()
            - (InvMx.Rpart() * Shift).newBaseFactor(InvMx.Tpart().BF()));
      }
  };

  inline ChOfBasisOp operator*(const ChOfBasisOp& lhs,
                               const ChOfBasisOp& rhs) {
    return ChOfBasisOp(
      (lhs.M() * rhs.M()).newBaseFactors(lhs.M()),
      (rhs.InvM() * lhs.InvM()).newBaseFactors(lhs.InvM()));
  }

} // namespace sgtbx

#endif // CCTBX_SGTBX_CHANGE_BASIS_H
