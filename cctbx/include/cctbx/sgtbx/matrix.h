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

#ifndef CCTBX_SGTBX_MATRIX_H
#define CCTBX_SGTBX_MATRIX_H

#include <iostream>
#include <string>
#include <cctbx/sgtbx/basic.h>
#include <cctbx/sgtbx/utils.h>
#include <cctbx/sgtbx/math.h>

namespace cctbx { namespace sgtbx {

#ifndef DOXYGEN_SHOULD_SKIP_THIS
  class TrVec : public af::int3 {
    private:
      int m_BF;
    public:
      typedef af::int3 base_type;
      explicit TrVec(int TranslationBaseFactor = STBF)
        : m_BF(TranslationBaseFactor) { for(int i=0;i<3;i++) elems[i] = 0; }
      explicit TrVec(const af::int3& v, int TranslationBaseFactor = STBF)
        : base_type(v), m_BF(TranslationBaseFactor) {}
      TrVec(const int& v0, const int& v1, const int& v2,
                   int TranslationBaseFactor = STBF)
        : m_BF(TranslationBaseFactor) {
        elems[0] = v0; elems[1] = v1; elems[2] = v2;
      }
      const int& BF() const { return m_BF; }
      int& BF() { return m_BF; }
      bool isValid() const { return m_BF != 0; }
      TrVec Null() const {
        return TrVec(BF());
      }
      bool isNull() const {
        return *this == Null();
      }
      TrVec newBaseFactor(int NewBF) const;
      TrVec scale(int factor) const {
        if (factor == 1) return *this;
        TrVec result(m_BF * factor);
        for(int i=0;i<3;i++) result[i] = elems[i] * factor;
        return result;
      }
      TrVec cancel() const;
      TrVec modPositive() const {
        TrVec result(m_BF);
        for(int i=0;i<3;i++) result[i] = sgtbx::modPositive(elems[i], m_BF);
        return result;
      }
      TrVec modShort() const {
        TrVec result(m_BF);
        for(int i=0;i<3;i++) result[i] = sgtbx::modShort(elems[i], m_BF);
        return result;
      }
      friend bool operator==(const TrVec& lhs, const TrVec& rhs) {
        if (lhs.m_BF != rhs.m_BF) return false;
        for(int i=0;i<3;i++) if (lhs.elems[i] != rhs.elems[i]) return false;
        return true;
      }
      TrVec operator-() const {
        TrVec result(m_BF);
        for(int i=0;i<3;i++) result[i] = -elems[i];
        return result;
      }
      friend TrVec operator+(const TrVec& lhs, const TrVec& rhs) {
        cctbx_assert(lhs.m_BF == rhs.m_BF);
        TrVec result(lhs.m_BF);
        for(int i=0;i<3;i++) result[i] = lhs[i] + rhs[i];
        return result;
      }
      TrVec plus(const TrVec& rhs) const;
      TrVec minus(const TrVec& rhs) const;
      TrVec& operator+=(const TrVec& rhs) {
        cctbx_assert(m_BF == rhs.m_BF);
        for(int i=0;i<3;i++) elems[i] += rhs.elems[i];
        return *this;
      }
      friend TrVec operator-(const TrVec& lhs, const TrVec& rhs) {
        cctbx_assert(lhs.m_BF == rhs.m_BF);
        TrVec result(lhs.m_BF);
        for(int i=0;i<3;i++) result[i] = lhs[i] - rhs[i];
        return result;
      }
      friend TrVec operator*(const TrVec& lhs, const int& rhs) {
        TrVec result(lhs.m_BF);
        for(int i=0;i<3;i++) result[i] = lhs[i] * rhs;
        return result;
      }
      friend TrVec operator/(const TrVec& lhs, int rhs);
      friend TrVec operator*(const int& lhs, const TrVec& rhs) {
        return rhs * lhs;
      }
      friend std::ostream& operator<<(std::ostream& os, const TrVec& T);
      template <class FloatType>
      af::tiny<FloatType, 3> as_array(const FloatType&) const {
        af::tiny<FloatType, 3> result;
        for(int i=0;i<3;i++) {
          result[i] = FloatType(elems[i]) / FloatType(m_BF);
        }
        return result;
      }
  };

  // Constructor for initialization of constants.
  class TrVec12 : public TrVec {
    public:
      TrVec12(int v0, int v1, int v2)
        : TrVec(TrVec(v0, v1, v2) * (STBF / 12)) {}
  };
#endif // DOXYGEN_SHOULD_SKIP_THIS

  //! Helper class for passing information about rotation matrices.
  /*! Determing the sense of rotation requires the rotation type
      and Eigenvector. Therefore it is most efficient to group
      all these properties so that the intermediate results
      can also be used.
   */
  class RotMxInfo {
    public:
      //! For internal use only.
      RotMxInfo() : m_Rtype(0), m_SenseOfRotation(0) {
        for(int i=0;i<3;i++) m_EV[i] = 0;
      }
      //! Rotation-part type (1, 2, 3, 4, 6, -1, -2=m, -3, -4, -6)
      const int& Rtype() const { return m_Rtype; }
      //! Axis direction (Eigenvector) of the corresponding proper rotation.
      /*! Only defined if abs(Rtype) > 1.<br>
          For Rtype > 0, the proper rotation is defined as R.<br>
          For Rtype < 0, the proper rotation is defined as -R.
       */
      const af::int3& EV() const { return m_EV; }
      //! Sense of rotation with respect to the axis direction.
      /*! Only defined if abs(Rtype) > 1.
       */
      const int& SenseOfRotation() const { return m_SenseOfRotation; }
    private:
#ifndef DOXYGEN_SHOULD_SKIP_THIS
      friend class RotMx;
#endif // DOXYGEN_SHOULD_SKIP_THIS
      int  m_Rtype;
      af::int3 m_EV;
      int  m_SenseOfRotation;
  };

#ifndef DOXYGEN_SHOULD_SKIP_THIS
  class RotMx : public af::int9 {
    private:
      int m_BF;
    public:
      typedef af::int9 base_type;
      explicit RotMx(int RotationBaseFactor = 1, int Diagonal = 1)
        : m_BF(RotationBaseFactor) {
        for(int i=0;i<9;i++) elems[i] = 0;
        if (Diagonal * m_BF) {
          for(int i=0;i<3;i++) elems[i * 4] = Diagonal * m_BF;
        }
      }
      explicit RotMx(const af::int9& m, int RotationBaseFactor = 1)
        : base_type(m), m_BF(RotationBaseFactor) {}
      RotMx(int m00, int m01, int m02,
                   int m10, int m11, int m12,
                   int m20, int m21, int m22,
                   const int RotationBaseFactor = 1)
        : m_BF(RotationBaseFactor) {
        elems[0] = m00; elems[1] = m01; elems[2] = m02;
        elems[3] = m10; elems[4] = m11; elems[5] = m12;
        elems[6] = m20; elems[7] = m21; elems[8] = m22;
      }
      const int& BF() const { return m_BF; }
      int& BF() { return m_BF; }
      bool isValid() const { return m_BF != 0; }
      RotMx Unit() const {
        return RotMx(BF());
      }
      bool isUnit() const {
        return *this == Unit();
      }
      RotMx minusUnit() const {
        RotMx result(*this);
        for (int i = 0; i < 9; i += 4) result[i] -= BF();
        return result;
      }
      RotMx newBaseFactor(int NewBF) const;
      RotMx scale(int factor) const {
        if (factor == 1) return *this;
        RotMx result(m_BF * factor);
        for(int i=0;i<9;i++) result[i] = elems[i] * factor;
        return result;
      }
      RotMx cancel() const;
      void setColumn(int c, const af::int3& v) {
        for(int i=0;i<3;i++) elems[i * 3 + c] = v[i];
      }
      af::int3 getColumn(int c) const {
        af::int3 result;
        for(int i=0;i<3;i++) result[i] = elems[i * 3 + c];
        return result;
      }
      void swapColumns(int c1, int c2) {
        for(int i=0;i<3;i++) {
          int e = elems[i * 3 + c1];
          elems[i * 3 + c1] = elems[i * 3 + c2];
          elems[i * 3 + c2] = e;
        }
      }
      int& operator()(int r, int c) {
        return elems[r * 3 + c];
      }
      const int& operator()(int r, int c) const {
        return elems[r * 3 + c];
      }
      int det() const {
        return MatrixLite::Determinant(*this);
      }
      int trace() const {
        return elems[0] + elems[4] + elems[8];
      }
      RotMx CoFactorMxTp() const;
      RotMx inverse() const;
      RotMx inverse_with_cancel() const;
      RotMx divide(int rhs) const;
      RotMx operator-() const {
        RotMx result(m_BF);
        for(int i=0;i<9;i++) result[i] = -elems[i];
        return result;
      }
      friend RotMx operator-(const RotMx& lhs, const RotMx& rhs) {
        cctbx_assert(lhs.m_BF == rhs.m_BF);
        RotMx result(lhs.m_BF);
        for(int i=0;i<9;i++) result[i] = lhs[i] - rhs[i];
        return result;
      }
      friend RotMx operator+(const RotMx& lhs, const RotMx& rhs) {
        cctbx_assert(lhs.m_BF == rhs.m_BF);
        RotMx result(lhs.m_BF);
        for(int i=0;i<9;i++) result[i] = lhs[i] + rhs[i];
        return result;
      }
      RotMx& operator+=(const RotMx& rhs) {
        cctbx_assert(m_BF == rhs.m_BF);
        for(int i=0;i<9;i++) elems[i] += rhs.elems[i];
        return *this;
      }
      friend bool operator==(const RotMx& lhs, const RotMx& rhs) {
        if (lhs.m_BF != rhs.m_BF) return false;
        for(int i=0;i<9;i++) if (lhs.elems[i] != rhs.elems[i]) return false;
        return true;
      }
      friend RotMx operator*(const RotMx& lhs, int rhs);
      friend RotMx operator*(const int lhs, RotMx& rhs) {
        return rhs * lhs;
      }
      RotMx& operator*=(int rhs) {
        for(int i=0;i<9;i++) elems[i] *= rhs;
        return *this;
      }
      friend RotMx operator*(const RotMx& lhs, const RotMx& rhs);
      friend RotMx operator/(const RotMx& lhs, int rhs);
      friend af::int3 operator*(const RotMx& lhs, const af::int3& rhs);
      friend TrVec operator*(const RotMx& lhs, const TrVec& rhs);
      friend TrVec operator*(const TrVec& lhs, const RotMx& rhs);
      RotMx multiply(const RotMx& rhs) const {
        return ((*this) * rhs).cancel();
      }
      TrVec multiply(const TrVec& rhs) const {
        return ((*this) * rhs).cancel();
      }
      friend std::ostream& operator<<(std::ostream& os, const RotMx& R);
      template <class FloatType>
      af::tiny<FloatType, 9> as_array(const FloatType&) const {
        af::tiny<FloatType, 9> result;
        for(int i=0;i<9;i++) {
          result[i] = FloatType(elems[i]) / FloatType(m_BF);
        }
        return result;
      }
      int getRtype() const;
      RotMxInfo getInfo() const;
      RotMx accumulate(int Rtype = 0) const;
  };

  int OrderOfRtype(int Rtype);

#endif // DOXYGEN_SHOULD_SKIP_THIS

  //! Result type for RTMx::analyzeTpart.
  /*! See RTMx::analyzeTpart.
   */
  class TranslationComponents {
    public:
      //! Default constructor. Some data members are not initialized!
      TranslationComponents() {}
      //! Constructor. For internal use only.
      TranslationComponents(const TrVec& IP,
                            const TrVec& LP,
                            const TrVec& OS)
        : m_IP(IP), m_LP(LP), m_OS(OS) {}
      //! Intrinsic (srew or glide) part.
      /*! See RTMx::getIntrinsicPart.
       */
      const TrVec& IntrinsicPart() const { return m_IP; }
      //! Location part.
      /*! See RTMx::getLocationPart.
       */
      const TrVec& LocationPart() const { return m_LP; }
      //! Origin shift.
      /*! See RTMx::getOriginShift.
       */
      const TrVec& OriginShift() const { return m_OS; }
    private:
      TrVec m_IP;
      TrVec m_LP;
      TrVec m_OS;
  };

  //! Rotation-Translation matrix.
  /*! Rpart() is the (3x3) rotation part or rotation matrix and Tpart()
      the (3x1) translation part or translation vector of the symmetry
      operation.  For efficiency, the elements of the matrix and the
      vector are stored as integer numbers and base factors. The actual
      value of an element is obtained by dividing the integer number by
      the corresponding base factor:<br>
      Rpart()[i] / RBF(), Tpart()[i] / TBF().
      <p>
      Python only:
      <dl>
      <dt>sgtbx.RTMx_from_tuple((R, T), RBF, TBF)
        <dd>R is a list or tuple with the 9 elements of the rotation part.
            T is a list or tuple with the 3 elements of the translation part.
            RBF and TBF are the rotation and translation base factors.
      <dt>sgtbx.RTMx().as_tuple()
        <dd>Return a tuple of two tuples: (R, T). Both the elements of
            R and T are Python floating point numbers.
      <dt>sgtbx.RTMx().as_tuple(RBF, TBF)
        <dd>Return a tuple of two tuples: (R, T). Both the elements of
            R and T are Python integers. The rotation and translation
            base factors are RBF and TBF, respectively. An exception
            is raised if the requested base factors are not compatible
            with the RTMx instance.
      </dl>
   */
  class RTMx {
    public:
      //! Initialize identity rotation matrix and zero translation vector.
      /*! The base factors used are RBF and TBF for the rotation part
          and the translation part, respectively.
       */
      explicit RTMx(int RBF = 1, int TBF = STBF) : m_R(RBF), m_T(TBF) {}
      //! Initialize with the given rotation matrix and translation vector.
      RTMx(const RotMx& rmx, const TrVec& trv) : m_R(rmx), m_T(trv) {}
      //! Initialize with the given rotation matrix.
      /*! The translation part is initialized as the zero translation
          with the supplied base factor TBF.
       */
      explicit RTMx(const RotMx& rmx, int TBF = STBF)
        : m_R(rmx), m_T(TBF) {}
      //! Initialize with the given translation vector.
      /*! The rotation part is initialized as the identity matrix
          with the supplied base factor RBF.
       */
      explicit RTMx(const TrVec& trv, int RBF = 1)
        : m_R(RBF), m_T(trv) {}
      //! Initialize by parsing a symbolic expression, e.g. "x+1/2,y,z".
      /*! Parsing will stop at any of the characters listed in StopChars.<br>
          The base factors of the new RTMx are RBF and TBF.<br>
          If an error occurs, an exception is thrown. parse_string can
          be investigated to locate the input character that triggered
          the error.
       */
      RTMx(parse_string& StrXYZ,
           const char* StopChars = "",
           int RBF = 1, int TBF = STBF);
      //! Initialize by parsing a symbolic expression, e.g. "x+1/2,y,z".
      /*! Identical to the constructor that takes a parse_string
          as the first argument. However, if an exception is thrown
          there is no way to locate the input character that triggered
          the error.
       */
      RTMx(const std::string& StrXYZ,
           const char* StopChars = "",
           int RBF = 1, int TBF = STBF);
      //! Access to rotation part.
      const RotMx& Rpart() const { return m_R; };
      //! Access to rotation part.
      RotMx& Rpart() { return m_R; };
      //! Access to translation part.
      const TrVec& Tpart() const { return m_T; };
      //! Access to translation part.
      TrVec& Tpart() { return m_T; };
      //! Convert to a symbolic expression, e.g. "x+1/2,y,z".
      /*! Constants can be formatted as fractional or decimal numbers.<br>
          E.g. "x+1/2,y,z" or "x+0.5,y,z".<br>
          The translation component can appear last or first.<br>
          E.g. "x+1/2,y,z" or "1/2+x,y,z".<br>
          LettersXYZ must contain three characters that are used to
          represent x, y, and z, respectively. Typically one only
          uses either LettersXYZ = "xyz" or LettersXYZ = "XYZ".
          Other letters could be used, but the resulting symbolic
          expression cannot be translated back with the constructors
          above.<br>
          Separator is inserted after the first and the second term
          of the symbolic expression. Typical strings used are
          Separator = "," and Separator = ", ".
       */
      std::string as_xyz(bool Decimal = false, bool TrFirst = false,
                         const char* LettersXYZ = "xyz",
                         const char* Separator = ",") const;
      //! Export matrix elements to a plain array of integers.
      /*! The 9 elements of the rotation part are copied to
          elements 0..8 of the plain array.
          The 3 elements of the translation part are copied to
          elements 9..11 of the plain array.
          Use RBF() and TBF() to access the base factors
          for conversion to floating point representations.
          <p>
          See also: as_array()
       */
      af::tiny<int, 12> as_int_array() const;
      //! Export matrix elements to a floating point array.
      /*! The 9 elements of the rotation part are divided by RBF()
          and copied to elements 0..8 of the array.
          The 3 elements of the translation part are divided by TBF()
          and copied to elements 9..11 of the array.
          <p>
          See also: as_int_array()
       */
      template <class FloatType>
      af::tiny<FloatType, 12> as_array(const FloatType&) const {
        af::tiny<FloatType, 12> result;
        int i;
        for(i=0;i<9;i++) {
          result[i    ] = FloatType(Rpart()[i]) / FloatType(RBF());
        }
        for(i=0;i<3;i++) {
          result[i + 9] = FloatType(Tpart()[i]) / FloatType(TBF());
        }
        return result;
      }
      //! Test if the matrix is valid.
      /*! A RTMx is valid only if both the rotation base factor and the
          translation base factor are not zero.
       */
      bool isValid() const {
        return m_R.isValid() && m_T.isValid();
      }
      //! Rotation base factor.
      /*! This is the factor by which the elements of the rotation
          part are multiplied.
       */
      int RBF() const { return m_R.BF(); }
      //! Translation base factor.
      /*! This is the factor by which the elements of the translation
          part are multiplied.
       */
      int TBF() const { return m_T.BF(); }
      //! Return a new unit matrix.
      /*! The new matrix inherits the rotation and translation base factors.
       */
      RTMx Unit() const {
        return RTMx(RBF(), TBF());
      }
      //! Test if the matrix is the unit matrix.
      bool isUnit() const {
        return *this == Unit();
      }
      //! Return a new copy with the base factors RBF and TBF.
      /*! An exception is thrown if the elements cannot be scaled to
          the new base factors.<br>
          RBF or TBF == 0 indicates that the corresponding old base
          factor is retained.
       */
      RTMx newBaseFactors(int RBF, int TBF = 0) const;
      //! Return a new copy of the matrix, but with the base factors of M.
      /*! An exception is thrown if the elements cannot be scaled to
          the new base factors.
       */
      RTMx newBaseFactors(const RTMx& M) const {
        return newBaseFactors(M.Rpart().BF(), M.Tpart().BF());
      }
      //! Multiply the elements and the base factors by the given factors.
      /*! if (factorT == 0) factorT = factorR;
       */
      RTMx scale(int factorR, int factorT = 0) const;
      //! Reduce base factors by cancellation.
      /*! The rotation base factor and the elements of the rotation
          part are divided by their greatest common denominator.
          The same procedure is applied to the translation part.
       */
      RTMx cancel() const;
      //! Apply modulus operation such that 0 <= x < TBF().
      /*! The operation is applied to the elements of the
          translation vector. The vector is modified in place.
       */
      void modPositiveInPlace() { m_T = m_T.modPositive(); }
      //! Apply modulus operation such that 0 <= x < TBF().
      /*! The operation is applied to the elements of the
          translation vector. A new instance of RTMx is created.
       */
      RTMx modPositive() const { return RTMx(m_R, m_T.modPositive()); }
      //! Apply modulus operation such that -TBF()/2+1 < x <= TBF()/2.
      /*! The operation is applied to the elements of the
          translation vector. The vector is modified in place.
       */
      void modShortInPlace() { m_T = m_T.modShort(); }
      //! Apply modulus operation such that -TBF()/2+1 < x <= TBF()/2.
      /*! The operation is applied to the elements of the
          translation vector. A new instance of RTMx is created.
       */
      RTMx modShort() const { return RTMx(m_R, m_T.modShort()); }
      //! Compute information about the rotation part.
      /*! See information about class RotMxInfo.
       */
      RotMxInfo getRotMxInfo() const { return m_R.getInfo(); }
      /*! \brief Test if the vector v is perpendicular to the axis
          direction of the rotation part.
       */
      bool isPerpendicular(const af::int3& v) const;
      //! Compute the intrinsic (screw or glide) part of the translation part.
      /*! Let N be the rotation-part type of Rpart().
          Following the procedure given by Fischer & Koch
          (International Tables for Crystallography, Volume A, 1983,
          Chapter 11), the first step in the analysis of the
          translation part Tpart() is the decomposition into the
          intrinsic part wi and the location part wl. For this,
          (Rpart()|Tpart())^n = (I|t) has to be computed, where the
          rotational order n = abs(N), except for N = -1
          and N = -3. For those two cases, n = -2*N. The
          intrinsic part is obtained as wi = (1/n)*t.<br>
          See also: analyzeTpart()
       */
      TrVec getIntrinsicPart() const;
      //! Compute the location part given the intrinsic part wi.
      /*! wi is the result of getIntrinsicPart(). The location part
          is simply the difference wl = Tpart() - wi.<br>
          See also: analyzeTpart()
       */
      TrVec getLocationPart(const TrVec& wi) const;
      //! Compute the origin shift given the location part wl.
      /*! wl is the result of getLocationPart(). The origin shift
          is obtained by solving the equation (Rpart()|wl)*x = x
          for x (see
          <A HREF="http://journals.iucr.org/a/issues/1999/02/02/au0146/"
          ><I>Acta Cryst.</I> 1999, <B>A55</B>:383-395</A>).
          For rotation-part types N > 1 and N = -2, the combination
          with the axis direction of Rpart() is a convenient
          description of all fixed points. For the other rotation-part
          types, the fixed point is unique.<br>
          See also: analyzeTpart()
       */
      TrVec getOriginShift(const TrVec& wl) const;
      //! Analyse the translation part.
      /*! The results of getIntrinsicPart(), getLocationPart(),
          and getOriginShift() are collected in class
          TranslationComponents.
       */
      TranslationComponents analyzeTpart() const;
      //! Efficient computation of (-I|InvT) * (R|T).
      /*! I is the identidy matrix. InvT is the translation part
          of a centre of inversion.
       */
      RTMx pre_multiply_InvT(const TrVec& InvT) {
        return RTMx(-m_R, -m_T + InvT);
      }
      //! Compute the inverse matrix.
      /*! An exception is thrown if the matrix cannot be inverted
          or if the result cannot be scaled to the rotation base factor
          or the translation base factor.
       */
      RTMx inverse() const;
      //! Compute the inverse matrix.
      /*! An exception is thrown if the matrix cannot be inverted.
       */
      RTMx inverse_with_cancel() const;
      /*! \brief Similar to /=, but multiply base factors instead
          of dividing elements.
       */
      void pseudo_divide(int rhs) {
        m_R.BF() *= rhs;
        m_T.BF() *= rhs;
      }
      //! Unary minus.
      RTMx operator-() const { return RTMx(-m_R, -m_T); }
      //! Addition operator.
      friend RTMx operator+(const RTMx& lhs, const RTMx& rhs);
      //! += operator.
      RTMx& operator+=(const RTMx& rhs) {
        m_R += rhs.m_R;
        m_T += rhs.m_T;
        return *this;
      }
      //! Multiplication of homogeneous RTMx with RBF()=1.
      /*! The rotation base factor of both matrices must be equal to 1.
          The translation base factors of the two matrices must be equal.
          <p>
          The base factors of the result are equal to the base factors
          of the operands.
          <p>
          operator*() is faster than multiply().
       */
      friend RTMx operator*(const RTMx& lhs, const RTMx& rhs);
      //! Multiplication with cancellation for general RTMx.
      /*! Similar to opertor*(). However, the operands may have any
          rotation base factor or translation base factor.
          <p>
          The base factors of the result are made as small as possible.
          <p>
          See also: RTMx::cancel()
       */
      RTMx multiply(const RTMx& rhs) const;
      //! Addition of translation vector to translation part.
      friend RTMx operator+(const RTMx& lhs, const TrVec& rhs);
      //! Test for equality.
      /*! Essentially this is a byte-wise comparison. Note that the
          translation vectors are not modified as part of this
          equality test.
       */
      friend bool operator==(const RTMx& lhs, const RTMx& rhs) {
        return (lhs.m_R == rhs.m_R) && (lhs.m_T == rhs.m_T);
      }
      //! Negation of test for equality.
      friend bool operator!=(const RTMx& lhs, const RTMx& rhs) {
        return !(lhs == rhs);
      }
      /*! \brief Refine gridding such that each grid point is
          mapped onto another grid point by the symmetry operation.
       */
      template <typename GridTupleType>
      GridTupleType refine_gridding(const GridTupleType& grid) {
        GridTupleType result;
        for(int ir=0;ir<3;ir++) {
          result[ir] = lcm(grid[ir], norm_denominator(m_T[ir], m_T.BF()));
          for(int ic=0;ic<3;ic++) {
            result[ir] = lcm(
              result[ir], norm_denominator(m_R(ir, ic), grid[ic]));
          }
        }
        return result;
      }
      //! Output operator, mainly for debugging.
      /*! Use the member function as_xyz() for high-level formatting.
       */
      friend std::ostream& operator<<(std::ostream& os, const RTMx& M);
    private:
      RotMx m_R;
      TrVec m_T;
  };

#ifndef DOXYGEN_SHOULD_SKIP_THIS

  inline af::double3
  operator+(const af::double3& lhs, const TrVec& rhs) {
    af::double3 result;
    for(int i=0;i<3;i++) result[i] = lhs[i] + double(rhs[i]) / rhs.BF();
    return result;
  }

  inline af::double3
  operator*(const RTMx& lhs, const af::double3& rhs) {
    const RotMx& R = lhs.Rpart();
    const TrVec& T = lhs.Tpart();
    const double RBF = R.BF();
    const double TBF = T.BF();
    af::double3 result;
    for(int i=0;i<3;i++) {
      result[i] = (  R[i * 3 + 0] * rhs[0]
                   + R[i * 3 + 1] * rhs[1]
                   + R[i * 3 + 2] * rhs[2]) / RBF + T[i] / TBF;
    }
    return result;
  }

#endif // DOXYGEN_SHOULD_SKIP_THIS

}} // namespace cctbx::sgtbx

#endif // CCTBX_SGTBX_MATRIX_H
