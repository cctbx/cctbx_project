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

#include <ctype.h> // cannot use cctype b/o non-conforming compilers
#include <string.h>
#include <stdio.h>
#include <cctbx/sgtbx/utils.h>
#include <cctbx/sgtbx/matrix.h>
#include <cctbx/basic/define_range.h>

namespace sgtbx {

  TrVec TrVec::newBaseFactor(int NewBF) const
  {
    TrVec result(NewBF);
    if (ChangeBaseFactor(elems, BF(), result.elems, NewBF, 3) != 0) {
      throw error_base_factor(
        __FILE__, __LINE__, "out of translation-base-factor range");
    }
    return result;
  }

  TrVec operator/(const TrVec& lhs, int rhs) {
    TrVec result(lhs.BF());
    for(int i=0;i<3;i++) {
      if (lhs[i] % rhs) throw cctbx_internal_error();
      result[i] = lhs[i] / rhs;
    }
    return result;
  }

  RotMx RotMx::newBaseFactor(int NewBF) const
  {
    RotMx result(NewBF);
    if (ChangeBaseFactor(elems, BF(), result.elems, NewBF, 9) != 0) {
      throw error_base_factor(
        __FILE__, __LINE__, "out of rotation-base-factor range");
    }
    return result;
  }

  RTMx RTMx::newBaseFactors(int RBF, int TBF) const
  {
    RTMx result(*this);
    if (RBF) result.m_R = result.m_R.newBaseFactor(RBF);
    if (TBF) result.m_T = result.m_T.newBaseFactor(TBF);
    return result;
  }

  RTMx RTMx::scale(int factorR, int factorT) const
  {
    if (factorT == 0) factorT = factorR;
    return RTMx(m_R.scale(factorR), m_T.scale(factorT));
  }

  Vec3 operator*(const RotMx& lhs, const Vec3& rhs)
  {
    Vec3 result;
    MatrixLite::multiply(lhs.elems, rhs.elems, 3, 3, 1, result.elems);
    return result;
  }

  TrVec operator*(const RotMx& lhs, const TrVec& rhs)
  {
    TrVec result(lhs.BF() * rhs.BF());
    MatrixLite::multiply(lhs.elems, rhs.elems, 3, 3, 1, result.elems);
    return result;
  }

  TrVec operator*(const TrVec& lhs, const RotMx& rhs)
  {
    TrVec result(lhs.BF() * rhs.BF());
    MatrixLite::multiply(lhs.elems, rhs.elems, 1, 3, 3, result.elems);
    return result;
  }

  std::ostream& operator<<(std::ostream& os, const TrVec& T)
  {
    os << "(";
    rangei(3) {
      if (i) os << ",";
      os << T[i];
    }
    os << ")/";
    os << T.BF();
    return os;
  }

  RotMx RotMx::CoFactorMxTp() const
  {
    RotMx result(BF() * BF());
    result[0] =  elems[4] * elems[8] - elems[5] * elems[7];
    result[1] = -elems[1] * elems[8] + elems[2] * elems[7];
    result[2] =  elems[1] * elems[5] - elems[2] * elems[4];
    result[3] = -elems[3] * elems[8] + elems[5] * elems[6];
    result[4] =  elems[0] * elems[8] - elems[2] * elems[6];
    result[5] = -elems[0] * elems[5] + elems[2] * elems[3];
    result[6] =  elems[3] * elems[7] - elems[4] * elems[6];
    result[7] = -elems[0] * elems[7] + elems[1] * elems[6];
    result[8] =  elems[0] * elems[4] - elems[1] * elems[3];
    return result;
  }

  RotMx RotMx::inverse() const
  {
    int detBF3 = det();
    if (detBF3 == 0) throw error("Rotation matrix is not invertible.");
    RotMx C(CoFactorMxTp());
    C.m_BF = BF();
    return (C * (BF() * BF())) / detBF3;
  }

  RotMx RotMx::inverse_with_cancel() const
  {
    boost::rational<int> d(det(), BF()*BF()*BF());
    if (d == 0) throw error("Rotation matrix is not invertible.");
    return (CoFactorMxTp() * d.denominator()).divide(d.numerator());
  }

  RotMx RotMx::divide(int rhs) const
  {
    RotMx result;
    if (rhs < 0) {
      result = -(*this);
      rhs = -rhs;
    }
    else {
      result = *this;
    }
    result.m_BF *= rhs;
    return result.cancel();
  }

  RotMx operator*(const RotMx& lhs, const RotMx& rhs)
  {
    RotMx result(lhs.BF() * rhs.BF());
    MatrixLite::multiply(lhs.elems, rhs.elems, 3, 3, 3, result.elems);
    return result;
  }

  RotMx operator*(const RotMx& lhs, int rhs)
  {
    RotMx result(lhs);
    rangei(9) result[i] *= rhs;
    return result;
  }

  RotMx operator/(const RotMx& lhs, int rhs)
  {
    RotMx result(lhs);
    if (ChangeBaseFactor(lhs.elems, rhs, result.elems, 1, 9) != 0) {
      throw error_base_factor(
        __FILE__, __LINE__, "out of rotation-base-factor range");
    }
    return result;
  }

  RTMx RTMx::inverse() const
  {
    RotMx InvR = m_R.inverse();
    TrVec InvT = (-InvR * m_T).newBaseFactor(m_T.BF());
    return RTMx(InvR, InvT);
  }

  RTMx RTMx::inverse_with_cancel() const
  {
    RotMx InvR = m_R.inverse_with_cancel();
    TrVec InvT = (-InvR).multiply(m_T);
    return RTMx(InvR, InvT);
  }

  std::ostream& operator<<(std::ostream& os, const RotMx& R)
  {
    os << "(";
    rangei(9) {
      if (i) os << ",";
      if (i == 3 || i == 6) os << " ";
      os << R[i];
    }
    os << ")/";
    os << R.BF();
    return os;
  }

  RTMx operator+(const RTMx& lhs, const RTMx& rhs)
  {
    return RTMx(lhs.m_R + rhs.m_R, lhs.m_T + rhs.m_T);
  }

  RTMx operator*(const RTMx& lhs, const RTMx& rhs)
  {
    cctbx_assert(lhs.m_T.BF() == rhs.m_T.BF());
    return RTMx(lhs.m_R * rhs.m_R,
                lhs.m_R * rhs.m_T + lhs.m_T.scale(lhs.m_R.BF()));
  }

  RTMx operator+(const RTMx& lhs, const TrVec& rhs)
  {
    return RTMx(lhs.m_R, lhs.m_T + rhs);
  }

  std::ostream& operator<<(std::ostream& os, const RTMx& M)
  {
    os << M.m_R << ", " << M.m_T;
    return os;
  }

  std::string RTMx::as_xyz(bool Decimal, bool TrFirst,
                           const char* LettersXYZ,
                           const char* Separator) const
  {
    std::string result;
    for (int i = 0; i < 3; i++) {
      std::string R_term;
      for (int j = 0; j < 3; j++) {
        rational R_frac(m_R(i, j), m_R.BF());
        if (R_frac != 0) {
          if (R_frac > 0) {
            if (!R_term.empty()) {
              R_term += "+";
            }
          }
          else {
            R_term += "-";
            R_frac *= -1;
          }
          if (R_frac != 1) {
            R_term += R_frac.format(Decimal) + "*";
          }
          R_term += LettersXYZ[j];
        }
      }
      if (i != 0) result += Separator;
      rational T_frac(m_T[i], m_T.BF());
      if (T_frac == 0) {
        if (R_term.empty()) result += "0";
        else                result += R_term;
      }
      else if (R_term.empty()) {
        result += T_frac.format(Decimal);
      }
      else if (TrFirst) {
        result += T_frac.format(Decimal);
        if (R_term[0] != '-') result += "+";
        result += R_term;
      }
      else {
        result += R_term;
        if (T_frac > 0) result += "+";
        result += T_frac.format(Decimal);
      }
    }
    return result;
  }

  boost::array<int, 12> RTMx::as_int_array() const
  {
    boost::array<int, 12> result;
    int i;
    for(i=0;i<9;i++) result[i    ] = Rpart()[i];
    for(i=0;i<3;i++) result[i + 9] = Tpart()[i];
    return result;
  }

  RTMx::RTMx(parse_string& StrXYZ, const char* StopChars, int RBF, int TBF)
    : m_R(0), m_T(0)
  {
    static const error parse_error("Parse error.");
    RTMx result(RBF, TBF);
    int Row    = 0;
    int Column = -1;
    int Sign   = 1;
    int Mult   = 0;
    double ValR[3];
    int i;
    for(i=0;i<3;i++) ValR[i] = 0.;
    double ValT   = 0.;
    double Value  = 0.;
    const unsigned int P_Add   = 0x01u;
    const unsigned int P_Mult  = 0x02u;
    const unsigned int P_Value = 0x04u;
    const unsigned int P_XYZ   = 0x08u;
    const unsigned int P_Comma = 0x10u;
    unsigned int P_mode = P_Add | P_Value | P_XYZ;
    for (;; StrXYZ.skip())
    {
      if (strchr(StopChars, StrXYZ()) || !isspace(StrXYZ()))
      {
        switch (strchr(StopChars, StrXYZ()) ? '\0' : StrXYZ())
        {
          case '_':
            break;
          case '+': Sign =  1; goto ProcessAdd;
          case '-': Sign = -1;
           ProcessAdd:
            if ((P_mode & P_Add) == 0) throw parse_error;
            if (Column >= 0) ValR[Column] += Value;
            else             ValT         += Value;
            Value = 0.;
            Column = -1;
            Mult = 0;
            P_mode = P_Value | P_XYZ;
            break;
          case '*':
            if ((P_mode & P_Mult) == 0) throw parse_error;
            Mult = 1;
            P_mode = P_Value | P_XYZ;
            break;
          case '/':
          case ':':
            if ((P_mode & P_Mult) == 0) throw parse_error;
            Mult = -1;
            P_mode = P_Value;
            break;
          case '.':
          case '0':
          case '1':
          case '2':
          case '3':
          case '4':
          case '5':
          case '6':
          case '7':
          case '8':
          case '9':
            if ((P_mode & P_Value) == 0) throw parse_error;
            {
            double V;
            int i = 1;
            int n = sscanf(StrXYZ.peek(), "%lf%n", &V, &i);
            StrXYZ.skip(i - 1);
            if (n != 1) throw parse_error;
            if (Sign == -1) { V = -V; Sign = 1; }
            if      (Mult ==  1)
              Value *= V;
            else if (Mult == -1) {
              if      (V     != 0.) Value /= V;
              else if (Value != 0.) throw parse_error;
            }
            else
              Value = V;
            }
            P_mode = P_Comma | P_Add | P_Mult | P_XYZ;
            break;
          case 'X':
          case 'x': Column = 0; goto Process_XYZ;
          case 'Y':
          case 'y': Column = 1; goto Process_XYZ;
          case 'Z':
          case 'z': Column = 2;
           Process_XYZ:
            if ((P_mode & P_XYZ) == 0) throw parse_error;
            if (Value == 0.) { Value = Sign; Sign = 1; }
            P_mode = P_Comma | P_Add | P_Mult;
            break;
          case ',':
          case ';':
            if (Row == 2) throw parse_error;
          case '\0':
            if ((P_mode & P_Comma) == 0) throw parse_error;
            if (Column >= 0) ValR[Column] += Value;
            else             ValT         += Value;
            for(i=0;i<3;i++) {
              if (rationalize(ValR[i], result.m_R(Row, i), RBF) != 0)
                throw error_base_factor(
                  __FILE__, __LINE__, "out of rotation-base-factor range");
            }
            if (rationalize(ValT, result.m_T[Row], TBF) != 0)
              throw error_base_factor(
                __FILE__, __LINE__, "out of translation-base-factor range");
            Row++;
            Column = -1;
            Sign = 1;
            Mult = 0;
            for(i=0;i<3;i++) ValR[i] = 0.;
            ValT = 0.;
            Value = 0.;
            P_mode = P_Add | P_Value | P_XYZ;
            break;
          default:
            throw parse_error;
        }
      }
      if (strchr(StopChars, StrXYZ()))
        break;
    }
    if (Row != 3) throw parse_error;
    m_R = result.m_R;
    m_T = result.m_T;
  }

  RTMx::RTMx(const std::string& StrXYZ, const char* StopChars,
             int RBF, int TBF)
    : m_R(0), m_T(0)
  {
    parse_string StrXYZ_PS(StrXYZ);
    *this = RTMx(StrXYZ_PS, StopChars, RBF, TBF);
  }

  TrVec RTMx::getIntrinsicPart() const
  {
    int Rtype = Rpart().getRtype();
    RotMx CumMx = Rpart().accumulate(Rtype);
    TrVec wi = CumMx * Tpart();
    wi = wi / OrderOfRtype(Rtype);
    return wi;
  }

  TrVec RTMx::getLocationPart(const TrVec& wi) const
  {
    return wi - Tpart();
  }

  TrVec RTMx::getOriginShift(const TrVec& wl) const
  {
    RotMx RmI = Rpart() - Rpart().Unit();
    RotMx P(1);
    (void) iRowEchelonFormT(RmI.elems, 3, 3, P.elems, 3);
    TrVec Pwl = P * wl;
    TrVec sh(0);
    sh.BF() = iREBacksubst(RmI.elems, Pwl.elems, 3, 3, sh.elems, 0);
    cctbx_assert(sh.BF() > 0);
    sh.BF() *= Pwl.BF();
    return sh;
  }

  TranslationComponents RTMx::analyzeTpart() const
  {
    TrVec wi = getIntrinsicPart();
    TrVec wl = getLocationPart(wi);
    TrVec sh = getOriginShift(wl);
    return TranslationComponents(wi, wl, sh);
  }

  TrVec TrVec::cancel() const
  {
    int g = BF();
    int i;
    for(i=0;i<3;i++) g = gcd(g, elems[i]);
    if (g == 0) g = 1;
    Vec3 result;
    for(i=0;i<3;i++) result[i] = elems[i] / g;
    return TrVec(result, BF() / g);
  }

  RotMx RotMx::cancel() const
  {
    int g = BF();
    int i;
    for(i=0;i<9;i++) g = gcd(g, elems[i]);
    if (g == 0) g = 1;
    Mx33 result;
    for(i=0;i<9;i++) result[i] = elems[i] / g;
    return RotMx(result, BF() / g);
  }

  RTMx RTMx::cancel() const
  {
    return RTMx(Rpart().cancel(), Tpart().cancel());
  }

  TrVec TrVec::plus(const TrVec& rhs) const
  {
    TrVec result(lcm(BF(), rhs.BF()));
    int l = result.BF() / BF();
    int r = result.BF() / rhs.BF();
    for(int i=0;i<3;i++) result[i] = elems[i] * l + rhs[i] * r;
    return result.cancel();
  }

  TrVec TrVec::minus(const TrVec& rhs) const
  {
    TrVec result(lcm(BF(), rhs.BF()));
    int l = result.BF() / BF();
    int r = result.BF() / rhs.BF();
    for(int i=0;i<3;i++) result[i] = elems[i] * l - rhs[i] * r;
    return result.cancel();
  }

  RTMx RTMx::multiply(const RTMx& rhs) const
  {
    if (TBF() == rhs.TBF()) {
      return RTMx(
        (Rpart() * rhs.Rpart()),
        (Rpart() * rhs.Tpart() + Tpart().scale(RBF()))).cancel();
    }
    else {
      int f = lcm(TBF(), rhs.TBF());
      int l = f / TBF() * RBF();
      int r = f / rhs.TBF();
      return RTMx(
        (Rpart() * rhs.Rpart()),
        (Rpart() * rhs.Tpart().scale(r) + Tpart().scale(l))).cancel();
    }
  }

} // namespace sgtbx
