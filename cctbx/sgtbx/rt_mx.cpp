/* Copyright (c) 2001-2002 The Regents of the University of California
   through E.O. Lawrence Berkeley National Laboratory, subject to
   approval by the U.S. Department of Energy.
   See files COPYRIGHT.txt and LICENSE.txt for further details.

   Revision history:
     2002 Sep: Refactored parts of sgtbx/matrix.cpp (rwgk)
     2001 Jul: Merged from CVS branch sgtbx_special_pos (rwgk)
     2001 May: merged from CVS branch sgtbx_type (R.W. Grosse-Kunstleve)
     2001 Apr: SourceForge release (R.W. Grosse-Kunstleve)
 */

#include <cctbx/sgtbx/rt_mx.h>
#include <cctbx/sgtbx/row_echelon.h>
#include <cctbx/rational.h>
#include <ctype.h>
#include <string.h>

namespace cctbx { namespace sgtbx {

  rt_mx rt_mx::new_denominators(int r_den, int t_den) const
  {
    rt_mx result(*this);
    if (r_den) result.r_ = result.r_.new_denominator(r_den);
    if (t_den) result.t_ = result.t_.new_denominator(t_den);
    return result;
  }

  rt_mx rt_mx::scale(int factor_r, int factor_t) const
  {
    if (factor_t == 0) factor_t = factor_r;
    return rt_mx(r_.scale(factor_r), t_.scale(factor_t));
  }

  rt_mx rt_mx::inverse() const
  {
    rot_mx r_inv = r_.inverse();
    tr_vec t_inv = (-r_inv * t_).new_denominator(t_.den());
    return rt_mx(r_inv, t_inv);
  }

  rt_mx rt_mx::inverse_cancel() const
  {
    rot_mx r_inv = r_.inverse_cancel();
    tr_vec t_inv = -r_inv.multiply(t_);
    return rt_mx(r_inv, t_inv);
  }

  rt_mx operator+(rt_mx const& lhs, rt_mx const& rhs)
  {
    return rt_mx(lhs.r_ + rhs.r_, lhs.t_ + rhs.t_);
  }

  rt_mx operator*(rt_mx const& lhs, rt_mx const& rhs)
  {
    CCTBX_ASSERT(lhs.t_.den() == rhs.t_.den());
    return rt_mx(lhs.r_ * rhs.r_,
                lhs.r_ * rhs.t_ + lhs.t_.scale(lhs.r_.den()));
  }

  rt_mx operator+(rt_mx const& lhs, tr_vec const& rhs)
  {
    return rt_mx(lhs.r_, lhs.t_ + rhs);
  }

  std::string rt_mx::as_xyz(bool decimal, bool t_first,
                            const char* letters_xyz,
                            const char* separator) const
  {
    CCTBX_ASSERT(strlen(letters_xyz) == 3);
    std::string result;
    for (int i = 0; i < 3; i++) {
      std::string R_term;
      for (int j = 0; j < 3; j++) {
        boost::rational<int> R_frac(r_(i, j), r_.den());
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
            R_term += format(R_frac, decimal) + "*";
          }
          R_term += letters_xyz[j];
        }
      }
      if (i != 0) result += separator;
      boost::rational<int> T_frac(t_[i], t_.den());
      if (T_frac == 0) {
        if (R_term.empty()) result += "0";
        else                result += R_term;
      }
      else if (R_term.empty()) {
        result += format(T_frac, decimal);
      }
      else if (t_first) {
        result += format(T_frac, decimal);
        if (R_term[0] != '-') result += "+";
        result += R_term;
      }
      else {
        result += R_term;
        if (T_frac > 0) result += "+";
        result += format(T_frac, decimal);
      }
    }
    return result;
  }

  namespace {

    void throw_parse_error() { throw error("Parse error."); }

    int rationalize(double fVal, int& iVal, int den)
    {
      if (den == 0) return -1;
          fVal *= den;
      if (fVal < 0.) iVal = int(fVal - .5);
      else           iVal = int(fVal + .5);
          fVal -= iVal;
          fVal /= den;
      if (fVal < 0.) fVal = -fVal;
      if (fVal > .0001) return -1;
      return 0;
    }

  } // namespace <anonymous>

  rt_mx::rt_mx(parse_string& str_xyz, const char* stop_chars,
               int r_den, int t_den)
    : r_(0), t_(0)
  {
    rt_mx result(r_den, t_den);
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
    for (;; str_xyz.skip())
    {
      if (strchr(stop_chars, str_xyz()) || !isspace(str_xyz()))
      {
        switch (strchr(stop_chars, str_xyz()) ? '\0' : str_xyz())
        {
          case '_':
            break;
          case '+': Sign =  1; goto ProcessAdd;
          case '-': Sign = -1;
           ProcessAdd:
            if ((P_mode & P_Add) == 0) throw_parse_error();
            if (Column >= 0) ValR[Column] += Value;
            else             ValT         += Value;
            Value = 0.;
            Column = -1;
            Mult = 0;
            P_mode = P_Value | P_XYZ;
            break;
          case '*':
            if ((P_mode & P_Mult) == 0) throw_parse_error();
            Mult = 1;
            P_mode = P_Value | P_XYZ;
            break;
          case '/':
          case ':':
            if ((P_mode & P_Mult) == 0) throw_parse_error();
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
            if ((P_mode & P_Value) == 0) throw_parse_error();
            {
            double V;
            int i = 1;
            int n = sscanf(str_xyz.peek(), "%lf%n", &V, &i);
            str_xyz.skip(i - 1);
            if (n != 1) throw_parse_error();
            if (Sign == -1) { V = -V; Sign = 1; }
            if      (Mult ==  1)
              Value *= V;
            else if (Mult == -1) {
              if      (V     != 0.) Value /= V;
              else if (Value != 0.) throw_parse_error();
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
            if ((P_mode & P_XYZ) == 0) throw_parse_error();
            if (Value == 0.) { Value = Sign; Sign = 1; }
            P_mode = P_Comma | P_Add | P_Mult;
            break;
          case ',':
          case ';':
            if (Row == 2) throw_parse_error();
          case '\0':
            if ((P_mode & P_Comma) == 0) throw_parse_error();
            if (Column >= 0) ValR[Column] += Value;
            else             ValT         += Value;
            for(i=0;i<3;i++) {
              if (rationalize(ValR[i], result.r_(Row, i), r_den) != 0) {
                throw_unsuitable_rot_mx(__FILE__, __LINE__);
              }
            }
            if (rationalize(ValT, result.t_[Row], t_den) != 0) {
              throw_unsuitable_tr_vec(__FILE__, __LINE__);
            }
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
            throw_parse_error();
        }
      }
      if (strchr(stop_chars, str_xyz()))
        break;
    }
    if (Row != 3) throw_parse_error();
    r_ = result.r_;
    t_ = result.t_;
  }

  rt_mx::rt_mx(std::string const& str_xyz, const char* stop_chars,
             int r_den, int t_den)
  : r_(0), t_(0)
  {
    parse_string parse_str_xyz(str_xyz);
    *this = rt_mx(parse_str_xyz, stop_chars, r_den, t_den);
  }

  af::tiny<int, 12> rt_mx::as_int_array() const
  {
    af::tiny<int, 12> result;
    for(std::size_t i=0;i<9;i++) result[i    ] = r_[i];
    for(std::size_t i=0;i<3;i++) result[i + 9] = t_[i];
    return result;
  }

  bool rt_mx::is_perpendicular(sg_vec3 const& v) const
  {
    return (r().accumulate() * v).is_zero();
  }

  tr_vec rt_mx::t_intrinsic_part() const
  {
    int type = r_.type();
    return r_.accumulate(type) * t_ / r_.order(type);
  }

  tr_vec rt_mx::t_location_part(tr_vec const& wi) const
  {
    return wi - t();
  }

  tr_vec rt_mx::t_origin_shift(tr_vec const& wl) const
  {
    rot_mx rmi = r_.minus_unit_mx();
    rot_mx p(1);
    scitbx::mat_ref<int> ref_rmi(rmi.num().begin(), 3, 3);
    scitbx::mat_ref<int> ref_p(p.num().begin(), 3, 3);
    row_echelon::form_t(ref_rmi, ref_p);
    tr_vec pwl = p * wl;
    tr_vec sh(0);
    sh.den() = row_echelon::back_substitution(
      ref_rmi, pwl.num().begin(), sh.num().begin());
    CCTBX_ASSERT(sh.den() > 0);
    sh.den() *= pwl.den();
    return sh;
  }

  rt_mx rt_mx::cancel() const
  {
    return rt_mx(r_.cancel(), t_.cancel());
  }

  rt_mx rt_mx::multiply(rt_mx const& rhs) const
  {
    if (t().den() == rhs.t().den()) {
      return rt_mx(
        (r_ * rhs.r()),
        (r_ * rhs.t() + t_.scale(r().den()))).cancel();
    }
    else {
      int f = boost::lcm(t().den(), rhs.t().den());
      int l = f / t().den() * r().den();
      int r = f / rhs.t().den();
      return rt_mx(
        (r_ * rhs.r()),
        (r_ * rhs.t().scale(r) + t_.scale(l))).cancel();
    }
  }

  translation_part_info::translation_part_info(rt_mx const& s)
  {
    ip_ = s.t_intrinsic_part();
    lp_ = s.t_location_part(ip_);
    os_ = s.t_origin_shift(lp_);
  }

}} // namespace cctbx::sgtbx
