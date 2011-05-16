#include <cctbx/sgtbx/rt_mx.h>
#include <scitbx/matrix/row_echelon.h>
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

  rt_mx rt_mx::operator+(rt_mx const& rhs) const
  {
    return rt_mx(r_ + rhs.r_, t_ + rhs.t_);
  }

  rt_mx rt_mx::operator*(rt_mx const& rhs) const
  {
    CCTBX_ASSERT(t_.den() == rhs.t_.den());
    return rt_mx(r_ * rhs.r_, r_ * rhs.t_ + t_.scale(r_.den()));
  }

  rt_mx rt_mx::operator+(tr_vec const& rhs) const
  {
    return rt_mx(r_, t_ + rhs);
  }

  namespace {

    void
    throw_parse_error(
      parse_string const& input,
      std::string const& info=": unexpected character")
    {
      throw error(
        "Parse error" + info + ":\n"
        + input.format_where_message(/* prefix */ "  "));
    }

    int rationalize(double fVal, int& iVal, int den)
    {
      if (den == 0) return -1;
          fVal *= den;
      if (fVal < 0.) iVal = int(fVal - .5);
      else           iVal = int(fVal + .5);
          fVal -= iVal;
          fVal /= den;
      if (fVal < 0.) fVal = -fVal;
      if (fVal > .0005) return -1;
      return 0;
    }

  } // namespace <anonymous>

  rt_mx_from_string::rt_mx_from_string(
    parse_string& input,
    const char* stop_chars,
    int r_den,
    int t_den,
    bool enable_xyz,
    bool enable_hkl,
    bool enable_abc)
  :
    rt_mx(r_den, t_den),
    have_xyz(false),
    have_hkl(false),
    have_abc(false)
  {
    int Row    = 0;
    int Column = -1;
    int Sign   = 1;
    int Mult   = 0;
    double ValR[3];
    int i;
    for(i=0;i<3;i++) ValR[i] = 0.;
    double ValT   = 0.;
    double Value  = 0.;
    bool have_value = false;
    const unsigned int P_Add   = 0x01u;
    const unsigned int P_Mult  = 0x02u;
    const unsigned int P_Value = 0x04u;
    const unsigned int P_XYZ   = 0x08u;
    const unsigned int P_Comma = 0x10u;
    unsigned int P_mode = P_Add | P_Value | P_XYZ;
    for (;; input.skip())
    {
      if (strchr(stop_chars, input()) || !isspace(input()))
      {
        switch (strchr(stop_chars, input()) ? '\0' : input())
        {
          case '_':
            break;
          case '+': Sign =  1; goto ProcessAdd;
          case '-': Sign = -1;
           ProcessAdd:
            if ((P_mode & P_Add) == 0) throw_parse_error(input);
            if (Column >= 0) ValR[Column] += Value;
            else             ValT         += Value;
            Value = 0.;
            have_value = false;
            Column = -1;
            Mult = 0;
            P_mode = P_Value | P_XYZ;
            break;
          case '*':
            if ((P_mode & P_Mult) == 0) throw_parse_error(input);
            Mult = 1;
            P_mode = P_Value | P_XYZ;
            break;
          case '/':
          case ':':
            if ((P_mode & P_Mult) == 0) throw_parse_error(input);
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
            if ((P_mode & P_Value) == 0) throw_parse_error(input);
            {
            const char *beginptr = input.peek();
            char *endptr;
            double V = std::strtod(beginptr, &endptr);
            if (endptr == beginptr) {
              throw_parse_error(input);
            }
            input.skip((endptr-beginptr)-1);
            if (Sign == -1) { V = -V; Sign = 1; }
            if      (Mult ==  1)
              Value *= V;
            else if (Mult == -1) {
              if      (V     != 0.) Value /= V;
              else if (Value != 0.) throw_parse_error(input);
            }
            else
              Value = V;
            }
            have_value = true;
            P_mode = P_Comma | P_Add | P_Mult | P_XYZ;
            break;
          case 'X':
          case 'x': Column = 0; have_xyz = true; goto Process_XYZ;
          case 'Y':
          case 'y': Column = 1; have_xyz = true; goto Process_XYZ;
          case 'Z':
          case 'z': Column = 2; have_xyz = true; goto Process_XYZ;
          case 'H':
          case 'h': Column = 0; have_hkl = true; goto Process_XYZ;
          case 'K':
          case 'k': Column = 1; have_hkl = true; goto Process_XYZ;
          case 'L':
          case 'l': Column = 2; have_hkl = true; goto Process_XYZ;
          case 'A':
          case 'a': Column = 0; have_abc = true; goto Process_XYZ;
          case 'B':
          case 'b': Column = 1; have_abc = true; goto Process_XYZ;
          case 'C':
          case 'c': Column = 2; have_abc = true; goto Process_XYZ;
           Process_XYZ:
            if (have_xyz && !enable_xyz) {
              throw_parse_error(
                input, ": x,y,z notation not supported in this context");
            }
            if (have_hkl && !enable_hkl) {
              throw_parse_error(
                input, ": h,k,l notation not supported in this context");
            }
            if (have_abc && !enable_abc) {
              throw_parse_error(
                input, ": a,b,c notation not supported in this context");
            }
            if (have_xyz && have_hkl) {
              throw_parse_error(input, ": mix of x,y,z and h,k,l notation");
            }
            if (have_xyz && have_abc) {
              throw_parse_error(input, ": mix of x,y,z and a,b,c notation");
            }
            if (have_hkl && have_abc) {
              throw_parse_error(input, ": mix of h,k,l and a,b,c notation");
            }
            if ((P_mode & P_XYZ) == 0) throw_parse_error(input);
            if (!have_value) { Value = Sign; Sign = 1; }
            P_mode = P_Comma | P_Add | P_Mult;
            break;
          case ',':
          case ';':
            if (Row == 2) {
              throw_parse_error(input, ": too many row expressions");
            }
          case '\0':
            if ((P_mode & P_Comma) == 0) {
              throw_parse_error(input, ": unexpected end of input");
            }
            if (Column >= 0) ValR[Column] += Value;
            else             ValT         += Value;
            for(i=0;i<3;i++) {
              if (rationalize(ValR[i], r()(Row, i), r_den) != 0) {
                throw_unsuitable_rot_mx(__FILE__, __LINE__);
              }
            }
            if (rationalize(ValT, t()[Row], t_den) != 0) {
              throw_unsuitable_tr_vec(__FILE__, __LINE__);
            }
            Row++;
            Column = -1;
            Sign = 1;
            Mult = 0;
            for(i=0;i<3;i++) ValR[i] = 0.;
            ValT = 0.;
            Value = 0.;
            have_value = false;
            P_mode = P_Add | P_Value | P_XYZ;
            break;
          default:
            throw_parse_error(input);
        }
      }
      if (strchr(stop_chars, input())) {
        break;
      }
    }
    if (Row != 3) throw_parse_error(input, ": not enough row expressions");
  }

  rt_mx::rt_mx(parse_string& symbol, const char* stop_chars,
               int r_den, int t_den)
    : r_(0), t_(0)
  {
    rt_mx_from_string result(
      symbol, stop_chars, r_den, t_den,
      /* enable_xyz */ true,
      /* enable_hkl */ true,
      /* enable_abc */ false);
    if (result.have_hkl) {
      CCTBX_ASSERT(result.t().is_zero());
      r_ = result.r().transpose();
    }
    else r_ = result.r();
    t_ = result.t();
  }

  rt_mx::rt_mx(std::string const& symbol, const char* stop_chars,
             int r_den, int t_den)
  : r_(0), t_(0)
  {
    parse_string parse_symbol(symbol);
    *this = rt_mx(parse_symbol, stop_chars, r_den, t_den);
  }

  rt_mx::rt_mx(scitbx::mat3<double> const& r,
               scitbx::vec3<double> const& t,
               int r_den, int t_den)
    : r_(0), t_(0)
  {
    rt_mx result(r_den, t_den);
    for(std::size_t i=0;i<9;i++) {
      if (rationalize(r[i], result.r_[i], r_den) != 0) {
        throw_unsuitable_rot_mx(__FILE__, __LINE__);
      }
    }
    for(std::size_t i=0;i<3;i++) {
      if (rationalize(t[i], result.t_[i], t_den) != 0) {
        throw_unsuitable_tr_vec(__FILE__, __LINE__);
      }
    }
    r_ = result.r_;
    t_ = result.t_;
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
    af::ref<int, af::mat_grid> ref_rmi(rmi.num().begin(), 3, 3);
    af::ref<int, af::mat_grid> ref_p(p.num().begin(), 3, 3);
    scitbx::matrix::row_echelon::form_t(ref_rmi, ref_p);
    tr_vec pwl = p * wl;
    tr_vec sh(0);
    sh.den() = scitbx::matrix::row_echelon::back_substitution_int(
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
