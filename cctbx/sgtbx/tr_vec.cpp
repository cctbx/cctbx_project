#include <cctbx/sgtbx/tr_vec.h>
#include <cctbx/sgtbx/utils.h>
#include <scitbx/math/gcd.h>

namespace cctbx { namespace sgtbx {

  void throw_unsuitable_tr_vec(const char* file, long line)
  {
    throw error_rational_vector(file, line,
      "Unsuitable value for rational translation vector.");
  }

  tr_vec tr_vec::new_denominator(int new_den) const
  {
    tr_vec result(new_den);
    if (utils::change_denominator(num_.begin(), den(),
                                  result.num_.begin(), new_den,
                                  num_.size()) != 0) {
      throw_unsuitable_tr_vec(__FILE__, __LINE__);
    }
    return result;
  }

  tr_vec operator/(tr_vec const& lhs, int rhs)
  {
    sg_vec3 new_num;
    for(std::size_t i=0;i<3;i++) {
      if (lhs.num_[i] % rhs) throw_unsuitable_tr_vec(__FILE__, __LINE__);
      new_num[i] = lhs.num_[i] / rhs;
    }
    return tr_vec(new_num, lhs.den_);
  }

  tr_vec tr_vec::cancel() const
  {
    int g = den();
    for(std::size_t i=0;i<3;i++) g = scitbx::math::gcd_int(g, num_[i]);
    if (g == 0) return *this;
    return tr_vec(num_ / g, den() / g);
  }

  tr_vec tr_vec::plus(tr_vec const& rhs) const
  {
    tr_vec result(boost::lcm(den(), rhs.den()));
    int l = result.den() / den();
    int r = result.den() / rhs.den();
    for(std::size_t i=0;i<3;i++) result[i] = num_[i] * l + rhs[i] * r;
    return result.cancel();
  }

  tr_vec tr_vec::minus(tr_vec const& rhs) const
  {
    tr_vec result(boost::lcm(den(), rhs.den()));
    int l = result.den() / den();
    int r = result.den() / rhs.den();
    for(std::size_t i=0;i<3;i++) result[i] = num_[i] * l - rhs[i] * r;
    return result.cancel();
  }

  std::string
  tr_vec::as_string(bool decimal, const char* separator) const
  {
    std::string result;
    for(int i=0;i<3;i++) {
      if (i != 0) result += separator;
      rat t_frac((*this)[i], den());
      result += scitbx::format(t_frac, decimal);
    }
    return result;
  }

}} // namespace cctbx::sgtbx
