#include <cctbx/sgtbx/rot_mx.h>
#include <cctbx/sgtbx/utils.h>
#include <scitbx/math/gcd.h>

namespace cctbx { namespace sgtbx {

  void throw_unsuitable_rot_mx(const char* file, long line)
  {
    throw error_rational_vector(file, line,
      "Unsuitable value for rational rotation matrix.");
  }

  rot_mx rot_mx::new_denominator(int new_den) const
  {
    rot_mx result(new_den);
    if (utils::change_denominator(num_.begin(), den(),
                                  result.num_.begin(), new_den,
                                  num_.size()) != 0) {
      throw_unsuitable_rot_mx(__FILE__, __LINE__);
    }
    return result;
  }

  rot_mx rot_mx::inverse() const
  {
    int det_den3 = num_.determinant();
    if (det_den3 == 0) throw error("Rotation matrix is not invertible.");
    return rot_mx(
      num_.co_factor_matrix_transposed() * (den_*den_), den_) / det_den3;
  }

  rot_mx operator/(rot_mx const& lhs, int rhs)
  {
    sg_mat3 new_num;
    for(std::size_t i=0;i<9;i++) {
      if (lhs.num_[i] % rhs) throw_unsuitable_rot_mx(__FILE__, __LINE__);
      new_num[i] = lhs.num_[i] / rhs;
    }
    return rot_mx(new_num, lhs.den_);
  }

  rot_mx rot_mx::cancel() const
  {
    int g = den();
    for(std::size_t i=0;i<9;i++) g = scitbx::math::gcd_int(g, num_[i]);
    if (g == 0) return *this;
    return rot_mx(num_ / g, den() / g);
  }

  rot_mx rot_mx::inverse_cancel() const
  {
    int det_den3 = num_.determinant();
    if (det_den3 == 0) throw error("Rotation matrix is not invertible.");
    rat d(det_den3, den_);
    return rot_mx(num_.co_factor_matrix_transposed() * d.denominator(), 1)
      .divide(d.numerator());
  }

  rot_mx rot_mx::divide(int rhs) const
  {
    sg_mat3 new_num;
    if (rhs < 0) {
      new_num = -num_;
      rhs = -rhs;
    }
    else {
      new_num = num_;
    }
    return rot_mx(new_num, den_ * rhs).cancel();
  }

  int rot_mx::type() const
  {
    int det = num_.determinant();
    if (det == -1 || det == 1) {
      switch (num_.trace()) {
        case -3:                return -1;
        case -2:                return -6;
        case -1: if (det == -1) return -4;
                 else           return  2;
        case  0: if (det == -1) return -3;
                 else           return  3;
        case  1: if (det == -1) return -2;
                 else           return  4;
        case  2:                return  6;
        case  3:                return  1;
      }
    }
    return 0;
  }

  int rot_mx::order(int type) const
  {
    if (type == 0) type = rot_mx::type();
    if (type > 0) return  type;
    if (type % 2) return -type * 2;
                  return -type;
  }

  rot_mx rot_mx::accumulate(int type) const
  {
    CCTBX_ASSERT(den_ == 1);
    int ord = order(type);
    if (ord == 1) return *this;
    CCTBX_ASSERT(ord != 0);
    sg_mat3 result(1);
    sg_mat3 a(num_);
    result += a;
    for(int i=2; i < ord; ++i) {
      a = a*num_;
      result += a;
    }
    return rot_mx(result, 1);
  }

}} // namespace cctbx::sgtbx
