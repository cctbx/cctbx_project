/* Copyright (c) 2001-2002 The Regents of the University of California
   through E.O. Lawrence Berkeley National Laboratory, subject to
   approval by the U.S. Department of Energy.
   See files COPYRIGHT.txt and LICENSE.txt for further details.

   Revision history:
     2002 Sep: Refactored parts of cctbx/sgtbx/matrix.h (rwgk)
     2001 Jul: Merged from CVS branch sgtbx_special_pos (rwgk)
     2001 May: merged from CVS branch sgtbx_type (R.W. Grosse-Kunstleve)
     2001 Apr: SourceForge release (R.W. Grosse-Kunstleve)
 */

#ifndef CCTBX_SGTBX_ROT_MX_H
#define CCTBX_SGTBX_ROT_MX_H

#include <cctbx/sgtbx/tr_vec.h>

namespace cctbx { namespace sgtbx {

  void throw_unsuitable_rot_mx(const char* file, long line);

  class rot_mx_info; // forward declaration

  class rot_mx
  {
    public:
      explicit
      rot_mx(int denominator = 1, int diagonal = 1)
      : num_(diagonal * denominator), den_(denominator)
      {}

      explicit
      rot_mx(sg_mat3 const& m, int denominator = 1)
      : num_(m), den_(denominator)
      {}

      rot_mx(int m00, int m01, int m02,
             int m10, int m11, int m12,
             int m20, int m21, int m22,
             int denominator = 1)
        : num_(m00,m01,m02, m10,m11,m12, m20,m21, m22),
          den_(denominator)
      {}

      sg_mat3 const& num() const { return num_; }
      sg_mat3&       num()       { return num_; }

      int const& operator[](std::size_t i) const { return num_[i]; }
      int&       operator[](std::size_t i)       { return num_[i]; }

      int const& operator()(int r, int c) const { return num_(r, c); }
      int&       operator()(int r, int c)       { return num_(r, c); }

      int const& den() const { return den_; }
      int&       den()       { return den_; }

      bool operator==(rot_mx const& rhs) const
      {
        if (den_ != rhs.den_) return false;
        return num_.const_ref().all_eq(rhs.num_.const_ref());
      }

      bool operator!=(rot_mx const& rhs) const
      {
        return !((*this) == rhs);
      }

      bool is_valid() const { return den_ != 0; }

      bool is_unit_mx() const { return num_ == sg_mat3(den_); }

      rot_mx minus_unit_mx() const
      {
        rot_mx result(*this);
        for (std::size_t i=0;i<9;i+=4) result[i] -= den_;
        return result;
      }

      rot_mx new_denominator(int new_den) const;

      rot_mx scale(int factor) const
      {
        if (factor == 1) return *this;
        return rot_mx(num_ * factor, den_ * factor);
      }

      rot_mx inverse() const;

      rot_mx operator-() const { return rot_mx(-num_, den_); }

      friend rot_mx operator+(rot_mx const& lhs, rot_mx const& rhs)
      {
        CCTBX_ASSERT(lhs.den_ == rhs.den_);
        return rot_mx(lhs.num_ + rhs.num_, lhs.den_);
      }

      friend rot_mx operator-(rot_mx const& lhs, rot_mx const& rhs)
      {
        CCTBX_ASSERT(lhs.den_ == rhs.den_);
        return rot_mx(lhs.num_ - rhs.num_, lhs.den_);
      }

      rot_mx& operator+=(rot_mx const& rhs)
      {
        CCTBX_ASSERT(den_ == rhs.den_);
        num_ += rhs.num_;
        return *this;
      }

      friend rot_mx operator*(rot_mx const& lhs, rot_mx const& rhs)
      {
        return rot_mx(lhs.num_ * rhs.num_, lhs.den_ * rhs.den_);
      }

      friend tr_vec operator*(rot_mx const& lhs, tr_vec const& rhs)
      {
        return tr_vec(lhs.num_ * rhs.num(), lhs.den_ * rhs.den());
      }

      friend tr_vec operator*(tr_vec const& lhs, rot_mx const& rhs)
      {
        return tr_vec(lhs.num() * rhs.num_, lhs.den() * rhs.den_);
      }

      friend sg_vec3 operator*(rot_mx const& lhs, sg_vec3 const& rhs)
      {
        return sg_vec3(lhs.num_ * rhs);
      }

      friend rot_mx operator*(rot_mx const& lhs, int rhs)
      {
        return rot_mx(lhs.num_ * rhs, lhs.den_);
      }

      friend rot_mx operator*(int lhs, rot_mx const& rhs)
      {
        return rhs * lhs;
      }

      rot_mx& operator*=(int rhs)
      {
        num_ *= rhs;
        return *this;
      }

      friend rot_mx operator/(rot_mx const& lhs, int rhs);

      rot_mx cancel() const;

      rot_mx inverse_cancel() const;

      rot_mx multiply(rot_mx const& rhs) const
      {
        return ((*this) * rhs).cancel();
      }

      tr_vec multiply(tr_vec const& rhs) const
      {
        return ((*this) * rhs).cancel();
      }

      rot_mx divide(int rhs) const;

      int type() const;

      int order(int type=0) const;

      rot_mx accumulate(int type=0) const;

      rot_mx_info info() const;

      template <typename FloatType>
      scitbx::mat3<FloatType>
      as_floating_point(scitbx::type_holder<FloatType>) const
      {
        return scitbx::mat3<FloatType>(num_) / FloatType(den_);
      }

      scitbx::mat3<double>
      as_double() const
      {
        return as_floating_point(scitbx::type_holder<double>());
      }

    private:
      sg_mat3 num_;
      int den_;
  };

}} // namespace cctbx::sgtbx

#endif // CCTBX_SGTBX_ROT_MX_H
