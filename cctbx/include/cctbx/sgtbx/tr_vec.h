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

#ifndef CCTBX_SGTBX_TR_VEC_H
#define CCTBX_SGTBX_TR_VEC_H

#include <cctbx/sgtbx/basic.h>
#include <cctbx/math/mod.h>
#include <scitbx/array_family/tiny.h>
#include <scitbx/type_holder.h>
#include <cctbx/import_scitbx_af.h>

namespace cctbx { namespace sgtbx {

  void throw_unsuitable_tr_vec(const char* file, long line);

  class tr_vec
  {
    public:
      explicit
      tr_vec(int tr_den=sg_t_den)
      : num_(0,0,0), den_(tr_den)
      {}

      explicit
      tr_vec(sg_vec3 const& v, int tr_den=sg_t_den)
      : num_(v), den_(tr_den)
      {}

      tr_vec(int v0, int v1, int v2, int tr_den=sg_t_den)
      : num_(v0,v1,v2), den_(tr_den)
      {}

      sg_vec3 const& num() const { return num_; }
      sg_vec3&       num()       { return num_; }

      int const& operator[](std::size_t i) const { return num_[i]; }
      int&       operator[](std::size_t i)       { return num_[i]; }

      int const& den() const { return den_; }
      int&       den()       { return den_; }

      bool operator==(tr_vec const& rhs) const
      {
        if (den_ != rhs.den_) return false;
        return num_.const_ref().all_eq(rhs.num_.const_ref());
      }

      bool operator!=(tr_vec const& rhs) const
      {
        return !((*this) == rhs);
      }

      bool is_valid() const { return den_ != 0; }

      bool is_zero() const { return num_.is_zero(); }

      tr_vec new_denominator(int new_bf) const;

      tr_vec scale(int factor) const
      {
        if (factor == 1) return *this;
        return tr_vec(num_ * factor, den_ * factor);
      }

      tr_vec mod_positive() const
      {
        tr_vec result(num_, den_);
        for(std::size_t i=0;i<3;i++) {
          result.num_[i] = math::mod_positive(result.num_[i], den_);
        }
        return result;
      }

      tr_vec mod_short() const
      {
        tr_vec result(num_, den_);
        for(std::size_t i=0;i<3;i++) {
          result.num_[i] = math::mod_short(result.num_[i], den_);
        }
        return result;
      }

      tr_vec operator-() const { return tr_vec(-num_, den_); }

      friend tr_vec operator+(tr_vec const& lhs, tr_vec const& rhs)
      {
        CCTBX_ASSERT(lhs.den_ == rhs.den_);
        return tr_vec(lhs.num_ + rhs.num_, lhs.den_);
      }

      friend tr_vec operator-(tr_vec const& lhs, tr_vec const& rhs)
      {
        CCTBX_ASSERT(lhs.den_ == rhs.den_);
        return tr_vec(lhs.num_ - rhs.num_, lhs.den_);
      }

      tr_vec& operator+=(tr_vec const& rhs)
      {
        CCTBX_ASSERT(den_ == rhs.den_);
        num_ += rhs.num_;
        return *this;
      }

      friend tr_vec operator*(tr_vec const& lhs, int rhs)
      {
        return tr_vec(lhs.num_ * rhs, lhs.den_);
      }

      friend tr_vec operator/(tr_vec const& lhs, int rhs);

      friend tr_vec operator*(int const& lhs, tr_vec const& rhs)
      {
        return rhs * lhs;
      }

      tr_vec cancel() const;

      tr_vec plus(tr_vec const& rhs) const;

      tr_vec minus(tr_vec const& rhs) const;

      template <typename FloatType>
      scitbx::vec3<FloatType>
      as_floating_point(scitbx::type_holder<FloatType>) const
      {
        return scitbx::vec3<FloatType>(num_) / FloatType(den_);
      }

      scitbx::vec3<double>
      as_double() const
      {
        return as_floating_point(scitbx::type_holder<double>());
      }

    private:
      sg_vec3 num_;
      int den_;
  };

  //! Constructor for initialization of constants.
  struct tr_vec_12 : tr_vec
  {
    tr_vec_12(int v0, int v1, int v2)
    : tr_vec(tr_vec(v0,v1,v2) * (sg_t_den / 12))
    {}
  };

}} // namespace cctbx::sgtbx

#endif // CCTBX_SGTBX_TR_VEC_H
