// $Id$
/* Copyright (c) 2001 The Regents of the University of California through
   E.O. Lawrence Berkeley National Laboratory, subject to approval by the
   U.S. Department of Energy. See files COPYRIGHT.txt and
   cctbx/LICENSE.txt for further details.

   Revision history:
     2002 Jul: Created from fragments of cctbx/miller.h (R.W. Grosse-Kunstleve)
 */

#ifndef CCTBX_MILLER_SELECTION_H
#define CCTBX_MILLER_SELECTION_H

#include <cctbx/miller.h>
#include <cctbx/array_family/shared.h>

namespace cctbx { namespace miller {

  class selection
  {
    public:
      selection() {}

      explicit
      selection(af::shared<Index> miller_indices, bool flag = true)
      : miller_indices_(miller_indices),
        flags_(miller_indices.size(), flag)
      {}

      // copy constructor: deep-copy flags
      selection(selection const& other)
      : miller_indices_(other.miller_indices_),
        flags_(other.flags_.deep_copy())
      {
        size_assert_intrinsic();
      }

      // assignment operator: deep-copy flags
      selection operator=(selection const& other)
      {
        other.size_assert_intrinsic();
        miller_indices_ = other.miller_indices_;
        flags_ = other.flags_.deep_copy();
        return *this;
      }

      std::size_t size_processed() const
      {
        return flags_.size();
      }

      void size_assert_intrinsic() const
      {
        cctbx_assert(miller_indices_.size() == flags_.size());
      }

      void size_assert(std::size_t sz) const
      {
        size_assert_intrinsic();
        cctbx_assert(sz == size_processed());
      }

      std::size_t n_selected() const;

      af::shared<Index> selected_miller_indices() const;

      template <typename DataType>
      af::shared<DataType>
      selected_data(af::shared<DataType> data) const;

      void negate();

      template <typename FloatType>
      void sigma_filter(
        af::shared<FloatType> data,
        af::shared<FloatType> sigmas,
        FloatType const& cutoff_factor);

    protected:
      af::shared<Index> miller_indices_;
      af::shared<bool> flags_;
  };

  template <typename DataType>
  af::shared<DataType>
  selection::selected_data(af::shared<DataType> data) const
  {
    af::shared<DataType> result;
    size_assert(data.size());
    for(std::size_t i=0;i<flags_.size();i++) {
      if (flags_[i] == true) result.push_back(data[i]);
    }
    return result;
  }

  template <typename FloatType>
  void selection::sigma_filter(
    af::shared<FloatType> data,
    af::shared<FloatType> sigmas,
    FloatType const& cutoff_factor)
  {
    size_assert(data.size());
    size_assert(sigmas.size());
    for(std::size_t i=0;i<flags_.size();i++) {
      if (math::abs(data[i]) < cutoff_factor * sigmas[i]) {
        flags_[i] = false;
      }
    }
  }

}} // namespace cctbx::miller

#endif // CCTBX_MILLER_SELECTION_H
