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

      selection(af::shared<Index> miller_indices, af::shared<bool> flags)
      : miller_indices_(miller_indices),
        flags_(flags.deep_copy())
      {
        size_assert_intrinsic();
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
      selected_data(af::shared<DataType> data) const
      {
        size_assert(data.size());
        af::shared<DataType> result;
        for(std::size_t i=0;i<flags_.size();i++) {
          if (flags_[i] == true) result.push_back(data[i]);
        }
        return result;
      }

      void negate();

      void operator()(af::shared<bool> flags);

    protected:
      af::shared<Index> miller_indices_;
      af::shared<bool> flags_;
  };

}} // namespace cctbx::miller

#endif // CCTBX_MILLER_SELECTION_H
