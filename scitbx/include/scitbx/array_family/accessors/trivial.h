/* Copyright (c) 2001-2002 The Regents of the University of California
   through E.O. Lawrence Berkeley National Laboratory, subject to
   approval by the U.S. Department of Energy.
   See files COPYRIGHT.txt and LICENSE.txt for further details.

   Revision history:
     2002 Oct: Created (R.W. Grosse-Kunstleve)
 */

#ifndef SCITBX_ARRAY_FAMILY_ACCESSORS_TRIVIAL_H
#define SCITBX_ARRAY_FAMILY_ACCESSORS_TRIVIAL_H

#include <cstddef>

namespace scitbx { namespace af {

  class trivial_accessor
  {
    public:
      typedef std::size_t index_type;
      struct index_value_type {};

      trivial_accessor() : size_(0) {}

      trivial_accessor(std::size_t const& n) : size_(n) {}

      std::size_t size_1d() const { return size_; }

    protected:
      std::size_t size_;
  };

}} // namespace scitbx::af

#endif // SCITBX_ARRAY_FAMILY_ACCESSORS_TRIVIAL_H
