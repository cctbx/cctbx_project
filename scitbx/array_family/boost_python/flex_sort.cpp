/* Copyright (c) 2001-2002 The Regents of the University of California
   through E.O. Lawrence Berkeley National Laboratory, subject to
   approval by the U.S. Department of Energy.
   See files COPYRIGHT.txt and LICENSE.txt for further details.

   Revision history:
     2003 Jan: Created (rwgk)
 */

#include <scitbx/array_family/boost_python/flex_fwd.h>

#include <scitbx/array_family/sort.h>
#include <boost/python/def.hpp>
#include <boost/python/overloads.hpp>

namespace scitbx { namespace af { namespace boost_python {

namespace {

  BOOST_PYTHON_FUNCTION_OVERLOADS(
    sort_permutation_overloads, sort_permutation, 1, 2)

} // namespace <anonymous>

  void wrap_flex_sort()
  {
    using namespace boost::python;

    def("sort_permutation",
      (shared<std::size_t>(*)(const_ref<int> const&, bool)) 0,
      sort_permutation_overloads());
    def("sort_permutation",
      (shared<std::size_t>(*)(const_ref<std::size_t> const&, bool)) 0,
      sort_permutation_overloads());
    def("sort_permutation",
      (shared<std::size_t>(*)(const_ref<double> const&, bool)) 0,
      sort_permutation_overloads());
  }

}}} // namespace scitbx::af::boost_python
