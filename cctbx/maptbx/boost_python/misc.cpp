/* Copyright (c) 2001-2002 The Regents of the University of California
   through E.O. Lawrence Berkeley National Laboratory, subject to
   approval by the U.S. Department of Energy.
   See files COPYRIGHT.txt and LICENSE.txt for further details.

   Revision history:
     2002 Dec: Created (R.W. Grosse-Kunstleve)
 */

#include <cctbx/boost_python/flex_fwd.h>

#include <boost/python/def.hpp>

namespace cctbx { namespace maptbx { namespace boost_python {

namespace {

  template <typename ElementType>
  struct set_if_less_than
  {
    static void
    run(
      af::versa<ElementType, af::flex_grid<> >& a,
      ElementType const& threshold_value,
      ElementType const& imposed_value)
    {
      for(std::size_t i=0;i<a.size();i++) {
        if (a[i] < threshold_value) a[i] = imposed_value;
      }
    }
  };

} // namespace <anoymous>

  void wrap_misc()
  {
    using namespace boost::python;
    def("set_if_less_than", set_if_less_than<float>::run);
    def("set_if_less_than", set_if_less_than<double>::run);
  }

}}} // namespace cctbx::maptbx::boost_python
