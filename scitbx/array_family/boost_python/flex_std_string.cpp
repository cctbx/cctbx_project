/* Copyright (c) 2001-2002 The Regents of the University of California
   through E.O. Lawrence Berkeley National Laboratory, subject to
   approval by the U.S. Department of Energy.
   See files COPYRIGHT.txt and LICENSE.txt for further details.

   Revision history:
     2002 Aug: Created (R.W. Grosse-Kunstleve)
 */

#include <scitbx/array_family/boost_python/flex_wrapper.h>
#include <scitbx/array_family/boost_python/flex_pickle_double_buffered.h>

namespace scitbx { namespace af { namespace boost_python {

  void add_flex_std_string(boost::python::module& m)
  {
    flex_wrapper<std::string>::ordered(m, "std_string")
      .def_pickle(flex_pickle_double_buffered<std::string>());
  }

}}} // namespace scitbx::af::boost_python
