/* Copyright (c) 2001-2002 The Regents of the University of California
   through E.O. Lawrence Berkeley National Laboratory, subject to
   approval by the U.S. Department of Energy.
   See files COPYRIGHT.txt and LICENSE.txt for further details.

   Revision history:
     2002 Sep: Created (R.W. Grosse-Kunstleve)
 */

#include <scitbx/array_family/boost_python/flex_wrapper.h>

namespace scitbx { namespace af { namespace boost_python {

  shared<size_t>
  column(
    const_ref<tiny<std::size_t, 2> > const& a,
    std::size_t i_column)
  {
    SCITBX_ASSERT(i_column < 2);
    shared<size_t> result((reserve(a.size())));
    for(std::size_t i=0;i<a.size();i++) {
      result.push_back(a[i][i_column]);
    }
    return result;
  }

  void wrap_flex_tiny_size_t_2()
  {
    flex_wrapper<tiny<std::size_t, 2> >::plain("tiny_size_t_2")
      .def("column", column);
    ;
  }

}}} // namespace scitbx::af::boost_python
