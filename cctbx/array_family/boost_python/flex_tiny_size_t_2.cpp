/* Copyright (c) 2001-2002 The Regents of the University of California
   through E.O. Lawrence Berkeley National Laboratory, subject to
   approval by the U.S. Department of Energy.
   See files COPYRIGHT.txt and LICENSE.txt for further details.

   Revision history:
     2002 Sep: Created (R.W. Grosse-Kunstleve)
 */

#include <scitbx/array_family/boost_python/flex_wrapper.h>

namespace scitbx { namespace af { namespace boost_python {

  void wrap_flex_tiny_size_t_2()
  {
    flex_wrapper<af::tiny<std::size_t, 2> >::plain("tiny_size_t_2");
  }

}}} // namespace scitbx::af::boost_python
