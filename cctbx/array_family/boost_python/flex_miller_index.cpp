/* Copyright (c) 2001-2002 The Regents of the University of California
   through E.O. Lawrence Berkeley National Laboratory, subject to
   approval by the U.S. Department of Energy.
   See files COPYRIGHT.txt and LICENSE.txt for further details.

   Revision history:
     2002 Aug: Created (R.W. Grosse-Kunstleve)
 */

#include <scitbx/array_family/boost_python/flex_wrapper.h>
#include <scitbx/boost_python/pickle_single_buffered.h>
#include <cctbx/miller.h>

namespace scitbx { namespace boost_python { namespace pickle_single_buffered {

  inline
  char* to_string(char* start, cctbx::miller::index<> const& value)
  {
    return
      to_string(to_string(to_string(start, value[0]), value[1]), value[2]);
  }

  template <>
  struct from_string<cctbx::miller::index<> >
  {
    from_string(const char* start)
    {
      end = start;
      for(std::size_t i=0;i<3;i++) {
        from_string<int> proxy(end);
        value[i] = proxy.value;
        end = proxy.end;
      }
    }

    cctbx::miller::index<> value;
    const char* end;
  };

}}} // namespace scitbx::boost_python::pickle_single_buffered

#include <scitbx/array_family/boost_python/flex_pickle_single_buffered.h>

namespace scitbx { namespace af { namespace boost_python {

  void wrap_flex_miller_index(boost::python::object const& flex_root_scope)
  {
    using namespace cctbx;

    flex_wrapper<miller::index<> >::ordered("miller_index", flex_root_scope)
      .def_pickle(flex_pickle_single_buffered<miller::index<>,
        3*pickle_size_per_element<miller::index<>::value_type>::value>());
  }

}}} // namespace scitbx::af::boost_python
