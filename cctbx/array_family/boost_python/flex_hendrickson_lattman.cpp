/* Copyright (c) 2001-2002 The Regents of the University of California
   through E.O. Lawrence Berkeley National Laboratory, subject to
   approval by the U.S. Department of Energy.
   See files COPYRIGHT.txt and LICENSE.txt for further details.

   Revision history:
     2002 Aug: Created (R.W. Grosse-Kunstleve)
 */

#include <cctbx/boost_python/flex_fwd.h>

#include <scitbx/array_family/boost_python/flex_wrapper.h>
#include <scitbx/boost_python/pickle_single_buffered.h>
#include <cctbx/hendrickson_lattman.h>

namespace scitbx { namespace boost_python { namespace pickle_single_buffered {

  inline
  char* to_string(char* start, cctbx::hendrickson_lattman<> const& value)
  {
    return to_string(to_string(to_string(to_string(start,
      value[0]), value[1]), value[2]), value[3]);
  }

  template <>
  struct from_string<cctbx::hendrickson_lattman<> >
  {
    from_string(const char* start)
    {
      end = start;
      for(std::size_t i=0;i<4;i++) {
        from_string<double> proxy(end);
        value[i] = proxy.value;
        end = proxy.end;
      }
    }

    cctbx::hendrickson_lattman<> value;
    const char* end;
  };

}}} // namespace scitbx::boost_python::pickle_single_buffered

#include <scitbx/array_family/boost_python/flex_pickle_single_buffered.h>

namespace scitbx { namespace af { namespace boost_python {

  void wrap_flex_hendrickson_lattman()
  {
    using namespace cctbx;

    flex_wrapper<hendrickson_lattman<> >::plain("hendrickson_lattman")
      .def_pickle(flex_pickle_single_buffered<hendrickson_lattman<>,
        4*pickle_size_per_element<
          hendrickson_lattman<>::base_type::value_type>::value>());
  }

}}} // namespace scitbx::af::boost_python
