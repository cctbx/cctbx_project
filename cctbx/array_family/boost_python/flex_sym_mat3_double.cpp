/* Copyright (c) 2001-2002 The Regents of the University of California
   through E.O. Lawrence Berkeley National Laboratory, subject to
   approval by the U.S. Department of Energy.
   See files COPYRIGHT.txt and LICENSE.txt for further details.

   Revision history:
     2002 Sep: Created (R.W. Grosse-Kunstleve)
 */

#include <scitbx/array_family/boost_python/flex_wrapper.h>
#include <scitbx/boost_python/pickle_single_buffered.h>
#include <scitbx/sym_mat3.h>

namespace scitbx { namespace boost_python { namespace pickle_single_buffered {

  inline
  char* to_string(char* start, sym_mat3<double> const& value)
  {
    for(std::size_t i=0;i<6;i++) start = to_string(start, value[i]);
    return start;
  }

  template <>
  struct from_string<sym_mat3<double> >
  {
    from_string(const char* start)
    {
      end = start;
      for(std::size_t i=0;i<6;i++) {
        from_string<double> proxy(end);
        value[i] = proxy.value;
        end = proxy.end;
      }
    }

    sym_mat3<double> value;
    const char* end;
  };

}}} // namespace scitbx::boost_python::pickle_single_buffered

#include <scitbx/array_family/boost_python/flex_pickle_single_buffered.h>

namespace scitbx { namespace af { namespace boost_python {

  void wrap_flex_sym_mat3_double()
  {
    flex_wrapper<sym_mat3<double> >::plain("sym_mat3_double")
      .def_pickle(flex_pickle_single_buffered<sym_mat3<double>,
        6*pickle_size_per_element<double>::value>());
  }

}}} // namespace scitbx::af::boost_python
