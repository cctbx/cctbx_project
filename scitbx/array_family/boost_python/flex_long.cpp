#include <scitbx/array_family/boost_python/flex_wrapper.h>
#include <scitbx/array_family/boost_python/flex_pickle_single_buffered.h>
#include <scitbx/array_family/boost_python/range_wrappers.h>
#include <scitbx/array_family/counts.h>
#include <scitbx/stl/map_fwd.h>
#include <boost/python/args.hpp>
#include <map>

namespace scitbx { namespace af { namespace boost_python {

  void wrap_flex_long()
  {
    using namespace boost::python;
    flex_wrapper<long>::signed_integer("long", boost::python::scope())
      .def_pickle(flex_pickle_single_buffered<long>())
      .def("counts", counts<long, std::map<long, long> >::unlimited)
      .def("counts", counts<long, std::map<long, long> >::limited, (
        arg_("max_keys")))
    ;
    range_wrappers<long, long>::wrap("long_range");
  }

}}} // namespace scitbx::af::boost_python
