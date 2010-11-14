#include <scitbx/array_family/boost_python/flex_wrapper.h>
#include <scitbx/array_family/boost_python/flex_pickle_single_buffered.h>
#include <scitbx/array_family/boost_python/range_wrappers.h>
#include <scitbx/array_family/counts.h>
#include <scitbx/stl/map_fwd.h>
#include <boost/python/args.hpp>
#include <map>

namespace scitbx { namespace af { namespace boost_python {

  boost::python::object
  ref_flex_as_numpy_array(
    ref<long, flex_grid<> > const& O);

  void wrap_flex_long()
  {
    using namespace boost::python;
    using boost::python::arg;
    flex_wrapper<long>::signed_integer("long", boost::python::scope())
      .def_pickle(flex_pickle_single_buffered<long>())
      .def("counts", counts<long, std::map<long, long> >::unlimited)
      .def("counts", counts<long, std::map<long, long> >::limited, (
        arg("max_keys")))
      .def("as_numpy_array", ref_flex_as_numpy_array)
    ;
    range_wrappers<long, long>::wrap("long_range");
  }

}}} // namespace scitbx::af::boost_python
