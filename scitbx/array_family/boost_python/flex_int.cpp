#include <scitbx/array_family/boost_python/flex_wrapper.h>
#include <scitbx/array_family/boost_python/flex_pickle_single_buffered.h>
#include <scitbx/array_family/counts.h>
#include <scitbx/stl/map_fwd.h>
#include <boost/python/args.hpp>
#include <map>

namespace scitbx { namespace af { namespace boost_python {

  void wrap_flex_int()
  {
    using namespace boost::python;
    flex_wrapper<int>::signed_integer("int", scope())
      .def_pickle(flex_pickle_single_buffered<int>())
      .def("counts", counts<int, std::map<long, long> >::unlimited)
      .def("counts", counts<int, std::map<long, long> >::limited, (
        arg_("max_keys")))
    ;
  }

}}} // namespace scitbx::af::boost_python
