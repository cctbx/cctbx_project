#include <scitbx/array_family/boost_python/flex_wrapper.h>
#include <scitbx/array_family/boost_python/flex_pickle_single_buffered.h>
#include <scitbx/array_family/boost_python/range_wrappers.h>
#include <scitbx/array_family/counts.h>
#include <scitbx/array_family/versa_matrix.h>
#include <scitbx/matrix/move.h>
#include <scitbx/stl/map_fwd.h>
#include <boost/python/args.hpp>
#include <boost/python/make_constructor.hpp>
#include <map>

namespace scitbx { namespace af { namespace boost_python {

  void wrap_flex_long()
  {
    using namespace boost::python;
    using boost::python::arg;
    flex_wrapper<long>::signed_integer("long", boost::python::scope())
      .def_pickle(flex_pickle_single_buffered<long>())
      .def("counts", counts<long, std::map<long, long> >::unlimited)
      .def("counts", counts<long, std::map<long, long> >::limited, (
        arg("max_keys")))
      .def("matrix_copy_block",
        (versa<long, c_grid<2> >(*)(
          const_ref<long, c_grid<2> > const&,
          unsigned, unsigned, unsigned, unsigned))
            matrix::copy_block, (
              arg("i_row"),
              arg("i_column"),
              arg("n_rows"),
              arg("n_columns")))
      .def("matrix_paste_block_in_place",
        (void(*)(
          ref<long, c_grid<2> > const&,
          const_ref<long, c_grid<2> > const&,
          unsigned, unsigned))
            matrix::paste_block_in_place, (
              arg("block"),
              arg("i_row"),
              arg("i_column")))
    ;
    range_wrappers<long, long>::wrap("long_range");
  }

}}} // namespace scitbx::af::boost_python
