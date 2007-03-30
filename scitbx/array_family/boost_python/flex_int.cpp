#include <scitbx/array_family/boost_python/flex_wrapper.h>
#include <scitbx/array_family/boost_python/flex_pickle_single_buffered.h>
#include <scitbx/array_family/boost_python/range_wrappers.h>
#include <scitbx/array_family/counts.h>
#include <scitbx/matrix/move.h>
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
      .def("matrix_copy_block",
        (versa<int, c_grid<2> >(*)(
          const_ref<int, c_grid<2> > const&,
          unsigned, unsigned, unsigned, unsigned))
            matrix::copy_block, (
              arg_("i_row"),
              arg_("i_column"),
              arg_("n_rows"),
              arg_("n_columns")))
      .def("matrix_paste_block_in_place",
        (void(*)(
          ref<int, c_grid<2> > const&,
          const_ref<int, c_grid<2> > const&,
          unsigned, unsigned))
            matrix::paste_block_in_place, (
              arg_("block"),
              arg_("i_row"),
              arg_("i_column")))
    ;
    range_wrappers<int, int>::wrap("int_range");
  }

}}} // namespace scitbx::af::boost_python
