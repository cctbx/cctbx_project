#include <scitbx/array_family/boost_python/flex_wrapper.h>
#include <scitbx/array_family/boost_python/flex_pickle_single_buffered.h>
#include <scitbx/array_family/boost_python/byte_str.h>
#include <scitbx/array_family/boost_python/range_wrappers.h>
#include <scitbx/array_family/boost_python/numpy_bridge.hpp>
#include <scitbx/array_family/counts.h>
#include <scitbx/array_family/versa_matrix.h>
#include <scitbx/matrix/packed.h>
#include <scitbx/matrix/move.h>
#include <scitbx/stl/map_fwd.h>
#include <boost/python/args.hpp>
#include <boost/format.hpp>
#include <boost/python/make_constructor.hpp>
#include <boost/lexical_cast.hpp>
#include <map>

namespace scitbx { namespace af { namespace boost_python {

  flex<int>::type*
  from_std_string(const_ref<std::string> const& s)
  {
    shared<int> result(reserve(s.size()));
    for(std::size_t i=0;i<s.size();i++) {
      result.push_back(boost::lexical_cast<int>(s[i]));
    }
    return new flex<int>::type(result, result.size());
  }

  shared<bool>
  as_bool(const_ref<int> const& self, bool strict=true)
  {
    shared<bool> result((reserve(self.size())));
    for(std::size_t i=0;i<self.size();i++) {
      int v = self[i];
      if (v == 0) {
        result.push_back(false);
      }
      else if (v == 1 || !strict) {
        result.push_back(true);
      }
      else {
        throw std::invalid_argument((
          boost::format(
            "scitbx.array_family.flex.int.as_bool(strict=True):"
            " all array elements must be 0 or 1,"
            " but value=%d at array index=%lu.") % v % i).str());
      }
    }
    return result;
  }

  /* For allowed syntax for the optional format_string argument see:
       http://www.boost.org/libs/format/doc/format.html#syntax
   */
  af::shared<std::string>
  as_string(af::const_ref<int, af::flex_grid<> > const& O,
            std::string format_string="%d")
  {
    af::shared<std::string> result((reserve(O.size())));
    std::size_t n = O.accessor().size_1d();
    for(std::size_t i=0;i<n;i++) {
      result.push_back((boost::format(format_string) %O[i]).str());
    }
    return result;
  }

  void wrap_flex_int()
  {
    using namespace boost::python;
    using boost::python::arg;
    flex_wrapper<int>::signed_integer("int", scope())
      .def_pickle(flex_pickle_single_buffered<int>())
      .def("__init__", make_constructor(
        from_std_string, default_call_policies()))
      .def("__init__", make_constructor(
        flex_int_from_numpy_array, default_call_policies()))
      .def("copy_to_byte_str", copy_to_byte_str<versa<int, flex_grid<> > >)
      .def("slice_to_byte_str",
        slice_to_byte_str<versa<int, flex_grid<> > >)
      .def("as_bool", as_bool, (arg("strict")=true))
      .def("as_string", as_string, (
          arg("other"),
          arg("format_string")="%d"))
      .def("counts", counts<int, std::map<long, long> >::unlimited)
      .def("counts", counts<int, std::map<long, long> >::limited, (
        arg("max_keys")))
      .def("matrix_is_symmetric",
        (bool(*)(
          const_ref<int, c_grid<2> > const&))
            matrix::is_symmetric)
      .def("matrix_copy_block",
        (versa<int, c_grid<2> >(*)(
          const_ref<int, c_grid<2> > const&,
          unsigned, unsigned, unsigned, unsigned))
            matrix::copy_block, (
              arg("i_row"),
              arg("i_column"),
              arg("n_rows"),
              arg("n_columns")))
      .def("matrix_transpose_in_place",
        (void(*)(versa<int, flex_grid<> >&)) matrix_transpose_in_place)
      .def("matrix_swap_rows_in_place",
        (void(*)(
          ref<int, c_grid<2> > const&, unsigned, unsigned))
            matrix::swap_rows_in_place, (
              arg("i"),
              arg("j")))
      .def("matrix_swap_columns_in_place",
        (void(*)(
          ref<int, c_grid<2> > const&, unsigned, unsigned))
            matrix::swap_columns_in_place, (
              arg("i"),
              arg("j")))
      .def("matrix_paste_block_in_place",
        (void(*)(
          ref<int, c_grid<2> > const&,
          const_ref<int, c_grid<2> > const&,
          unsigned, unsigned))
            matrix::paste_block_in_place, (
              arg("block"),
              arg("i_row"),
              arg("i_column")))
      .def("as_numpy_array", flex_int_as_numpy_array)
    ;
    def(
      "int_from_byte_str",
      shared_from_byte_str<int>,
      (arg("byte_str")));
    range_wrappers<int, int>::wrap("int_range");
  }

}}} // namespace scitbx::af::boost_python
