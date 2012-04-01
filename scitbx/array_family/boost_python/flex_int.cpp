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
      if (s[i].size() == 0) {
        throw std::invalid_argument(
          "Empty string (integer value expected).");
      }
      int value = 0;
      try {
        value = boost::lexical_cast<int>(s[i]);
      }
      catch (boost::bad_lexical_cast const&) {
        throw std::invalid_argument(
          "Invalid integer value: \"" + s[i] + "\"");
      }
      result.push_back(value);
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

  af::versa<long, af::flex_grid<> >
  as_long(
    af::const_ref<int, af::flex_grid<> > const& O)
  {
    af::versa<long, af::flex_grid<> > result(
      O.accessor(), af::init_functor_null<long>());
    std::size_t n = O.accessor().size_1d();
    long* r = result.begin();
    for(std::size_t i=0;i<n;i++) {
      r[i] = static_cast<long>(O[i]);
    }
    return result;
  }

  /* For allowed syntax for the optional format_string argument see:
       http://www.boost.org/libs/format/doc/format.html#syntax
   */
  af::shared<std::string>
  as_string(
    af::const_ref<int, af::flex_grid<> > const& O,
    std::string format_string="%d")
  {
    std::size_t n = O.accessor().size_1d();
    af::shared<std::string> result((af::reserve(n)));
    for(std::size_t i=0;i<n;i++) {
      result.push_back((boost::format(format_string) % O[i]).str());
    }
    return result;
  }

  std::string
  as_rgb_scale_string(
    af::const_ref<int, af::flex_grid<> > const& O,
    af::tiny<double, 3> const& rgb_scales,
    int saturation)
  {
    SCITBX_ASSERT(rgb_scales.const_ref().all_ge(0));
    SCITBX_ASSERT(rgb_scales.const_ref().all_le(1));
    SCITBX_ASSERT(saturation != 0);
    double scale = 1. / saturation;
    std::size_t n = O.accessor().size_1d();
    std::string result(n*3, '\0');
    std::size_t j = 0;
    for(std::size_t i=0;i<n;i++) {
      double f = O[i] * scale;
      if      (f < 0) f = 0;
      else if (f > 1) f = 1;
      for(unsigned k=0;k<3;k++) {
        double fs = f * rgb_scales[k];
        result[j++] = static_cast<char>(static_cast<int>((1-fs) * 255 + 0.5));
      }
    }
    return result;
  }

  std::string
  as_rgb_gray_scale_string(
    af::const_ref<int, af::flex_grid<> > const& O,
    int saturation)
  {
    return as_rgb_scale_string(O, af::tiny<double, 3>(1,1,1), saturation);
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
      .def("as_long", as_long)
      .def("as_string", as_string, (arg("format_string")="%d"))
      .def("as_rgb_scale_string", as_rgb_scale_string, (
        arg("rgb_scales"),
        arg("saturation")))
      .def("as_rgb_gray_scale_string", as_rgb_gray_scale_string, (
        arg("saturation")))
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
      .def("matrix_rot90",
        (versa<int, c_grid<2> >(*)(
           const_ref<int, c_grid<2> > const&, int)) matrix_rot90)
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
      .def("as_numpy_array", flex_int_as_numpy_array, (
        arg("optional")=false))
    ;
    def(
      "int_from_byte_str",
      shared_from_byte_str<int>,
      (arg("byte_str")));
    range_wrappers<int, int>::wrap("int_range");
  }

}}} // namespace scitbx::af::boost_python
