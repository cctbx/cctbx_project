#include <scitbx/array_family/boost_python/flex_wrapper.h>
#include <scitbx/array_family/boost_python/flex_pickle_single_buffered.h>
#include <scitbx/array_family/boost_python/flex_wrapper_complex.h>
#include <scitbx/array_family/versa_matrix.h>
#include <scitbx/matrix/move.h>
#include <boost/python/make_constructor.hpp>

namespace scitbx { namespace af { namespace boost_python {

namespace {

  versa<std::complex<double>, flex_grid<> >*
  from_pair_of_flex_double(
    versa<double, flex_grid<> > const& reals,
    versa<double, flex_grid<> > const& imags)
  {
    SCITBX_ASSERT(reals.size() == imags.size());
    typedef versa<std::complex<double>, flex_grid<> > result_type;
    result_type result(
      reals.accessor(), init_functor_null<std::complex<double> >());
    std::complex<double>*c = result.begin();
    std::complex<double>*c_end = result.end();
    const double* r = reals.begin();
    const double* i = imags.begin();
    while (c != c_end) {
      *c++ = std::complex<double>(*r++, *i++);
    }
    return new result_type(result);
  }

  versa<std::complex<double>, flex_grid<> >
  mul_ac_ar(
    versa<std::complex<double>, flex_grid<> > const& a1,
    versa<             double , flex_grid<> > const& a2)
  {
    return a1 * a2;
  }

  versa<std::complex<double>, c_grid<2> >
  matrix_multiply_complex_matrix_complex_matrix(
    const_ref<std::complex<double>, c_grid<2> > const& a,
    const_ref<std::complex<double>, c_grid<2> > const& b)
  {
    return matrix_multiply(a, b);
  }

  versa<std::complex<double>, c_grid<2> >
  matrix_multiply_complex_matrix_real_matrix(
    const_ref<std::complex<double>, c_grid<2> > const& a,
    const_ref<double, c_grid<2> > const& b)
  {
    return matrix_multiply(a, b);
  }

} // namespace <anonymous>

  void wrap_flex_complex_double()
  {
    using namespace boost::python;
    flex_wrapper<std::complex<double> >::numeric_common(
      "complex_double", boost::python::scope())
      .def("__init__", make_constructor(
        from_pair_of_flex_double, default_call_policies(), (
          arg_("reals"), arg_("imags"))))
      .def_pickle(flex_pickle_single_buffered<std::complex<double> >())
      .def("__mul__", mul_ac_ar)
      .def("__rmul__", mul_ac_ar)
      .def("matrix_multiply", matrix_multiply_complex_matrix_complex_matrix)
      .def("matrix_multiply", matrix_multiply_complex_matrix_real_matrix)
      .def("matrix_transpose",
        (versa<std::complex<double>, c_grid<2> >(*)(
           const_ref<std::complex<double>, c_grid<2> > const&))
             matrix_transpose)
      .def("matrix_packed_u_as_symmetric",
        (versa<std::complex<double>, c_grid<2> >(*)(
          const_ref<std::complex<double> > const&))
            matrix::packed_u_as_symmetric)
      .def("matrix_copy_block",
        (versa<std::complex<double>, c_grid<2> >(*)(
          const_ref<std::complex<double>, c_grid<2> > const&,
          unsigned, unsigned, unsigned, unsigned))
            matrix::copy_block, (
              arg_("i_row"),
              arg_("i_column"),
              arg_("n_rows"),
              arg_("n_columns")))
      .def("matrix_paste_block_in_place",
        (void(*)(
          ref<std::complex<double>, c_grid<2> > const&,
          const_ref<std::complex<double>, c_grid<2> > const&,
          unsigned, unsigned))
            matrix::paste_block_in_place, (
              arg_("block"),
              arg_("i_row"),
              arg_("i_column")))
    ;
    flex_wrapper_complex_functions<double>::wrap(boost::python::scope());
  }

}}} // namespace scitbx::af::boost_python
