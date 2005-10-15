#include <scitbx/array_family/boost_python/flex_wrapper.h>
#include <scitbx/array_family/versa_matrix.h>
#include <scitbx/matrix/outer_product.h>
#include <scitbx/matrix/cholesky.h>
#include <boost/python/args.hpp>
#include <boost/python/overloads.hpp>
#include <boost/python/return_value_policy.hpp>
#include <boost/python/return_by_value.hpp>

namespace scitbx { namespace af {

namespace {

  // test here in lack of a better place
  void
  exercise_packed_u_accessor()
  {
    af::versa<double, matrix::packed_u_accessor>
      a(matrix::packed_u_accessor(5));
    af::ref<double, matrix::packed_u_accessor> r = a.ref();
    SCITBX_ASSERT(a.size() == 5*(5+1)/2);
    SCITBX_ASSERT(r.size() == 5*(5+1)/2);
    SCITBX_ASSERT(a.accessor().n == 5);
    SCITBX_ASSERT(r.accessor().n == 5);
    for(unsigned i=0;i<r.size();i++) r[i] = i+1;
    unsigned v = 1;
    for(unsigned i=0;i<r.accessor().n;i++) {
      for(unsigned j=i;j<r.accessor().n;j++,v++) {
        SCITBX_ASSERT(a(i,j) == v);
        SCITBX_ASSERT(r(i,j) == v);
      }
    }
  }

  struct matrix_cholesky_gill_murray_wright_wrappers
  {
    typedef matrix::cholesky::gill_murray_wright w_t;

    static void
    wrap()
    {
      using namespace boost::python;
      typedef return_value_policy<return_by_value> rbv;
      class_<w_t>("matrix_cholesky_gill_murray_wright_wrappers", no_init)
        .def(init<ref<double> const&>(arg_("u")))
        .add_property("e", make_getter(&w_t::e, rbv()))
        .add_property("pivots", make_getter(&w_t::pivots, rbv()))
      ;
    }

    static w_t
    factory(ref<double> const& u) { return w_t(u); }
  };

  bool
  is_square_matrix(
    af::versa<double, af::flex_grid<> > const& self)
  {
    return self.accessor().is_square_matrix();
  }

} // namespace <anonymous>

namespace boost_python {

  BOOST_PYTHON_FUNCTION_OVERLOADS(
    matrix_symmetric_as_packed_u_overloads,
    matrix::symmetric_as_packed_u, 1, 2)
  BOOST_PYTHON_FUNCTION_OVERLOADS(
    matrix_symmetric_as_packed_l_overloads,
    matrix::symmetric_as_packed_l, 1, 2)

  BOOST_PYTHON_FUNCTION_OVERLOADS(
    matrix_cholesky_decomposition_overloads,
    matrix::cholesky::decomposition, 1, 2)

  void
  wrap_flex_double_matrix(
    flex_wrapper<double>::class_f_t& class_f_t)
  {
    exercise_packed_u_accessor();
    matrix_cholesky_gill_murray_wright_wrappers::wrap();

    using namespace boost::python;

    class_f_t
      .def("is_square_matrix", is_square_matrix)
      .def("matrix_diagonal",
        (shared<double>(*)(
          const_ref<double, c_grid<2> > const&)) matrix_diagonal)
      .def("matrix_diagonal_set_in_place",
        (void(*)(
          ref<double, c_grid<2> > const&,
          double const&)) matrix_diagonal_set_in_place,
            arg_("value"))
      .def("matrix_diagonal_add_in_place",
        (void(*)(
          ref<double, c_grid<2> > const&,
          double const&)) matrix_diagonal_add_in_place,
            arg_("value"))
      .def("matrix_diagonal_sum",
        (double(*)(
          const_ref<double, c_grid<2> > const&)) matrix_diagonal_sum)
      .def("matrix_trace",
        (double(*)(
          const_ref<double, c_grid<2> > const&)) matrix_diagonal_sum)
      .def("matrix_diagonal_product",
        (double(*)(
          const_ref<double, c_grid<2> > const&)) matrix_diagonal_product)
      .def("matrix_multiply",
        (versa<double, c_grid<2> >(*)(
          const_ref<double, c_grid<2> > const&,
          const_ref<double, c_grid<2> > const&)) matrix_multiply)
      .def("matrix_multiply",
        (shared<double>(*)(
          const_ref<double, c_grid<2> > const&,
          const_ref<double> const&)) matrix_multiply)
      .def("matrix_multiply",
        (shared<double>(*)(
          const_ref<double> const&,
          const_ref<double, c_grid<2> > const&)) matrix_multiply)
      .def("matrix_multiply",
        (double(*)(
          const_ref<double> const&,
          const_ref<double> const&)) matrix_multiply)
      .def("dot",
        (double(*)(
          const_ref<double> const&,
          const_ref<double> const&)) matrix_multiply)
      .def("matrix_transpose_multiply_as_packed_u",
        (shared<double>(*)(
          const_ref<double, c_grid<2> > const&))
            matrix_transpose_multiply_as_packed_u)
      .def("matrix_transpose",
        (versa<double, c_grid<2> >(*)(
           const_ref<double, c_grid<2> > const&)) matrix_transpose)
      .def("matrix_transpose_in_place",
        (void(*)(versa<double, flex_grid<> >&)) matrix_transpose_in_place)
      .def("matrix_outer_product",
        (versa<double, c_grid<2> >(*)(
           const_ref<double> const&,
           const_ref<double> const&)) matrix::outer_product, (arg_("rhs")))
      .def("matrix_lu_decomposition_in_place",
        (shared<std::size_t>(*)(
          ref<double, c_grid<2> > const&)) matrix_lu_decomposition_in_place)
      .def("matrix_lu_back_substitution",
        (shared<double>(*)(
          const_ref<double, c_grid<2> > const&,
          const_ref<std::size_t> const&,
          const_ref<double> const&)) matrix_lu_back_substitution, (
        arg_("pivot_indices"), arg_("b")))
      .def("matrix_determinant_via_lu",
        (double(*)(
          const_ref<double, c_grid<2> > const&,
          const_ref<std::size_t> const&)) matrix_determinant_via_lu, (
        arg_("pivot_indices")))
      .def("matrix_determinant_via_lu",
        (double(*)(
          const_ref<double, c_grid<2> > const&)) matrix_determinant_via_lu)
      .def("matrix_inversion_in_place",
        (void(*)(
          ref<double, c_grid<2> > const&,
          ref<double, c_grid<2> > const&)) matrix_inversion_in_place, (
        arg_("b")))
      .def("matrix_inversion_in_place",
        (void(*)(ref<double, c_grid<2> > const&)) matrix_inversion_in_place)
      .def("matrix_upper_triangle_as_packed_u",
        (af::shared<double>(*)(
          const_ref<double, c_grid<2> > const&))
            matrix::upper_triangle_as_packed_u)
      .def("matrix_packed_u_as_upper_triangle",
        (af::versa<double, c_grid<2> >(*)(
          const_ref<double> const&))
            matrix::packed_u_as_upper_triangle)
      .def("matrix_lower_triangle_as_packed_l",
        (af::shared<double>(*)(
          const_ref<double, c_grid<2> > const&))
            matrix::lower_triangle_as_packed_l)
      .def("matrix_packed_l_as_lower_triangle",
        (af::versa<double, c_grid<2> >(*)(
          const_ref<double> const&))
            matrix::packed_l_as_lower_triangle)
      .def("matrix_symmetric_as_packed_u",
        (af::shared<double>(*)(
          const_ref<double, af::c_grid<2> > const&, double const&))
            matrix::symmetric_as_packed_u,
              matrix_symmetric_as_packed_u_overloads((
                arg_("self"),
                arg_("relative_epsilon")=1.e-15)))
      .def("matrix_symmetric_as_packed_l",
        (af::shared<double>(*)(
          const_ref<double, af::c_grid<2> > const&, double const&))
            matrix::symmetric_as_packed_l,
              matrix_symmetric_as_packed_l_overloads((
                arg_("self"),
                arg_("relative_epsilon")=1.e-15)))
      .def("matrix_packed_u_as_symmetric",
        (af::versa<double, c_grid<2> >(*)(
          const_ref<double> const&))
            matrix::packed_u_as_symmetric)
      .def("matrix_packed_l_as_symmetric",
        (af::versa<double, c_grid<2> >(*)(
          const_ref<double> const&))
            matrix::packed_l_as_symmetric)
      .def("matrix_cholesky_decomposition",
        (af::shared<double>(*)(
          const_ref<double> const&, double const&))
            matrix::cholesky::decomposition,
              matrix_cholesky_decomposition_overloads((
                arg_("self"),
                arg_("relative_epsilon")=1.e-15)))
      .def("matrix_cholesky_gill_murray_wright",
        matrix_cholesky_gill_murray_wright_wrappers::factory)
      .def("matrix_copy_upper_to_lower_triangle_in_place",
        (void(*)(
          ref<double, c_grid<2> > const&))
            matrix::copy_upper_to_lower_triangle_in_place)
      .def("matrix_copy_lower_to_upper_triangle_in_place",
        (void(*)(
          ref<double, c_grid<2> > const&))
            matrix::copy_lower_to_upper_triangle_in_place)
      .def("matrix_swap_rows_in_place",
        (void(*)(
          ref<double, c_grid<2> > const&, unsigned, unsigned))
            matrix::swap_rows_in_place, (
                arg_("self"),
                arg_("i"),
                arg_("j")))
      .def("matrix_swap_columns_in_place",
        (void(*)(
          ref<double, c_grid<2> > const&, unsigned, unsigned))
            matrix::swap_columns_in_place, (
                arg_("self"),
                arg_("i"),
                arg_("j")))
      .def("matrix_symmetric_upper_triangle_swap_rows_and_columns_in_place",
        (void(*)(
          ref<double, c_grid<2> > const&, unsigned, unsigned))
            matrix::symmetric_upper_triangle_swap_rows_and_columns_in_place, (
                arg_("self"),
                arg_("i"),
                arg_("j")))
      .def("matrix_packed_u_swap_rows_and_columns_in_place",
        (void(*)(
          ref<double> const&, unsigned, unsigned))
            matrix::packed_u_swap_rows_and_columns_in_place, (
                arg_("self"),
                arg_("i"),
                arg_("j")))
      .def("cos_angle",
        (boost::optional<double>(*)(
          const_ref<double> const&,
          const_ref<double> const&)) cos_angle, (
        arg_("b")))
      .def("cos_angle",
        (double(*)(
          const_ref<double> const&,
          const_ref<double> const&,
          const double&)) cos_angle, (
        arg_("b"), arg_("value_if_undefined")))
      .def("angle",
        (boost::optional<double>(*)(
          const_ref<double> const&,
          const_ref<double> const&)) angle, (
        arg_("b")))
      .def("angle",
        (boost::optional<double>(*)(
          const_ref<double> const&,
          const_ref<double> const&,
          bool)) angle, (
        arg_("b"), arg_("deg")=false))
    ;
  }

}}} // namespace scitbx::af::boost_python
