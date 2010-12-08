#include <scitbx/array_family/boost_python/flex_wrapper.h>
#include <scitbx/array_family/versa_matrix.h>
#include <scitbx/matrix/outer_product.h>
#include <scitbx/matrix/norms.h>
#include <scitbx/matrix/move.h>
#include <scitbx/matrix/matrix_vector_operations.h>
#include <scitbx/array_family/boost_python/packed_to_flex_conversions.h>
#include <boost/python/args.hpp>
#include <boost/python/return_value_policy.hpp>
#include <boost/python/return_by_value.hpp>
#include <boost_adaptbx/std_pair_conversion.h>

namespace scitbx { namespace af {

namespace {

  void
  exercise_packed_u_accessor()
  {
    versa<double, packed_u_accessor>
      a(packed_u_accessor(5));
    ref<double, packed_u_accessor> r = a.ref();
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

  versa<double, packed_u_accessor> exercise_versa_packed_u_to_flex() {
    versa<double, packed_u_accessor> result(3);
    for (int i=0; i<3; ++i) for (int j=i; j<3; ++j) result(i,j) = 10*(i+1) + j+1;
    return result;
  }


  bool
  is_square_matrix(
    versa<double, flex_grid<> > const& self)
  {
    return self.accessor().is_square_matrix();
  }

} // namespace <anonymous>

namespace boost_python {

  versa<double, c_grid<2> >
  matrix_multiply_real_matrix_real_matrix(
    const_ref<double, c_grid<2> > const& a,
    const_ref<double, c_grid<2> > const& b)
  {
    return matrix_multiply(a, b);
  }

  versa<std::complex<double>, c_grid<2> >
  matrix_multiply_real_matrix_complex_matrix(
    const_ref<double, c_grid<2> > const& a,
    const_ref<std::complex<double>, c_grid<2> > const& b)
  {
    return matrix_multiply(a, b);
  }

  versa<double, c_grid<2> >
  matrix_transpose_multiply_real_matrix_real_matrix(
    const_ref<double, c_grid<2> > const& a,
    const_ref<double, c_grid<2> > const& b)
  {
    return matrix_transpose_multiply(a, b);
  }

  versa<double, c_grid<2> >
  matrix_multiply_transpose_real_matrix_real_matrix(
    const_ref<double, c_grid<2> > const& a,
    const_ref<double, c_grid<2> > const& b)
  {
    return matrix_multiply_transpose(a, b);
  }

  versa<double, c_grid<2> >
  matrix_multiply_packed_u_real_matrix_real_u(
    const_ref<double, c_grid<2> > const& a,
    const_ref<double> const& b)
  {
    return matrix_multiply_packed_u(a, b);
  }

  versa<std::complex<double>, c_grid<2> >
  matrix_multiply_packed_u_real_matrix_complex_u(
    const_ref<double, c_grid<2> > const& a,
    const_ref<std::complex<double> > const& b)
  {
    return matrix_multiply_packed_u(a, b);
  }

  shared<double>
  matrix_multiply_packed_u_multiply_lhs_transpose_real_matrix_real_u(
    const_ref<double, c_grid<2> > const& a,
    const_ref<double> const& b)
  {
    return matrix_multiply_packed_u_multiply_lhs_transpose(a, b);
  }

  shared<std::complex<double> >
  matrix_multiply_packed_u_multiply_lhs_transpose_real_matrix_complex_u(
    const_ref<double, c_grid<2> > const& a,
    const_ref<std::complex<double> > const& b)
  {
    return matrix_multiply_packed_u_multiply_lhs_transpose(a, b);
  }

  double (*matrix_norm_1)(const_ref<double, mat_grid> const &) = matrix::norm_1;

  double matrix_symmetric_upper_triangle_quadratic_form(
    const_ref<double, packed_u_accessor> const &q,
    const_ref<double> const &x)
  {
    SCITBX_ASSERT(q.n_columns() == x.size());
    return matrix::quadratic_form_packed_u(x.size(), q.begin(), x.begin());
  }


  void
  wrap_flex_double_matrix(
    flex_wrapper<double>::class_f_t& class_f_t)
  {
    exercise_packed_u_accessor();
    default_packed_flex_conversions<double>();
    boost::python::def("exercise_versa_packed_u_to_flex",
                       exercise_versa_packed_u_to_flex);

    using namespace boost::python;
    using boost::python::arg;

    boost_adaptbx::std_pair_conversions::to_tuple<shared<double>,
                                                  shared<double> >();

    class_f_t
      .def("is_square_matrix", is_square_matrix)
      .def("matrix_diagonal",
        (shared<double>(*)(
          const_ref<double, c_grid<2> > const&)) matrix_diagonal)
      .def("matrix_upper_bidiagonal", matrix_upper_bidiagonal<double>)
      .def("matrix_lower_bidiagonal", matrix_lower_bidiagonal<double>)
      .def("matrix_diagonal_set_in_place",
        (void(*)(
          ref<double, c_grid<2> > const&,
          double const&)) matrix_diagonal_set_in_place,
            arg("value"))
      .def("matrix_diagonal_set_in_place",
        (void(*)(
          ref<double, c_grid<2> > const&,
          const_ref<double> const&)) matrix_diagonal_set_in_place,
            arg("diagonal"))
      .def("matrix_diagonal_add_in_place",
        (void(*)(
          ref<double, c_grid<2> > const&,
          double const&)) matrix_diagonal_add_in_place,
            arg("value"))
      .def("matrix_diagonal_sum",
        (double(*)(
          const_ref<double, c_grid<2> > const&)) matrix_diagonal_sum)
      .def("matrix_trace",
        (double(*)(
          const_ref<double, c_grid<2> > const&)) matrix_diagonal_sum)
      .def("matrix_diagonal_product",
        (double(*)(
          const_ref<double, c_grid<2> > const&)) matrix_diagonal_product)
      .def("matrix_norm_1", matrix_norm_1)
      .def("matrix_norm_inf", matrix::norm_inf<double>)
      .def("matrix_norm_frobenius", matrix::norm_frobenius<double>)
      .def("matrix_multiply", matrix_multiply_real_matrix_real_matrix)
      .def("matrix_multiply", matrix_multiply_real_matrix_complex_matrix)
      .def("matrix_transpose_multiply",
        matrix_transpose_multiply_real_matrix_real_matrix)
      .def("matrix_multiply_transpose",
        matrix_multiply_transpose_real_matrix_real_matrix)
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
      .def("matrix_multiply_packed_u",
        matrix_multiply_packed_u_real_matrix_real_u)
      .def("matrix_multiply_packed_u",
        matrix_multiply_packed_u_real_matrix_complex_u)
      .def("matrix_multiply_packed_u_multiply_lhs_transpose",
        matrix_multiply_packed_u_multiply_lhs_transpose_real_matrix_real_u, (
          arg("packed_u")))
      .def("matrix_multiply_packed_u_multiply_lhs_transpose",
        matrix_multiply_packed_u_multiply_lhs_transpose_real_matrix_complex_u,(
          arg("packed_u")))
      .def("matrix_transpose_multiply_as_packed_u",
        (shared<double>(*)(
          const_ref<double, c_grid<2> > const&))
            matrix_transpose_multiply_as_packed_u)
      .def("matrix_transpose_multiply_diagonal_multiply_as_packed_u",
        (shared<double>(*)(
          const_ref<double, c_grid<2> > const&,
          const_ref<double> const&))
            matrix_transpose_multiply_diagonal_multiply_as_packed_u, (
              arg("diagonal_elements")))
      .def("matrix_transpose",
        (versa<double, c_grid<2> >(*)(
           const_ref<double, c_grid<2> > const&)) matrix_transpose)
      .def("matrix_transpose_in_place",
        (void(*)(versa<double, flex_grid<> >&)) matrix_transpose_in_place)
      .def("matrix_outer_product",
        (versa<double, c_grid<2> >(*)(
           const_ref<double> const&,
           const_ref<double> const&)) matrix::outer_product, (arg("rhs")))
      .def("matrix_lu_decomposition_in_place",
        (shared<std::size_t>(*)(
          ref<double, c_grid<2> > const&)) matrix_lu_decomposition_in_place)
      .def("matrix_lu_back_substitution",
        (shared<double>(*)(
          const_ref<double, c_grid<2> > const&,
          const_ref<std::size_t> const&,
          const_ref<double> const&)) matrix_lu_back_substitution, (
        arg("pivot_indices"), arg("b")))
      .def("matrix_forward_substitution", matrix_forward_substitution<double>,
           (arg("l"), arg("b"), arg("unit_diag")=false))
      .def("matrix_back_substitution", matrix_back_substitution<double>,
           (arg("u"), arg("b"), arg("unit_diag")=false))
      .def("matrix_forward_substitution_given_transpose",
           matrix_forward_substitution_given_transpose<double>,
           (arg("u"), arg("b"), arg("unit_diag")=false))
      .def("matrix_back_substitution_given_transpose",
           matrix_back_substitution_given_transpose<double>,
           (arg("l"), arg("b"), arg("unit_diag")=false))
      .def("matrix_determinant_via_lu",
        (double(*)(
          const_ref<double, c_grid<2> > const&,
          const_ref<std::size_t> const&)) matrix_determinant_via_lu, (
        arg("pivot_indices")))
      .def("matrix_determinant_via_lu",
        (double(*)(
          const_ref<double, c_grid<2> > const&)) matrix_determinant_via_lu)
      .def("matrix_inversion_in_place",
        (void(*)(
          ref<double, c_grid<2> > const&,
          ref<double, c_grid<2> > const&)) matrix_inversion_in_place, (
        arg("b")))
      .def("matrix_inversion_in_place",
        (void(*)(ref<double, c_grid<2> > const&)) matrix_inversion_in_place)
      .def("matrix_upper_triangle_as_packed_u",
        (shared<double>(*)(
          const_ref<double, c_grid<2> > const&))
            matrix::upper_triangle_as_packed_u)
      .def("matrix_packed_u_as_upper_triangle",
        (versa<double, c_grid<2> >(*)(
          const_ref<double> const&))
            matrix::packed_u_as_upper_triangle)
      .def("matrix_lower_triangle_as_packed_l",
        (shared<double>(*)(
          const_ref<double, c_grid<2> > const&))
            matrix::lower_triangle_as_packed_l)
      .def("matrix_packed_l_as_lower_triangle",
        (versa<double, c_grid<2> >(*)(
          const_ref<double> const&))
            matrix::packed_l_as_lower_triangle)
      .def("matrix_symmetric_as_packed_u",
        (shared<double>(*)(
          const_ref<double, c_grid<2> > const&, double const&))
            matrix::symmetric_as_packed_u, (
              arg("relative_epsilon")=1e-12))
      .def("matrix_symmetric_as_packed_l",
        (shared<double>(*)(
          const_ref<double, c_grid<2> > const&, double const&))
            matrix::symmetric_as_packed_l, (
              arg("relative_epsilon")=1e-12))
      .def("matrix_is_symmetric",
        (bool(*)(
          const_ref<double, c_grid<2> > const&, double const&))
            matrix::is_symmetric, ((
              arg("relative_epsilon"))))
      .def("matrix_packed_u_as_symmetric",
        (versa<double, c_grid<2> >(*)(
          const_ref<double> const&))
            matrix::packed_u_as_symmetric)
      .def("matrix_packed_l_as_symmetric",
        (versa<double, c_grid<2> >(*)(
          const_ref<double> const&))
            matrix::packed_l_as_symmetric)
      .def("matrix_packed_u_diagonal",
        (shared<double>(*)(
          const_ref<double> const&))
            matrix::packed_u_diagonal)
      .def("matrix_packed_u_diagonal_add_in_place",
           (void (*)(ref<double> const &, double))
            matrix::packed_u_diagonal_add_in_place)
      .def("matrix_copy_upper_to_lower_triangle_in_place",
        (void(*)(
          ref<double, c_grid<2> > const&))
            matrix::copy_upper_to_lower_triangle_in_place)
      .def("matrix_copy_lower_to_upper_triangle_in_place",
        (void(*)(
          ref<double, c_grid<2> > const&))
            matrix::copy_lower_to_upper_triangle_in_place)
      .def("matrix_copy_column",
        (shared<double>(*)(
          const_ref<double, c_grid<2> > const&,
          unsigned))
            matrix::copy_column, (
              arg("i_column")))
      .def("matrix_copy_block",
        (versa<double, c_grid<2> >(*)(
          const_ref<double, c_grid<2> > const&,
          unsigned, unsigned, unsigned, unsigned))
            matrix::copy_block, (
              arg("i_row"),
              arg("i_column"),
              arg("n_rows"),
              arg("n_columns")))
      .def("matrix_paste_block_in_place",
        (void(*)(
          ref<double, c_grid<2> > const&,
          const_ref<double, c_grid<2> > const&,
          unsigned, unsigned))
            matrix::paste_block_in_place, (
              arg("block"),
              arg("i_row"),
              arg("i_column")))
      .def("matrix_paste_column_in_place",
           matrix::paste_column_in_place<double>,
           (arg("column"), "i_column"))
      .def("matrix_copy_upper_triangle", matrix::copy_upper_triangle<double>)
      .def("matrix_copy_lower_triangle", matrix::copy_lower_triangle<double>)
      .def("matrix_swap_rows_in_place",
        (void(*)(
          ref<double, c_grid<2> > const&, unsigned, unsigned))
            matrix::swap_rows_in_place, (
              arg("i"),
              arg("j")))
      .def("matrix_swap_columns_in_place",
        (void(*)(
          ref<double, c_grid<2> > const&, unsigned, unsigned))
            matrix::swap_columns_in_place, (
              arg("i"),
              arg("j")))
      .def("matrix_symmetric_upper_triangle_swap_rows_and_columns_in_place",
        (void(*)(
          ref<double, c_grid<2> > const&, unsigned, unsigned))
            matrix::symmetric_upper_triangle_swap_rows_and_columns_in_place, (
              arg("i"),
              arg("j")))
      .def("matrix_symmetric_upper_triangle_quadratic_form",
           matrix_symmetric_upper_triangle_quadratic_form)
      .def("matrix_packed_u_swap_rows_and_columns_in_place",
        (void(*)(
          ref<double> const&, unsigned, unsigned))
            matrix::packed_u_swap_rows_and_columns_in_place, (
              arg("i"),
              arg("j")))
      .def("cos_angle",
        (boost::optional<double>(*)(
          const_ref<double> const&,
          const_ref<double> const&)) cos_angle, (
        arg("b")))
      .def("cos_angle",
        (double(*)(
          const_ref<double> const&,
          const_ref<double> const&,
          const double&)) cos_angle, (
        arg("b"), arg("value_if_undefined")))
      .def("angle",
        (boost::optional<double>(*)(
          const_ref<double> const&,
          const_ref<double> const&)) angle, (
        arg("b")))
      .def("angle",
        (boost::optional<double>(*)(
          const_ref<double> const&,
          const_ref<double> const&,
          bool)) angle, (
        arg("b"), arg("deg")=false))
    ;
  }

}}} // namespace scitbx::af::boost_python
