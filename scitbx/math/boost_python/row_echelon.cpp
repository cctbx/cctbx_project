#include <scitbx/array_family/boost_python/flex_fwd.h>

#include <scitbx/matrix/row_echelon.h>
#include <scitbx/matrix/row_echelon_full_pivoting.h>
#include <scitbx/array_family/versa.h>
#include <boost/python/class.hpp>
#include <boost/python/def.hpp>
#include <boost/python/return_value_policy.hpp>
#include <boost/python/return_by_value.hpp>

namespace scitbx { namespace matrix { namespace boost_python {

namespace {

  template <typename ElementType>
  af::ref<ElementType, af::mat_grid>
  flex_as_mat_ref(af::versa<ElementType, af::flex_grid<> >& a)
  {
    SCITBX_ASSERT(a.accessor().nd() == 2);
    SCITBX_ASSERT(a.accessor().is_0_based());
    SCITBX_ASSERT(!a.accessor().is_padded());
    return af::ref<ElementType, af::mat_grid>(
      a.begin(),
      a.accessor().all()[0],
      a.accessor().all()[1]);
  }

  std::size_t
  row_echelon_form_t(
    af::versa<int, af::flex_grid<> >& m,
    af::versa<int, af::flex_grid<> >& t)
  {
    af::ref<int, af::mat_grid> m_ref = flex_as_mat_ref(m);
    af::ref<int, af::mat_grid> t_ref = flex_as_mat_ref(t);
    std::size_t rank = row_echelon::form_t(m_ref, t_ref);
    m.resize(af::flex_grid<>(m_ref.n_rows(), m_ref.n_columns()));
    return rank;
  }

  std::size_t
  row_echelon_form(
    af::versa<int, af::flex_grid<> >& m)
  {
    af::ref<int, af::mat_grid> m_ref = flex_as_mat_ref(m);
    std::size_t rank = row_echelon::form(m_ref);
    m.resize(af::flex_grid<>(m_ref.n_rows(), m_ref.n_columns()));
    return rank;
  }

  int
  row_echelon_back_substitution_int(
    af::versa<int, af::flex_grid<> >& re_mx,
    af::const_ref<int> const& v,
    af::ref<int> const& sol,
    af::ref<bool> const& indep)
  {
    af::ref<int, af::mat_grid> re_mx_ref = flex_as_mat_ref(re_mx);
    const int* v_ptr = 0;
    int* sol_ptr = 0;
    bool* indep_ptr = 0;
    if (v.size()) {
      SCITBX_ASSERT(v.size() == re_mx_ref.n_rows());
      v_ptr = v.begin();
    }
    if (sol.size()) {
      SCITBX_ASSERT(sol.size() == re_mx_ref.n_columns());
      sol_ptr = sol.begin();
    }
    if (indep.size()) {
      SCITBX_ASSERT(indep.size() == re_mx_ref.n_columns());
      indep_ptr = indep.begin();
    }
    return row_echelon::back_substitution_int(
      re_mx_ref, v_ptr, sol_ptr, indep_ptr);
  }

  bool
  row_echelon_back_substitution_float(
    af::versa<int, af::flex_grid<> >& re_mx,
    af::const_ref<double> const& v,
    af::ref<double> const& sol)
  {
    af::ref<int, af::mat_grid> re_mx_ref = flex_as_mat_ref(re_mx);
    const double* v_ptr = 0;
    double* sol_ptr = 0;
    if (v.size()) {
      SCITBX_ASSERT(v.size() == re_mx_ref.n_rows());
      v_ptr = v.begin();
    }
    if (sol.size()) {
      SCITBX_ASSERT(sol.size() == re_mx_ref.n_columns());
      sol_ptr = sol.begin();
    }
    return row_echelon::back_substitution_float(
      re_mx_ref, v_ptr, sol_ptr);
  }

  void wrap_row_echelon()
  {
    using namespace boost::python;
    def("row_echelon_form_t", row_echelon_form_t);
    def("row_echelon_form", row_echelon_form);
    def("row_echelon_back_substitution_int",
      row_echelon_back_substitution_int);
    def("row_echelon_back_substitution_float",
      row_echelon_back_substitution_float);
  }

  struct full_pivoting_wrapper
  {
    typedef row_echelon::full_pivoting<double> wt;

    static void wrap()
    {
      using namespace boost::python;
      typedef return_value_policy<return_by_value> rbv;
      class_<wt>("row_echelon_full_pivoting", no_init)
        .def(init<
          af::versa<double, af::flex_grid<> >,
          double const&,
          int>((
            arg("a_work"),
            arg("min_abs_pivot")=0,
            arg("max_rank")=-1)))
        .def(init<
          af::versa<double, af::flex_grid<> >,
          af::shared<double>,
          double const&,
          int>((
            arg("a_work"),
            arg("b_work"),
            arg("min_abs_pivot")=0,
            arg("max_rank")=-1)))
        .add_property("col_perm", make_getter(&wt::col_perm, rbv()))
        .def_readonly("rank", &wt::rank)
        .def_readonly("nullity", &wt::nullity)
        .def("is_in_row_space", &wt::is_in_row_space, (
          arg("x"),
          arg("epsilon")))
        .def("back_substitution", &wt::back_substitution, (
          arg("free_values"),
          arg("epsilon")=0))
      ;
    }
  };

}}} // namespace matrix::boost_python::<anonymous>

namespace math { namespace boost_python {

  void wrap_row_echelon()
  {
    matrix::boost_python::wrap_row_echelon();
    matrix::boost_python::full_pivoting_wrapper::wrap();
  }

}}} // namespace scitbx::math::boost_python
