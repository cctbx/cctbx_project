#include <scitbx/array_family/boost_python/flex_fwd.h>

#include <scitbx/matrix/row_echelon.h>
#include <boost/python/def.hpp>

namespace scitbx { namespace matrix { namespace boost_python {

namespace {

  template <typename ElementType>
  scitbx::mat_ref<ElementType>
  flex_as_mat_ref(af::versa<ElementType, af::flex_grid<> >& a)
  {
    SCITBX_ASSERT(a.accessor().nd() == 2);
    SCITBX_ASSERT(a.accessor().is_0_based());
    SCITBX_ASSERT(!a.accessor().is_padded());
    return scitbx::mat_ref<ElementType>(
      a.begin(),
      a.accessor().all()[0],
      a.accessor().all()[1]);
  }

  std::size_t
  row_echelon_form_t(
    af::versa<int, af::flex_grid<> >& m,
    af::versa<int, af::flex_grid<> >& t)
  {
    scitbx::mat_ref<int> m_ref = flex_as_mat_ref(m);
    scitbx::mat_ref<int> t_ref = flex_as_mat_ref(t);
    std::size_t rank = row_echelon::form_t(m_ref, t_ref);
    m.resize(af::flex_grid<>(m_ref.n_rows(), m_ref.n_columns()));
    return rank;
  }

  std::size_t
  row_echelon_form(
    af::versa<int, af::flex_grid<> >& m)
  {
    scitbx::mat_ref<int> m_ref = flex_as_mat_ref(m);
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
    scitbx::mat_ref<int> re_mx_ref = flex_as_mat_ref(re_mx);
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
    scitbx::mat_ref<int> re_mx_ref = flex_as_mat_ref(re_mx);
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

}}} // namespace matrix::boost_python::<anonymous>

namespace math { namespace boost_python {

  void wrap_row_echelon()
  {
    matrix::boost_python::wrap_row_echelon();
  }

}}} // namespace scitbx::math::boost_python
