/* Copyright (c) 2001-2002 The Regents of the University of California
   through E.O. Lawrence Berkeley National Laboratory, subject to
   approval by the U.S. Department of Energy.
   See files COPYRIGHT.txt and LICENSE.txt for further details.

   Revision history:
     2002 Dec: Created (rwgk)
 */

#include <cctbx/boost_python/flex_fwd.h>

#include <cctbx/sgtbx/row_echelon.h>
#include <boost/python/def.hpp>
#include <boost/python/class.hpp> // work around gcc-3.3-darwin bug

namespace cctbx { namespace sgtbx { namespace boost_python {

namespace {

  template <typename ElementType>
  scitbx::mat_ref<ElementType>
  flex_as_mat_ref(af::versa<ElementType, af::flex_grid<> >& a)
  {
    CCTBX_ASSERT(a.accessor().nd() == 2);
    CCTBX_ASSERT(a.accessor().is_0_based());
    CCTBX_ASSERT(!a.accessor().is_padded());
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
  row_echelon_back_substitution(
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
      CCTBX_ASSERT(v.size() == re_mx_ref.n_rows());
      v_ptr = v.begin();
    }
    if (sol.size()) {
      CCTBX_ASSERT(sol.size() == re_mx_ref.n_columns());
      sol_ptr = sol.begin();
    }
    if (indep.size()) {
      CCTBX_ASSERT(indep.size() == re_mx_ref.n_columns());
      indep_ptr = indep.begin();
    }
    return row_echelon::back_substitution(
      re_mx_ref, v_ptr, sol_ptr, indep_ptr);
  }

  struct dummy {}; // work around gcc-3.3-darwin bug

} // namespace <anoymous>

  void wrap_row_echelon()
  {
    using namespace boost::python;
#if defined(__APPLE__) && defined(__MACH__) \
 && defined(__GNUC__) && __GNUC__ == 3 && __GNUC_MINOR__ == 3
    class_<dummy>("_dummy", no_init);
#endif
    def("row_echelon_form_t", row_echelon_form_t);
    def("row_echelon_form", row_echelon_form);
    def("row_echelon_back_substitution", row_echelon_back_substitution);
  }

}}} // namespace cctbx::sgtbx::boost_python
