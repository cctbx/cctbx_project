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

} // namespace <anoymous>

  void wrap_row_echelon()
  {
    using namespace boost::python;
    def("row_echelon_form_t", row_echelon_form_t);
    def("row_echelon_form", row_echelon_form);
  }

}}} // namespace cctbx::sgtbx::boost_python
