/* Copyright (c) 2001-2002 The Regents of the University of California
   through E.O. Lawrence Berkeley National Laboratory, subject to
   approval by the U.S. Department of Energy.
   See files COPYRIGHT.txt and LICENSE.txt for further details.

   Revision history:
     2002 Aug: Copied from cctbx/arraytbx/flexmodule.cpp (rwgk)
     2002 Aug: Created, based on sharedmodule.cpp, shared_bpl.h (rwgk)
 */

#include <scitbx/array_family/accessors/flex_grid.h>
#include <boost/python/tuple.hpp>
#include <boost/python/class.hpp>
#include <boost/python/overloads.hpp>
#include <boost/python/return_value_policy.hpp>
#include <boost/python/copy_const_reference.hpp>

namespace scitbx { namespace af { namespace boost_python {

namespace {

  struct flex_grid_wrappers : boost::python::pickle_suite
  {
    typedef flex_grid_default_index_type df_i_t;
    typedef flex_grid<> w_t;
    typedef flex_grid_default_index_type::value_type ivt;
    typedef boost::python::class_<w_t> c_w_t;

    BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(
      set_layout_tuple_overloads, set_layout, 1, 2)

    BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(
      set_layout_convenience_overloads, set_layout, 1, 6)

    BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(
      last_overloads, last, 0, 1)

    BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(
      layout_overloads, layout, 0, 1)

    static boost::python::tuple
    getinitargs(w_t const& fg)
    {
      bool open_range = true;
      return boost::python::make_tuple(
        fg.origin(),
        fg.last(open_range),
        open_range);
    }

    static df_i_t
    getstate(w_t const& fg)
    {
      return fg.layout();
    }

    static void
    setstate(w_t& fg, df_i_t const& state)
    {
      fg.set_layout(state);
    }

    static void
    wrap()
    {
      using namespace boost::python;
      typedef return_value_policy<copy_const_reference> copy_const_reference;
      c_w_t("grid")
        .def(init<df_i_t const&>())
        .def(init<ivt const&, optional<ivt const&, ivt const&,
                  ivt const&, ivt const&, ivt const&> >())
        .def(init<df_i_t const&, df_i_t const&, optional<bool> >())
        .def("set_layout",
          (w_t(w_t::*)(df_i_t const&, bool)) 0,
            set_layout_tuple_overloads())
        .def("set_layout",
          (w_t(w_t::*)(ivt const&, ivt const&, ivt const&,
                       ivt const&, ivt const&, ivt const&)) 0,
            set_layout_convenience_overloads())
        .def("nd", &w_t::nd)
        .def("size_1d", &w_t::size_1d)
        .def("has_origin", &w_t::has_origin)
        .def("origin", &w_t::origin)
        .def("grid", &w_t::grid, copy_const_reference())
        .def("last", (df_i_t(w_t::*)(bool)) 0, last_overloads())
        .def("has_layout", &w_t::has_layout)
        .def("layout", (df_i_t(w_t::*)(bool)) 0, layout_overloads())
        .def("layout_size_1d", &w_t::layout_size_1d)
        .def("is_0_based", &w_t::is_0_based)
        .def("is_padded", &w_t::is_padded)
        .def("shift_origin", &w_t::shift_origin)
        .def("is_valid_index", &w_t::is_valid_index)
        .def("__call__", &w_t::operator())
        .def("__eq__", &w_t::operator==)
        .def("__ne__", &w_t::operator!=)
        .def_pickle(flex_grid_wrappers())
      ;
    }
  };

} // namespace <anoymous>

  void wrap_flex_grid()
  {
    flex_grid_wrappers::wrap();
  }

}}} // namespace scitbx::af::boost_python
