/* Copyright (c) 2001-2002 The Regents of the University of California
   through E.O. Lawrence Berkeley National Laboratory, subject to
   approval by the U.S. Department of Energy.
   See files COPYRIGHT.txt and LICENSE.txt for further details.

   Revision history:
     2002 Aug: Copied from cctbx/arraytbx/flexmodule.cpp (rwgk)
     2002 Aug: Created, based on sharedmodule.cpp, shared_bpl.h (rwgk)
 */

#include <scitbx/array_family/flex_grid_accessor.h>
#include <scitbx/boost_python/utils.h>
#include <boost/python/module.hpp>
#include <boost/python/class.hpp>
#include <scitbx/array_family/boost_python/small_conversions.h>
#include <boost/python/return_value_policy.hpp>
#include <boost/python/copy_const_reference.hpp>

namespace scitbx { namespace af { namespace boost_python {

  //BOOST_PYTHON_MEM_FUN_GENERATOR(flex_grid_last_stubs, last, 1, 1)

  struct flex_grid_wrappers : boost::python::pickle_suite
  {
    static flex_grid<>
    set_layout(
      flex_grid<>& fg,
      flex_grid_default_index_type const& layout)
    {
      return fg.set_layout(layout);
    }

    static flex_grid_default_index_type
    last_0(flex_grid<> const& fg)
    {
      return fg.last();
    }

    static flex_grid_default_index_type
    last_1(flex_grid<> const& fg, bool open_range)
    {
      return fg.last(open_range);
    }

    static boost::python::tuple
    getinitargs(flex_grid<> const& fg)
    {
      bool open_range = true;
      return boost::python::make_tuple(
        fg.origin(),
        fg.last(open_range),
        open_range);
    }

    static flex_grid_default_index_type
    getstate(flex_grid<> const& fg)
    {
      return fg.layout();
    }

    static void
    setstate(flex_grid<>& fg, flex_grid_default_index_type const& state)
    {
      fg.set_layout(state);
    }

  };

  void add_flex_bool(boost::python::module& m);
  void add_flex_size_t(boost::python::module& m);
  void add_flex_int(boost::python::module& m);
  void add_flex_long(boost::python::module& m);
  void add_flex_float(boost::python::module& m);
  void add_flex_double(boost::python::module& m);
  void add_flex_complex_double(boost::python::module& m);
  void add_flex_std_string(boost::python::module& m);

  void init_module(boost::python::module& m)
  {
    using namespace boost::python;

    typedef boost::python::return_value_policy<
      boost::python::copy_const_reference>
        copy_const_reference;

    register_flex_grid_default_index_type_conversions();

    m
      .setattr("__version__",
        scitbx::boost_python::cvs_revision("$Revision$"))
    ;

    m.add(
      class_<flex_grid<> >("grid",
                  args<>())
        .def_init(args<
          flex_grid_default_index_type const&>())
        .def_init(args<
          flex_grid_default_index_type const&,
          flex_grid_default_index_type const&>())
        .def_init(args<
          flex_grid_default_index_type const&,
          flex_grid_default_index_type const&,
          bool>())
        .def("set_layout", flex_grid_wrappers::set_layout)
        .def("nd", &flex_grid<>::nd)
        .def("size_1d", &flex_grid<>::size_1d)
        .def("origin", &flex_grid<>::origin, copy_const_reference())
        .def("grid", &flex_grid<>::grid, copy_const_reference())
        .def("last", flex_grid_wrappers::last_0)
        .def("last", flex_grid_wrappers::last_1)
        //.def("last", &flex_grid<>::last, flex_grid_last_stubs())
        .def("layout", &flex_grid<>::layout, copy_const_reference())
        .def("is_0_based", &flex_grid<>::is_0_based)
        .def("is_padded", &flex_grid<>::is_padded)
        .def("__call__", &flex_grid<>::operator())
        .def("is_valid_index", &flex_grid<>::is_valid_index)
        .def("__eq__", &flex_grid<>::operator==)
        .def("__ne__", &flex_grid<>::operator!=)
        .def_pickle(flex_grid_wrappers())
    );

    add_flex_bool(m);
    add_flex_size_t(m);
    add_flex_int(m);
    add_flex_long(m);
    add_flex_float(m);
    add_flex_double(m);
    add_flex_complex_double(m);
    add_flex_std_string(m);
  }

}}} // namespace scitbx::af::boost_python

BOOST_PYTHON_MODULE_INIT(flex)
{
  boost::python::module this_module("flex");
  scitbx::af::boost_python::init_module(this_module);
}
