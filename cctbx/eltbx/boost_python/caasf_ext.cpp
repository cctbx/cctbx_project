/* Copyright (c) 2001-2002 The Regents of the University of California
   through E.O. Lawrence Berkeley National Laboratory, subject to
   approval by the U.S. Department of Energy.
   See files COPYRIGHT.txt and LICENSE.txt for further details.

   Revision history:
     2002 Sep: Created (rwgk)
 */

#include <cctbx/eltbx/caasf.h>
#include <scitbx/array_family/tiny.h>
#include <cctbx/import_scitbx_af.h>
#include <scitbx/boost_python/utils.h>
#include <scitbx/boost_python/iterator_wrappers.h>
#include <boost/python/module.hpp>
#include <boost/python/scope.hpp>
#include <boost/python/class.hpp>

namespace cctbx { namespace eltbx { namespace caasf { namespace boost_python {

namespace {

  template <std::size_t N>
  struct base_wrappers
  {
    typedef base<N> w_t;

    static af::tiny<float, N>
    a(w_t const& o)
    {
      af::tiny<float, N> result;
      for(std::size_t i=0;i<N;i++) result[i] = o.a(i);
      return result;
    }

    static af::tiny<float, N>
    b(w_t const& o)
    {
      af::tiny<float, N> result;
      for(std::size_t i=0;i<N;i++) result[i] = o.b(i);
      return result;
    }

    static void
    wrap(const char* python_name)
    {
      using namespace boost::python;
      class_<w_t>(python_name, no_init)
        .def("table", &w_t::table)
        .def("label", &w_t::label)
        .def("a", a)
        .def("b", b)
        .def("c", &w_t::c)
        .def("at_stol_sq", &w_t::at_stol_sq)
        .def("at_stol", &w_t::at_stol)
        .def("at_d_star_sq", &w_t::at_d_star_sq)
      ;
    }
  };

  struct it1992_wrappers
  {
    typedef it1992 w_t;

    static void
    wrap()
    {
      using namespace boost::python;
      base_wrappers<4>::wrap("base_4");
      class_<w_t, bases<base<4> > >("it1992", no_init)
        .def(init<std::string const&, optional<bool> >())
      ;
    }
  };

  struct wk1995_wrappers
  {
    typedef wk1995 w_t;

    static void
    wrap()
    {
      using namespace boost::python;
      base_wrappers<5>::wrap("base_5");
      class_<w_t, bases<base<5> > >("wk1995", no_init)
        .def(init<std::string const&, optional<bool> >())
      ;
    }
  };

  void init_module()
  {
    using namespace boost::python;

    scope().attr("__version__") = scitbx::boost_python::cvs_revision(
      "$Revision$");

    it1992_wrappers::wrap();
    scitbx::boost_python::iterator_wrappers<
      it1992, it1992_iterator>::wrap("it1992_iterator");

    wk1995_wrappers::wrap();
    scitbx::boost_python::iterator_wrappers<
      wk1995, wk1995_iterator>::wrap("wk1995_iterator");
  }

} // namespace <anonymous>
}}}} // namespace cctbx::eltbx::caasf::boost_python

BOOST_PYTHON_MODULE(caasf_ext)
{
  cctbx::eltbx::caasf::boost_python::init_module();
}
