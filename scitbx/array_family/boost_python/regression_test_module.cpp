/* Copyright (c) 2001-2002 The Regents of the University of California
   through E.O. Lawrence Berkeley National Laboratory, subject to
   approval by the U.S. Department of Energy.
   See files COPYRIGHT.txt and LICENSE.txt for further details.

   Revision history:
     2002 Aug: Created (R.W. Grosse-Kunstleve)
 */

#include <vector>
#include <list>
#include <boost/array.hpp>
#include <scitbx/error.h>
#include <scitbx/array_family/small.h>
#include <scitbx/boost_python/utils.h>
#include <scitbx/boost_python/container_conversions.h>
#include <boost/python/module.hpp>
#include <boost/python/class.hpp>

namespace {

  using namespace scitbx;

  int std_vector(std::vector<double> const& a)
  {
    double s = 0;
    for(std::size_t i=0;i<a.size();i++) s += a[i];
    return (int) s + .5;
  }

  int std_list(std::list<double> const& lst)
  {
    double s = 0;
    for(std::list<double>::const_iterator i=lst.begin();i!=lst.end();i++) {
      s += *i;
    }
    return (int) s + .5;
  }

  int boost_array_3(boost::array<double, 3> const& a)
  {
    double s = 0;
    for(std::size_t i=0;i<a.size();i++) s += a[i];
    return (int) s + .5;
  }

  int boost_array_4(boost::array<double, 4> const& a)
  {
    double s = 0;
    for(std::size_t i=0;i<a.size();i++) s += a[i];
    return (int) s + .5;
  }

  int small_6(af::small<double, 6> const& a)
  {
    double s = 0;
    for(std::size_t i=0;i<a.size();i++) s += a[i];
    return (int) s + .5;
  }

  void init_module(boost::python::module& this_module)
  {
    this_module
      .setattr("__version__",
        scitbx::boost_python::cvs_revision("$Revision$"))
    ;

    scitbx::boost_python::import_module("scitbx_boost.array_family.flex");

    this_module
      .def("std_vector", std_vector)
      .def("std_list", std_list)
      .def("boost_array", boost_array_3)
      .def("boost_array", boost_array_4)
      .def("small", small_6)
    ;

    scitbx::boost_python::container_conversions::from_python_sequence<
      std::vector<double>,
      scitbx::boost_python::container_conversions::variable_size_policy>();

    scitbx::boost_python::container_conversions::from_python_sequence<
      std::list<double>,
      scitbx::boost_python::container_conversions::linked_list_policy>();

    scitbx::boost_python::container_conversions::from_python_sequence<
      boost::array<double, 3>,
      scitbx::boost_python::container_conversions::fixed_size_policy>();

    scitbx::boost_python::container_conversions::from_python_sequence<
      boost::array<double, 4>,
      scitbx::boost_python::container_conversions::fixed_size_policy>();

    scitbx::boost_python::container_conversions::from_python_sequence<
      af::small<double, 6>,
      scitbx::boost_python::container_conversions::fixed_capacity_policy>();
  }

}

BOOST_PYTHON_MODULE_INIT(regression_test_ext)
{
  boost::python::module this_module("regression_test_ext");
  init_module(this_module);
}
