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
#include <scitbx/boost_python/utils.h>
#include <scitbx/boost_python/container_conversions.h>
#include <scitbx/array_family/small.h>
#include <scitbx/array_family/boost_python/shared_flex_conversions.h>
#include <boost/python/module_init.hpp>
#include <boost/python/scope.hpp>
#include <boost/python/def.hpp>

namespace scitbx { namespace {

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

  af::shared<double> make_shared()
  {
    af::shared<double> a;
    a.push_back(3);
    a.push_back(1);
    a.push_back(2);
    return a;
  }

  int use_shared(af::shared<double> const& a)
  {
    double s = 0;
    for(std::size_t i=0;i<a.size();i++) s += a[i];
    return (int) s + .5;
  }

  void modify_shared(af::shared<double>& a)
  {
    for(std::size_t i=0;i<a.size();i++) a[i] *= 2;
  }

  int use_const_ref(af::const_ref<double> const& a)
  {
    double s = 0;
    for(std::size_t i=0;i<a.size();i++) s += a[i];
    return (int) s + .5;
  }

  void modify_ref(af::ref<double> a)
  {
    for(std::size_t i=0;i<a.size();i++) a[i] *= 2;
  }

  boost::array<int, 2>
  make_boost_int_2(int x0=7, int x1=2)
  {
    boost::array<int, 2> result = {x0, x1};
    return result;
  }

  BOOST_PYTHON_FUNCTION_OVERLOADS(
    make_boost_int_2_stubs, make_boost_int_2, 0, 2)

  void init_module()
  {
    using namespace boost::python;

    scope().attr("__version__") = scitbx::boost_python::cvs_revision(
      "$Revision$");

    scitbx::boost_python::import_module("scitbx_boost.array_family.flex");

    def("std_vector", std_vector);
    def("std_list", std_list);
    def("boost_array", boost_array_3);
    def("boost_array", boost_array_4);
    def("small", small_6);
    def("make_shared", make_shared);
    def("use_shared", use_shared);
    def("modify_shared", modify_shared);
    def("use_const_ref", use_const_ref);
    def("modify_ref", modify_ref);
    def("make_boost_int_2", (boost::array<int, 2>(*)(int, int)) 0,
      make_boost_int_2_stubs());

    boost::python::to_python_converter<
      boost::array<int, 2>,
      scitbx::boost_python::container_conversions::to_tuple<
        boost::array<int, 2> > >();

    scitbx::boost_python::container_conversions::from_python_sequence<
      std::vector<double>,
      scitbx::boost_python::container_conversions::variable_capacity_policy>();

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

}} // namespace::scitbx::<anonymous>

BOOST_PYTHON_MODULE_INIT(regression_test_ext)
{
  scitbx::init_module();
}
