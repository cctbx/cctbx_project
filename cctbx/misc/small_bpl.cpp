// $Id$
/* Copyright (c) 2001 The Regents of the University of California through
   E.O. Lawrence Berkeley National Laboratory, subject to approval by the
   U.S. Department of Energy. See files COPYRIGHT.txt and
   cctbx/LICENSE.txt for further details.

   Revision history:
     2002 Aug: Created (R.W. Grosse-Kunstleve)
 */

#include <cctbx/array_family/small_bpl.h>

BOOST_PYTHON_BEGIN_CONVERSION_NAMESPACE

  cctbx::af::flex_grid_default_index_type
  from_python(PyObject* p,
    boost::python::type<cctbx::af::flex_grid_default_index_type const&>)

  {
    typedef
      cctbx::af::flex_grid_default_index_type::value_type
        index_value_type;
    boost::python::tuple
    tup = cctbx::bpl_utils::tuple_from_python_list_or_tuple(p);
    cctbx::af::flex_grid_default_index_type result;
    if (tup.size() >= result.capacity()) {
      PyErr_SetString(PyExc_ValueError, "too many values in tuple.");
      boost::python::throw_error_already_set();
    }
    for(int i=0;i<tup.size();i++) {
      result.push_back(
        from_python(tup[i].get(), boost::python::type<index_value_type>()));
    }
    return result;
  }

  PyObject* to_python(cctbx::af::flex_grid_default_index_type const& tobj)
  {
    boost::python::tuple result(tobj.size());
    for(int i=0;i<tobj.size();i++) {
      result.set_item(i, boost::python::ref(to_python(tobj[i])));
    }
    return result.reference().release();
  }

BOOST_PYTHON_END_CONVERSION_NAMESPACE
