// This is an automatically generated file. Do not edit.
/* Copyright (c) 2001 The Regents of the University of California through
   E.O. Lawrence Berkeley National Laboratory, subject to approval by the
   U.S. Department of Energy. See files COPYRIGHT.txt and
   cctbx/LICENSE.txt for further details.

   Revision history:
     Apr 2001: SourceForge release (R.W. Grosse-Kunstleve)
 */

#include <boost/array.hpp>
#include <cctbx/bpl_utils.h>

BOOST_PYTHON_BEGIN_CONVERSION_NAMESPACE


  boost::array<int, 3> from_python(PyObject* p,
    boost::python::type<const boost::array<int, 3>&>)

  {
    boost::python::tuple
    tup = cctbx::bpl_utils::tuple_from_python_list_or_tuple(p);
    boost::array<int, 3> result;
    if (tup.size() != result.size()) {
      PyErr_SetString(PyExc_ValueError,
        "incorrect number of values in tuple.");
      throw boost::python::error_already_set();
    }
    for(int i=0;i<result.size();i++) {
      result[i] = from_python(tup[i].get(), boost::python::type<int>());
    }
    return result;
  }

  PyObject* to_python(const boost::array<int, 3>& tobj)

  {
    boost::python::tuple result(tobj.size());
    for(int i=0;i<tobj.size();i++) {
      result.set_item(i, boost::python::ref(to_python(tobj[i])));
    }
    return result.reference().release();
  }


  boost::array<int, 9> from_python(PyObject* p,
    boost::python::type<const boost::array<int, 9>&>)

  {
    boost::python::tuple
    tup = cctbx::bpl_utils::tuple_from_python_list_or_tuple(p);
    boost::array<int, 9> result;
    if (tup.size() != result.size()) {
      PyErr_SetString(PyExc_ValueError,
        "incorrect number of values in tuple.");
      throw boost::python::error_already_set();
    }
    for(int i=0;i<result.size();i++) {
      result[i] = from_python(tup[i].get(), boost::python::type<int>());
    }
    return result;
  }

  PyObject* to_python(const boost::array<int, 9>& tobj)

  {
    boost::python::tuple result(tobj.size());
    for(int i=0;i<tobj.size();i++) {
      result.set_item(i, boost::python::ref(to_python(tobj[i])));
    }
    return result.reference().release();
  }


  boost::array<std::size_t, 3> from_python(PyObject* p,
    boost::python::type<const boost::array<std::size_t, 3>&>)

  {
    boost::python::tuple
    tup = cctbx::bpl_utils::tuple_from_python_list_or_tuple(p);
    boost::array<std::size_t, 3> result;
    if (tup.size() != result.size()) {
      PyErr_SetString(PyExc_ValueError,
        "incorrect number of values in tuple.");
      throw boost::python::error_already_set();
    }
    for(int i=0;i<result.size();i++) {
      result[i] = from_python(tup[i].get(), boost::python::type<std::size_t>());
    }
    return result;
  }

  PyObject* to_python(const boost::array<std::size_t, 3>& tobj)

  {
    boost::python::tuple result(tobj.size());
    for(int i=0;i<tobj.size();i++) {
      result.set_item(i, boost::python::ref(to_python(tobj[i])));
    }
    return result.reference().release();
  }


  boost::array<double, 3> from_python(PyObject* p,
    boost::python::type<const boost::array<double, 3>&>)

  {
    boost::python::tuple
    tup = cctbx::bpl_utils::tuple_from_python_list_or_tuple(p);
    boost::array<double, 3> result;
    if (tup.size() != result.size()) {
      PyErr_SetString(PyExc_ValueError,
        "incorrect number of values in tuple.");
      throw boost::python::error_already_set();
    }
    for(int i=0;i<result.size();i++) {
      result[i] = from_python(tup[i].get(), boost::python::type<double>());
    }
    return result;
  }

  PyObject* to_python(const boost::array<double, 3>& tobj)

  {
    boost::python::tuple result(tobj.size());
    for(int i=0;i<tobj.size();i++) {
      result.set_item(i, boost::python::ref(to_python(tobj[i])));
    }
    return result.reference().release();
  }


  boost::array<double, 6> from_python(PyObject* p,
    boost::python::type<const boost::array<double, 6>&>)

  {
    boost::python::tuple
    tup = cctbx::bpl_utils::tuple_from_python_list_or_tuple(p);
    boost::array<double, 6> result;
    if (tup.size() != result.size()) {
      PyErr_SetString(PyExc_ValueError,
        "incorrect number of values in tuple.");
      throw boost::python::error_already_set();
    }
    for(int i=0;i<result.size();i++) {
      result[i] = from_python(tup[i].get(), boost::python::type<double>());
    }
    return result;
  }

  PyObject* to_python(const boost::array<double, 6>& tobj)

  {
    boost::python::tuple result(tobj.size());
    for(int i=0;i<tobj.size();i++) {
      result.set_item(i, boost::python::ref(to_python(tobj[i])));
    }
    return result.reference().release();
  }


  boost::array<double, 9> from_python(PyObject* p,
    boost::python::type<const boost::array<double, 9>&>)

  {
    boost::python::tuple
    tup = cctbx::bpl_utils::tuple_from_python_list_or_tuple(p);
    boost::array<double, 9> result;
    if (tup.size() != result.size()) {
      PyErr_SetString(PyExc_ValueError,
        "incorrect number of values in tuple.");
      throw boost::python::error_already_set();
    }
    for(int i=0;i<result.size();i++) {
      result[i] = from_python(tup[i].get(), boost::python::type<double>());
    }
    return result;
  }

  PyObject* to_python(const boost::array<double, 9>& tobj)

  {
    boost::python::tuple result(tobj.size());
    for(int i=0;i<tobj.size();i++) {
      result.set_item(i, boost::python::ref(to_python(tobj[i])));
    }
    return result.reference().release();
  }


BOOST_PYTHON_END_CONVERSION_NAMESPACE

