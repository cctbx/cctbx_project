// This is an automatically generated file. Do not edit.
/* Copyright (c) 2001 The Regents of the University of California through
   E.O. Lawrence Berkeley National Laboratory, subject to approval by the
   U.S. Department of Energy. See files COPYRIGHT.txt and
   cctbx/LICENSE.txt for further details.

   Revision history:
     Apr 2001: SourceForge release (R.W. Grosse-Kunstleve)
 */

#ifndef CCTBX_BASIC_BOOST_ARRAY_BPL_H
#define CCTBX_BASIC_BOOST_ARRAY_BPL_H

#include <boost/array.hpp>

BOOST_PYTHON_BEGIN_CONVERSION_NAMESPACE


  boost::array<int, 3> from_python(PyObject* p,
    boost::python::type<const boost::array<int, 3>&>);

  PyObject* to_python(const boost::array<int, 3>& tobj);

  boost::array<int, 9> from_python(PyObject* p,
    boost::python::type<const boost::array<int, 9>&>);

  PyObject* to_python(const boost::array<int, 9>& tobj);

  boost::array<double, 3> from_python(PyObject* p,
    boost::python::type<const boost::array<double, 3>&>);

  PyObject* to_python(const boost::array<double, 3>& tobj);

  boost::array<double, 6> from_python(PyObject* p,
    boost::python::type<const boost::array<double, 6>&>);

  PyObject* to_python(const boost::array<double, 6>& tobj);

  boost::array<double, 9> from_python(PyObject* p,
    boost::python::type<const boost::array<double, 9>&>);

  PyObject* to_python(const boost::array<double, 9>& tobj);

BOOST_PYTHON_END_CONVERSION_NAMESPACE

#endif // CCTBX_BASIC_BOOST_ARRAY_BPL_H

