// $Id$
/* Copyright (c) 2001 The Regents of the University of California through
   E.O. Lawrence Berkeley National Laboratory, subject to approval by the
   U.S. Department of Energy. See files COPYRIGHT.txt and
   cctbx/LICENSE.txt for further details.

   Revision history:
     Created 2001 Jul 03 (R.W. Grosse-Kunstleve)
 */

#ifndef CCTBX_COORDINATES_BPL_H
#define CCTBX_COORDINATES_BPL_H

#include <cctbx/coordinates.h>

BOOST_PYTHON_BEGIN_CONVERSION_NAMESPACE

  cctbx::coordinates::cartesian<double> from_python(PyObject* p,
    boost::python::type<const cctbx::coordinates::cartesian<double>&>)
  {
    return cctbx::coordinates::cartesian<double>(from_python(p,
      boost::python::type<const boost::array<double, 3>&>()));
  }

  PyObject* to_python(const cctbx::coordinates::cartesian<double>& Xc) {
    return to_python(static_cast<const boost::array<double, 3>&>(Xc));
  }

  cctbx::coordinates::fractional<double> from_python(PyObject* p,
    boost::python::type<const cctbx::coordinates::fractional<double>&>)
  {
    return cctbx::coordinates::fractional<double>(from_python(p,
      boost::python::type<const boost::array<double, 3>&>()));
  }

  PyObject* to_python(const cctbx::coordinates::fractional<double>& Xf) {
    return to_python(static_cast<const boost::array<double, 3>&>(Xf));
  }

BOOST_PYTHON_END_CONVERSION_NAMESPACE

#endif // CCTBX_COORDINATES_BPL_H
