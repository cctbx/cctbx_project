// $Id$
/* Copyright (c) 2001 The Regents of the University of California through
   E.O. Lawrence Berkeley National Laboratory, subject to approval by the
   U.S. Department of Energy. See files COPYRIGHT.txt and
   cctbx/LICENSE.txt for further details.

   Revision history:
     2002 Jan 14: Created (rwgk)
 */

#ifndef CCTBX_SHARED_NDIM_BPL_H
#define CCTBX_SHARED_NDIM_BPL_H

#include <cctbx/carray_bpl.h>
#include <cctbx/ndim.h>

BOOST_PYTHON_BEGIN_CONVERSION_NAMESPACE

  cctbx::dimension<3> from_python(PyObject* p,
    boost::python::type<const cctbx::dimension<3>&>)
  {
    return cctbx::dimension<3>(from_python(p,
      boost::python::type<const cctbx::carray<int, 3>&>()));
  }

  cctbx::dimension<3> from_python(PyObject* p,
    boost::python::type<cctbx::dimension<3>&>)
  {
    return cctbx::dimension<3>(from_python(p,
      boost::python::type<const cctbx::carray<int, 3>&>()));
  }

  PyObject* to_python(const cctbx::dimension<3>& dim) {
    return to_python(static_cast<const cctbx::carray<int, 3>&>(dim));
  }

BOOST_PYTHON_END_CONVERSION_NAMESPACE

#endif // CCTBX_SHARED_NDIM_BPL_H
