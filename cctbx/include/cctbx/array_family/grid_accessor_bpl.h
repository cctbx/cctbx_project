// $Id$
/* Copyright (c) 2001 The Regents of the University of California through
   E.O. Lawrence Berkeley National Laboratory, subject to approval by the
   U.S. Department of Energy. See files COPYRIGHT.txt and
   cctbx/LICENSE.txt for further details.

   Revision history:
     Mar 2002: Modified copy of miller_bpl.h (rwgk)
     Apr 2001: SourceForge release (R.W. Grosse-Kunstleve)
 */

#ifndef CCTBX_GRID_ACCESSOR_BPL_H
#define CCTBX_GRID_ACCESSOR_BPL_H

#include <cctbx/array_family/tiny_bpl.h>
#include <cctbx/array_family/grid_accessor.h>

BOOST_PYTHON_BEGIN_CONVERSION_NAMESPACE

  inline
  cctbx::af::grid<3> from_python(PyObject* p,
    boost::python::type<const cctbx::af::grid<3>&>)
  {
    return cctbx::af::grid<3>(from_python(p,
      boost::python::type<const cctbx::af::grid<3>::index_type&>()));
  }

  inline
  cctbx::af::grid<3> from_python(PyObject* p,
    boost::python::type<cctbx::af::grid<3> >)
  {
    return cctbx::af::grid<3>(from_python(p,
      boost::python::type<const cctbx::af::grid<3>::index_type&>()));
  }

  inline
  PyObject* to_python(const cctbx::af::grid<3>& g) {
    return to_python(static_cast<const cctbx::af::grid<3>::index_type&>(g));
  }

BOOST_PYTHON_END_CONVERSION_NAMESPACE

#endif // CCTBX_GRID_ACCESSOR_BPL_H
