// $Id$
/* Copyright (c) 2001 The Regents of the University of California through
   E.O. Lawrence Berkeley National Laboratory, subject to approval by the
   U.S. Department of Energy. See files COPYRIGHT.txt and
   cctbx/LICENSE.txt for further details.

   Revision history:
     Apr 2001: SourceForge release (R.W. Grosse-Kunstleve)
 */

#ifndef CCTBX_MILLER_BPL_H
#define CCTBX_MILLER_BPL_H

#include <cctbx/miller.h>

BOOST_PYTHON_BEGIN_CONVERSION_NAMESPACE

  cctbx::Miller::Index from_python(PyObject* p,
    boost::python::type<const cctbx::Miller::Index&>)
  {
    return cctbx::Miller::Index(from_python(p,
      boost::python::type<const cctbx::int3&>()));
  }

  cctbx::Miller::Index from_python(PyObject* p,
    boost::python::type<cctbx::Miller::Index>)
  {
    return cctbx::Miller::Index(from_python(p,
      boost::python::type<const cctbx::int3&>()));
  }

  PyObject* to_python(const cctbx::Miller::Index& MIx) {
    return to_python(static_cast<const cctbx::int3&>(MIx));
  }

BOOST_PYTHON_END_CONVERSION_NAMESPACE

#endif // CCTBX_MILLER_BPL_H
