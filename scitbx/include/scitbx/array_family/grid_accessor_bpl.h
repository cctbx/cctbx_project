/* Copyright (c) 2001-2002 The Regents of the University of California
   through E.O. Lawrence Berkeley National Laboratory, subject to
   approval by the U.S. Department of Energy.
   See files COPYRIGHT.txt and LICENSE.txt for further details.

   Revision history:
     2002 Aug: Copied from cctbx/array_family (R.W. Grosse-Kunstleve)
     2002 Mar: Modified copy of miller_bpl.h (rwgk)
     2001 Apr: SourceForge release (R.W. Grosse-Kunstleve)
 */

#ifndef SCITBX_GRID_ACCESSOR_BPL_H
#define SCITBX_GRID_ACCESSOR_BPL_H

#include <scitbx/array_family/tiny_bpl.h>
#include <scitbx/array_family/grid_accessor.h>

BOOST_PYTHON_BEGIN_CONVERSION_NAMESPACE

  inline
  scitbx::af::grid<3> from_python(PyObject* p,
    boost::python::type<scitbx::af::grid<3> const&>)
  {
    return scitbx::af::grid<3>(from_python(p,
      boost::python::type<scitbx::af::grid<3>::index_type const&>()));
  }

  inline
  scitbx::af::grid<3> from_python(PyObject* p,
    boost::python::type<scitbx::af::grid<3> >)
  {
    return scitbx::af::grid<3>(from_python(p,
      boost::python::type<scitbx::af::grid<3>::index_type const&>()));
  }

  inline
  PyObject* to_python(scitbx::af::grid<3> const& g) {
    return to_python(static_cast<scitbx::af::grid<3>::index_type const&>(g));
  }

BOOST_PYTHON_END_CONVERSION_NAMESPACE

#endif // SCITBX_GRID_ACCESSOR_BPL_H
