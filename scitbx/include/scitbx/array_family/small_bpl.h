/* Copyright (c) 2001-2002 The Regents of the University of California
   through E.O. Lawrence Berkeley National Laboratory, subject to
   approval by the U.S. Department of Energy.
   See files COPYRIGHT.txt and LICENSE.txt for further details.

   Revision history:
     2002 Aug: Copied from cctbx/array_family (R.W. Grosse-Kunstleve)
     2002 Aug: Created (R.W. Grosse-Kunstleve)
 */

#ifndef SCITBX_ARRAY_FAMILY_SMALL_BPL_H
#define SCITBX_ARRAY_FAMILY_SMALL_BPL_H

#include <scitbx/bpl_utils.h>
#include <scitbx/array_family/flex_grid_accessor.h>

BOOST_PYTHON_BEGIN_CONVERSION_NAMESPACE

  scitbx::af::flex_grid_default_index_type
  from_python(PyObject* p,
    boost::python::type<scitbx::af::flex_grid_default_index_type const&>);

  PyObject* to_python(scitbx::af::flex_grid_default_index_type const& tobj);

BOOST_PYTHON_END_CONVERSION_NAMESPACE

#endif // SCITBX_ARRAY_FAMILY_SMALL_BPL_H
