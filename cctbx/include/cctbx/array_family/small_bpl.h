// $Id$
/* Copyright (c) 2001 The Regents of the University of California through
   E.O. Lawrence Berkeley National Laboratory, subject to approval by the
   U.S. Department of Energy. See files COPYRIGHT.txt and
   cctbx/LICENSE.txt for further details.

   Revision history:
     2002 Aug: Created (R.W. Grosse-Kunstleve)
 */

#ifndef CCTBX_ARRAY_FAMILY_SMALL_BPL_H
#define CCTBX_ARRAY_FAMILY_SMALL_BPL_H

#include <cctbx/bpl_utils.h>
#include <cctbx/array_family/flex_grid_accessor.h>

BOOST_PYTHON_BEGIN_CONVERSION_NAMESPACE

  cctbx::af::flex_grid_default_index_type
  from_python(PyObject* p,
    boost::python::type<cctbx::af::flex_grid_default_index_type const&>);

  PyObject* to_python(cctbx::af::flex_grid_default_index_type const& tobj);

BOOST_PYTHON_END_CONVERSION_NAMESPACE

#endif // CCTBX_ARRAY_FAMILY_SMALL_BPL_H
