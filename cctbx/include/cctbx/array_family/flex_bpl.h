// $Id$
/* Copyright (c) 2001 The Regents of the University of California through
   E.O. Lawrence Berkeley National Laboratory, subject to approval by the
   U.S. Department of Energy. See files COPYRIGHT.txt and
   cctbx/LICENSE.txt for further details.

   Revision history:
     2002 Aug: Created (R.W. Grosse-Kunstleve)
 */

#ifndef CCTBX_ARRAY_FAMILY_FLEX_BPL_H
#define CCTBX_ARRAY_FAMILY_FLEX_BPL_H

#include <cctbx/error.h>
#include <cctbx/array_family/flex_grid_accessor.h>
#include <cctbx/array_family/versa.h>

#define CCTBX_ARRAY_FAMILY_FLEX_IMPORT(ElementType, flex_name) \
    { \
      boost::python::import_converters< \
        cctbx::af::versa<ElementType, cctbx::af::flex_grid<> > > \
      py_flex("cctbx_boost.arraytbx.flex", flex_name); \
    }

#endif // CCTBX_ARRAY_FAMILY_FLEX_BPL_H
