// $Id$
/* Copyright (c) 2001 The Regents of the University of California through
   E.O. Lawrence Berkeley National Laboratory, subject to approval by the
   U.S. Department of Energy. See files COPYRIGHT.txt and
   cctbx/LICENSE.txt for further details.

   Revision history:
     2002 Aug: Created (R.W. Grosse-Kunstleve)
 */

#ifndef CCTBX_ARRAY_FAMILY_SHARED_BPL__H
#define CCTBX_ARRAY_FAMILY_SHARED_BPL__H

#include <cctbx/error.h>
#include <cctbx/array_family/flex_grid_accessor.h>
#include <cctbx/array_family/versa.h>
#include <cctbx/array_family/shared.h>

#define CCTBX_ARRAY_FAMILY_FLEX_IMPORT(ElementType, flex_name) \
    { \
      boost::python::import_converters< \
        cctbx::af::versa<ElementType, cctbx::af::flex_grid<> > > \
      py_flex("cctbx_boost.arraytbx.flex", flex_name); \
    }

#define CCTBX_ARRAY_FAMILY_IMPLICIT_SHARED_CONVERTERS(ElementType) \
\
BOOST_PYTHON_BEGIN_CONVERSION_NAMESPACE \
\
  cctbx::af::shared<ElementType > \
  from_python( \
    PyObject* p, \
    boost::python::type<cctbx::af::shared<ElementType > const&>) \
  { \
    typedef cctbx::af::versa<ElementType, cctbx::af::flex_grid<> > flex_type; \
    flex_type flex = from_python(p, boost::python::type<flex_type const&>()); \
    cctbx_assert(flex.accessor().nd() == 1); \
    cctbx_assert(flex.accessor().origin()[0] == 0); \
    cctbx::af::shared<ElementType > b = flex.as_base_array(); \
    cctbx_assert(flex.accessor().grid()[0] == b.size()); \
    return b; \
  } \
\
  cctbx::af::shared<ElementType > \
  from_python( \
    PyObject* p, \
    boost::python::type<cctbx::af::shared<ElementType > >) \
  { \
    return from_python(p, \
      boost::python::type<cctbx::af::shared<ElementType > const&>()); \
  } \
\
  PyObject* to_python(cctbx::af::shared<ElementType > const& a) \
  { \
    typedef cctbx::af::versa<ElementType, cctbx::af::flex_grid<> > flex_type; \
    cctbx::af::flex_grid_default_index_type grid; \
    grid.append(a.size()); \
    return to_python(flex_type(a, cctbx::af::flex_grid<>(grid))); \
  } \
\
BOOST_PYTHON_END_CONVERSION_NAMESPACE

#endif // CCTBX_ARRAY_FAMILY_SHARED_BPL__H
