/* Copyright (c) 2001-2002 The Regents of the University of California
   through E.O. Lawrence Berkeley National Laboratory, subject to
   approval by the U.S. Department of Energy.
   See files COPYRIGHT.txt and LICENSE.txt for further details.

   Revision history:
     2002 Aug: Copied from cctbx/array_family (R.W. Grosse-Kunstleve)
     2002 Aug: Created (R.W. Grosse-Kunstleve)
 */

#ifndef SCITBX_ARRAY_FAMILY_FLEX_BPL_H
#define SCITBX_ARRAY_FAMILY_FLEX_BPL_H

#include <scitbx/error.h>
#include <scitbx/array_family/flex_grid_accessor.h>
#include <scitbx/array_family/versa.h>

#define SCITBX_ARRAY_FAMILY_FLEX_IMPORT(ElementType, flex_name) \
    { \
      boost::python::import_converters< \
        scitbx::af::versa<ElementType, scitbx::af::flex_grid<> > > \
      py_flex("scitbx_boost.arraytbx.flex", flex_name); \
    }

#define SCITBX_ARRAY_FAMILY_IMPLICIT_VERSA_GRID_CONVERTERS(ElementType, Nd) \
\
BOOST_PYTHON_BEGIN_CONVERSION_NAMESPACE \
\
  scitbx::af::versa<ElementType, scitbx::af::grid<Nd> > \
  from_python( \
    PyObject* p, \
    boost::python::type< \
      scitbx::af::versa<ElementType, scitbx::af::grid<Nd> > const&>) \
  { \
    typedef scitbx::af::versa<ElementType, scitbx::af::flex_grid<> > flex_type; \
    flex_type flex = from_python(p, boost::python::type<flex_type const&>()); \
    scitbx_assert(flex.accessor().nd() == Nd); \
    scitbx_assert(flex.accessor().is_0_based()); \
    scitbx_assert(!flex.accessor().is_padded()); \
    return scitbx::af::versa<ElementType, scitbx::af::grid<Nd> >( \
      flex, scitbx::af::adapt(flex.accessor().grid())); \
  } \
\
  scitbx::af::versa<ElementType, scitbx::af::grid<Nd> > \
  from_python( \
    PyObject* p, \
    boost::python::type< \
      scitbx::af::versa<ElementType, scitbx::af::grid<Nd> > >) \
  { \
    return from_python(p, \
      boost::python::type< \
        scitbx::af::versa<ElementType, scitbx::af::grid<Nd> > const&>()); \
  } \
 \
  PyObject* to_python( \
    scitbx::af::versa<ElementType, scitbx::af::grid<Nd> > const& a) \
  { \
    typedef scitbx::af::versa<ElementType, scitbx::af::flex_grid<> > flex_type; \
    return to_python(flex_type(a, scitbx::af::adapt(a.accessor()))); \
  } \
\
BOOST_PYTHON_END_CONVERSION_NAMESPACE

#define SCITBX_ARRAY_FAMILY_IMPLICIT_SHARED_CONVERTERS(ElementType) \
\
BOOST_PYTHON_BEGIN_CONVERSION_NAMESPACE \
\
  scitbx::af::shared<ElementType > \
  from_python( \
    PyObject* p, \
    boost::python::type<scitbx::af::shared<ElementType > const&>) \
  { \
    typedef scitbx::af::versa<ElementType, scitbx::af::flex_grid<> > flex_type; \
    flex_type flex = from_python(p, boost::python::type<flex_type const&>()); \
    scitbx_assert(flex.accessor().nd() == 1); \
    scitbx_assert(flex.accessor().is_0_based()); \
    scitbx_assert(!flex.accessor().is_padded()); \
    scitbx::af::shared<ElementType > b = flex.as_base_array(); \
    scitbx_assert(flex.accessor().grid()[0] == b.size()); \
    return b; \
  } \
\
  scitbx::af::shared<ElementType > \
  from_python( \
    PyObject* p, \
    boost::python::type<scitbx::af::shared<ElementType > >) \
  { \
    return from_python(p, \
      boost::python::type<scitbx::af::shared<ElementType > const&>()); \
  } \
\
  PyObject* to_python(scitbx::af::shared<ElementType > const& a) \
  { \
    typedef scitbx::af::versa<ElementType, scitbx::af::flex_grid<> > flex_type; \
    scitbx::af::flex_grid_default_index_type grid; \
    grid.append(a.size()); \
    return to_python(flex_type(a, scitbx::af::flex_grid<>(grid))); \
  } \
\
BOOST_PYTHON_END_CONVERSION_NAMESPACE

#endif // SCITBX_ARRAY_FAMILY_FLEX_BPL_H
