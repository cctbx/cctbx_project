// $Id$
/* Copyright (c) 2001 The Regents of the University of California through
   E.O. Lawrence Berkeley National Laboratory, subject to approval by the
   U.S. Department of Energy. See files COPYRIGHT.txt and
   cctbx/LICENSE.txt for further details.

   Revision history:
     Apr 2001: SourceForge release (R.W. Grosse-Kunstleve)
 */

#include <cctbx/bpl_utils.h>

namespace cctbx { namespace bpl_utils {

  void raise_index_error()
  {
    PyErr_SetString(PyExc_IndexError, "Index out of range.");
    boost::python::throw_error_already_set();
  }

  void raise_must_be_0_based_1d()
  {
    PyErr_SetString(PyExc_RuntimeError,
      "Array must be 0-based 1-dimensional.");
    boost::python::throw_error_already_set();
  }

  void raise_must_be_0_based_3d()
  {
    PyErr_SetString(PyExc_RuntimeError,
      "Array must be 0-based 3-dimensional.");
    boost::python::throw_error_already_set();
  }

  void assert_0_based_1d(af::flex_grid<> const& grid)
  {
    if (grid.nd() != 1 || !grid.is_0_based()) raise_must_be_0_based_1d();
  }

  void assert_0_based_3d(af::flex_grid<> const& grid)
  {
    if (grid.nd() != 3 || !grid.is_0_based()) raise_must_be_0_based_3d();
  }

  void raise_shared_size_mismatch()
  {
    PyErr_SetString(PyExc_RuntimeError, "Shared size mismatch.");
    boost::python::throw_error_already_set();
  }

  void raise_incompatible_arrays()
  {
    PyErr_SetString(PyExc_RuntimeError, "Incompatible arrays.");
    boost::python::throw_error_already_set();
  }

  boost::python::tuple tuple_from_python_list_or_tuple(PyObject* p) {
    using namespace boost::python;
    if (PyList_Check(p))
      return tuple(ref(PyList_AsTuple(p)));
    else
      return tuple(ref(p, ref::increment_count));
  }

}}
