/* Copyright (c) 2001-2002 The Regents of the University of California
   through E.O. Lawrence Berkeley National Laboratory, subject to
   approval by the U.S. Department of Energy.
   See files COPYRIGHT.txt and LICENSE.txt for further details.

   Revision history:
     2002 Aug: Fragments from cctbx/misc/bpl_utils.cpp (R.W. Grosse-Kunstleve)
     2001 Apr: SourceForge release (R.W. Grosse-Kunstleve)
 */

#include <scitbx/array_family/boost_python/utils.h>
#include <boost/python/errors.hpp>

namespace scitbx { namespace af { namespace boost_python {

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

}}} // namespace scitbx::af::boost_python
