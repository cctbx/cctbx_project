// $Id$
/* Copyright (c) 2001 The Regents of the University of California through
   E.O. Lawrence Berkeley National Laboratory, subject to approval by the
   U.S. Department of Energy. See files COPYRIGHT.txt and
   cctbx/LICENSE.txt for further details.

   Revision history:
     Jul 2002: created, copy of miller_bpl.h (R.W. Grosse-Kunstleve)
 */

#ifndef CCTBX_HENDRICKSON_LATTMAN_BPL_H
#define CCTBX_HENDRICKSON_LATTMAN_BPL_H

#include <cctbx/array_family/tiny_bpl.h>
#include <cctbx/hendrickson_lattman.h>

BOOST_PYTHON_BEGIN_CONVERSION_NAMESPACE

  inline
  cctbx::hendrickson_lattman<double> from_python(PyObject* p,
    boost::python::type<cctbx::hendrickson_lattman<double> const&>)
  {
    return cctbx::hendrickson_lattman<double>(from_python(p,
      boost::python::type<cctbx::af::tiny<double, 4> const&>()));
  }

  inline
  cctbx::hendrickson_lattman<double> from_python(PyObject* p,
    boost::python::type<cctbx::hendrickson_lattman<double> >)
  {
    return cctbx::hendrickson_lattman<double>(from_python(p,
      boost::python::type<cctbx::af::tiny<double, 4> const&>()));
  }

  inline
  PyObject* to_python(cctbx::hendrickson_lattman<double> const& coeff)
  {
    return to_python(coeff.array());
  }

BOOST_PYTHON_END_CONVERSION_NAMESPACE

#endif // CCTBX_HENDRICKSON_LATTMAN_BPL_H
