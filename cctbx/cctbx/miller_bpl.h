// $Id$

#ifndef CCTBX_MILLER_BPL_H
#define CCTBX_MILLER_BPL_H

#include <cctbx/miller.h>

BOOST_PYTHON_BEGIN_CONVERSION_NAMESPACE

  Miller::Index from_python(PyObject* p,
    boost::python::type<const Miller::Index&>)
  {
    return Miller::Index(from_python(p,
      boost::python::type<const Miller::Vec3&>()).data());
  }

  PyObject* to_python(const Miller::Index& MIx) {
    return to_python(static_cast<const Miller::Vec3&>(MIx));
  }

BOOST_PYTHON_END_CONVERSION_NAMESPACE

#endif // CCTBX_MILLER_BPL_H
