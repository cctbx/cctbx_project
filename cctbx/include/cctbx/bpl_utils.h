// $Id$
/* Copyright (c) 2001 The Regents of the University of California through
   E.O. Lawrence Berkeley National Laboratory, subject to approval by the
   U.S. Department of Energy. See files COPYRIGHT.txt and
   cctbx/LICENSE.txt for further details.

   Revision history:
     Apr 2001: SourceForge release (R.W. Grosse-Kunstleve)
 */

#ifndef CCTBX_BPL_UTILS_H
#define CCTBX_BPL_UTILS_H

#include <boost/python/class_builder.hpp>

namespace cctbx { namespace bpl_utils {

  void throw_index_out_of_range();
  boost::python::tuple tuple_from_python_list_or_tuple(PyObject* p);
}}

#endif // CCTBX_BPL_UTILS_H
