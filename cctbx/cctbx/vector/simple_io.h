// $Id$
/* Copyright (c) 2001 The Regents of the University of California through
   E.O. Lawrence Berkeley National Laboratory, subject to approval by the
   U.S. Department of Energy. See files COPYRIGHT.txt and
   cctbx/LICENSE.txt for further details.

   Revision history:
     Jan 2002: Created (R.W. Grosse-Kunstleve)
 */

#ifndef CCTBX_VECTOR_SIMPLE_IO_H
#define CCTBX_VECTOR_SIMPLE_IO_H

#include <iostream>

namespace cctbx { namespace vector {

  template <typename VectorType>
  std::ostream& simple_output(std::ostream& os, const VectorType& v) {
    os << "(";
    if (v.size() > 0) {
      for (std::size_t i = 0;;) {
        os << v[i];
        i++;
        if (i == v.size()) break;
        os << ",";
      }
    }
    os << ")";
    return os;
  }

}} // namespace cctbx::vector

#endif // CCTBX_VECTOR_SIMPLE_IO_H
