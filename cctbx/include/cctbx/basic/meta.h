// $Id$
/* Copyright (c) 2001 The Regents of the University of California through
   E.O. Lawrence Berkeley National Laboratory, subject to approval by the
   U.S. Department of Energy. See files COPYRIGHT.txt and
   cctbx/LICENSE.txt for further details.

   Revision history:
     Jan 2002: Created (R.W. Grosse-Kunstleve)
 */

#ifndef CCTBX_BASIC_META_H
#define CCTBX_BASIC_META_H

namespace cctbx {

    template <class T> struct type_holder { typedef T type; };

} // namespace cctbx

#endif // CCTBX_BASIC_META_H
