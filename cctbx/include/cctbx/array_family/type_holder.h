// $Id$
/* Copyright (c) 2001 The Regents of the University of California through
   E.O. Lawrence Berkeley National Laboratory, subject to approval by the
   U.S. Department of Energy. See files COPYRIGHT.txt and
   cctbx/LICENSE.txt for further details.

   Revision history:
     2002 Aug: Renamed: cctbx/basic/meta.h -> cctbx/array_family/type_holder.h
     Jan 2002: Created (R.W. Grosse-Kunstleve)
 */

#ifndef CCTBX_ARRAY_FAMILY_TYPE_HOLDER_H
#define CCTBX_ARRAY_FAMILY_TYPE_HOLDER_H

namespace cctbx { namespace af {

    template <class T> struct type_holder { typedef T type; };

}} // namespace cctbx::af

#endif // CCTBX_ARRAY_FAMILY_TYPE_HOLDER_H
