/* Copyright (c) 2001-2002 The Regents of the University of California
   through E.O. Lawrence Berkeley National Laboratory, subject to
   approval by the U.S. Department of Energy.
   See files COPYRIGHT.txt and LICENSE.txt for further details.

   Revision history:
     2002 Aug: Copied from cctbx/array_family (R.W. Grosse-Kunstleve)
     2002 Aug: Renamed: scitbx/basic/meta.h -> scitbx/array_family/type_holder.h
     2002 Jan: Created (R.W. Grosse-Kunstleve)
 */

#ifndef SCITBX_ARRAY_FAMILY_TYPE_HOLDER_H
#define SCITBX_ARRAY_FAMILY_TYPE_HOLDER_H

namespace scitbx { namespace af {

    template <class T> struct type_holder { typedef T type; };

}} // namespace scitbx::af

#endif // SCITBX_ARRAY_FAMILY_TYPE_HOLDER_H
