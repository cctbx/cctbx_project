// $Id$
/* Copyright (c) 2001 The Regents of the University of California through
   E.O. Lawrence Berkeley National Laboratory, subject to approval by the
   U.S. Department of Energy. See files COPYRIGHT.txt and
   cctbx/LICENSE.txt for further details.

   Revision history:
     2002 Aug: Created (R.W. Grosse-Kunstleve)
 */

#ifndef CCTBX_ARRAY_FAMILY_ARRAY_ADAPTOR_H
#define CCTBX_ARRAY_FAMILY_ARRAY_ADAPTOR_H

namespace cctbx { namespace af {

  template <typename T>
  struct array_adaptor
  {
    const T* pointee;
    array_adaptor(T const& a) : pointee(&a) {}
  };

  template <typename T>
  inline
  array_adaptor<T>
  adapt(T const& a)
  {
    return array_adaptor<T>(a);
  }

}} // namespace cctbx::af

#endif // CCTBX_ARRAY_FAMILY_ARRAY_ADAPTOR_H
