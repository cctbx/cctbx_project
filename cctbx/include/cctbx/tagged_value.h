/* Copyright (c) 2001-2002 The Regents of the University of California
   through E.O. Lawrence Berkeley National Laboratory, subject to
   approval by the U.S. Department of Energy.
   See files COPYRIGHT.txt and LICENSE.txt for further details.

   Revision history:
     2002 Oct: Fragment from cctbx/maps/sym_tags.h (rwgk)
     2002 Jan: Created (R.W. Grosse-Kunstleve)
 */

#ifndef CCTBX_TAGGED_VALUE_H
#define CCTBX_TAGGED_VALUE_H

namespace cctbx {

  template <typename ValueType, typename TagType = bool>
  struct tagged_value
  {
    typedef ValueType value_type;
    typedef TagType tag_type;

    tagged_value() {}

    tagged_value(ValueType const& v)
    : value(v)
    {}

    tagged_value(ValueType const& v, TagType const& t)
    : value(v), tag(t)
    {}

    ValueType value;
    TagType tag;
  };

} // namespace cctbx

#endif // CCTBX_TAGGED_VALUE_H
