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
