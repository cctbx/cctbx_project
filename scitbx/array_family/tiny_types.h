#ifndef SCITBX_ARRAY_FAMILY_TINY_TYPES_H
#define SCITBX_ARRAY_FAMILY_TINY_TYPES_H

#include <scitbx/array_family/tiny.h>

namespace scitbx { namespace af {

  typedef tiny<int, 2> int2;
  typedef tiny<int, 3> int3;
  typedef tiny<int, 4> int4;
  typedef tiny<int, 6> int6;
  typedef tiny<int, 9> int9;
  typedef tiny<long, 3> long3;
  typedef tiny<double, 2> double2;
  typedef tiny<double, 3> double3;
  typedef tiny<double, 4> double4;
  typedef tiny<double, 6> double6;
  typedef tiny<double, 9> double9;
  typedef tiny<float, 2> float2;
  typedef tiny<float, 3> float3;
  typedef tiny<float, 4> float4;
  typedef tiny<float, 6> float6;
  typedef tiny<float, 9> float9;

}} // namespace scitbx::af

#if !defined(BOOST_NO_TEMPLATE_PARTIAL_SPECIALIZATION)

namespace boost {

  template<>
  struct has_trivial_destructor<scitbx::af::int2> {
    static const bool value = true;
  };

  template<>
  struct has_trivial_destructor<scitbx::af::int3> {
    static const bool value = true;
  };

  template<>
  struct has_trivial_destructor<scitbx::af::int4> {
    static const bool value = true;
  };

  template<>
  struct has_trivial_destructor<scitbx::af::int6> {
    static const bool value = true;
  };

  template<>
  struct has_trivial_destructor<scitbx::af::int9> {
    static const bool value = true;
  };

  template<>
  struct has_trivial_destructor<scitbx::af::long3> {
    static const bool value = true;
  };

  template<>
  struct has_trivial_destructor<scitbx::af::double2> {
    static const bool value = true;
  };

  template<>
  struct has_trivial_destructor<scitbx::af::double3> {
    static const bool value = true;
  };

  template<>
  struct has_trivial_destructor<scitbx::af::double4> {
    static const bool value = true;
  };

  template<>
  struct has_trivial_destructor<scitbx::af::double6> {
    static const bool value = true;
  };

  template<>
  struct has_trivial_destructor<scitbx::af::double9> {
    static const bool value = true;
  };

  template<>
  struct has_trivial_destructor<scitbx::af::float2> {
    static const bool value = true;
  };

  template<>
  struct has_trivial_destructor<scitbx::af::float3> {
    static const bool value = true;
  };

  template<>
  struct has_trivial_destructor<scitbx::af::float4> {
    static const bool value = true;
  };

  template<>
  struct has_trivial_destructor<scitbx::af::float6> {
    static const bool value = true;
  };

  template<>
  struct has_trivial_destructor<scitbx::af::float9> {
    static const bool value = true;
  };

}

#endif // !defined(BOOST_NO_TEMPLATE_PARTIAL_SPECIALIZATION)

#endif // SCITBX_ARRAY_FAMILY_TINY_TYPES_H
