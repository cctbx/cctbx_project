#include <scitbx/array_family/boost_python/flex_fwd.h>

#include <boost/python/def.hpp>
#include <boost/python/args.hpp>

#include <scitbx/math/jacks_expf.h>
#include <scitbx/array_family/shared.h>
#include <tbxx/error_utils.hpp>

namespace scitbx { namespace math {

namespace {

  af::shared<float>
  jacks_expf_a(
    af::const_ref<float> const& array_of_float)
  {
    af::shared<float> result(
      array_of_float.size(),
      af::init_functor_null<float>());
    for(std::size_t i=0;i<array_of_float.size();i++) {
      result[i] = jacks_expf(array_of_float[i]);
    }
    return result;
  }

  union float_unsigned
  {
    float f;
    unsigned u;
  };

  namespace float_bits {
    static const unsigned sign = 0x80000000U;
    static const unsigned expo = 0x7F800000U;
    static const unsigned mant = 0x007FFFFFU;
    static const unsigned expo_size = 8U;
    static const unsigned mant_size = 23U;
    static const int expo_min = -127;
    static const int expo_max = 127;
    static const int expo_bias = 127;
  }

  af::shared<float>
  exercise_jacks_expf(
    bool negative_sign,
    int exponent,
    unsigned mantissa_step_size,
    unsigned j_sample)
  {
    TBXX_ASSERT(sizeof(float) == 4);
    TBXX_ASSERT(sizeof(int) == 4);
    TBXX_ASSERT(exponent >= float_bits::expo_min);
    TBXX_ASSERT(exponent <= float_bits::expo_max);
    float_unsigned value;
    value.f = -3.9572185e-13f;
    TBXX_ASSERT(value.u == 2866726279U);
      // expected to be true for IEEE 754 binary32
    af::shared<float> result;
    value.u = 0U;
    if (negative_sign) {
      value.u |= float_bits::sign;
    }
    unsigned uexp = static_cast<unsigned>(exponent + float_bits::expo_bias);
    value.u |= uexp << float_bits::mant_size;
    unsigned j_end = 1U << float_bits::mant_size;
    float je;
    for(unsigned j=0;j<j_end;j+=mantissa_step_size) {
      value.u &= float_bits::sign | float_bits::expo;
      value.u |= j;
      je = jacks_expf(value.f);
      if (j == 0 || j == j_sample) {
        result.push_back(value.f);
        result.push_back(je);
      }
    }
    result.push_back(value.f);
    result.push_back(je);
    return result;
  }

} // namespace <anonymous>

namespace boost_python {

  void wrap_exp_functions()
  {
    using namespace boost::python;
    def("jacks_expf", jacks_expf_a, (arg("array_of_float")));
    def("exercise_jacks_expf", exercise_jacks_expf, (
      arg("negative_sign"),
      arg("exponent"),
      arg("mantissa_step_size"),
      arg("j_sample")));
  }

}}} // namespace scitbx::math::boost_python
