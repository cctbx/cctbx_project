#include <scitbx/array_family/simple_io.h>
#include <scitbx/array_family/tiny_algebra.h>
#include <scitbx/array_family/small_algebra.h>
#include <scitbx/array_family/shared_algebra.h>
#include <scitbx/array_family/versa_algebra.h>
#include <scitbx/array_family/ref_algebra.h>

using namespace scitbx;

namespace {

# include "tst_af_helpers.cpp"

  template <typename ArrayType1,
            typename ArrayType2>
  void
  exercise_ew_bool(ArrayType1 const& a1, ArrayType2 const& a2)
  {
    typedef typename ArrayType1::value_type element_type1;
    verify(__LINE__, a1 == a2,
      af::tiny<bool, 3>(false, false, false));
    verify(__LINE__, a1 == element_type1(1),
      af::tiny<bool, 3>(false, true, false));
    verify(__LINE__, element_type1(1) == a1,
      af::tiny<bool, 3>(false, true, false));

    verify(__LINE__, a1 != a2,
      af::tiny<bool, 3>(true, true, true));
    verify(__LINE__, a1 != element_type1(1),
      af::tiny<bool, 3>(true, false, true));
    verify(__LINE__, element_type1(1) != a1,
      af::tiny<bool, 3>(true, false, true));

    verify(__LINE__, a1 > a2,
      af::tiny<bool, 3>(false, false, false));
    verify(__LINE__, a1 > element_type1(1),
      af::tiny<bool, 3>(false, false, true));
    verify(__LINE__, element_type1(1) > a1,
      af::tiny<bool, 3>(true, false, false));

    verify(__LINE__, a1 < a2,
      af::tiny<bool, 3>(true, true, true));
    verify(__LINE__, a1 < element_type1(1),
      af::tiny<bool, 3>(true, false, false));
    verify(__LINE__, element_type1(1) < a1,
      af::tiny<bool, 3>(false, false, true));

    verify(__LINE__, a1 >= a2,
      af::tiny<bool, 3>(false, false, false));
    verify(__LINE__, a1 >= element_type1(1),
      af::tiny<bool, 3>(false, true, true));
    verify(__LINE__, element_type1(1) >= a1,
      af::tiny<bool, 3>(true, true, false));

    verify(__LINE__, a1 <= a2,
      af::tiny<bool, 3>(true, true, true));
    verify(__LINE__, a1 <= element_type1(1),
      af::tiny<bool, 3>(true, true, false));
    verify(__LINE__, element_type1(1) <= a1,
      af::tiny<bool, 3>(false, true, true));
  }

  template <typename ArrayType1,
            typename ArrayType2>
  void
  exercise_logical(ArrayType1 const& a1, ArrayType2& a2)
  {
    typedef typename ArrayType1::value_type element_type1;
    verify(__LINE__, !a1,
      af::tiny<bool, 3>(true,false,false));
    verify(__LINE__, a1 && a2,
      af::tiny<bool, 3>(false,true,true));
    verify(__LINE__, a1 && element_type1(1),
      af::tiny<bool, 3>(false,true,true));
    verify(__LINE__, element_type1(2) && a1,
      af::tiny<bool, 3>(false,true,true));
  }

  template <typename ArrayType1,
            typename ArrayType2>
  void
  exercise_arithmetic(ArrayType1 const& a1, ArrayType2& a2)
  {
    typedef typename ArrayType1::value_type element_type1;
    typedef typename ArrayType2::value_type element_type2;
    verify(__LINE__, -a1, af::tiny<element_type1, 3>(0,-1,-2));
    verify(__LINE__, a1 + a2, af::tiny<element_type2, 3>(3,5,7));
    verify(__LINE__, a1 + element_type1(1), af::tiny<element_type1, 3>(1,2,3));
    verify(__LINE__, element_type1(2) + a1, af::tiny<element_type1, 3>(2,3,4));
    a2 += a1;
    verify(__LINE__, a2, af::tiny<element_type2, 3>(3,5,7));
    a2 += element_type2(1);
    verify(__LINE__, a2, af::tiny<element_type2, 3>(4,6,8));
  }

  template <typename ArrayType,
            typename BoolArrayType,
            typename ResultArrayType>
  void
  exercise_functions(
    ArrayType& a, ResultArrayType const&, BoolArrayType const&)
  {
    typedef typename ArrayType::value_type element_type;
    typedef typename ResultArrayType::value_type result_element_type;
    a[0] = 0.1;
    a[1] = 0.2;
    a[2] = 0.3;
    { ResultArrayType r = af::absolute(a);
      check_true(__LINE__, fn::absolute(r[2] - fn::absolute(a[2])) < 1.e-6); }
    { ResultArrayType r = af::pow2(a);
      check_true(__LINE__, fn::absolute(r[2] - fn::pow2(a[2])) < 1.e-6); }
    { ResultArrayType r = af::acos(a);
      check_true(__LINE__, fn::absolute(r[2] - std::acos(a[2])) < 1.e-6); }
    { ResultArrayType r = af::pow(a, a);
      check_true(__LINE__, r[2] == std::pow(a[2], a[2])); }
    { ResultArrayType r = af::pow(a, a[0]);
      check_true(__LINE__, r[2] == std::pow(a[2], a[0])); }
    { ResultArrayType r = af::pow(a[0], a);
      check_true(__LINE__, r[2] == std::pow(a[0], a[2])); }
    { BoolArrayType r = af::approx_equal(a, a, element_type(1));
      check_true(__LINE__, r[2] == fn::approx_equal(
        a[2], a[2], element_type(1))); }
    { BoolArrayType r = af::approx_equal(a, a[0], element_type(1));
      check_true(__LINE__, r[2] == fn::approx_equal(
        a[2], a[0], element_type(1))); }
    { BoolArrayType r = af::approx_equal(a[0], a, element_type(1));
      check_true(__LINE__, r[2] == fn::approx_equal(
        a[0], a[2], element_type(1))); }
  }

  template <typename ArrayType1,
            typename ArrayType2,
            typename ArrayType3,
            typename ArrayType4,
            typename ArrayType5>
  void
  exercise_all(ArrayType1& a1,
               ArrayType2& a2,
               ArrayType3& a3,
               ArrayType4& a4,
               ArrayType5& a5)
  {
    exercise_ew_bool(a1, a2);
    exercise_ew_bool(a1, a3);
    exercise_logical(a1, a2);
    exercise_logical(a1, a3);
    exercise_arithmetic(a1, a2);
    exercise_arithmetic(a1, a3);
    exercise_functions(a3, a4, a5);
  }

  template <typename IntType, typename FloatType>
  struct exercise_main
  {
    static void run()
    {
      af::tiny<IntType, 3> t1(0,1,2);
      af::tiny<IntType, 3> t2(3,4,5);
      {
        if (verbose) std::cout << __LINE__ << std::endl;
        af::tiny<IntType, 3> a1 = t1;
        af::tiny<IntType, 3> a2 = t2;
        af::tiny<FloatType, 3> a3 = t2;
        af::tiny<bool, 3> a4;
        exercise_all(a1, a2, a3, a3, a4);
      }
      {
        if (verbose) std::cout << __LINE__ << std::endl;
        af::small<IntType, 3> a1;
        af::small<IntType, 3> a2;
        af::small<FloatType, 3> a3;
        af::small<bool, 3> a4;
        a1.assign(t1);
        a2.assign(t2);
        a3.assign(t2);
        exercise_all(a1, a2, a3, a3, a4);
      }
      {
        if (verbose) std::cout << __LINE__ << std::endl;
        af::shared<IntType> a1;
        af::shared<IntType> a2;
        af::shared<FloatType> a3;
        af::shared<bool> a4;
        a1.assign(t1);
        a2.assign(t2);
        a3.assign(t2);
        exercise_all(a1, a2, a3, a3, a4);
      }
      {
        if (verbose) std::cout << __LINE__ << std::endl;
        af::versa<IntType> a1(t1.size());
        af::versa<IntType> a2(t2.size());
        af::versa<FloatType> a3(t2.size());
        af::versa<bool> a4;
        a1.as_base_array().assign(t1);
        a2.as_base_array().assign(t2);
        a3.as_base_array().assign(t2);
        exercise_all(a1, a2, a3, a3, a4);
      }
      {
        if (verbose) std::cout << __LINE__ << std::endl;
        af::versa<IntType> a1(t1.size());
        af::versa<IntType> a2(t2.size());
        af::versa<FloatType> a3(t2.size());
        af::versa<bool> a4;
        a1.as_base_array().assign(t1);
        a2.as_base_array().assign(t2);
        a3.as_base_array().assign(t2);
        af::ref<IntType> r1 = a1.ref();
        af::ref<IntType> r2 = a2.ref();
        af::ref<FloatType> r3 = a3.ref();
        exercise_all(r1, r2, r3, a3, a4);
      }
    }
  };

  static long a_value_allocation = 0;

  // a type with a non-trivial destructor
  template <typename ValueType>
  struct a_value {

    a_value() { m_init(ValueType()); }
    template <typename OtherValueType>
    a_value(OtherValueType const& v) { m_init(static_cast<ValueType>(v)); }
    a_value(a_value<ValueType> const& other) { m_init(*(other.m_value)); }
    a_value<ValueType>&
    operator=(a_value<ValueType> const& other) {
      m_destory();
      m_init(*(other.m_value));
      return *this;
    }
    ~a_value() { m_destory(); }

    a_value<ValueType>&
    operator++() {
      (*m_value)++;
      return *this;
    }

    template <typename OtherValueType>
    a_value<ValueType>&
    operator+=(a_value<OtherValueType> const& other) {
      (*m_value) += *(other.m_value);
      return *this;
    }

    template <typename OtherValueType>
    a_value<ValueType>&
    operator*=(a_value<OtherValueType> const& other) {
      (*m_value) *= *(other.m_value);
      return *this;
    }

    template <typename OtherValueType>
    a_value<ValueType>&
    operator+=(OtherValueType const& other) {
      (*m_value) += other;
      return *this;
    }

    operator double() const { return *m_value; }

    void m_init(ValueType const& v)
    {
      m_value = new ValueType(v);
      a_value_allocation++;
    }

    void m_destory()
    {
// This destructor is called too often under
// Mac OS 10.2 g++ (GCC) 3.1 20020420 (prerelease)
#if !(defined(__APPLE__) && defined(__MACH__) && __APPLE_CC__ <= 1161)
      delete m_value;
#endif
      a_value_allocation--;
    }

    ValueType* m_value;
  };

} // namespace <anonymous>

int main(int argc, char* argv[])
{
  for(;;) {
#if !(defined(__GNUC__) && __GNUC__ < 3) \
 && !(defined(BOOST_MSVC) && BOOST_MSVC == 1310)
    exercise_main<int, double>::run();
    exercise_main<a_value<int>, double>::run();
    exercise_main<int, a_value<double> >::run();
    exercise_main<a_value<int>, a_value<double> >::run();
#endif
    if (argc == 1) break;
  }
  std::cout << "Total OK: " << ok_counter << std::endl;
  if (error_counter || verbose) {
    std::cout << "Total Errors: " << error_counter << std::endl;
  }
  if (a_value_allocation || verbose) {
    std::cout << "a_value_allocation: " << a_value_allocation << std::endl;
  }
  return 0;
}
