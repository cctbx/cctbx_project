#include <cctbx/array_family/simple_io.h>

#include <cctbx/array_family/tiny.h>
#include <cctbx/array_family/tiny_algebra.h>
#include <cctbx/array_family/small.h>
#include <cctbx/array_family/small_algebra.h>
#include <cctbx/array_family/shared.h>
#include <cctbx/array_family/shared_algebra.h>
#include <cctbx/array_family/versa.h>
#include <cctbx/array_family/versa_algebra.h>

using namespace cctbx;

namespace {

# include "tst_af_helpers.cpp"

  template <typename ArrayType1,
            typename ArrayType2>
  void
  exercise_reducing_bool(const ArrayType1& a1, const ArrayType2& a2)
  {
    typedef typename ArrayType1::value_type element_type1;
    check_false(__LINE__, a1 == a2);
    check_true(__LINE__, a1 == a1);
    check_false(__LINE__, a1 == element_type1(0));
    check_false(__LINE__, element_type1(0) == a1);
    check_true(__LINE__, a1 != a2);
    check_false(__LINE__, a1 != a1);
    check_true(__LINE__, a1 != element_type1(0));
    check_true(__LINE__, element_type1(0) != a1);
    check_false(__LINE__, a1 > a2);
    check_true(__LINE__, a1 > element_type1(0));
    check_false(__LINE__, element_type1(0) > a1);
    check_true(__LINE__, a1 < a2);
    check_false(__LINE__, a1 < element_type1(0));
    check_true(__LINE__, element_type1(0) < a1);
    check_false(__LINE__, a1 >= a2);
    check_true(__LINE__, a1 >= element_type1(0));
    check_false(__LINE__, element_type1(0) >= a1);
    check_true(__LINE__, a1 <= a2);
    check_false(__LINE__, a1 <= element_type1(0));
    check_true(__LINE__, element_type1(0) <= a1);
  }

  template <typename ArrayType1,
            typename ArrayType2>
  void
  exercise_ew_bool(const ArrayType1& a1, const ArrayType2& a2)
  {
    typedef typename ArrayType1::value_type element_type1;
    verify(__LINE__, af::equal_to(a1, a2),
      af::tiny<bool, 3>(false, false, false));
    verify(__LINE__, af::equal_to(a1, element_type1(1)),
      af::tiny<bool, 3>(false, true, false));
    verify(__LINE__, af::equal_to(element_type1(1), a1),
      af::tiny<bool, 3>(false, true, false));

    verify(__LINE__, af::not_equal_to(a1, a2),
      af::tiny<bool, 3>(true, true, true));
    verify(__LINE__, af::not_equal_to(a1, element_type1(1)),
      af::tiny<bool, 3>(true, false, true));
    verify(__LINE__, af::not_equal_to(element_type1(1), a1),
      af::tiny<bool, 3>(true, false, true));

    verify(__LINE__, af::greater(a1, a2),
      af::tiny<bool, 3>(false, false, false));
    verify(__LINE__, af::greater(a1, element_type1(1)),
      af::tiny<bool, 3>(false, false, true));
    verify(__LINE__, af::greater(element_type1(1), a1),
      af::tiny<bool, 3>(true, false, false));

    verify(__LINE__, af::less(a1, a2),
      af::tiny<bool, 3>(true, true, true));
    verify(__LINE__, af::less(a1, element_type1(1)),
      af::tiny<bool, 3>(true, false, false));
    verify(__LINE__, af::less(element_type1(1), a1),
      af::tiny<bool, 3>(false, false, true));

    verify(__LINE__, af::greater_equal(a1, a2),
      af::tiny<bool, 3>(false, false, false));
    verify(__LINE__, af::greater_equal(a1, element_type1(1)),
      af::tiny<bool, 3>(false, true, true));
    verify(__LINE__, af::greater_equal(element_type1(1), a1),
      af::tiny<bool, 3>(true, true, false));

    verify(__LINE__, af::less_equal(a1, a2),
      af::tiny<bool, 3>(true, true, true));
    verify(__LINE__, af::less_equal(a1, element_type1(1)),
      af::tiny<bool, 3>(true, true, false));
    verify(__LINE__, af::less_equal(element_type1(1), a1),
      af::tiny<bool, 3>(false, true, true));
  }

  template <typename ArrayType1,
            typename ArrayType2>
  void
  exercise_arithmetic(const ArrayType1& a1, ArrayType2& a2)
  {
    typedef typename ArrayType1::value_type element_type1;
    typedef typename ArrayType2::value_type element_type2;
    verify(__LINE__, -a1, af::tiny<element_type1, 3>(0,-1,-2));
    verify(__LINE__, !a1, af::tiny<bool, 3>(true,false,false));
    verify(__LINE__, a1 + a2, af::tiny<element_type2, 3>(3,5,7));
    verify(__LINE__, a1 + element_type1(1), af::tiny<element_type1, 3>(1,2,3));
    verify(__LINE__, element_type1(2) + a1, af::tiny<element_type1, 3>(2,3,4));
    a2 += a1;
    verify(__LINE__, a2, af::tiny<element_type2, 3>(3,5,7));
    a2 += element_type2(1);
    verify(__LINE__, a2, af::tiny<element_type2, 3>(4,6,8));
  }

  template <typename ArrayType1,
            typename ArrayType2>
  void
  exercise_functions(const ArrayType1& a1, const ArrayType2& a2)
  {
    typedef typename ArrayType1::value_type element_type1;
    typedef typename ArrayType2::value_type element_type2;
    af::acos(a1);
    af::pow(a1, a2);
    af::pow(a1, element_type1(2));
    af::pow(element_type2(2), a2);
    af::approx_equal_scaled(a1, a2, element_type1(1));
  }

  template <typename ArrayType1,
            typename ArrayType2,
            typename ArrayType3>
  void
  exercise_all(ArrayType1& a1,
               ArrayType2& a2,
               ArrayType3& a3)
  {
    exercise_reducing_bool(a1, a2);
    exercise_reducing_bool(a1, a3);
    exercise_ew_bool(a1, a2);
    exercise_ew_bool(a1, a3);
    exercise_arithmetic(a1, a2);
    exercise_arithmetic(a1, a3);
    exercise_functions(a3, a3);
  }

// XXX exercise_logical
// XXX exercise_complex_special
// XXX exercise_1arg_reductions
// XXX exercise_2arg_reductions

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
        exercise_all(a1, a2, a3);
      }
      {
        if (verbose) std::cout << __LINE__ << std::endl;
        af::small<IntType, 3> a1;
        af::small<IntType, 3> a2;
        af::small<FloatType, 3> a3;
        a1.assign(t1);
        a2.assign(t2);
        a3.assign(t2);
        exercise_all(a1, a2, a3);
      }
      {
        if (verbose) std::cout << __LINE__ << std::endl;
        af::shared<IntType> a1;
        af::shared<IntType> a2;
        af::shared<FloatType> a3;
        a1.assign(t1);
        a2.assign(t2);
        a3.assign(t2);
        exercise_all(a1, a2, a3);
      }
      {
        if (verbose) std::cout << __LINE__ << std::endl;
        af::versa<IntType> a1(af::grid<1>(t1.size()));
        af::versa<IntType> a2(af::grid<1>(t2.size()));
        af::versa<FloatType> a3(af::grid<1>(t1.size()));
        a1.as_shared().assign(t1);
        a2.as_shared().assign(t2);
        a3.as_shared().assign(t2);
        exercise_all(a1, a2, a3);
      }
    }
  };

  // a type with a non-trivial destructor
  template <typename ValueType>
  struct a_value {

    a_value() { m_init(ValueType()); }
    template <typename OtherValueType>
    a_value(const OtherValueType& v) { m_init(v); }
    a_value(const a_value<ValueType>& other) { m_init(*(other.m_value)); }
    a_value<ValueType>&
    operator=(const a_value<ValueType>& other) {
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
    operator+=(const a_value<OtherValueType>& other) {
      (*m_value) += *(other.m_value);
      return *this;
    }

    template <typename OtherValueType>
    a_value<ValueType>&
    operator+=(const OtherValueType& other) {
      (*m_value) += other;
      return *this;
    }

    operator double() const { return *m_value; }

    void m_init(const ValueType& v) { m_value = new ValueType(v); }
    void m_destory() { delete m_value; }
    ValueType* m_value;
  };
}

int main()
{
  exercise_main<int, double>::run();
  exercise_main<a_value<int>, double>::run();
  exercise_main<int, a_value<double> >::run();
  exercise_main<a_value<int>, a_value<double> >::run();
  std::cout << "Total OK: " << ok_counter << std::endl;
  if (error_counter || verbose) {
    std::cout << "Total Errors: " << error_counter << std::endl;
  }

  return 0;
}
