#include <cctbx/array_family/simple_io.h>
#if !(defined(BOOST_MSVC) && BOOST_MSVC <= 1300) // VC++ 7.0
#include <cctbx/array_family_ew.h>
#endif

using namespace cctbx;

namespace {

# include "tst_af_helpers.cpp"

#if !(defined(BOOST_MSVC) && BOOST_MSVC <= 1300) // VC++ 7.0

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
  exercise_logical(const ArrayType1& a1, ArrayType2& a2)
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
  exercise_arithmetic(const ArrayType1& a1, ArrayType2& a2)
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
    ArrayType& a, const ResultArrayType&, const BoolArrayType&)
  {
    typedef typename ArrayType::value_type element_type;
    a[0] = 0.1;
    a[1] = 0.2;
    a[2] = 0.3;
    { ResultArrayType r = af::acos(a);
      check_true(__LINE__, r[2] == std::acos(a[2])); }
    { ResultArrayType r = af::pow(a, a);
      check_true(__LINE__, r[2] == std::pow(a[2], a[2])); }
    { ResultArrayType r = af::pow(a, a[0]);
      check_true(__LINE__, r[2] == std::pow(a[2], a[0])); }
    { ResultArrayType r = af::pow(a[0], a);
      check_true(__LINE__, r[2] == std::pow(a[0], a[2])); }
    { BoolArrayType r = af::approx_equal_scaled(a, a, element_type(1));
      check_true(__LINE__, r[2] == af::approx_equal_scaled(
        a[2], a[2], element_type(1))); }
    { BoolArrayType r = af::approx_equal_scaled(a, a[0], element_type(1));
      check_true(__LINE__, r[2] == af::approx_equal_scaled(
        a[2], a[0], element_type(1))); }
    { BoolArrayType r = af::approx_equal_scaled(a[0], a, element_type(1));
      check_true(__LINE__, r[2] == af::approx_equal_scaled(
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

  template <typename ElementType>
  struct exercise_complex_special
  {
    static void run()
    {
      af::shared<int> i;
      i.assign(af::tiny<int, 3>(0, 1, 2));
      af::shared<ElementType> e;
      e.assign(af::tiny<ElementType, 3>(0.1, 0.2, 0.3));
      af::shared<std::complex<ElementType> > c;
      c.assign(af::tiny<ElementType, 3>(0.4, 0.5, 2.0));
      { af::shared<ElementType> r = af::real(c);
        check_true(__LINE__, r[2] == std::real(c[2])); }
      { af::shared<ElementType> r = af::imag(c);
        check_true(__LINE__, r[2] == std::imag(c[2])); }
      { af::shared<ElementType> r = af::abs(c);
        check_true(__LINE__, r[2] == std::abs(c[2])); }
      { af::shared<ElementType> r = af::arg(c);
        check_true(__LINE__, r[2] == std::arg(c[2])); }
      { af::shared<ElementType> r = af::norm(c);
        check_true(__LINE__, r[2] == std::norm(c[2])); }
      { af::shared<std::complex<ElementType> > r = af::pow(c, i);
        check_true(__LINE__, r[2] == std::pow(c[2], i[2])); }
      { af::shared<std::complex<ElementType> > r = af::pow(c, i[0]);
        check_true(__LINE__, r[2] == std::pow(c[2], i[0])); }
      { af::shared<std::complex<ElementType> > r = af::pow(c[0], i);
        check_true(__LINE__, r[2] == std::pow(c[0], i[2])); }
      { af::shared<std::complex<ElementType> > r = af::pow(c, e);
        check_true(__LINE__, r[2] == std::pow(c[2], e[2])); }
      { af::shared<std::complex<ElementType> > r = af::pow(c, e[0]);
        check_true(__LINE__, r[2] == std::pow(c[2], e[0])); }
      { af::shared<std::complex<ElementType> > r = af::pow(c[0], e);
        check_true(__LINE__, r[2] == std::pow(c[0], e[2])); }
      { af::shared<std::complex<ElementType> > r = af::pow(c, c);
        check_true(__LINE__, r[2] == std::pow(c[2], c[2])); }
      { af::shared<std::complex<ElementType> > r = af::pow(c, c[0]);
        check_true(__LINE__, r[2] == std::pow(c[2], c[0])); }
      { af::shared<std::complex<ElementType> > r = af::pow(c[0], c);
        check_true(__LINE__, r[2] == std::pow(c[0], c[2])); }
      { af::shared<std::complex<ElementType> > r = af::pow(e, c);
        check_true(__LINE__, r[2] == std::pow(e[2], c[2])); }
      { af::shared<std::complex<ElementType> > r = af::pow(e, c[0]);
        check_true(__LINE__, r[2] == std::pow(e[2], c[0])); }
      { af::shared<std::complex<ElementType> > r = af::pow(e[0], c);
        check_true(__LINE__, r[2] == std::pow(e[0], c[2])); }
      { af::shared<std::complex<ElementType> > r = af::polar(e, e);
        check_true(__LINE__, r[2] == std::polar(e[2], e[2])); }
      { af::shared<std::complex<ElementType> > r = af::polar(e, e[0]);
        check_true(__LINE__, r[2] == std::polar(e[2], e[0])); }
      { af::shared<std::complex<ElementType> > r = af::polar(e[0], e);
        check_true(__LINE__, r[2] == std::polar(e[0], e[2])); }
    }
  };

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
        af::versa<IntType> a1(af::grid<1>(t1.size()));
        af::versa<IntType> a2(af::grid<1>(t2.size()));
        af::versa<FloatType> a3(af::grid<1>(t2.size()));
        af::versa<bool> a4;
        a1.as_base_array().assign(t1);
        a2.as_base_array().assign(t2);
        a3.as_base_array().assign(t2);
        exercise_all(a1, a2, a3, a3, a4);
      }
      {
        if (verbose) std::cout << __LINE__ << std::endl;
        af::versa<IntType> a1(af::grid<1>(t1.size()));
        af::versa<IntType> a2(af::grid<1>(t2.size()));
        af::versa<FloatType> a3(af::grid<1>(t2.size()));
        af::versa<bool> a4;
        a1.as_base_array().assign(t1);
        a2.as_base_array().assign(t2);
        a3.as_base_array().assign(t2);
        af::ref<IntType> r1 = a1.ref();
        af::ref<IntType> r2 = a2.ref();
        af::ref<FloatType> r3 = a3.ref();
        exercise_all(r1, r2, r3, a3, a4);
      }
      exercise_complex_special<double>::run();
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
    operator*=(const a_value<OtherValueType>& other) {
      (*m_value) *= *(other.m_value);
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
#endif // ! VC++ 7.0

} // namespace <anonymous>

int main(int argc, char* argv[])
{
  for(;;) {
#if !(defined(__GNUC__) && __GNUC__ < 3)
#if !(defined(BOOST_MSVC) && BOOST_MSVC <= 1300) // VC++ 7.0
    exercise_main<int, double>::run();
    exercise_main<a_value<int>, double>::run();
    exercise_main<int, a_value<double> >::run();
    exercise_main<a_value<int>, a_value<double> >::run();
#endif
#endif
    if (argc == 1) break;
  }
  std::cout << "Total OK: " << ok_counter << std::endl;
  if (error_counter || verbose) {
    std::cout << "Total Errors: " << error_counter << std::endl;
  }
  return 0;
}
