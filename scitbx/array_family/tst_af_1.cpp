#include <scitbx/array_family/tiny_plain_apply.h>
#include <scitbx/array_family/tiny_apply.h>
#include <scitbx/array_family/small_plain_apply.h>
#include <scitbx/array_family/small_apply.h>
#include <scitbx/array_family/shared_plain_apply.h>
#include <scitbx/array_family/shared_apply.h>
#include <scitbx/array_family/versa_plain_apply.h>
#include <scitbx/array_family/versa_apply.h>
#include <scitbx/array_family/ref_apply.h>
#include <scitbx/array_family/misc_functions.h>
#include <scitbx/array_family/ref_reductions.h>
#include <scitbx/array_family/misc_functions.h>
#include <scitbx/array_family/simple_io.h>
#include <boost/bind.hpp>
#include <vector>

using namespace scitbx;

namespace {

# include "tst_af_helpers.cpp"

  template <typename ArrayType>
  struct array_exercise {
    typedef typename ArrayType::value_type element_type;
    static void run() {
      run_1();
      run_2();
    }
    static void run_1() {
      ArrayType a1;
      check_true(__LINE__, a1.size() == 0);
      ArrayType a2(10);
      check_true(__LINE__, a2.size() == 10);
      ArrayType a4(10, element_type(123));
      check_true(__LINE__, a4.size() == 10);
      ArrayType a5(a4.begin(), a4.begin());
      check_true(__LINE__, a5.size() == 0);
      ArrayType a6(a4.begin(), a4.begin() + 3);
      check_true(__LINE__, a6.size() == 3);
    }
    static void run_2() {
      std::vector<element_type> u;
      for(int i=0;i<10;i++) u.push_back(element_type(i+1));

      verify(__LINE__, u, af::make_const_ref(u));
      verify(__LINE__, u, af::make_ref(u));

      ArrayType a1;
      a1.assign(12, element_type(3));
      check_true(__LINE__, a1.size() == 12);
      a1.assign(0, element_type(4));
      check_true(__LINE__, a1.size() == 0);
      a1.assign(24, element_type(5));
      check_true(__LINE__, a1.size() == 24);
      while (a1.size()) a1.pop_back();
      check_true(__LINE__, a1.size() == 0);
      {
        ArrayType a2;
        a2.insert(a2.begin(), element_type(3));
        a2.insert(a2.begin(), element_type(2));
        a2.insert(a2.begin(), element_type(1));
        a2.insert(a2.begin(), element_type(0));
        check_true(__LINE__, a2.size() == 4);
        bool ok = true;
        for(int i=0;i<a2.size();i++) {
          if (a2[i] != i) {
            ok = false;
            break;
          }
        }
        check_true(__LINE__, ok);
        a1.insert(a1.end(), element_type(0));
        a1.insert(a1.end(), element_type(1));
        a1.insert(a1.end(), element_type(2));
        a1.insert(a1.end(), element_type(3));
        verify(__LINE__, a1, a2);
        a1.clear();
        a1.insert(a1.end(), element_type(3));
        a1.insert(a1.begin(), element_type(0));
        a1.insert(a1.begin()+1, element_type(2));
        a1.insert(a1.begin()+1, element_type(1));
        verify(__LINE__, a1, a2);
        check_true(__LINE__, ok && a1.size() == a2.size());
      }

      ArrayType a(a1.begin(), a1.end());
      std::vector<element_type> v(a1.begin(), a1.end());
      verify(__LINE__, v, a);

      a = ArrayType();
      v = std::vector<element_type>();
      verify(__LINE__, v, a);

      v.insert(v.begin(), &u[0], &u[2]);
      a.insert(a.begin(), &u[0], &u[2]);
      verify(__LINE__, v, a);

      v.insert(v.begin()+2, &u[0], &u[2]);
      a.insert(a.begin()+2, &u[0], &u[2]);
      verify(__LINE__, v, a);

      v.insert(v.end(), &u[0], &u[2]);
      a.insert(a.end(), &u[0], &u[2]);

      v.insert(v.begin(), 3, element_type(13));
      a.insert(a.begin(), 3, element_type(13));
      verify(__LINE__, v, a);

      v.insert(v.end(), 3, element_type(14));
      a.insert(a.end(), 3, element_type(14));
      verify(__LINE__, v, a);

      {
        std::size_t n = a.size() + a.capacity();
        if (n > a.max_size()) n = 20;
        v.insert(v.begin()+5, n, element_type(15));
        a.insert(a.begin()+5, n, element_type(15));
        verify(__LINE__, v, a);
      }

      u = v;

      v.erase(v.begin());
      a.erase(a.begin());
      verify(__LINE__, v, a);

      v.erase(v.begin()+8);
      a.erase(a.begin()+8);
      verify(__LINE__, v, a);

      v.erase(v.end()-2);
      a.erase(a.end()-2);
      verify(__LINE__, v, a);

      v.erase(v.end()-1);
      a.erase(a.end()-1);
      verify(__LINE__, v, a);

      while (v.size()) v.erase(v.begin());
      while (a.size()) a.erase(a.begin());
      verify(__LINE__, v, a);

      v = u;
      a.assign(&*u.begin(), &*u.end());

      v.erase(v.begin(), v.begin());
      a.erase(a.begin(), a.begin());
      verify(__LINE__, v, a);

      v.erase(v.begin(), v.begin()+4);
      a.erase(a.begin(), a.begin()+4);
      verify(__LINE__, v, a);

      v.erase(v.end()-4, v.end());
      a.erase(a.end()-4, a.end());
      verify(__LINE__, v, a);

      v.erase(v.begin()+4, v.end()-4);
      a.erase(a.begin()+4, a.end()-4);
      verify(__LINE__, v, a);

      a = ArrayType();
      v = std::vector<element_type>();

      v.resize(50, element_type(3));
      a.resize(50, element_type(3));
      verify(__LINE__, v, a);

      if (a.max_size() >= 500) {
        v.resize(500);
        a.resize(500);
        verify(__LINE__, v, a);
      }

      v.resize(5, element_type(3));
      a.resize(5, element_type(3));
      verify(__LINE__, v, a);

      a = ArrayType();
      a.assign(v);
      verify(__LINE__, v, a);

      std::vector<int> vi(20);
      a.assign(vi);
      verify(__LINE__, vi, a);
    }
  };

  template <typename ArrayType>
  struct shared_exercise {
    typedef typename ArrayType::value_type element_type;
    static void run() {
      run_1();
      run_2();
      run_3();
      run_4();
    }
    static void run_1() {
      ArrayType a1;
      check_true(__LINE__, a1.size() == 0);
      check_true(__LINE__, a1.capacity() == 0);
      check_true(__LINE__, a1.begin() == 0);
      ArrayType a2(10, af::reserve_flag());
      check_true(__LINE__, a2.size() == 0);
      check_true(__LINE__, a2.capacity() == 10);
      ArrayType a4(10, element_type(123));
      check_true(__LINE__, a4.size() == 10);
      check_true(__LINE__, a4.capacity() == 10);
      ArrayType a5(a4.begin(), a4.begin());
      check_true(__LINE__, a5.capacity() == 0);
      ArrayType a6(a4.begin(), a4.begin() + 3);
      check_true(__LINE__, a6.capacity() == 3);
      ArrayType a7(10, af::init_functor_null<element_type>());
      check_true(__LINE__, a7.size() == 10);
      check_true(__LINE__, a7.capacity() == 10);
    }
    static void run_2() {
      ArrayType a1;
      {
        ArrayType a2 = a1.deep_copy();
        verify(__LINE__, a1, a2);
      }
      check_true(__LINE__, a1.use_count() == 1);
      {
        ArrayType a2(10, element_type(123));
        check_true(__LINE__, a2.use_count() == 1);
        check_true(__LINE__, a2.size() == 10);
        check_true(__LINE__, a2.capacity() == 10);
      }
      check_true(__LINE__, a1.use_count() == 1);
      {
        ArrayType a2(a1);
        check_true(__LINE__, a1.use_count() == 2);
        check_true(__LINE__, a2.use_count() == 2);
      }
      check_true(__LINE__, a1.use_count() == 1);
      {
        ArrayType a2(a1);
        ArrayType a3;
        a3 = a2;
        check_true(__LINE__, a3.use_count() == 3);
      }
      check_true(__LINE__, a1.use_count() == 1);
      {
        ArrayType a2(a1);
        ArrayType a3;
        ArrayType a4(a3);
        a3 = a2;
        check_true(__LINE__, a3.use_count() == 3);
        check_true(__LINE__, a4.use_count() == 1);
      }
    }
    static void run_3() {
      ArrayType a1;
      check_true(__LINE__, a1.use_count() == 1);
      {
        ArrayType a2(a1);
        check_true(__LINE__, a1.id() == a2.id());
        a2 = ArrayType();
        check_true(__LINE__, a1.id() != a2.id());
        a2 = a1;
        check_true(__LINE__, a1.id() == a2.id());
      }
      {
        ArrayType a2(a1);
        bool ok = true;
        for(int i=0;i<256;i++) {
          a2.push_back(element_type(i));
          if (a1.size() != a2.size()) {
            ok = false;
            break;
          }
        }
        check_true(__LINE__, ok);
      }
      check_true(__LINE__, a1.use_count() == 1);
      check_true(__LINE__, a1.size() == 256);
      {
        ArrayType a2(a1);
        a1.reserve(1024);
        check_true(__LINE__, a1.capacity() == 1024);
        check_true(__LINE__, a2.capacity() == 1024);
      }
      check_true(__LINE__, a1.size() == 256);
      {
        ArrayType a2;
        a1 = a2;
      }
    }
    static void run_4() {
      ArrayType w1;
      {
        ArrayType a1(3);
        w1 = ArrayType(a1, af::weak_ref_flag());
        check_false(__LINE__, a1.is_weak_ref());
        check_true(__LINE__, w1.is_weak_ref());
        check_true(__LINE__, a1.use_count() == 1);
        check_true(__LINE__, a1.weak_count() == 1);
        {
          ArrayType a2(a1);
          check_true(__LINE__, a2.use_count() == 2);
          check_true(__LINE__, a2.weak_count() == 1);
          ArrayType w2(a1, af::weak_ref_flag());
          check_true(__LINE__, w2.use_count() == 2);
          check_true(__LINE__, w2.weak_count() == 2);
        }
        check_true(__LINE__, a1.use_count() == 1);
        check_true(__LINE__, a1.weak_count() == 1);
        check_true(__LINE__, w1.begin() != 0);
        check_true(__LINE__, w1.size() == 3);
      }
      check_true(__LINE__, w1.use_count() == 0);
      check_true(__LINE__, w1.weak_count() == 1);
      check_true(__LINE__, w1.begin() == 0);
      check_true(__LINE__, w1.size() == 0);
      check_true(__LINE__, w1.capacity() == 0);
      {
        ArrayType w2 = w1;
        check_true(__LINE__, w1.use_count() == 0);
        check_true(__LINE__, w1.weak_count() == 2);
      }
      check_true(__LINE__, w1.use_count() == 0);
      check_true(__LINE__, w1.weak_count() == 1);
      {
        ArrayType w2;
        check_true(__LINE__, w2.use_count() == 1);
        check_true(__LINE__, w2.weak_count() == 0);
        w2 = w1;
        check_true(__LINE__, w1.use_count() == 0);
        check_true(__LINE__, w1.weak_count() == 2);
      }
      check_true(__LINE__, w1.use_count() == 0);
      check_true(__LINE__, w1.weak_count() == 1);
      {
        ArrayType w2(w1, af::weak_ref_flag());
        check_true(__LINE__, w1.use_count() == 0);
        check_true(__LINE__, w1.weak_count() == 2);
      }
      check_true(__LINE__, w1.use_count() == 0);
      check_true(__LINE__, w1.weak_count() == 1);
      {
        ArrayType w2 = w1.weak_ref();
        check_true(__LINE__, w1.use_count() == 0);
        check_true(__LINE__, w1.weak_count() == 2);
      }
      check_true(__LINE__, w1.use_count() == 0);
      check_true(__LINE__, w1.weak_count() == 1);
    }
  };

  template <typename ArrayType, typename AltArrayType>
  struct versa_exercise {
    typedef typename ArrayType::value_type element_type;
    static void run() {
      shared_exercise<ArrayType>::run_2();
      run_1();
      run_2();
    }
    static void run_1() {
      ArrayType a1;
      a1 = ArrayType(af::grid<1>(3));
      a1 = ArrayType(3);
      a1 = ArrayType(af::grid<1>(3), element_type(123));
      a1 = ArrayType(3, element_type(123));
      ArrayType w1(a1, af::weak_ref_flag());
      check_true(__LINE__, w1.use_count() == 1);
      check_true(__LINE__, w1.weak_count() == 1);
      {
        AltArrayType a2(af::grid<2>(3, 4));
        ArrayType a3(a2, af::grid<1>(12));
        ArrayType a4(a2, 12);
        ArrayType a5(a2, af::grid<1>(14), element_type(1));
        ArrayType a6(a2, 16, element_type(2));
        check_true(__LINE__, a2.use_count() == 5);
        check_true(__LINE__, a2.size() == 12);
        check_true(__LINE__, a3.size() == 12);
        check_true(__LINE__, a4.size() == 12);
        check_true(__LINE__, a5.size() == 14);
        check_true(__LINE__, a6.size() == 16);
        check_true(__LINE__, a2.end() - a2.begin() == a2.size());
        check_true(__LINE__, a3.end() - a3.begin() == a3.size());
        check_true(__LINE__, a4.end() - a4.begin() == a4.size());
        check_true(__LINE__, a5.end() - a5.begin() == a5.size());
        check_true(__LINE__, a6.end() - a6.begin() == a6.size());
        a2.resize(af::grid<2>(4, 5), element_type(3));
        ArrayType a2_1d = a2.as_1d();
        check_true(__LINE__, a2.use_count() == 6);
        af::small<element_type, 20> v;
        v.insert(v.end(), 12, element_type());
        v.insert(v.end(), 2, element_type(1));
        v.insert(v.end(), 2, element_type(2));
        v.insert(v.end(), 4, element_type(3));
        verify(__LINE__, a2_1d, v);
        check_true(__LINE__, a2(0,0) == element_type(0));
        check_true(__LINE__, a2(1,1) == element_type(0));
        check_true(__LINE__, a2(3,4) == element_type(3));
        w1 = a2_1d.weak_ref();
        check_true(__LINE__, w1.use_count() == 6);
      }
      check_true(__LINE__, w1.use_count() == 0);
    }
    static void run_2() {
      ArrayType a(10);
      typename ArrayType::base_array_type b = a.as_base_array();
      check_true(__LINE__, a.check_shared_size());
      b.resize(12);
      check_true(__LINE__, a.check_shared_size());
      b.resize(8);
      check_false(__LINE__, a.check_shared_size());
      b.resize(10);
      check_true(__LINE__, a.check_shared_size());
    }
  };

  double foo(int x) { return .1 + x; }

  template <typename ArrayType1, typename ArrayType2>
  void exercise_apply(ArrayType1 const& a1, ArrayType2 const&)
  {
#if !(defined(BOOST_MSVC) && BOOST_MSVC <= 1200) // VC++ 6.0
    ArrayType2 r = af::apply(boost::bind(foo, _1), a1);
    for(std::size_t i=0;i<a1.size();i++) {
      check_true(__LINE__, fn::absolute(r[i] - foo(a1[i])) < 1.e-6);
    }
#endif
  }

  template <typename IntType, typename FloatType>
  struct exercise_apply_all
  {
    static void run()
    {
      af::tiny<IntType, 3> t1(0,1,2);
      {
        if (verbose) std::cout << __LINE__ << std::endl;
        af::tiny<IntType, 3> a1 = t1;
        exercise_apply(a1, af::tiny<FloatType, 3>());
        af::tiny_plain<IntType, 3> a2 = t1;
        exercise_apply(a2, af::tiny_plain<FloatType, 3>());
      }
      {
        if (verbose) std::cout << __LINE__ << std::endl;
        af::small<IntType, 3> a1;
        a1.assign(t1);
        exercise_apply(a1, af::small<FloatType, 3>());
        af::small_plain<IntType, 3> a2;
        a2.assign(t1);
        exercise_apply(a2, af::small_plain<FloatType, 3>());
      }
      {
        if (verbose) std::cout << __LINE__ << std::endl;
        af::shared<IntType> a1;
        a1.assign(t1);
        exercise_apply(a1, af::shared<FloatType>());
        af::shared_plain<IntType> a2;
        a2.assign(t1);
        exercise_apply(a2, af::shared_plain<FloatType>());
      }
      {
        if (verbose) std::cout << __LINE__ << std::endl;
        af::versa<IntType> a1(af::grid<1>(t1.size()));
        a1.as_base_array().assign(t1);
        exercise_apply(a1, af::versa<FloatType>());
        af::versa_plain<IntType> a2(af::grid<1>(t1.size()));
        a2.as_base_array().assign(t1);
        exercise_apply(a2, af::versa_plain<FloatType>());
        exercise_apply(a2.const_ref(), af::versa_plain<FloatType>());
        exercise_apply(a2.ref(), af::versa_plain<FloatType>());
      }
    }
  };

  void exercise_adapt()
  {
    af::tiny<int, 3> a0(1,2,3);
    af::tiny_plain<int, 3> a1(af::adapt(a0));
    verify(__LINE__, a0, a1);
    af::tiny<int, 3> a2(af::adapt(a0));
    verify(__LINE__, a0, a2);
    af::small_plain<int, 5> a3(af::adapt(a0));
    verify(__LINE__, a0, a3);
    af::small<int, 5> a4(af::adapt(a0));
    verify(__LINE__, a0, a4);
    af::shared_plain<int> a5(af::adapt(a0));
    verify(__LINE__, a0, a5);
    af::shared<int> a6(af::adapt(a0));
    verify(__LINE__, a0, a6);
  }

  void exercise_reductions()
  {
    af::tiny<int, 3> a1(0,1,2);
    af::tiny<int, 3> a2(3,4,5);
    af::tiny<int, 3> a5(3,-5,3);
    af::const_ref<int> r1 = a1.const_ref();
    af::const_ref<int> r2 = a2.const_ref();
    af::const_ref<int> r5 = a5.const_ref();
    check_true(__LINE__, r1.all_eq(r1));
    check_false(__LINE__, r1.all_eq(r2));
    check_false(__LINE__, r1.all_eq(0));
    check_true(__LINE__, r1.all_ne(r2));
    check_false(__LINE__, r1.all_ne(r1));
    check_false(__LINE__, r1.all_ne(0));
    check_true(__LINE__, r1.all_ne(3));
    check_true(__LINE__, r1.all_lt(r2));
    check_false(__LINE__, r1.all_lt(r1));
    check_true(__LINE__, r1.all_lt(5));
    check_false(__LINE__, r1.all_lt(2));
    check_true(__LINE__, r2.all_gt(r1));
    check_false(__LINE__, r1.all_gt(r1));
    check_true(__LINE__, r1.all_gt(-1));
    check_false(__LINE__, r1.all_gt(0));
    check_true(__LINE__, r1.all_le(r1));
    check_true(__LINE__, r1.all_le(r2));
    check_true(__LINE__, r1.all_le(2));
    check_false(__LINE__, r1.all_le(0));
    check_true(__LINE__, r1.all_ge(r1));
    check_false(__LINE__, r1.all_ge(r2));
    check_true(__LINE__, r1.all_ge(0));
    check_false(__LINE__, r1.all_ge(2));
    check_true(__LINE__, r1.all_approx_equal(r1, 0));
    check_false(__LINE__, r1.all_approx_equal(r2, 0));
    check_false(__LINE__, r1.all_approx_equal(0, 0));
    check_true(__LINE__, af::order(r1, r2) == -1);
    check_true(__LINE__, af::order(r1, r1) == 0);
    check_true(__LINE__, af::order(r2, r2) == 0);
    check_true(__LINE__, af::order(r2, r1) == 1);
    check_true(__LINE__, af::max_index(r1) == 2);
    check_true(__LINE__, af::min_index(r1) == 0);
    check_true(__LINE__, af::max(r1) == r1[2]);
    check_true(__LINE__, af::min(r1) == r1[0]);
    check_true(__LINE__, af::max_absolute(r5) == -r5[1]);
    check_true(__LINE__, af::sum(r1) == r1[0] + r1[1] + r1[2]);
    check_true(__LINE__, af::product(r1) == r1[0] * r1[1] * r1[2]);
    check_true(__LINE__, af::mean(r1) == (r1[0] + r1[1] + r1[2]) / 3);
    af::tiny<double, 3> a3(3,4,5);
    af::tiny<double, 3> a4(4,5,6);
    af::const_ref<double> r3 = a3.const_ref();
    af::const_ref<double> r4 = a4.const_ref();
    check_true(__LINE__, fn::absolute(
      af::mean_sq(r3)
      - (r3[0]*r3[0] + r3[1]*r3[1] + r3[2]*r3[2]) / 3) < 1.e-6);
    check_true(__LINE__, fn::absolute(
      af::mean_weighted(r3, r4)
      - ((r3[0]*r4[0] + r3[1]*r4[1] + r3[2]*r4[2]) / af::sum(r4))) < 1.e-6);
    check_true(__LINE__, fn::absolute(
      af::mean_sq_weighted(r3, r4)
      - ((  r3[0]*r3[0]*r4[0]
          + r3[1]*r3[1]*r4[1]
          + r3[2]*r3[2]*r4[2]) / af::sum(r4))) < 1.e-6);
  }

}

int main(int argc, char* argv[])
{
  for(;;)
  {
    if (verbose) std::cout << __LINE__ << ":" << std::endl;
    array_exercise<af::small_plain<int, 128> >::run();
    if (verbose) std::cout << __LINE__ << ":" << std::endl;
    array_exercise<af::small<int, 128> >::run();

    if (verbose) std::cout << __LINE__ << ":" << std::endl;
    array_exercise<af::shared_plain<int> >::run();
    if (verbose) std::cout << __LINE__ << ":" << std::endl;
    shared_exercise<af::shared_plain<int> >::run();

    if (verbose) std::cout << __LINE__ << ":" << std::endl;
    array_exercise<af::shared<int> >::run();
    if (verbose) std::cout << __LINE__ << ":" << std::endl;
    shared_exercise<af::shared<int> >::run();

    if (verbose) std::cout << __LINE__ << ":" << std::endl;
    versa_exercise<af::versa_plain<int>,
                    af::versa_plain<int, af::grid<2> > >::run();
    versa_exercise<af::versa<int>,
                    af::versa<int, af::grid<2> > >::run();

    exercise_apply_all<int, double>::run();

    exercise_adapt();
    exercise_reductions();

    if (argc == 1) break;
  }

  std::cout << "Total OK: " << ok_counter << std::endl;
  if (error_counter || verbose) {
    std::cout << "Total Errors: " << error_counter << std::endl;
  }

  return 0;
}
