#include <cctbx/array_family/simple_io.h>

#include <cctbx/array_family.h>

#include <vector>

using namespace cctbx;

namespace {

  bool verbose = false;

  static std::size_t ok_counter = 0;
  static std::size_t error_counter = 0;

  template <typename VectorType1, typename VectorType2>
  void verify(long line, const VectorType1& v, const VectorType2& a)
  {
    if (v.size() != a.size()) {
      std::cout << line << ": size mismatch: "
                << v.size() << ", " << a.size() << std::endl;
      error_counter++;
      return;
    }
    else {
      for(std::size_t i=0;i<v.size();i++) {
        if (v[i] != a[i]) {
          std::cout << line << ": value mismatch, index " << i << std::endl;
          error_counter++;
          return;
        }
      }
    }
    if (verbose) std::cout << line << ": OK" << std::endl;
    ok_counter++;
  }

  void check_true(long line, bool stat)
  {
    if (!stat) {
      std::cout << line << ": Error" << std::endl;
      error_counter++;
    }
    else {
      if (verbose) std::cout << line << ": OK" << std::endl;
      ok_counter++;
    }
  }

  void check_false(long line, bool stat) {
    check_true(line, !stat);
  }

  template <typename ArrayType>
  struct array_excercise {
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
      ArrayType a3(10, af::no_initialization_flag());
      check_true(__LINE__, a3.size() == 10);
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
      a.assign(u.begin(), u.end());

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
  struct shared_excercise {
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
      ArrayType a3(10, af::no_initialization_flag());
      check_true(__LINE__, a3.size() == 10);
      check_true(__LINE__, a3.capacity() == 10);
      ArrayType a4(10, element_type(123));
      check_true(__LINE__, a4.size() == 10);
      check_true(__LINE__, a4.capacity() == 10);
      ArrayType a5(a4.begin(), a4.begin());
      check_true(__LINE__, a5.capacity() == 0);
      ArrayType a6(a4.begin(), a4.begin() + 3);
      check_true(__LINE__, a6.capacity() == 3);
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
  struct versa_excercise {
    typedef typename ArrayType::value_type element_type;
    static void run() {
      shared_excercise<ArrayType>::run_2();
      run_1();
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
  };

}

int main(void)
{
  if (verbose) std::cout << __LINE__ << ":" << std::endl;
  array_excercise<af::small_plain<int, 128> >::run();
  if (verbose) std::cout << __LINE__ << ":" << std::endl;
  array_excercise<af::small<int, 128> >::run();

  if (verbose) std::cout << __LINE__ << ":" << std::endl;
  array_excercise<af::shared_plain<int> >::run();
  if (verbose) std::cout << __LINE__ << ":" << std::endl;
  shared_excercise<af::shared_plain<int> >::run();

  if (verbose) std::cout << __LINE__ << ":" << std::endl;
  array_excercise<af::shared<int> >::run();
  if (verbose) std::cout << __LINE__ << ":" << std::endl;
  shared_excercise<af::shared<int> >::run();

  if (verbose) std::cout << __LINE__ << ":" << std::endl;
  versa_excercise<af::versa_plain<int>,
                  af::versa_plain<int, af::grid<2> > >::run();
  versa_excercise<af::versa<int>,
                  af::versa<int, af::grid<2> > >::run();

  std::cout << "Total OK: " << ok_counter << std::endl;
  if (error_counter || verbose) {
    std::cout << "Total Errors: " << error_counter << std::endl;
  }

  return 0;
}
