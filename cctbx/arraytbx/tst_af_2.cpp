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

  template <typename ArrayType>
  void
  exercise_reducing_bool(const ArrayType& a1, const ArrayType& a2)
  {
    typedef typename ArrayType::value_type element_type;
    check_false(__LINE__, a1 == a2);
    check_true(__LINE__, a1 == a1);
    check_false(__LINE__, a1 == int(0));
    check_false(__LINE__, int(0) == a1);
    check_true(__LINE__, a1 != a2);
    check_false(__LINE__, a1 != a1);
    check_true(__LINE__, a1 != int(0));
    check_true(__LINE__, int(0) != a1);
    check_false(__LINE__, a1 > a2);
    check_true(__LINE__, a1 > int(0));
    check_false(__LINE__, int(0) > a1);
    check_true(__LINE__, a1 < a2);
    check_false(__LINE__, a1 < int(0));
    check_true(__LINE__, int(0) < a1);
    check_false(__LINE__, a1 >= a2);
    check_true(__LINE__, a1 >= int(0));
    check_false(__LINE__, int(0) >= a1);
    check_true(__LINE__, a1 <= a2);
    check_false(__LINE__, a1 <= int(0));
    check_true(__LINE__, int(0) <= a1);
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
    verify(__LINE__, -a1, af::tiny<int, 3>(0,-1,-2));
    verify(__LINE__, !a1, af::tiny<bool, 3>(true,false,false));
    verify(__LINE__, a1 + a2, af::tiny<int, 3>(3,5,7));
    verify(__LINE__, a1 / a3, af::tiny<double, 3>(0.,1/4.,2/5.));
    a1 += a3;
    verify(__LINE__, a1, af::tiny<int, 3>(3,5,7));
    a1 += int(1);
    verify(__LINE__, a1, af::tiny<int, 3>(4,6,8));
    af::acos(a3);
    af::pow(a3, 2.);
  }

}

int main(void)
{
  af::tiny<int, 3> t1(0,1,2);
  af::tiny<int, 3> t2(3,4,5);
  {
    if (verbose) std::cout << __LINE__ << std::endl;
    af::tiny<int, 3> a1 = t1;
    af::tiny<int, 3> a2 = t2;
    af::tiny<double, 3> a3 = t2;
    exercise_all(a1, a2, a3);
  }
  {
    if (verbose) std::cout << __LINE__ << std::endl;
    af::small<int, 3> a1;
    af::small<int, 3> a2;
    af::small<double, 3> a3;
    a1.assign(t1);
    a2.assign(t2);
    a3.assign(t2);
    exercise_all(a1, a2, a3);
  }
  {
    if (verbose) std::cout << __LINE__ << std::endl;
    af::shared<int> a1;
    af::shared<int> a2;
    af::shared<double> a3;
    a1.assign(t1);
    a2.assign(t2);
    a3.assign(t2);
    exercise_all(a1, a2, a3);
  }
  {
    if (verbose) std::cout << __LINE__ << std::endl;
    af::versa<int> a1(af::grid<1>(t1.size()));
    af::versa<int> a2(af::grid<1>(t2.size()));
    af::versa<double> a3(af::grid<1>(t1.size()));
    a1.as_shared().assign(t1);
    a2.as_shared().assign(t2);
    a3.as_shared().assign(t2);
    exercise_all(a1, a2, a3);
  }

  std::cout << "Total OK: " << ok_counter << std::endl;
  if (error_counter || verbose) {
    std::cout << "Total Errors: " << error_counter << std::endl;
  }

  return 0;
}
