#include <scitbx/mat_ref.h>
#include <iostream>

using namespace scitbx;

namespace {

# include "tst_af_helpers.cpp"

}

int main(int argc, char* argv[])
{
  {
    mat_ref<int> va;
    mat_ref<int> vb(0, mat_grid(0,0));
    mat_ref<int> vc(0, 0, 0);
  }
  {
    std::vector<int> v(6);
    mat_ref<int> a(&*v.begin(), 2, 3);
    check_true(__LINE__, a.grid()[0] == 2);
    check_true(__LINE__, a.grid()[1] == 3);
    check_true(__LINE__, a.n_rows() == 2);
    check_true(__LINE__, a.n_columns() == 3);
    check_true(__LINE__, a.is_same_grid(a));
    mat_ref<int> b(&*v.begin(), 3, 2);
    check_false(__LINE__, a.is_same_grid(b));
    check_true(__LINE__, a == a);
    check_false(__LINE__, a == b);
    check_true(__LINE__, a == 0);
    check_false(__LINE__, a == 1);
    check_true(__LINE__, 0 == a);
    check_false(__LINE__, 1 == a);
    check_false(__LINE__, a != a);
    check_true(__LINE__, a != b);
    check_false(__LINE__, a != 0);
    check_true(__LINE__, a != 1);
    check_false(__LINE__, 0 != a);
    check_true(__LINE__, 1 != a);
  }
  {
    af::tiny<int, 6> t0(1,2,3,4,5,6);
    af::tiny<int, 6> t1(1,2,3,4,5,6);
    mat_ref<int> a(t1.begin(), 2, 3);
    mat_ref<int> b(t1.begin(), 3, 2);
    check_true(__LINE__, a(0,0) == 1);
    check_true(__LINE__, a(0,1) == 2);
    check_true(__LINE__, a(0,2) == 3);
    check_true(__LINE__, a(1,0) == 4);
    check_true(__LINE__, a(1,1) == 5);
    check_true(__LINE__, a(1,2) == 6);
    check_true(__LINE__, b(0,0) == 1);
    check_true(__LINE__, b(0,1) == 2);
    check_true(__LINE__, b(1,0) == 3);
    check_true(__LINE__, b(1,1) == 4);
    check_true(__LINE__, b(2,0) == 5);
    check_true(__LINE__, b(2,1) == 6);
    b = mat_ref<int>(t0.begin(), 2, 3);
    a += b;
    check_true(__LINE__,
      a == mat_ref<int>(af::tiny<int, 6>(2,4,6,8,10,12).begin(), 2, 3));
    a -= b;
    check_true(__LINE__, a == b);
    a += 3;
    check_true(__LINE__,
      a == mat_ref<int>(af::tiny<int, 6>(4,5,6,7,8,9).begin(), 2, 3));
    a -= 3;
    check_true(__LINE__, a == b);
    a *= 3;
    check_true(__LINE__,
      a == mat_ref<int>(af::tiny<int, 6>(3,6,9,12,15,18).begin(), 2, 3));
    a /= 3;
    check_true(__LINE__, a == b);
    a %= 3;
    check_true(__LINE__,
      a == mat_ref<int>(af::tiny<int, 6>(1,2,0,1,2,0).begin(), 2, 3));
  }
  {
    af::tiny<int, 6> ta(1,2,3,4,5,6);
    af::tiny<int, 12> tb;
    for(int i=0;i<12;i++) tb[i] = i+1;
    af::tiny<int, 8> tab;
    mat_ref<int> a(ta.begin(), 2, 3);
    mat_ref<int> b(tb.begin(), 3, 4);
    mat_ref<int> ab(tab.begin(), 2, 4);
    multiply(a, b, ab);
    check_true(__LINE__,
      ab == mat_ref<int>(
        af::tiny<int, 8>(38,44,50,56,83,98,113,128).begin(), 2, 4));
  }
  {
    af::tiny<int, 6> ta(1,2,3,4,5,6);
    af::tiny<int, 6> tb(1,2,3,4);
    af::tiny<int, 6> tc(2,0,0,2);
    mat_ref<int> a(ta.begin(), 2, 3);
    mat_ref<int> b(tb.begin(), 2, 2);
    mat_ref<int> c(tc.begin(), 2, 2);
    check_false(__LINE__, a.is_square());
    check_false(__LINE__, a.is_diagonal());
    check_true(__LINE__, b.is_square());
    check_false(__LINE__, b.is_diagonal());
    check_true(__LINE__, c.is_square());
    check_true(__LINE__, c.is_diagonal());
  }
  {
    af::tiny<int, 6> ta(1,2,3,4,5,6);
    mat_ref<int> a(ta.begin(), 3, 2);
    a.swap_rows(0, 1);
    check_true(__LINE__,
      a == mat_ref<int>(af::tiny<int, 6>(3,4,1,2,5,6).begin(), 3, 2));
    a.swap_rows(1, 2);
    check_true(__LINE__,
      a == mat_ref<int>(af::tiny<int, 6>(3,4,5,6,1,2).begin(), 3, 2));
    a.swap_rows(0, 2);
    check_true(__LINE__,
      a == mat_ref<int>(af::tiny<int, 6>(1,2,5,6,3,4).begin(), 3, 2));
    a.swap_rows(2, 1);
    check_true(__LINE__,
      a == mat_ref<int>(af::tiny<int, 6>(1,2,3,4,5,6).begin(), 3, 2));
    a.swap_columns(0, 1);
    check_true(__LINE__,
      a == mat_ref<int>(af::tiny<int, 6>(2,1,4,3,6,5).begin(), 3, 2));
    a.swap_columns(1, 0);
    check_true(__LINE__,
      a == mat_ref<int>(af::tiny<int, 6>(1,2,3,4,5,6).begin(), 3, 2));
    a = mat_ref<int>(a.begin(), 2, 3);
    a.swap_columns(2, 1);
    check_true(__LINE__,
      a == mat_ref<int>(af::tiny<int, 6>(1,3,2,4,6,5).begin(), 2, 3));
    a.swap_columns(0, 2);
    check_true(__LINE__,
      a == mat_ref<int>(af::tiny<int, 6>(2,3,1,5,6,4).begin(), 2, 3));
  }
  {
    af::tiny<int, 9> ta;
    mat_ref<int> a(ta.begin(), 3, 3);
    a.set_diagonal(3);
    check_true(__LINE__,
      a == mat_ref<int>(af::tiny<int, 9>(3,0,0,0,3,0,0,0,3).begin(), 3, 3));
    check_true(__LINE__, a.is_diagonal());
    a.set_identity();
    check_true(__LINE__,
      a == mat_ref<int>(af::tiny<int, 9>(1,0,0,0,1,0,0,0,1).begin(), 3, 3));
    check_true(__LINE__, a.is_diagonal());
  }
  {
    af::tiny<int, 9> ta(1,2,3,4,5,6,7,8,9);
    mat_ref<int> a(ta.begin(), 3, 3);
    a.transpose_in_place();
    check_true(__LINE__,
      a == mat_ref<int>(af::tiny<int, 9>(1,4,7,2,5,8,3,6,9).begin(), 3, 3));
    a.transpose_in_place();
    check_true(__LINE__,
      a == mat_ref<int>(af::tiny<int, 9>(1,2,3,4,5,6,7,8,9).begin(), 3, 3));
  }
  {
    af::tiny<int, 6> ta(1,2,3,4,5,6);
    mat_ref<int> a(ta.begin(), 2, 3);
    a.transpose_in_place();
    check_true(__LINE__,
      a == mat_ref<int>(af::tiny<int, 6>(1,4,2,5,3,6).begin(), 3, 2));
    a.transpose_in_place();
    check_true(__LINE__,
      a == mat_ref<int>(af::tiny<int, 6>(1,2,3,4,5,6).begin(), 2, 3));
  }

  std::cout << "Total OK: " << ok_counter << std::endl;
  if (error_counter || verbose) {
    std::cout << "Total Errors: " << error_counter << std::endl;
  }

  return 0;
}
