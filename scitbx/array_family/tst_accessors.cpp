#include <scitbx/array_family/flex_types.h>
#include <scitbx/array_family/accessors/c_grid.h>
#include <iostream>

using namespace scitbx;

namespace {

# include "tst_af_helpers.cpp"

}

int main(int argc, char* argv[])
{
  {
    // many more tests via Python
    af::flex_int a;
    check_true(__LINE__, a.size() == 0);
    a = af::flex_int(af::flex_grid<>(1));
    check_true(__LINE__, a.size() == 1);
    check_true(__LINE__, a.accessor().nd() == 1);
    a = af::flex_int(af::flex_grid<>(1, 2));
    check_true(__LINE__, a.size() == 2);
    check_true(__LINE__, a.accessor().nd() == 2);
    a = af::flex_int(af::flex_grid<>(1, 2, 3));
    check_true(__LINE__, a.size() == 6);
    check_true(__LINE__, a.accessor().nd() == 3);
    a = af::flex_int(af::flex_grid<>(1, 2, 3, 4));
    check_true(__LINE__, a.size() == 24);
    check_true(__LINE__, a.accessor().nd() == 4);
    a = af::flex_int(af::flex_grid<>(1, 2, 3, 4, 5));
    check_true(__LINE__, a.size() == 120);
    check_true(__LINE__, a.accessor().nd() == 5);
    a = af::flex_int(af::flex_grid<>(1, 2, 3, 4, 5, 6));
    check_true(__LINE__, a.size() == 720);
    check_true(__LINE__, a.accessor().nd() == 6);
    af::flex_int_const_ref cr = a.const_ref();
    check_true(__LINE__, cr.size() == 720);
    check_true(__LINE__, cr.accessor().nd() == 6);
    af::flex_int_ref r = a.ref();
    check_true(__LINE__, r.size() == 720);
    check_true(__LINE__, r.accessor().nd() == 6);
  }
  {
    af::c_grid<1> a;
    check_true(__LINE__, a.size() == 1);
    check_true(__LINE__, a.size_1d() == 0);
    a = af::c_grid<1>(af::tiny<std::size_t, 1>(3));
    check_true(__LINE__, a.size() == 1);
    check_true(__LINE__, a.size_1d() == 3);
    a = af::c_grid<1>(af::flex_grid<>(3));
    check_true(__LINE__, a.size() == 1);
    check_true(__LINE__, a.size_1d() == 3);
    try {
      a = af::c_grid<1>(af::flex_grid<>(3, 4));
      check_true(__LINE__, false);
    }
    catch (...) {
      check_true(__LINE__, true);
    }
    a = af::c_grid<1>(af::adapt(af::tiny<std::size_t, 1>(5)));
    check_true(__LINE__, a.size() == 1);
    check_true(__LINE__, a.size_1d() == 5);
    verify(__LINE__, a.index_nd(3), af::tiny<std::size_t, 1>(3));
    check_true(__LINE__, a.is_valid_index(af::tiny<std::size_t, 1>(4)));
    check_false(__LINE__, a.is_valid_index(af::tiny<std::size_t, 1>(5)));
    for(std::size_t i=0;i<5;i++) {
      check_true(__LINE__, a(af::tiny<std::size_t, 1>(i)) == i);
    }
  }
  {
    af::c_grid<4> a;
    check_true(__LINE__, a.size() == 4);
    check_true(__LINE__, a.size_1d() == 0);
    a = af::c_grid<4>(af::adapt(af::tiny<std::size_t, 4>(3,2,5,4)));
    check_true(__LINE__, a.size() == 4);
    check_true(__LINE__, a.size_1d() == 120);
    a = af::c_grid<4>(af::flex_grid<>(3,2,5,4));
    check_true(__LINE__, a.size() == 4);
    check_true(__LINE__, a.size_1d() == 120);
    check_false(__LINE__, a.is_valid_index(af::tiny<std::size_t, 4>(3,1,4,3)));
    check_false(__LINE__, a.is_valid_index(af::tiny<std::size_t, 4>(2,1,4,4)));
    std::size_t i = 0;
    af::tiny<std::size_t, 4> j;
    for(j[0]=0;j[0]<3;j[0]++)
    for(j[1]=0;j[1]<2;j[1]++)
    for(j[2]=0;j[2]<5;j[2]++)
    for(j[3]=0;j[3]<4;j[3]++, i++) {
      verify(__LINE__, a.index_nd(i), j);
      check_true(__LINE__, a.is_valid_index(j));
      check_true(__LINE__, a(j) == i);
    }
  }
  {
    af::c_grid<2> a;
    check_true(__LINE__, a.size() == 2);
    check_true(__LINE__, a.size_1d() == 0);
    a = af::c_grid<2>(af::adapt(af::tiny<std::size_t, 2>(3,2)));
    check_true(__LINE__, a.size() == 2);
    check_true(__LINE__, a.size_1d() == 6);
    a = af::c_grid<2>(af::flex_grid<>(3,2));
    check_true(__LINE__, a.size() == 2);
    check_true(__LINE__, a.size_1d() == 6);
    check_false(__LINE__, a.is_valid_index(af::tiny<std::size_t, 2>(3,1)));
    check_false(__LINE__, a.is_valid_index(af::tiny<std::size_t, 2>(2,2)));
    std::size_t i = 0;
    af::tiny<std::size_t, 2> j;
    for(j[0]=0;j[0]<3;j[0]++)
    for(j[1]=0;j[1]<2;j[1]++, i++) {
      verify(__LINE__, a.index_nd(i), j);
      check_true(__LINE__, a.is_valid_index(j));
      check_true(__LINE__, a(j) == i);
    }
  }
  {
    af::c_grid<3> a;
    check_true(__LINE__, a.size() == 3);
    check_true(__LINE__, a.size_1d() == 0);
    a = af::c_grid<3>(af::adapt(af::tiny<std::size_t, 3>(3,2,5)));
    check_true(__LINE__, a.size() == 3);
    check_true(__LINE__, a.size_1d() == 30);
    a = af::c_grid<3>(af::flex_grid<>(3,2,5));
    check_true(__LINE__, a.as_flex_grid().size_1d() == 30);
    check_true(__LINE__, a.size() == 3);
    check_true(__LINE__, a.size_1d() == 30);
    check_false(__LINE__, a.is_valid_index(af::tiny<std::size_t, 3>(3,1,4)));
    check_false(__LINE__, a.is_valid_index(af::tiny<std::size_t, 3>(2,2,4)));
    check_false(__LINE__, a.is_valid_index(af::tiny<std::size_t, 3>(2,1,5)));
    std::size_t i = 0;
    af::tiny<std::size_t, 3> j;
    for(j[0]=0;j[0]<3;j[0]++)
    for(j[1]=0;j[1]<2;j[1]++)
    for(j[2]=0;j[2]<5;j[2]++, i++) {
      verify(__LINE__, a.index_nd(i), j);
      check_true(__LINE__, a.is_valid_index(j));
      check_true(__LINE__, a(j) == i);
    }
  }

  std::cout << "Total OK: " << ok_counter << std::endl;
  if (error_counter || verbose) {
    std::cout << "Total Errors: " << error_counter << std::endl;
  }

  return 0;
}
