#include <scitbx/array_family/flex_types.h>

using namespace scitbx;

namespace {

# include "tst_af_helpers.cpp"

}

int main(int argc, char* argv[])
{
  {
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

  std::cout << "Total OK: " << ok_counter << std::endl;
  if (error_counter || verbose) {
    std::cout << "Total Errors: " << error_counter << std::endl;
  }

  return 0;
}
