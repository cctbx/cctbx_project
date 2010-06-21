#include <fem.hpp> // Fortran EMulation library of fable module

namespace placeholder_please_replace {

using namespace fem::major_types;

struct common_fields
{
  static const int n1 = 110;
  static const int n2 = 120;
  static const int n3 = 130;

  fem::arr_ref_dims<3> field1_dims;
#ifdef USE_CMN_ARR
# define ARR fem::cmn_arr
#else
# define ARR arr
#endif
  ARR<float, 3> field1;
  ARR<float, 3> field2;
  ARR<float, 3> field3;
#undef ARR

  common_fields() :
    field1_dims(dimension(n1, n2, n3)),
    field1(field1_dims, fem::fill0),
    field2(field1_dims, fem::fill0),
    field3(field1_dims, fem::fill0)
  {}
};

const int common_fields::n1;
const int common_fields::n2;
const int common_fields::n3;

struct common :
  fem::common,
  common_fields
{};

void
set_random(
  int& jr,
  arr_ref<float> a,
  int const& n)
{
  a(dimension(star));
  int i = fem::int0;
  FEM_DO_SAFE(i, 1, n) {
    jr = fem::mod(jr * 1366 + 150889, 714025);
    a(i) = (fem::mod(jr, 20000) - 10000) / 10000.0f;
  }
}

void
find_max_sq(
  common& cmn,
  arr_ref<int> max_indices)
{
  max_indices(dimension(3));
  int max_f = 0;
  int k = fem::int0;
  const int n3 = 130;
  int j = fem::int0;
  const int n2 = 120;
  int i = fem::int0;
  const int n1 = 110;
  float f = fem::float0;
  FEM_DO_SAFE(k, 1, n3) {
    FEM_DO_SAFE(j, 1, n2) {
      FEM_DO_SAFE(i, 1, n1) {
        f = fem::pow2(((cmn.field1(i, j, k) - cmn.field2(i, j, k)) * cmn.field3(i, j, k)));
        if (max_f < f) {
          max_f = f;
          max_indices(1) = i;
          max_indices(2) = j;
          max_indices(3) = k;
        }
      }
    }
  }
}

void
program_prog(
  int argc,
  char const* argv[])
{
  if (argc != 1) {
    throw std::runtime_error("Unexpected command-line arguments.");
  }
  common cmn;
  common_write write(cmn);
  int jr = 0;
  const int n1 = 110;
  const int n2 = 120;
  const int n3 = 130;
  set_random(jr, cmn.field1, n1 * n2 * n3);
  set_random(jr, cmn.field2, n1 * n2 * n3);
  set_random(jr, cmn.field3, n1 * n2 * n3);
  int i = 0;
  int j = 0;
  int k = 0;
  int iter = fem::int0;
  arr_1d<3, int> max_indices(fem::fill0);
  FEM_DO_SAFE(iter, 1, 300) {
    i = fem::mod(i + iter * 13, n1) + 1;
    j = fem::mod(j + iter * 17, n2) + 1;
    k = fem::mod(k + iter * 19, n3) + 1;
    cmn.field1(i, j, k) = cmn.field2(i, j, k) *  - 1.3f;
    find_max_sq(cmn, max_indices);
  }
  write(6, star), max_indices;
}

} // namespace placeholder_please_replace

int
main(
  int argc,
  char const* argv[])
{
  return fem::main_with_catch(
    argc, argv,
    placeholder_please_replace::program_prog);
}
