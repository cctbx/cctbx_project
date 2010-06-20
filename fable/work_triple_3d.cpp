#include <fem.hpp> // Fortran EMulation library of fable module

namespace placeholder_please_replace {

using namespace fem::major_types;

struct dynamic_parameters_t
{
  int n3;

  dynamic_parameters_t(
    int argc,
    char const* argv[])
  :
    n3(130)
  {
    fem::dynamic_parameters_from_argv(argc, argv, 1)
      .reset_if_given(n3)
    ;
  }
};

struct common_fields
{
  static const int n1 = 110;
  static const int n2 = 120;
  const int n3;

  arr<float, 3> field1;
  arr<float, 3> field2;
  arr<float, 3> field3;

  common_fields(
    dynamic_parameters_t const& dynamic_parameters)
  :
    n3(dynamic_parameters.n3),
    field1(dimension(n1, n2, n3), fem::fill0),
    field2(dimension(n1, n2, n3), fem::fill0),
    field3(dimension(n1, n2, n3), fem::fill0)
  {}
};

const int common_fields::n1;
const int common_fields::n2;

struct common :
  fem::common,
  common_fields
{
  dynamic_parameters_t dynamic_parameters;

  common(
    dynamic_parameters_t const& dynamic_parameters_)
  :
    common_fields(dynamic_parameters_),
    dynamic_parameters(dynamic_parameters_)
  {}
};

void
set_random(
  common& cmn,
  int& jr,
  arr_ref<float, 3> field,
  int const& n1,
  int const& n2,
  int const& n3)
{
  field(dimension(n1, n2, n3));
  int k = fem::int0;
  int j = fem::int0;
  int i = fem::int0;
  FEM_DO(k, 1, n3) {
    FEM_DO(j, 1, n2) {
      FEM_DO(i, 1, n1) {
        jr = fem::mod(jr * 1366 + 150889, 714025);
        field(i, j, k) = (fem::mod(jr, 20000) - 10000) / 10000.0f;
      }
    }
  }
}

void
find_max_sq(
  common& cmn,
  float& result)
{
  int k = fem::int0;
  const int n3 = cmn.dynamic_parameters.n3;
  int j = fem::int0;
  const int n2 = 120;
  int i = fem::int0;
  const int n1 = 110;
  float f = fem::float0;
  FEM_DO(k, 1, n3) {
    FEM_DO(j, 1, n2) {
      FEM_DO(i, 1, n1) {
        f = fem::pow2(((cmn.field1(i, j, k) - cmn.field2(i, j, k)) * cmn.field3(i, j, k)));
        if (result < f) {
          result = f;
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
  common cmn(dynamic_parameters_t(argc, argv));
  common_write write(cmn);
  int jr = 0;
  const int n1 = 110;
  const int n2 = 120;
  const int n3 = cmn.dynamic_parameters.n3;
  set_random(cmn, jr, cmn.field1, n1, n2, n3);
  set_random(cmn, jr, cmn.field2, n1, n2, n3);
  set_random(cmn, jr, cmn.field3, n1, n2, n3);
  float result = 0;
  int iter = fem::int0;
  FEM_DO(iter, 1, 300) {
    find_max_sq(cmn, result);
  }
  write(6, star), result;
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
