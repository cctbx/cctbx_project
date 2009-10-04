#include <scitbx/array_family/versa.h>
#include <scitbx/array_family/accessors/flex_grid.h>
#include <scitbx/array_family/shared.h>
#include <scitbx/array_family/small.h>
#include <scitbx/array_family/tiny.h>
#include <boost_adaptbx/tests/tst_optional_copy.h>
#include <iostream>
#include <cstdio>

int main(int argc, char* argv[])
{
  unsigned long n_iter = 1;
  if (argc == 2) {
    std::sscanf(argv[1], "%lu", &n_iter);
  }
  for(unsigned long i_iter=0;i_iter==0||i_iter<n_iter;)
  {
    namespace af = scitbx::af;
    using boost_adaptbx::tst_optional_copy::exercise;
    {
      af::tiny<int, 1> v1(1);
      af::tiny<int, 1> v4(4);
      exercise(v1, v4, /*value_is_shared*/ false);
    }
    {
      std::size_t n = 1;
      af::small<int, 1> v1(n, 1);
      af::small<int, 1> v4(n, 4);
      exercise(v1, v4, /*value_is_shared*/ false);
    }
    {
      std::size_t n = 1;
      af::shared<int> v1(n, 1);
      af::shared<int> v4(n, 4);
      exercise(v1, v4, /*value_is_shared*/ true);
    }
    {
      af::flex_grid<> fg(1);
      af::versa<int, af::flex_grid<> > v1(fg, 1);
      af::versa<int, af::flex_grid<> > v4(fg, 4);
      exercise(v1, v4, /*value_is_shared*/ true);
    }
    if (n_iter != 0) i_iter++;
  }
  std::cout << "OK\n";
  return 0;
}
