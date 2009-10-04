#include <boost_adaptbx/tests/tst_optional_copy.h>
#include <vector>
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
    using boost_adaptbx::tst_optional_copy::exercise;
    {
      std::size_t n = 1;
      std::vector<int> v1(n, 1);
      std::vector<int> v4(n, 4);
      exercise(v1, v4, /*value_is_shared*/ false);
    }
    if (n_iter != 0) i_iter++;
  }
  std::cout << "OK\n";
  return 0;
}
