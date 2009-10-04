#include <scitbx/array_family/versa.h>
#include <scitbx/array_family/accessors/flex_grid.h>
#include <scitbx/array_family/shared.h>
#include <scitbx/array_family/small.h>
#include <scitbx/array_family/tiny.h>
#include <scitbx/optional_copy.h>
#include <vector>
#include <iostream>
#include <cstdio>

using namespace scitbx;

namespace {

  template <typename ValueType>
  optional_copy<ValueType> const&
  oc_const_ref(optional_copy<ValueType> const& o)
  {
    return o;
  }

  template <typename ValueType>
  void
  exercise(
    ValueType const& v1,
    ValueType const& v4,
    bool value_is_shared)
  {
    typedef optional_copy<ValueType> o;
    {
      o o1;
      SCITBX_ASSERT(!o1);
      SCITBX_ASSERT(o1.get() == 0);
      o o2(o1);
      SCITBX_ASSERT(o2.get() == 0);
      o2 = o1;
      SCITBX_ASSERT(o2.get() == 0);
      o2 = v1;
      SCITBX_ASSERT(o2);
      SCITBX_ASSERT((*o2.get())[0] == 1);
      SCITBX_ASSERT(o2->begin() == o2.get()->begin());
      (*o2)[0] = 2;
      SCITBX_ASSERT((*oc_const_ref(o2))[0] == 2);
      (*o2)[0] = 1;
      SCITBX_ASSERT(o2[0] == 1);
      o1 = o2;
      SCITBX_ASSERT(o1[0] == 1);
      o2.release();
      SCITBX_ASSERT(o2.get() == 0);
      SCITBX_ASSERT(o1[0] == 1);
      SCITBX_ASSERT(oc_const_ref(o2).get() == 0);
      SCITBX_ASSERT(o1->begin() == oc_const_ref(o1).get()->begin());
      SCITBX_ASSERT(oc_const_ref(o1)[0] == 1);
    }
    {
      o o1(v1);
      SCITBX_ASSERT(o1.get() != 0);
      SCITBX_ASSERT(o1[0] == 1);
      o o2(o1);
      SCITBX_ASSERT(o2.get() != 0);
      SCITBX_ASSERT(o2.get() != o1.get());
      SCITBX_ASSERT(o2[0] == 1);
      o2[0] = 2;
      SCITBX_ASSERT(o1[0] == (value_is_shared ? 2 : 1));
      SCITBX_ASSERT(o2[0] == 2);
      o1 = o2;
      SCITBX_ASSERT(o1[0] == 2);
      SCITBX_ASSERT(o2[0] == 2);
      o1[0] = 3;
      SCITBX_ASSERT(o1[0] == 3);
      SCITBX_ASSERT(o2[0] == (value_is_shared ? 3 : 2));
      o2 = v4;
      SCITBX_ASSERT(o1[0] == 3);
      SCITBX_ASSERT(o2[0] == 4);
    }
  };

} // namespace <anonymous>

int main(int argc, char* argv[])
{
  unsigned long n_iter = 1;
  if (argc == 2) {
    std::sscanf(argv[1], "%lu", &n_iter);
  }
  for(unsigned long i_iter=0;i_iter==0||i_iter<n_iter;)
  {
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
      std::vector<int> v1(n, 1);
      std::vector<int> v4(n, 4);
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
