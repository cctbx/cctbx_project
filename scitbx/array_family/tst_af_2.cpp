#include <scitbx/array_family/flex_grid_accessor.h>
#include <scitbx/array_family/versa.h>
#include <scitbx/array_family/loops.h>
#include <scitbx/array_family/simple_io.h>

using namespace scitbx;

namespace {

# include "tst_af_helpers.cpp"

  void exercise()
  {
    af::flex_grid<> g1;
    check_true(__LINE__, g1.size_1d() == 0);
    af::flex_grid_default_index_type grid;
    check_true(__LINE__, grid.size() == 0);
    af::flex_grid_default_index_type origin;
    check_true(__LINE__, origin.size() == 0);
    af::flex_grid_default_index_type last;
    check_true(__LINE__, last.size() == 0);
    {
      af::flex_grid<> g2(grid);
      check_true(__LINE__, g2.nd() == 0);
      check_true(__LINE__, g2.size_1d() == 0);
      af::flex_grid<> g3(origin, last);
      check_true(__LINE__, g3.nd() == 0);
      check_true(__LINE__, g3.size_1d() == 0);
      af::flex_grid<> g4(origin, last, true);
      check_true(__LINE__, g4.nd() == 0);
      check_true(__LINE__, g4.size_1d() == 0);
    }
    {
      af::tiny<long, 3> g(2,3,5);
      af::tiny<long, 3> o(1,2,3);
      af::tiny<long, 3> l(3,5,7);
      grid.insert(grid.begin(), g.begin(), g.end());
      origin.insert(origin.begin(), o.begin(), o.end());
      last.insert(last.begin(), l.begin(), l.end());
    }
    {
      af::flex_grid<> g2(grid);
      check_true(__LINE__, g2.nd() == 3);
      check_true(__LINE__, g2.size_1d() == 30);
      verify(__LINE__, g2.origin(), af::tiny<long, 3>(0,0,0));
      verify(__LINE__, g2.grid(), grid);
      verify(__LINE__, g2.last(), grid);
      verify(__LINE__, g2.last(true), grid);
      verify(__LINE__, g2.last(false), af::tiny<long, 3>(1,2,4));
      af::flex_grid<> g3(origin, last, false);
      check_true(__LINE__, g3.nd() == 3);
      check_true(__LINE__, g3.size_1d() == 60);
      verify(__LINE__, g3.origin(), origin);
      verify(__LINE__, g3.grid(), af::tiny<long, 3>(3,4,5));
      verify(__LINE__, g3.last(false), last);
      verify(__LINE__, g3.last(true), af::tiny<long, 3>(4,6,8));
      verify(__LINE__, g3.last(), af::tiny<long, 3>(4,6,8));
      af::flex_grid<> g4(origin, last, true);
      check_true(__LINE__, g4.nd() == 3);
      check_true(__LINE__, g4.size_1d() == 24);
      verify(__LINE__, g4.origin(), origin);
      verify(__LINE__, g4.grid(), af::tiny<long, 3>(2,3,4));
      verify(__LINE__, g4.last(false), af::tiny<long, 3>(2,4,6));
      verify(__LINE__, g4.last(true), last);
      verify(__LINE__, g4.last(), last);
      af::flex_grid<> g5(origin, last);
      check_true(__LINE__, g5.size_1d() == 24);
    }
    {
      af::flex_grid<> g2(grid);
      std::size_t i=0;
      af::nested_loop<af::flex_grid_default_index_type> loop(grid);
      for(af::flex_grid_default_index_type const& gpt = loop();
          !loop.over();
          loop.incr(), i++)
      {
        check_true(__LINE__, g2.is_valid_index(gpt));
        check_true(__LINE__, g2(gpt) == i);
      }
      af::flex_grid_default_index_type x(3, 9L);
      check_false(__LINE__, g2.is_valid_index(x));
    }
    {
      bool open_range[] = {true, false};
      for(std::size_t i_o_r=0;i_o_r<2;i_o_r++) {
        af::flex_grid<> g3(origin, last, open_range[i_o_r]);
        std::size_t i=0;
        af::nested_loop<af::flex_grid_default_index_type>
        loop(origin, last, open_range[i_o_r]);
        for(af::flex_grid_default_index_type const& gpt = loop();
            !loop.over();
            loop.incr(), i++)
        {
          check_true(__LINE__, g3.is_valid_index(gpt));
          check_true(__LINE__, g3(gpt) == i);
        }
        af::flex_grid_default_index_type x(3, 9L);
        check_false(__LINE__, g3.is_valid_index(x));
      }
    }
    {
      af::versa<double, af::flex_grid<> > f1;
      check_true(__LINE__, f1.size() == 0);
      af::versa<double, af::flex_grid<> > f2; // XXX have to split resize()!?
      f2.resize(af::flex_grid<>(grid));
      check_true(__LINE__, f2.size() == 30);
      af::versa<double, af::flex_grid<> > f3;
      f3.resize(af::flex_grid<>(origin, last, false));
      check_true(__LINE__, f3.size() == 60);
      af::versa<double, af::flex_grid<> > f4;
      f4.resize(af::flex_grid<>(origin, last, true));
      check_true(__LINE__, f4.size() == 24);
      f4.resize(af::flex_grid<>(grid));
      check_true(__LINE__, f4.size() == 30);
    }
  }
}

int main(int argc, char* argv[])
{
  for(;;)
  {
    exercise();
    if (argc == 1) break;
  }

  std::cout << "Total OK: " << ok_counter << std::endl;
  if (error_counter || verbose) {
    std::cout << "Total Errors: " << error_counter << std::endl;
  }

  return 0;
}
