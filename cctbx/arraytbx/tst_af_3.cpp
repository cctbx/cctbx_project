#include <cctbx/array_family/flex_grid_accessor.h>
#include <cctbx/array_family/versa.h>
#include <cctbx/array_family/loops.h>
#include <cctbx/array_family/simple_io.h>

using namespace cctbx;

namespace {

# include "tst_af_helpers.cpp"

  void exercise()
  {
    af::flex_grid<> g1;
    check_true(__LINE__, g1.size1d() == 0);
    af::flex_grid_default_index_type grid;
    check_true(__LINE__, grid.size() == 0);
    af::flex_grid_default_index_type origin;
    check_true(__LINE__, origin.size() == 0);
    af::flex_grid_default_index_type last;
    check_true(__LINE__, last.size() == 0);
    {
      af::flex_grid<> g2(grid);
      check_true(__LINE__, g2.nd() == 0);
      check_true(__LINE__, g2.size1d() == 0);
      af::flex_grid<> g3(origin, last);
      check_true(__LINE__, g3.nd() == 0);
      check_true(__LINE__, g3.size1d() == 0);
      af::flex_grid<> g4(origin, last, true);
      check_true(__LINE__, g4.nd() == 0);
      check_true(__LINE__, g4.size1d() == 0);
    }
    {
      af::tiny<long, 3> g(2,3,5);
      af::tiny<long, 3> o(1,2,3);
      af::tiny<long, 3> l(3,5,7);
      grid.insert(grid.begin(), g.begin(), g.end());
      origin.insert(origin.begin(), o.begin(), o.end());
      last.insert(last.begin(), l.begin(), l.end());
      af::flex_grid<> g2(grid);
      check_true(__LINE__, g2.nd() == 3);
      check_true(__LINE__, g2.size1d() == 30);
      af::flex_grid<> g3(origin, last, false);
      check_true(__LINE__, g3.nd() == 3);
      check_true(__LINE__, g3.size1d() == 60);
      af::flex_grid<> g4(origin, last, true);
      check_true(__LINE__, g4.nd() == 3);
      check_true(__LINE__, g4.size1d() == 24);
    }
    {
      af::tiny<long, 3> g(2,3,5);
      af::tiny<long, 3> o(1,2,3);
      af::tiny<long, 3> l(3,5,7);
      af::flex_grid<> g2(g);
      check_true(__LINE__, g2.nd() == 3);
      check_true(__LINE__, g2.size1d() == 30);
      af::flex_grid<> g3(o, l, false);
      check_true(__LINE__, g3.nd() == 3);
      check_true(__LINE__, g3.size1d() == 60);
      af::flex_grid<> g4(o, l, true);
      check_true(__LINE__, g4.nd() == 3);
      check_true(__LINE__, g4.size1d() == 24);
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
