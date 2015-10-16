#include <boost/format.hpp>
#include <iostream>
#include <cmath>

struct show_resolution_save
{
  bool first;
  float ass, bss, css;
  show_resolution_save() : first(true) {}
};

struct common
{
  float a, b, c;
  show_resolution_save show_resolution_sve;
  std::ostream & out_stream;
  common(std::ostream& out_stream_) : out_stream(out_stream_) {}
};

void
show_resolution(common& cmn, int h, int k, int l)
{
  show_resolution_save& sve = cmn.show_resolution_sve;
  if (sve.first) {
    sve.first = false;
    if (cmn.a <= 0 || cmn.b <= 0 || cmn.c <= 0) {
      throw std::runtime_error(
        "invalid unit cell constants.");
    }
    sve.ass = 1/(cmn.a*cmn.a);
    sve.bss = 1/(cmn.b*cmn.b);
    sve.css = 1/(cmn.c*cmn.c);
  }
  float dss = h*h*sve.ass + k*k*sve.bss + l*l*sve.css;
  std::ostream& cout = cmn.out_stream;
  if (dss == 0) {
    cout << boost::format(" %3d %3d %3d     infinity\n")
      % h % k % l;
  }
  else {
    cout << boost::format(" %3d %3d %3d %12.6f\n")
      % h % k % l % std::sqrt(1/dss);
  }
}

void
conv_recipe(common& cmn)
{
  cmn.a = 11.0;
  cmn.b = 12.0;
  cmn.c = 13.0;
  show_resolution(cmn, 0, 0, 0);
  show_resolution(cmn, 1, 2, 3);
}

int
main()
{
  common cmn(std::cout);
  conv_recipe(cmn);
  return 0;
}
