#include "iotbx/xplor/map_reader.h"

#include <fstream>
#include <cstdlib>
#include <cstdio>
#include <iotbx/error.h>
#include <scitbx/array_family/accessors/c_grid.h>


namespace iotbx { namespace xplor {

namespace af = scitbx::af;

map_reader::map_reader(
  std::string const& file_name,
  std::size_t n_header_lines,
  af::flex_grid<> const& grid)
:
  data(grid, 0)
{
  IOTBX_ASSERT(grid.nd() == 3);
  IOTBX_ASSERT(grid.all().all_gt(0));
  std::ifstream cin(file_name.c_str());
  std::string line;
  for (std::size_t i=0;i<n_header_lines;i++) {
    std::getline(cin, line);
  }
  af::ref<double, af::c_grid<3> > data_ref(
    data.begin(),
    af::c_grid<3>(af::adapt(data.accessor().all())));
  for(std::size_t iz=0;iz<data_ref.accessor()[2];iz++) {
    std::getline(cin, line); // reads section number
    std::size_t i_fld = 6;
    for(std::size_t iy=0;iy<data_ref.accessor()[1];iy++) {
      for(std::size_t ix=0;ix<data_ref.accessor()[0];ix++) {
        if (i_fld == 6) {
          std::getline(cin, line);
          i_fld = 0;
        }
        data_ref(ix,iy,iz) = std::atof(line.substr(i_fld*12,12).c_str());
        i_fld++;
      }
    }
  }
  std::getline(cin, line);
  if (line.size() == 0) {
    average = -1;
    standard_deviation = -1;
  }
  else {
    int expected_9999 = std::atoi(line.substr(0,8).c_str());
    IOTBX_ASSERT(expected_9999 == -9999);
    std::getline(cin, line);
    average = std::atof(line.substr(0,12).c_str());
    standard_deviation = std::atof(line.substr(12,12).c_str());
  }
  cin.close();
}

}}
