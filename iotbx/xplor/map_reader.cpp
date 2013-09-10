#include "iotbx/xplor/map_reader.h"

#include <fstream>
#include <cstdlib>
#include <cstdio>
#include <iotbx/error.h>
#include <scitbx/array_family/accessors/c_grid.h>
#include <boost/algorithm/string/trim.hpp>

namespace iotbx { namespace xplor {

namespace af = scitbx::af;

map_reader::map_reader(
  std::string const& file_name,
  std::size_t n_header_lines,
  af::flex_grid<> const& grid)
:
  data(grid, 0)
{
  this->load(file_name, n_header_lines, grid);
}

map_reader::map_reader(const std::string &file_name, bool header_only)
{
  std::ifstream f(file_name.c_str());
  this->read(f, header_only);
  f.close();
}

void map_reader::read(std::istream &f, bool header_only)
{
  f.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
  std::string line;
  std::getline(f, line);
  boost::algorithm::trim(line);
  int ntitle=std::atoi(line.substr(0, line.find_first_of('!')).c_str());
  while( ntitle-- )
  {
    std::getline(f,line);
    this->title_lines.push_back(line);
  }
  std::getline(f,line);
  std::list<int> values;
  for(short d=0; d<3; ++d)
  {
    int p=d*8*3;
    this->grid_size[d]     = std::atoi(line.substr(p,8).c_str());
    IOTBX_ASSERT( this->grid_size[d]>0 );
    this->grid_first[d] = std::atoi(line.substr(p+8,8).c_str());
    this->grid_last[d]  = std::atoi(line.substr(p+8*2,8).c_str());
  }
  std::getline(f,line);
  for(short p=0; p<6; ++p)
    this->unit_cell_parameters[p] = std::atof(line.substr(p*12,12).c_str());
  std::getline(f,line);
  boost::algorithm::trim(line);
  IOTBX_ASSERT( line=="ZYX" );
  af::flex_grid<> grid(af::adapt(this->grid_first), af::adapt(this->grid_last),
    false);
  if( !header_only )
  {
    data.resize(grid, 0.);
    this->read(f, 0, grid);
  }
}

void map_reader::load(std::string const& file_name, std::size_t n_header_lines,
  af::flex_grid<> const& grid)
{
  std::ifstream cin(file_name.c_str());
  this->read(cin, n_header_lines, grid);
  cin.close();
}

void map_reader::read(std::istream &cin, std::size_t n_header_lines,
  af::flex_grid<> const& grid)
{
  IOTBX_ASSERT(grid.nd() == 3);
  IOTBX_ASSERT(grid.all().all_gt(0));
  for (std::size_t i=0;i<n_header_lines;i++) {
    cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
  }
  std::string line;
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
}

}}
