#ifndef IOTBX_XPLOR_MAP_READER_H
#define IOTBX_XPLOR_MAP_READER_H

#include <string>
#include <iosfwd>
#include <list>
#include <scitbx/array_family/flex_types.h>
#include <scitbx/array_family/tiny_types.h>

namespace iotbx { namespace xplor {

//! Reads density maps stored in xplor format
class map_reader
{
public:
  map_reader() {}

  map_reader(std::string const& file_name, std::size_t n_header_lines,
    scitbx::af::flex_grid<> const& grid);

  explicit map_reader(const std::string &file_name, bool header_only=false);

  scitbx::af::versa<double, scitbx::af::flex_grid<> > data;
  double average;
  double standard_deviation;
  std::list<std::string> title_lines;

  scitbx::af::int3 grid_size, grid_first, grid_last;
  scitbx::af::double6 unit_cell_parameters;

private:

  void load(std::string const& file_name, std::size_t n_header_lines,
    scitbx::af::flex_grid<> const& grid);

  void read(std::istream &text_stream, std::size_t n_header_lines,
    scitbx::af::flex_grid<> const& grid);

  void read(std::istream &text_stream, bool header_only);
};

}}

#endif
