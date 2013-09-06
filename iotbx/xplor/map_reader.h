#ifndef IOTBX_XPLOR_MAP_READER_H
#define IOTBX_XPLOR_MAP_READER_H

#include <string>
#include <scitbx/array_family/flex_types.h>

namespace iotbx { namespace xplor {

//! Reads density maps stored in xplor format
class map_reader
{
public:
  map_reader() {}

  map_reader(
    std::string const& file_name,
    std::size_t n_header_lines,
    scitbx::af::flex_grid<> const& grid);

  scitbx::af::versa<double, scitbx::af::flex_grid<> > data;
  double average;
  double standard_deviation;
};

}}

#endif
