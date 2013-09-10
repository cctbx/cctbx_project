#ifndef IOTBX_XPLOR_MAP_WRITER_H
#define IOTBX_XPLOR_MAP_WRITER_H

#include <string>
#include <scitbx/array_family/flex_types.h>
#include <scitbx/array_family/accessors/c_grid_padded_periodic.h>
#include <scitbx/array_family/tiny_types.h>

//forward declaration
namespace cctbx { namespace uctbx { class unit_cell; } }

namespace iotbx { namespace xplor {

//! Saves density map in xplor format
void
map_writer_p1_cell(
  std::string const& file_name,
  cctbx::uctbx::unit_cell const& unit_cell,
  scitbx::af::int3 const& gridding_first,
  scitbx::af::int3 const& gridding_last,
  scitbx::af::const_ref<double, scitbx::af::c_grid_padded_periodic<3> > const&,
  double average,
  double standard_deviation);

//! Saves density map in xplor format
void
map_writer_box(
  std::string const& file_name,
  cctbx::uctbx::unit_cell const& unit_cell,
  scitbx::af::const_ref<double, scitbx::af::flex_grid<> > const& data,
  double average,
  double standard_deviation);

//! Saves density map in xplor format, computes mean and esd
void
map_writer(
  std::string const& file_name,
  cctbx::uctbx::unit_cell const& unit_cell,
  scitbx::af::const_ref<double, scitbx::af::flex_grid<> > const& data,
  const scitbx::af::tiny<unsigned,3> &whole_unit_cell_size,
  const std::string &remark="");

}}

#endif
