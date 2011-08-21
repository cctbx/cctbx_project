#ifndef IOTBX_CIF_BUILDER_H
#define IOTBX_CIF_BUILDER_H

#include <string>
#include <scitbx/array_family/shared.h>

namespace iotbx { namespace cif {

// Builders should inherit from builder_base and provide the functions
// defined below.
// See for example py_builder in iotbx/cif/boost_python/cif_ext.pp which in
// turn calls the Python cif_model_builder in iotbx/cif/builders.py.

struct builder_base
{
  virtual ~builder_base() {}
  virtual void start_save_frame(std::string const& save_frame_heading) = 0;
  virtual void end_save_frame() = 0;
  virtual void add_data_item(std::string const& tag, std::string const& value) = 0;
  virtual void add_loop(scitbx::af::shared<std::string> const& loop_headers,
                        scitbx::af::shared<std::string> const& values) = 0;
  virtual void add_data_block(std::string const& data_block_heading) = 0;
};

}} // namespace iotbx::cif

#endif // GUARD
