#ifndef IOTBX_CIF_BUILDER_H
#define IOTBX_CIF_BUILDER_H

#include <string>

namespace iotbx { namespace cif {


/*! The mininum functions required for an array interface.
    Arrays are used for storing lists of loop headers and looped data items,
    which are passed to the builders.
    Arrays are also used to store lists of any lexer and parser error messages
    that are encountered.
 */
struct array_wrapper_base
{
  virtual ~array_wrapper_base() {}
  virtual void push_back(std::string const&) = 0;
  virtual std::string operator[](unsigned const&) = 0;
  virtual unsigned size() = 0;
};


/*! Builders should inherit from builder_base and provide the functions
    defined below.
    See for example py_builder in iotbx/cif/boost_python/cif_ext.pp which in
    turn calls the Python cif_model_builder in iotbx/cif/builders.py.
 */
struct builder_base
{
  virtual ~builder_base() {}
  virtual void start_save_frame(std::string const& save_frame_heading) = 0;
  virtual void end_save_frame() = 0;
  virtual void add_data_item(std::string const& tag, std::string const& value) = 0;
  virtual void add_loop(array_wrapper_base const& loop_headers,
                        array_wrapper_base const& values) = 0;
  virtual void add_data_block(std::string const& data_block_heading) = 0;
  virtual array_wrapper_base* new_array() = 0;

};

}} // namespace iotbx::cif

#endif // GUARD
