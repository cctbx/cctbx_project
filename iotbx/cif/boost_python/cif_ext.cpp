#include <boost/python/module.hpp>
#include <boost/python/class.hpp>
#include <boost/python/def.hpp>
#include <boost/python/args.hpp>
#include <boost/python/list.hpp>
#include <boost/python/make_constructor.hpp>
#include <boost/python/object.hpp>
#include <boost/python/return_by_value.hpp>
#include <boost/python/return_value_policy.hpp>

#include <iotbx/cif/parser.h>
#include <iotbx/error.h>

namespace iotbx { namespace cif {

  struct shared_array_wrapper : array_wrapper_base
  {
    scitbx::af::shared<std::string> array;

    shared_array_wrapper()
    :
    array()
    {}

    virtual void push_back(std::string const& value)
    {
      array.push_back(value);
    }

    virtual std::string operator[](unsigned const& i)
    {
      return array[i];
    }

    virtual unsigned size()
    {
      return array.size();
    }

  };

  struct py_builder : builder_base
  {
    boost::python::object builder;

    py_builder(boost::python::object builder_)
      :
    builder(builder_) {}

    virtual void start_save_frame(std::string const& save_frame_heading)
    {
      builder.attr("start_save_frame")(save_frame_heading);
    }

    virtual void end_save_frame()
    {
      builder.attr("end_save_frame")();
    }

    virtual void add_data_item(std::string const& tag, std::string const& value)
    {
      builder.attr("add_data_item")(tag, value);
    }

    virtual void add_loop(array_wrapper_base const& loop_headers,
                          array_wrapper_base const& values)
    {
      builder.attr("add_loop")(
        dynamic_cast<shared_array_wrapper const&>(loop_headers).array,
        dynamic_cast<shared_array_wrapper const&>(values).array
      );
    }

    virtual void add_data_block(std::string const& data_block_heading)
    {
      builder.attr("add_data_block")(data_block_heading);
    }

    virtual array_wrapper_base* new_array()
    {
      return new shared_array_wrapper();
    }

  };

  struct parser_wrapper : parser
  {
    parser_wrapper(std::string filename, builder_base* builder,
                   bool strict=true)
    : parser(filename, builder, strict) {}

    scitbx::af::shared<std::string>& tree_walker_errors() {
      if (tree_psr != NULL) {
        return dynamic_cast<shared_array_wrapper*>(tree_psr->errors)->array;
      }
      else {
        return *(new scitbx::af::shared<std::string>);
      }
    }

    scitbx::af::shared<std::string>& parser_errors() {
      return dynamic_cast<shared_array_wrapper*>(psr->errors)->array;
    }

    scitbx::af::shared<std::string>& lexer_errors() {
      return dynamic_cast<shared_array_wrapper*>(lxr->errors)->array;
    }

  };

  static iotbx::cif::parser_wrapper* run_cif_parser(
    std::string filename, boost::python::object& builder_, bool strict)
  {
    iotbx::cif::py_builder builder(builder_);
    return new iotbx::cif::parser_wrapper(filename, &builder, strict);
  }

  // Convenience function for sorting a single array of
  // looped data into a given number of columns
  boost::python::list looped_data_as_columns(
      scitbx::af::shared<std::string> data, unsigned n_columns) {
    unsigned n_rows = data.size()/n_columns;
    IOTBX_ASSERT(n_rows * n_columns == data.size());
    boost::python::list result;
    for (std::size_t i=0; i<n_columns; i++) {
      result.append(scitbx::af::shared<std::string>(n_rows));
    }
    for (std::size_t i=0; i<n_columns; i++) {
      for (std::size_t j=0; j<n_rows; j++) {
        result[i][j] = data[i + j*n_columns];
      }
    }
    return result;
  }

namespace boost_python {

  struct cif_wrapper
  {
    typedef iotbx::cif::parser_wrapper wt;

    static void wrap(char const *name) {
      using namespace boost::python;
      typedef return_value_policy<return_by_value> rbv;
      class_<wt, boost::noncopyable>(name, no_init)
        .def("__init__", make_constructor(run_cif_parser,
          default_call_policies(),
          (arg("input"), arg("builder"), arg("strict")=true)))
        .def("tree_walker_errors", &wt::tree_walker_errors, rbv())
        .def("parser_errors", &wt::parser_errors, rbv())
        .def("lexer_errors", &wt::lexer_errors, rbv())
        ;
    }
  };

  void init_module() {
    using namespace boost::python;

    cif_wrapper::wrap("fast_reader");
    def("looped_data_as_columns", looped_data_as_columns);
  }

}}} //iotbx::cif::boost_python


BOOST_PYTHON_MODULE(iotbx_cif_ext)
{
        iotbx::cif::boost_python::init_module();
}
