#include <boost/python/module.hpp>
#include <boost/python/class.hpp>
#include <boost/python/def.hpp>
#include <boost/python/args.hpp>
#include <boost/python/list.hpp>
#include <boost/python/make_constructor.hpp>
#include <boost/python/object.hpp>
#include <boost/python/return_by_value.hpp>
#include <boost/python/return_value_policy.hpp>

#include <scitbx/array_family/shared.h>
#include <ucif/parser.h>
#include <iotbx/cif/loop_format.h>
#include <iotbx/error.h>

namespace iotbx { namespace cif {

  struct shared_array_wrapper : ucif::array_wrapper_base
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

    virtual std::string operator[](unsigned const& i) const
    {
      return array[i];
    }

    virtual unsigned size() const
    {
      return array.size();
    }

  };

  struct py_builder : ucif::builder_base
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

    virtual void add_loop(ucif::array_wrapper_base const& loop_headers,
                          std::vector<ucif::array_wrapper_base*> const& values)
    {
      boost::python::list result;
      for (std::size_t i=0; i<loop_headers.size(); i++) {
        result.append(
          dynamic_cast<shared_array_wrapper* const>(values[i])->array);
      }
      builder.attr("add_loop")(
        dynamic_cast<shared_array_wrapper const&>(loop_headers).array,
        result
      );
    }

    virtual void add_data_block(std::string const& data_block_heading)
    {
      builder.attr("add_data_block")(data_block_heading);
    }

    virtual ucif::array_wrapper_base* new_array()
    {
      return new shared_array_wrapper();
    }

  };

  struct parser_wrapper : ucif::parser
  {
    parser_wrapper(ucif::builder_base* builder,
                   std::string input_string,
                   std::string filename="memory",
                   bool strict=true)
    : ucif::parser(builder, input_string, filename, strict) {}

    scitbx::af::shared<std::string>& parser_errors() {
      return dynamic_cast<shared_array_wrapper*>(psr->errors)->array;
    }

    scitbx::af::shared<std::string>& lexer_errors() {
      return dynamic_cast<shared_array_wrapper*>(lxr->errors)->array;
    }

  };

  static iotbx::cif::parser_wrapper* run_cif_parser(
    boost::python::object& builder_,
    std::string input_string,
    std::string filename, bool strict)
  {
    iotbx::cif::py_builder builder(builder_);
    return new iotbx::cif::parser_wrapper(
      &builder, input_string, filename, strict);
  }

  // columns: a list of flex.std_string (the shared <-> flex converters
  // share storage); widths: a sequence of int
  static std::string format_loop_rows_from_python(
    boost::python::list const& columns,
    std::size_t start,
    std::size_t stop,
    std::string const& prefix,
    boost::python::object const& widths)
  {
    using namespace boost::python;
    scitbx::af::shared<scitbx::af::shared<std::string> > cols;
    std::size_t n = len(columns);
    cols.reserve(n);
    for (std::size_t i = 0; i < n; i++) {
      cols.push_back(
        extract<scitbx::af::shared<std::string> >(columns[i])());
    }
    scitbx::af::shared<int> w;
    std::size_t nw = len(widths);
    w.reserve(nw);
    for (std::size_t i = 0; i < nw; i++) {
      w.push_back(extract<int>(widths[i])());
    }
    return format_loop_rows(
      cols.const_ref(), start, stop, prefix, w.const_ref());
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
          (arg("builder"), arg("input_string"),
           arg("filename"), arg("strict")=true)))
        .def("parser_errors", &wt::parser_errors, rbv())
        .def("lexer_errors", &wt::lexer_errors, rbv())
        ;
    }
  };

  void init_module() {
    using namespace boost::python;

    cif_wrapper::wrap("fast_reader");
    def("format_loop_rows", format_loop_rows_from_python,
      (arg("columns"), arg("start"), arg("stop"), arg("prefix"),
       arg("widths")));
    def("column_needs_format_value", column_needs_format_value,
      (arg("column")));
    def("first_non_plain_number", first_non_plain_number,
      (arg("column"), arg("start")));
  }

}}} //iotbx::cif::boost_python


BOOST_PYTHON_MODULE(iotbx_cif_ext)
{
        iotbx::cif::boost_python::init_module();
}
