#include <boost/python/module.hpp>
#include <boost/python/class.hpp>
#include <boost/python/def.hpp>
#include <boost/python/args.hpp>
#include <boost/python/list.hpp>
#include <boost/python/return_by_value.hpp>
#include <boost/python/return_value_policy.hpp>

#include <iotbx/cif/parser.h>
#include <iotbx/error.h>

namespace iotbx { namespace cif {

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
    typedef iotbx::cif::parser wt;

    static void wrap(char const *name) {
      using namespace boost::python;
      typedef return_value_policy<return_by_value> rbv;
      class_<wt, boost::noncopyable>(name, no_init)
        .def(init<std::string, boost::python::object&>(
          (arg("input"), arg("builder"))))
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
