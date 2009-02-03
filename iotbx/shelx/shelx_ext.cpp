#include <boost/python/module.hpp>
#include <boost/python/class.hpp>
#include <boost/python/def.hpp>
#include <boost/python/args.hpp>

#include <iotbx/shelx/hklf_simple.h>
#include <iotbx/shelx/hklf.h>
#include <strstream>

namespace iotbx { namespace shelx { namespace boost_python {

  struct hklf_reader
  {
    hklf_reader(std::string const &content, bool strict) {
      std::istringstream input(content);
      iotbx::shelx::hklf_reader delegate(input, strict);
      indices_ = delegate.indices();
      data_ = delegate.data();
      sigmas_ = delegate.sigmas();
      extra_ = delegate.alphas();
    }

    scitbx::af::shared<cctbx::miller::index<> > indices_;
    scitbx::af::shared<double> data_, sigmas_;
    scitbx::af::shared<int> extra_; // batch numbers or phases

    scitbx::af::shared<cctbx::miller::index<> > indices() { return indices_; };

    scitbx::af::shared<double> data() { return data_; }

    scitbx::af::shared<double> sigmas() { return sigmas_; }

    scitbx::af::shared<int> alphas() { return extra_; }

    scitbx::af::shared<int> batch_numbers() { return extra_; }
  };

  template<class WrappedType>
  struct any_hklf_reader_wrapper
  {
    typedef WrappedType wt;

    static void wrap(char const *name) {
      using namespace boost::python;
      class_<wt>(name, no_init)
        .def(init<std::string const &, bool>((arg("content"),
                                              arg("strict")=true)))
        .def("indices", &wt::indices)
        .def("data", &wt::data)
        .def("sigmas", &wt::sigmas)
        .def("alphas", &wt::alphas)
        .def("batch_numbers", &wt::batch_numbers)
        ;
    }
  };

  void init_module() {
    any_hklf_reader_wrapper<simple_hklf_reader>::wrap("simple_hklf_reader");
    any_hklf_reader_wrapper<hklf_reader>::wrap("hklf_reader");
  }

}}} //iotbx::shelx::boost_python


BOOST_PYTHON_MODULE(iotbx_shelx_ext)
{
        iotbx::shelx::boost_python::init_module();
}
