#include <boost/python/module.hpp>
#include <boost/python/def.hpp>
#include <boost/python/class.hpp>
#include <boost/python/args.hpp>
#include <boost/python/return_value_policy.hpp>
#include <boost/python/return_by_value.hpp>
#include <boost/optional.hpp>

#include <sstream>
#include <vector>
#include <string>
#include <mmtbx/cablam/cablam_align_utils.h>
#include <boost/python.hpp>


namespace mmtbx { namespace cablam {
namespace boost_python {

  void wrap_cablam_align_utils_proxies ()
  {
    using namespace boost::python;
    def("get_similar_regions", get_similar_regions); 
    
    typedef index_mean i_m;
    typedef return_value_policy<return_by_value> rbv;
    class_<i_m>("index_mean", no_init)
      .def(init<int, int, double, int >((
          arg("i_1"),
          arg("i_2"),
          arg("mean"),
          arg("window_length"))))
      .add_property("i_1", make_getter(&i_m::i_1, rbv()))
      .add_property("i_2", make_getter(&i_m::i_2, rbv()))
      .add_property("mean", make_getter(&i_m::mean, rbv()))
      .add_property("window_length", make_getter(&i_m::window_length, rbv()))
    ;
  }
  
  void wrap_cablam_align_utils ()
  {
    wrap_cablam_align_utils_proxies();
  }

}}} 
BOOST_PYTHON_MODULE(mmtbx_cablam_align_utils_ext)
{
  mmtbx::cablam::boost_python::wrap_cablam_align_utils();
}

