#include <cctbx/xray/thickness.h>

#include <boost/python/class.hpp>
#include <boost/python/def.hpp>
#include <boost/python/args.hpp>


namespace cctbx { namespace xray { namespace boost_python {

namespace {

  template <typename FloatType>
  struct thickness_wrapper {
    static void wrap() {
      typedef thickness<FloatType> wt;
      using namespace boost::python;
      return_internal_reference<> rir;
      class_<wt>("thickness", no_init)
        .def(init<FloatType, optional<bool> >
             ((arg("value"), arg("grad"))))
        .def(init<FloatType, int, bool>
          ((arg("value"), arg("tag"), arg("grad"))))
        .def(init<thickness<FloatType> const&>
             ((arg("source"))))
        .def_readwrite("grad_index", &wt::grad_index)
        .def_readwrite("grad", &wt::grad)
        .def_readwrite("value", &wt::value)
        .def_readonly("tag", &wt::tag)
        .def("deep_copy", &wt::deep_copy)
        ;
    }
  };

} // namespace anonymous

  void wrap_thickness()
  {
    using namespace boost::python;
    thickness_wrapper<double>::wrap();
  }

}}} // namespace cctbx::xray::boost_python
