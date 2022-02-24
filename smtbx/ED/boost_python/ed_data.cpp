#include <boost/python/module.hpp>
#include <boost/python/class.hpp>
#include <boost/python/def.hpp>

#include <boost/python/return_value_policy.hpp>
#include <boost/python/return_by_value.hpp>
#include <boost/python/return_internal_reference.hpp>

#include <smtbx/ED/ed_data.h>

namespace cctbx { namespace smtbx { namespace ED {

namespace boost_python {
  using namespace cctbx::smtbx::ED;

  template <typename FloatType>
  struct ed_data_wrapper {
    static void wrap_frame() {
      using namespace boost::python;
      return_value_policy<return_by_value> rbv;
      typedef FrameInfo<FloatType> wt;

      class_<wt, std::auto_ptr<wt> >("frame_info", no_init)
        .def(init<typename wt::cart_t const&,
          FloatType, FloatType, FloatType, FloatType, FloatType>
             ((arg("normal"),
               arg("alpha"), arg("beta"), arg("omega"), arg("angle"), arg("scale"))))
        .add_property("normal", make_getter(&wt::normal, rbv))
        .add_property("alpha", &wt::alpha)
        .add_property("beta", &wt::alpha)
        .add_property("omega", &wt::alpha)
        .add_property("angle", &wt::alpha)
        .add_property("scale", &wt::alpha)
        ;
    }

    static void wrap_beam() {
      using namespace boost::python;
      return_value_policy<return_by_value> rbv;
      return_internal_reference<> rir;
      typedef BeamInfo<FloatType> wt;

      class_<wt, std::auto_ptr<wt> >("beam_info", no_init)
        .def(init<typename wt::parent_ptr_t, const miller::index<>&, FloatType, FloatType>
          ((arg("parent"), arg("index"), arg("I"), arg("sig"))))
        .add_property("parent", make_getter(&wt::parent, rir))
        .add_property("index", make_getter(&wt::index, rbv))
        .add_property("I", &wt::I)
        .add_property("sig", &wt::sig)
        ;
    }

    static void wrap() {
      wrap_frame();
      wrap_beam();
    }
  };


  namespace {
    void init_module() {
      ed_data_wrapper<double>::wrap();
    }
  }

}}}} // namespace cctbx::smtbx::ED::boost_python

BOOST_PYTHON_MODULE(smtbx_ed_data_ext)
{
  cctbx::smtbx::ED::boost_python::init_module();
}
