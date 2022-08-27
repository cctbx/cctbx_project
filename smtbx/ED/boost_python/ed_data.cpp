#include <boost/python/module.hpp>
#include <boost/python/class.hpp>
#include <boost/python/def.hpp>

#include <boost/python/return_value_policy.hpp>
#include <boost/python/return_by_value.hpp>
#include <boost/python/return_internal_reference.hpp>

#include <scitbx/array_family/boost_python/shared_wrapper.h>
#include <smtbx/ED/ed_data.h>

// what is going on here???
//#ifdef __WIN32__
#include <scitbx/boost_python/slice.cpp>
//#endif

namespace smtbx { namespace ED {

namespace boost_python {
  using namespace smtbx::ED;

  template <typename FloatType>
  struct ed_data_wrapper {
    static void wrap_frame() {
      using namespace boost::python;
      typedef return_internal_reference<> rir_t;
      //return_value_policy<return_by_value> rbv;
      typedef FrameInfo<FloatType> wt;

      class_<wt, std::auto_ptr<wt> >("frame_info", no_init)
        .def(init<int, typename wt::cart_t const&,
          FloatType, FloatType, FloatType, FloatType, FloatType,
          typename wt::mat3_t const&>
             ((arg("id"), arg("normal"),
               arg("alpha"), arg("beta"), arg("omega"), arg("angle"), arg("scale"),
               arg("UB"))))
        .def_readonly("id", &wt::id)
        .def_readwrite("tag", &wt::tag)
        .add_property("normal", make_getter(&wt::normal, rir_t()))
        .add_property("RM", make_getter(&wt::RM, rir_t()))
        .add_property("RMf", make_getter(&wt::RMf, rir_t()))
        .add_property("alpha", &wt::alpha)
        .add_property("beta", &wt::beta)
        .add_property("omega", &wt::omega)
        .add_property("angle", &wt::angle)
        .add_property("scale", &wt::scale)
        .def("is_excited", &wt::is_excited);
      scitbx::af::boost_python::shared_wrapper<wt, rir_t>::wrap("shared_frame_info");
    }

    static void wrap_beam() {
      using namespace boost::python;
      return_value_policy<return_by_value> rbv;
      typedef return_internal_reference<> rir_t;
      rir_t rir;
      typedef BeamInfo<FloatType> wt;

      class_<wt, std::auto_ptr<wt> >("beam_info", no_init)
        .def(init<typename wt::parent_t *,
          const miller::index<>&, FloatType, FloatType>
          ((arg("parent"), arg("index"), arg("I"), arg("sig"))))
        .add_property("parent", make_function(&wt::get_parent, rir))
        .add_property("index", make_getter(&wt::index, rbv))
        .add_property("I", &wt::I)
        .add_property("sig", &wt::sig);
      scitbx::af::boost_python::shared_wrapper<wt, rir_t>::wrap("shared_beam_info");
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

}}} // namespace cctbx::smtbx::ED::boost_python

BOOST_PYTHON_MODULE(smtbx_ed_data_ext)
{
  smtbx::ED::boost_python::init_module();
}
