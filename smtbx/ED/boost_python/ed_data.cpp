#include <boost/python/module.hpp>
#include <boost/python/class.hpp>
#include <boost/python/def.hpp>

#include <boost/python/return_value_policy.hpp>
#include <boost/python/return_by_value.hpp>
#include <boost/python/return_internal_reference.hpp>

#include <scitbx/array_family/boost_python/shared_wrapper.h>
#include <smtbx/ED/ed_data.h>
#include <smtbx/ED/utils.h>
#include <smtbx/ED/frame_profiler.h>

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
      return_value_policy<return_by_value> rbv;
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
        .add_property("normal", make_getter(&wt::normal, rbv))
        .add_property("RMf", make_getter(&wt::RMf, rbv))
        .add_property("alpha", &wt::alpha)
        .add_property("beta", &wt::beta)
        .add_property("omega", &wt::omega)
        .add_property("angle", &wt::angle)
        .add_property("scale", &wt::scale)
        .add_property("indices", make_getter(&wt::indices, rbv))
        .add_property("beams", make_getter(&wt::beams, rbv))
        .add_property("strong_beams", make_getter(&wt::strong_beams, rbv))
        .add_property("strong_measured_beams", make_getter(&wt::strong_measured_beams, rbv))
        .add_property("weak_beams", make_getter(&wt::weak_beams, rbv))
        .add_property("weak_measured_beams", make_getter(&wt::weak_measured_beams, rbv))
        .def("is_excited_index", &wt::is_excited_index)
        .def("is_excited_beam", &wt::is_excited_beam)
        .def("is_fully_covered", &wt::is_excited_beam)
        .def("add_beam", &wt::add_beam)
        .def("set_beams", &wt::set_beams)
        .def("top_up", &wt::top_up)
        .def("unify", &wt::unify)
        .def("add_indices", &wt::add_indices)
        .def("analyse_strength", &wt::analyse_strength)
        .def("update_alpha", &wt::update_alpha)
        .def("update_angles", &wt::update_angles)
        .def("Sg_to_angle", &wt::Sg_to_angle)
        .def("angle_to_Sg", &wt::angle_to_Sg)
        .def("get_int_angles", &wt::get_int_angles)
        .def("get_angles", &wt::get_angles)
        .staticmethod("get_angles")
        ;
      scitbx::af::boost_python::shared_wrapper<wt, rir_t>::wrap("shared_frame_info");
    }

    static void wrap_beam() {
      using namespace boost::python;
      return_value_policy<return_by_value> rbv;
      typedef return_internal_reference<> rir_t;
      typedef BeamInfo<FloatType> wt;

      class_<wt, std::auto_ptr<wt> >("beam_info", no_init)
        .add_property("h", make_getter(&wt::index, rbv))
        .add_property("I", &wt::I)
        .add_property("s", &wt::sig)
        .def_readwrite("diffraction_angle", &wt::diffraction_angle)
        ;
      scitbx::af::boost_python::shared_wrapper<wt, rir_t>::wrap("shared_beam_info");
    }

    static void wrap_peak_profile_point() {
      using namespace boost::python;
      typedef return_internal_reference<> rir_t;
      typedef PeakProfilePoint<FloatType> wt;

      class_<wt, std::auto_ptr<wt> >("peak_profile_point", no_init)
        .add_property("I", &wt::I)
        .add_property("Sg", &wt::Sg)
        .add_property("angle", &wt::angle)
        .add_property("g", &wt::g)
        ;
      scitbx::af::boost_python::shared_wrapper<wt, rir_t>::wrap("peak_profile_point");
    }

    static void wrap_frame_profiler() {
      using namespace boost::python;
      typedef frame_profiler<FloatType> wt;
      typedef refinement::least_squares::f_calc_function_base<FloatType> f_calc_f_t;

      class_<wt, std::auto_ptr<wt> >("frame_profiler", no_init)
        .def(init< const FrameInfo<FloatType> &,
          f_calc_f_t&,
          cctbx::xray::fc_correction<FloatType> const&,
          sgtbx::space_group const&,
          bool,
          cctbx::xray::thickness<FloatType> const&,
          RefinementParams<FloatType> const&>(
            (arg("frame"),
              arg("f_calc_function"), arg("fc_correction"),
              arg("space_group"), arg("anomalous_flag"),
              arg("thickness"),
              arg("params"))))
        .def("build_profile", &wt::build_profile)
        ;
    }

    static void wrap_refinement_params() {
      using namespace boost::python;
      typedef return_internal_reference<> rir_t;
      typedef RefinementParams<FloatType> wt;

      class_<wt, std::auto_ptr<wt> >("refinement_params", no_init)
        .def(init<const af::shared<FloatType> &>
          ((arg("values"))))
        .add_property("Kl_val", &wt::getKl_vac)
        .add_property("Kl", &wt::getKl)
        .add_property("Fc2Ug", &wt::getFc2Ug)
        .add_property("epsilon", &wt::getEpsilon)
        .add_property("matrix_type", &wt::getMatrixType, &wt::setMatrixType)
        .add_property("beam_n", &wt::getBeamN)
        .add_property("thread_n", &wt::getThreadN)
        .add_property("int_span", &wt::getIntSpan)
        .add_property("int_step", &wt::getIntStep)
        .add_property("int_points", &wt::getIntPoints)
        ;
    }

    static void wrap() {
      wrap_frame();
      wrap_beam();
      wrap_peak_profile_point();
      wrap_frame_profiler();
      wrap_refinement_params();
    }
  };

  template <typename FloatType>
  struct ed_utils_wrapper {
    static void wrap_excited_beam() {
      using namespace boost::python;
      typedef typename utils<FloatType>::ExcitedBeam wt;
      class_<wt>("ExcitedBeam", no_init)
        .add_property("weight", &wt::w)
        .add_property("Sg", &wt::Sg)
        .add_property("h", &wt::h)
        .add_property("g", &wt::g)
        ;
      typedef return_internal_reference<> rir_t;
      scitbx::af::boost_python::shared_wrapper<wt, rir_t>::wrap("shared_excited_beams");
    }

    static void wrap_utils() {
      using namespace boost::python;
      typedef utils<FloatType> wt;
      class_<wt>("utils", no_init)
        .def("build_eigen_matrix_recipro", &wt::build_eigen_matrix_recipro)
        .staticmethod("build_eigen_matrix_recipro")

        .def("calc_amps_recipro", &wt::calc_amps_recipro)
        .staticmethod("calc_amps_recipro")

        .def("is_excited_g", &wt::is_excited_g)
        .staticmethod("is_excited_g")
        .def("is_excited_h", &wt::is_excited_h)
        .staticmethod("is_excited_h")
        .def("generate_index_set", &wt::generate_index_set)
        .staticmethod("generate_index_set")
        ;
    }
    static void wrap() {
      wrap_excited_beam();
      wrap_utils();
    }
  };

  namespace {
    void init_module() {
      ed_data_wrapper<double>::wrap();
      ed_utils_wrapper<double>::wrap();
    }
  }

}}} // namespace cctbx::smtbx::ED::boost_python

BOOST_PYTHON_MODULE(smtbx_ed_data_ext)
{
  smtbx::ED::boost_python::init_module();
}
