#include <boost/python/module.hpp>
#include <boost/python/class.hpp>
#include <boost/python/def.hpp>

#include <boost/python/return_value_policy.hpp>
#include <boost/python/return_by_value.hpp>
#include <boost/python/return_internal_reference.hpp>
#include <boost/python/pure_virtual.hpp>

#include <scitbx/array_family/boost_python/shared_wrapper.h>
#include <scitbx/boost_python/std_pair.h>
#include <smtbx/ED/ed_data.h>
#include <smtbx/ED/utils.h>
#include <smtbx/ED/frame_profiler.h>
#include <smtbx/ED/dyn_calculator.h>
#include <smtbx/ED/n_beam.h>

// what is going on here???
//#ifdef __WIN32__
#include <scitbx/boost_python/slice.cpp>
//#endif
namespace bp = boost::python;

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
        .def("get_angles_Sg", &wt::get_angles_Sg)
        .staticmethod("get_angles")
        .def("PL_correctionROD", &wt::PL_correctionROD)
        .def("get_diffraction_angle", &wt::get_diffraction_angle,
          (arg("h"), arg("K"), arg("sweep_angle")=3.0))
        .def("compute_RMf_N", &wt::compute_RMf_N)
        ;
      scitbx::af::boost_python::shared_wrapper<wt, rir_t>::wrap("shared_frame_info");
    }

    static void wrap_beam() {
      using namespace boost::python;
      return_value_policy<return_by_value> rbv;
      typedef return_internal_reference<> rir_t;
      typedef BeamInfo<FloatType> wt;

      class_<wt, std::auto_ptr<wt> >("beam_info", no_init)
        .def(init<const miller::index<> &,
          FloatType, FloatType>
          ((arg("h"), arg("I"),
            arg("sig"))))
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
      return_internal_reference<> rir;
      return_value_policy<return_by_value> rbv;

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
        .def("build_incident_profile", &wt::build_incident_profile)
        .add_property("mi_lookup", make_getter(&wt::mi_lookup, rir))
        .add_property("Fcs_k", make_getter(&wt::Fcs_k, rbv))
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
    ED_UTIL_TYPEDEFS;

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

    typedef typename std::pair< af::shared<miller::index<> >, cmat_t> Ug_rt;
    static Ug_rt build_Ug_matrix_N(
      const af::shared<complex_t>& Fcs_k,
      const lookup_t& mi_lookup,
      const af::shared<miller::index<> >& index_selection,
      const cart_t& K,
      const miller::index<>& h,
      const mat3_t& RMf,
      size_t num, bool use_Sg, FloatType wght)
    {
      Ug_rt rv;
      rv.first = utils<FloatType>::build_Ug_matrix_N(rv.second, Fcs_k,
        mi_lookup, index_selection,
        K, h, RMf, num, use_Sg, wght);
      return rv;
    }

    static void wrap_utils() {
      using namespace boost::python;
      typedef utils<FloatType> wt;
      typedef FloatType(*calc_Sg_t1)(const mat3_t&,
        const miller::index<>&, const cart_t&);

      class_<wt>("utils", no_init)
        .def("build_Ug_matrix_N", &build_Ug_matrix_N)
        .staticmethod("build_Ug_matrix_N")
        .def("calc_Sg", (calc_Sg_t1) &wt::calc_Sg)
        .staticmethod("calc_Sg")
        .def("calc_g", &wt::calc_g)
        .staticmethod("calc_g")
        //.def("generate_index_set", &wt::generate_index_set)
        //.staticmethod("generate_index_set")
        ;
    }
    static void wrap() {
      wrap_excited_beam();
      wrap_utils();
    }
  };
  template <typename FloatType>
  struct dyn_calculator_wrapper {
    ED_UTIL_TYPEDEFS;

    static void wrap_base() {
      using namespace boost::python;
      typedef a_dyn_calculator<FloatType> wt;
      typedef wt& (wt::*reset_t1)(const cmat_t&, const mat3_t&, const cart_t&);
      typedef wt& (wt::*reset_t2)(const cmat_t&, const std::pair<mat3_t, cart_t>&);
      typedef wt& (wt::*reset_t3)(const af::shared<miller::index<> >&,
        const cmat_t&, const mat3_t&, const cart_t&);
      return_value_policy<reference_existing_object> reo;
      return_value_policy<return_by_value> rbv;
      class_<wt, boost::noncopyable>("a_dyn_calculator", no_init)
        .def("reset", (reset_t1)&wt::reset, (
          arg("m"), arg("RMf"), arg("N")), reo)
        .def("reset", (reset_t2) &wt::reset, (
          arg("m"), arg("fi")), reo)
        .def("reset", (reset_t3) &wt::reset, (
          arg("indices"), arg("m"), arg("RMf"), arg("N")),
          reo)
        .def("calc_amps", &wt::calc_amps)
        .def("calc_amps_ext", &wt::calc_amps_ext)
        .def("calc_amps_ext_1", &wt::calc_amps_ext_1)
        .def("build", pure_virtual(&wt::build), reo)
        .def("matrix", &wt::get_matrix, rbv)
        ;
    }

    static void wrap_factory() {
      using namespace boost::python;
      typedef dyn_calculator_factory<FloatType> wt;
      typedef boost::shared_ptr<a_dyn_calculator<FloatType> >(wt::*make_t1)(
        const af::shared<miller::index<> >&,
        const cmat_t&, const cart_t&, const mat3_t&,
        const cart_t&, FloatType) const;
      typedef boost::shared_ptr<a_dyn_calculator<FloatType> >(wt::*make_t2)(
        const af::shared<miller::index<> >&,
        const cart_t&, FloatType) const;

      class_<wt, boost::shared_ptr<wt>,
        boost::noncopyable>("dyn_calculator_factory", no_init)
        .def(init<int>((arg("type"))))
        .def("make", (make_t1) &wt::make,
          (arg("indices"), arg("mat_Ug"), arg("K"), arg("RMf"), arg("N"),
            arg("thickness")))
        .def("make", (make_t2) &wt::make,
          (arg("indices"), arg("K"), arg("thickness")))
        ;
    }

    static void wrap_n_beam() {
      using namespace boost::python;
      typedef dyn_calculator_n_beam<FloatType> wt;
      return_value_policy<reference_existing_object> reo;
      return_value_policy<return_by_value> rbv;

      typedef wt& (wt::*init_t1)(const miller::index<> &, FloatType,
        const af::shared<complex_t> &, const lookup_t &);
      typedef wt& (wt::*init_t2)(const miller::index<>&,
        const mat3_t&, const af::shared<complex_t>&, const lookup_t&);
      
      class_<wt, boost::shared_ptr<wt>,
        boost::noncopyable>("dyn_calculator_n_beam", no_init)
        .def(init<
          size_t, int, const FrameInfo<FloatType>&,
          const cart_t &, FloatType, bool, FloatType>(
            (arg("N"), arg("mat_type"),arg("frame"), arg("K"),
              arg("thickness"), arg("useSG"), arg("wght"))))
        .def("calc_amp", &wt::calc_amp,
          (arg("fi"), arg("idx")=1))
        .def("init", (init_t1) & wt::init, reo)
        .def("init", (init_t2)&wt::init, reo)
        .def("build", &wt::build)
        .def("indices", make_getter(&wt::indices, rbv))
        .def("matrix", &wt::get_matrix, rbv)
        .def("dc", &wt::get_dc, reo)
        ;
    }

    static void wrap() {
      wrap_base();
      wrap_factory();
      wrap_n_beam();
    }
  };

  namespace {
    void init_module() {
      ed_data_wrapper<double>::wrap();
      ed_utils_wrapper<double>::wrap();
      dyn_calculator_wrapper<double>::wrap();

      scitbx::boost_python::RegisterPyPair<scitbx::mat3<double>, scitbx::vec3<double> >();
    }
  }

}}} // namespace cctbx::smtbx::ED::boost_python

BOOST_PYTHON_MODULE(smtbx_ed_data_ext)
{
  smtbx::ED::boost_python::init_module();
}
