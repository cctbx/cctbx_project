#include <cctbx/boost_python/flex_fwd.h>

#include <boost/python/module.hpp>
#include <boost/python/class.hpp>
#include <boost/python/def.hpp>
#include <boost/python/dict.hpp>
#include <boost/python/list.hpp>
#include <scitbx/array_family/flex_types.h>
#include <scitbx/array_family/shared.h>
#include <scitbx/math/mean_and_variance.h>
#include <scitbx/array_family/boost_python/shared_wrapper.h>
#include <cctbx/examples/merging/xscale_prototype_core.h>
#include <Eigen/Sparse>
#include <boost/python/return_internal_reference.hpp>
#include <boost/math/tools/precision.hpp>

using namespace boost::python;
namespace cctbx{ namespace merging {

double d_product (double const&A, double const&B, double const&dA, double const&dB){
  return dA * B + A * dB;
}

static const double logmax = boost::math::tools::log_max_value<double>()/2.;
static const double logmin = boost::math::tools::log_min_value<double>()/2.;

static boost::python::tuple
task6c_f_g(intensity_data::vecd const& x,
           const int& N_I, const int& N_G,
           const intensity_data & FSIM,
           bool compute_curv) {
  double f = 0.;
  const double* Iptr = x.begin(); // indices into the array of parameters
  const double* Gptr = Iptr + N_I;
  const double* Bptr = Gptr + N_G;

  intensity_data::vecd g(x.size(),0);
  intensity_data::vecd c(x.size(),0);
  double* Iptrg = g.begin(); // indices into the array of derivatives
  double* Gptrg = Iptrg + N_I;
  double* Bptrg = Gptrg + N_G;
  double* Iptrc = c.begin(); // indices into the array of curvatures
  double* Gptrc = Iptrc + N_I;
  double* Bptrc = Gptrc + N_G;

  for (size_t i = 0; i<FSIM.raw_obs.size(); ++i){
      double weight = 1./FSIM.exp_var[i];
      double Gitem  = Gptr[ FSIM.frame[i] ];
      double Bargument = -2.* Bptr[ FSIM.frame[i] ] * FSIM.stol_sq[i];

      if (logmax < Bargument){
        throw SCITBX_ERROR("exp argument greater than logmax");
      }
      if (logmin > Bargument){
        throw SCITBX_ERROR("exp argument less than logmin");
      }

      double Bitem  = std::exp(Bargument); // :=exp(beta)
      double Iitem  = Iptr[ FSIM.miller[i] ];

      double residual = - FSIM.raw_obs[i] +
                        Gitem *
                        Bitem *
                        Iitem;

      f += weight * residual * residual;

      Gptrg[ FSIM.frame[i] ] += weight * residual * Bitem * Iitem;
      Iptrg[ FSIM.miller[i] ]+= weight * residual * Bitem * Gitem;

      double d_Bitem_d_B =    Bitem * ( -2. * FSIM.stol_sq[i] );
      double d_residual_d_B = Gitem * Iitem * d_Bitem_d_B;

      Bptrg[ FSIM.frame[i] ] += weight * residual * d_residual_d_B;

      if (compute_curv) {
        Gptrc[ FSIM.frame[i] ] +=  weight * Bitem * Bitem * Iitem * Iitem;
        Iptrc[ FSIM.miller[i] ]+=  weight * Bitem * Bitem * Gitem * Gitem;
        Bptrc[ FSIM.frame[i] ] +=  weight * d_residual_d_B * d_residual_d_B;
      }
  }

  return (boost::python::make_tuple(f/2., g, c));
}

void f_g_error(std::string const& message){
  throw SCITBX_ERROR(message);
}

namespace boost_python { namespace {

  void
  large_scale_merging_init_module() {
    using namespace boost::python;


    def("task6c_f_g", &cctbx::merging::task6c_f_g);
    def("f_g_error", &cctbx::merging::f_g_error);

    typedef return_value_policy<return_by_value> rbv;
    typedef default_call_policies dcp;

    typedef cctbx::merging::intensity_data cmid;
    class_<cmid>("intensity_data", init<>())
      .add_property("frame",
                    make_getter(&cmid::frame, rbv()),
                    make_setter(&cmid::frame, dcp()))
      .add_property("miller",
                    make_getter(&cmid::miller, rbv()),
                    make_setter(&cmid::miller, dcp()))
      .add_property("raw_obs",
                    make_getter(&cmid::raw_obs, rbv()),
                    make_setter(&cmid::raw_obs, dcp()))
      .add_property("exp_var",
                    make_getter(&cmid::exp_var, rbv()),
                    make_setter(&cmid::exp_var, dcp()))
      .add_property("stol_sq",
                    make_getter(&cmid::stol_sq, rbv()),
                    make_setter(&cmid::stol_sq, dcp()))
      .def("estimate_G", &cmid::estimate_G)
      .def("estimate_I", &cmid::estimate_I)
      .def("reset_mem", &cmid::reset_mem)
    ;

    typedef scaling_common_functions scf;
    class_<scf >(
      "scaling_common_functions", no_init)
      .def("set_cpp_data",&scf::set_cpp_data,
        (arg("fsim")))
    ;

    typedef scitbx::example::non_linear_ls_eigen_wrapper nllsew;
    typedef xscale6e wt6e;
    class_<wt6e,
           bases<nllsew, scf  > >(
      "xscale6e", no_init)
      .def(init<int>(arg("n_parameters")))
      .def("access_cpp_build_up_directly_eigen_eqn",&wt6e::access_cpp_build_up_directly_eigen_eqn,
        (arg("objective_only"),arg("current_values")))
      .def("reset_mem", &wt6e::reset_mem)
    ;

  }

}}
}} // namespace xfel::boost_python::<anonymous>

BOOST_PYTHON_MODULE(cctbx_large_scale_merging_ext)
{
  cctbx::merging::boost_python::large_scale_merging_init_module();

}
