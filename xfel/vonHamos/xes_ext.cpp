#include <cctbx/boost_python/flex_fwd.h>

#include <boost/python/module.hpp>
#include <boost/python/class.hpp>
#include <boost/python/def.hpp>
#include <scitbx/array_family/flex_types.h>
#include <scitbx/array_family/shared.h>
#include <scitbx/array_family/boost_python/shared_wrapper.h>

#include <vector>
#include <cmath>
#include <stdexcept>
#include <scitbx/lstbx/normal_equations.h>

using namespace boost::python;
namespace xfel{

namespace xes { namespace example {
  class gaussian_fit_inheriting_from_non_linear_ls: public scitbx::lstbx::normal_equations::non_linear_ls<double> {

    public:
      gaussian_fit_inheriting_from_non_linear_ls(int n_parameters):
        scitbx::lstbx::normal_equations::non_linear_ls<double>(n_parameters),
        jacobian_one_row(n_parameters)
        {}

      void set_cpp_data(scitbx::af::shared<double> x,
                        scitbx::af::shared<double> y,
                        double const& g,
                        double const& s) {
                        free_x=x, free_y=y, gain_to_sigma=g, sigmafac=s; sigmafac_sq=s*s;
                        residual_terms_0 = std::vector<double>(x.size());
                        residual_terms_1 = std::vector<double>(x.size());
                        residuals = scitbx::af::shared<double>(x.size());
      }

      void access_cpp_build_up_directly(bool objective_only, scitbx::af::shared<double> current_values) {
        fvec_callable(current_values);
        if (objective_only) {
          add_residuals(residuals.const_ref(), scitbx::af::shared<double>().const_ref());
        }else{
          // add one of the normal equations for each observation
          for (int ix = 0; ix < free_x.size(); ++ix) {
            /*
              z_mean = current_values[0]
              z_ampl = current_values[1]
              z_sigm = current_values[2]
              o_ampl = current_values[3]
              o_mean = z_mean + z_sigm * GAIN_TO_SIGMA
              o_sigm = z_sigm * SIGMAFAC
              # leaving out small cross terms where one-photon peak influences
              # derivatives with respect to z_mean and z_sigm.
              # XXX not accounting for charge sharing at all
             */
            double udiff = free_x[ix] - current_values[0];
            double Afactor = udiff/(current_values[2] * current_values[2]);
            double Sfactor = Afactor * udiff/current_values[2];
            jacobian_one_row[0] = ( residual_terms_0[ix] * Afactor );
            jacobian_one_row[1] = ( residual_terms_0[ix] / current_values[1] );
            jacobian_one_row[2] = ( residual_terms_0[ix] * Sfactor );
            jacobian_one_row[3] = ( residual_terms_1[ix] / current_values[3] );

            add_equation(residuals[ix], jacobian_one_row.const_ref(), 1.);
          }
        }
      }

      void
      fvec_callable(scitbx::af::shared<double> &current_values) {
        /*
          z_mean = current_values[0]
          z_ampl = current_values[1]
          z_sigm = current_values[2]
          o_ampl = current_values[3]
          o_mean = z_mean + z_sigm * GAIN_TO_SIGMA
          o_sigm = z_sigm * SIGMAFAC

          # model minus obs
          # sqrt2pi_inv = 1./math.sqrt(2.*math.pi)
          # the gaussian function is sqrt2pi_inv * exp( - (x-mean)**2 / (2.*(sigma**2)))/sigma
          # take off the coefficient sqrt2pi_inv / sigma, use ampl
        */

        for (int ix = 0; ix < free_x.size(); ++ix) {
          double model=0;
          //zero-photon Gaussian
          double diff = free_x[ix] - current_values[0];
          residual_terms_0[ix] = (
            current_values[1] * std::exp( - (diff*diff) / (2. * current_values[2] * current_values[2])));
          model += residual_terms_0[ix];

          //one-photon Gaussian
          diff = free_x[ix] - ( current_values[0] + current_values[2] * gain_to_sigma );
          residual_terms_1[ix] = (
            current_values[3] * std::exp( - (diff*diff) / (2. * current_values[2] * current_values[2] * sigmafac_sq)));
          model += residual_terms_1[ix];

          //terms = [
          //ampl * flex.exp(-flex.pow2(free_x - mean) / (2.*sigm*sigm))
          //for mean,ampl,sigm in [(z_mean,z_ampl,z_sigm),(o_mean,o_ampl,o_sigm)]]
          //#print "residual", math.sqrt(flex.sum((model-free_y)*(model-free_y)))
          residuals[ix] = ( model - free_y[ix] );
        }
      }

    private:
      scitbx::af::shared<double> free_x, free_y, jacobian_one_row, residuals;
      double gain_to_sigma, sigmafac, sigmafac_sq;
      std::vector<double> residual_terms_0, residual_terms_1;

  };

  class gaussian_3fit_inheriting_from_non_linear_ls: public scitbx::lstbx::normal_equations::non_linear_ls<double> {

    public:
      gaussian_3fit_inheriting_from_non_linear_ls(int n_parameters):
        scitbx::lstbx::normal_equations::non_linear_ls<double>(n_parameters),
        jacobian_one_row(n_parameters)
        {}

      void set_cpp_data(scitbx::af::shared<double> const& val,
                        scitbx::af::shared<double> const& x,
                        scitbx::af::shared<double> const& y) {
                        free_x=x, free_y=y;
                        residual_terms_0 = std::vector<double>(x.size());
                        residual_terms_1 = std::vector<double>(x.size());
                        residual_terms_2 = std::vector<double>(x.size());
                        residuals = scitbx::af::shared<double>(x.size());
                        constants = val;
      }

      void access_cpp_build_up_directly(bool objective_only, scitbx::af::shared<double> current_values) {

        fvec_callable(current_values);

        if (objective_only) {

          add_residuals(residuals.const_ref(), scitbx::af::shared<double>().const_ref());
        }else{
          // add one of the normal equations for each observation
          for (int ix = 0; ix < free_x.size(); ++ix) {
            /*
              z_mean = current_values[0]
              z_ampl = current_values[1]
              z_sigm = current_values[2]
              inelast_mean = zmean + zsigm * constants[0]
              inelast_ampl = current_values[3]
              inelast_sigm = (constants[0]/constants[1])*(elast_sigm[5] - zsigm[2])+zsigm[2]
              elast_mean = z_mean + z_sigm * constants[1]
              elast_ampl = current_values[4]
              elast_sigm = current_values[5]
              # leaving out small cross terms where one-photon peak influences
              # derivatives with respect to z_mean and z_sigm.
              # XXX not accounting for charge sharing or 2-photon peak
             */
            double udiff = free_x[ix] - current_values[0];
            double Afactor = udiff/(current_values[2] * current_values[2]);
            double Sfactor = Afactor * udiff/current_values[2];
            jacobian_one_row[0] = ( residual_terms_0[ix] * Afactor );
            jacobian_one_row[1] = ( residual_terms_0[ix] / current_values[1] );
            jacobian_one_row[2] = ( residual_terms_0[ix] * Sfactor );
            jacobian_one_row[3] = ( residual_terms_1[ix] / current_values[3] );
            jacobian_one_row[4] = ( residual_terms_2[ix] / current_values[4] );
            double udiff3 = free_x[ix] - (current_values[0] + current_values[2] * constants[1]);
            double Afactor3 = udiff3/(current_values[5] * current_values[5]);
            double Sfactor3 = Afactor3 * udiff3/current_values[5];
            jacobian_one_row[5] = ( residual_terms_2[ix] * Sfactor3 );

            add_equation(residuals[ix], jacobian_one_row.const_ref(), 1.);
          }
        }

      }

      void
      fvec_callable(scitbx::af::shared<double> &current_values) {
        /*
          # model minus obs
          # sqrt2pi_inv = 1./math.sqrt(2.*math.pi)
          # the gaussian function is sqrt2pi_inv * exp( - (x-mean)**2 / (2.*(sigma**2)))/sigma
          # take off the coefficient sqrt2pi_inv / sigma, use ampl
        */


        for (int ix = 0; ix < free_x.size(); ++ix) {
          double model=0;


          //zero-photon Gaussian
          double diff = free_x[ix] - current_values[0];
          residual_terms_0[ix] = (
            current_values[1] * std::exp( - (diff*diff) / (2. * current_values[2] * current_values[2])));
          model += residual_terms_0[ix];

          //inelastic-photon Gaussian
          {
          double pmean = current_values[0] + current_values[2] * constants[0];
          diff = free_x[ix] - pmean;
          double psigma = (constants[0]/constants[1])*(current_values[5] - current_values[2])+current_values[2];
          residual_terms_1[ix] = (
            current_values[3] * std::exp( - (diff*diff) / (2. * psigma * psigma)));
          model += residual_terms_1[ix];
          }

          //elastic-photon Gaussian
          {
          double pmean = current_values[0] + current_values[2] * constants[1];
          diff = free_x[ix] - pmean;
          double psigma = current_values[5];
          residual_terms_2[ix] = (
            current_values[4] * std::exp( - (diff*diff) / (2. * psigma * psigma)));
          model += residual_terms_2[ix];
          }

          //terms = [
          //ampl * flex.exp(-flex.pow2(free_x - mean) / (2.*sigm*sigm))
          //for mean,ampl,sigm in [(z_mean,z_ampl,z_sigm),(o_mean,o_ampl,o_sigm)]]
          //#print "residual", math.sqrt(flex.sum((model-free_y)*(model-free_y)))
          residuals[ix] = ( model - free_y[ix] );
        }
      }

    private:
      scitbx::af::shared<double> free_x, free_y, jacobian_one_row, residuals;
      scitbx::af::shared<double> constants;
      double gain_to_sigma, sigmafac, sigmafac_sq;
      std::vector<double> residual_terms_0, residual_terms_1, residual_terms_2;

  };

}} //xes::example

namespace boost_python { namespace {

  void
  xes_ext_init_module() {
    using namespace boost::python;

    typedef xes::example::gaussian_fit_inheriting_from_non_linear_ls wt;
    typedef xes::example::gaussian_3fit_inheriting_from_non_linear_ls wt3;

    class_<xes::example::gaussian_fit_inheriting_from_non_linear_ls,
           bases<scitbx::lstbx::normal_equations::non_linear_ls<double> > >(
      "gaussian_fit_inheriting_from_non_linear_ls", no_init)
      .def(init<int>(arg("n_parameters")))
      .def("access_cpp_build_up_directly",&wt::access_cpp_build_up_directly,
        (arg("objective_only"),arg("current_values")))
      .def("set_cpp_data",&wt::set_cpp_data,
        (arg("free_x"),arg("free_y"),arg("gain_to_sigma"),arg("sigmafac")))
    ;
    class_<xes::example::gaussian_3fit_inheriting_from_non_linear_ls,
           bases<scitbx::lstbx::normal_equations::non_linear_ls<double> > >(
      "gaussian_3fit_inheriting_from_non_linear_ls", no_init)
      .def(init<int>(arg("n_parameters")))
      .def("access_cpp_build_up_directly",&wt3::access_cpp_build_up_directly,
        (arg("objective_only"),arg("current_values")))
      .def("set_cpp_data",&wt3::set_cpp_data,
        (arg("constants"),arg("free_x"),arg("free_y")))
    ;
  }

}
}} // namespace xfel::boost_python::<anonymous>

BOOST_PYTHON_MODULE(xes_ext)
{
  xfel::boost_python::xes_ext_init_module();

}
