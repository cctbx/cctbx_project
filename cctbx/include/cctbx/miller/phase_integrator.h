#ifndef CCTBX_MILLER_PHASE_INTEGRATOR_H
#define CCTBX_MILLER_PHASE_INTEGRATOR_H

#include <cctbx/sgtbx/space_group.h>

namespace cctbx { namespace miller {

  /*! Conversion of Hendrickson-Lattman coefficients to
      centroid phases and figures of merit.
   */
  class phase_integrator
  {
    public:
      //! Initialization of internal cos and sin lookup table.
      phase_integrator(unsigned n_steps=360/5)
      :
        n_steps_(n_steps)
      {
        CCTBX_ASSERT(n_steps > 0);
        angular_step_ = scitbx::constants::two_pi / n_steps_;
        cos_sin_table_.reserve(n_steps_);
        for(unsigned i_step=0;i_step<n_steps_;i_step++) {
          double angle = i_step * angular_step_;
          cos_sin_table_.push_back(af::tiny<double, 4>(
            std::cos(angle),
            std::sin(angle),
            std::cos(angle+angle),
            std::sin(angle+angle)));
        }
      }

      //! Number of integration steps as passed to the constructor.
      unsigned
      n_steps() const { return n_steps_; }

      //! Computation of phase integral.
      template <typename FloatType>
      std::complex<FloatType>
      operator()(
        sgtbx::phase_info const& phase_info,
        cctbx::hendrickson_lattman<FloatType> const& hendrickson_lattman) const
      {
        typedef FloatType f_t;
        if (phase_info.is_centric()) {
          f_t angle = phase_info.ht_angle();
          f_t cos_angle = std::cos(angle);
          f_t sin_angle = std::sin(angle);
          f_t exp_arg = hendrickson_lattman.a() * cos_angle
                      + hendrickson_lattman.b() * sin_angle;
          f_t arg_term;
          if (exp_arg < 0) {
            arg_term = std::log(std::exp( 2*exp_arg) + 1) - exp_arg;
          }
          else {
            arg_term = std::log(std::exp(-2*exp_arg) + 1) + exp_arg;
          }
          f_t fom = std::exp( exp_arg - arg_term)
                  - std::exp(-exp_arg - arg_term);
          return std::complex<f_t>(fom*cos_angle, fom*sin_angle);
        }
        std::vector<f_t> exp_args;
        exp_args.reserve(n_steps_);
        f_t max_exp_arg = 0;
        for(unsigned i_step=0;i_step<n_steps_;i_step++) {
          f_t exp_arg = 0;
          for(unsigned i=0;i<4;i++) {
            exp_arg += hendrickson_lattman[i] * cos_sin_table_[i_step][i];
          }
          exp_args.push_back(exp_arg);
          max_exp_arg = std::max(max_exp_arg, exp_arg);
        }
        f_t sum_exp = 0;
        for(unsigned i_step=0;i_step<n_steps_;i_step++) {
          sum_exp += std::exp(exp_args[i_step] - max_exp_arg);
        }
        f_t arg_term = std::log(sum_exp * angular_step_) + max_exp_arg;
        std::complex<f_t> result;
        for(unsigned i_step=0;i_step<n_steps_;i_step++) {
          f_t fom = std::exp(exp_args[i_step] - arg_term);
          result += std::complex<f_t>(
            fom * cos_sin_table_[i_step][0],
            fom * cos_sin_table_[i_step][1]);
        }
        result *= angular_step_;
        return result;
      }

      //! Computation of phase integrals.
      template <typename FloatType>
      af::shared<std::complex<FloatType> >
      operator()(
        sgtbx::space_group const& space_group,
        af::const_ref<miller::index<> > const& miller_indices,
        af::const_ref<hendrickson_lattman<FloatType> > const&
          hendrickson_lattman_coefficients) const
      {
        CCTBX_ASSERT(hendrickson_lattman_coefficients.size()
                  == miller_indices.size());
        af::shared<std::complex<FloatType> >
          result((af::reserve(miller_indices.size())));
        for(std::size_t i=0;i<miller_indices.size();i++) {
          result.push_back((*this)(
            space_group.phase_restriction(miller_indices[i]),
            hendrickson_lattman_coefficients[i]));
        }
        return result;
      }

    protected:
      unsigned n_steps_;
      double angular_step_;
      std::vector<af::tiny<double, 4> > cos_sin_table_;
  };

}} // namespace cctbx::miller

#endif // CCTBX_MILLER_PHASE_INTEGRATOR_H
