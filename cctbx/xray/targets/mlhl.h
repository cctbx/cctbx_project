#ifndef CCTBX_XRAY_TARGETS_MLHL_H
#define CCTBX_XRAY_TARGETS_MLHL_H

#include <cctbx/hendrickson_lattman.h>

namespace cctbx { namespace xray { namespace targets { namespace mlhl {

  //namespace detail {

    inline
    double
    similar(double y)
    {
      double epsilon = 1.0e-15;
      double lowerlim = 20.0;
      int maxterms = 150;
      double x = std::abs(y);
      double tot0 = 1;
      double subtot0 = 1;
      double tot1 = 1;
      double subtot1 = 1;
      if (x < lowerlim) {
        int n = 1;
        while ((n <= maxterms) && (subtot0 >= epsilon)) {
          double dpn = static_cast<double>(n);
          subtot0 = x*x*subtot0/(4*dpn*dpn);
          subtot1 = x*x*subtot1/(4*dpn*(dpn+1));
          tot0 += subtot0;
          tot1 += subtot1;
          n++;
        }
        tot0 = tot1*x/(2*tot0);
      }
      else {
        int n = 1;
        while ((n <= maxterms) && (std::abs(subtot0) >= epsilon)) {
          double dpn = static_cast<double>(2*n);
          subtot0 = (dpn - 1)*(dpn - 1) / (4*x*dpn)*subtot0;
          tot0 += subtot0;
          tot1 += (2/(1 - dpn) - 1) * subtot0;
          n++;
        }
        tot0 = tot1/tot0;
      }
      if (y < 0) tot0 = -tot0;
      return tot0;
    }

  //} // namespace detail

    inline
    double
    target_one_h(
      double fo,
      double fc,
      double pc,
      double alpha,
      double beta,
      double k,
      double epsilon,
      bool cf,
      cctbx::hendrickson_lattman<double> const& abcd,
      const af::tiny_plain<double, 4>* cos_sin_table,
      int n_steps,
      double integration_step_size,
      double* workspace)
    {
      CCTBX_ASSERT(fo >= 0);
      CCTBX_ASSERT(fc >= 0);
      const double small = 1.e-9;
      CCTBX_ASSERT(std::abs(k) > small);
      if(alpha <= 0 || beta <= 0) return 0;
      double target = 0;
      alpha *= k;
      beta *= k*k;
      double hl_a = abcd.a();
      double hl_b = abcd.b();
      double hl_c = abcd.c();
      double hl_d = abcd.d();
      // acentric reflection
      if(!cf) {
        double arg = 2*alpha*fo*fc/(beta*epsilon);
        double a_prime = arg * std::cos(pc) + hl_a;
        double b_prime = arg * std::sin(pc) + hl_b;
        // calculate target analytically
        if(std::abs(hl_c) < small && std::abs(hl_d) < small) {
          double val = std::sqrt(a_prime*a_prime + b_prime*b_prime);
          if(std::abs(hl_a) < small && std::abs(hl_b) < small) {
            val = arg;
          }
          target = scitbx::math::bessel::ln_of_i0(val);
        }
        // calculate target numerically
        else {
          double maxv = 0;
          for(int i=0;i<n_steps;i++) {
            const double* tab = cos_sin_table[i].begin();
            double term = a_prime * tab[0]
                        + b_prime * tab[1]
                        + hl_c    * tab[2]
                        + hl_d    * tab[3];
            if (maxv < term) maxv = term;
            workspace[i] = term;
          }
          for(int i=0;i<n_steps;i++) {
            target += std::exp(-maxv+workspace[i]);
          }
          target *= integration_step_size;
          target = std::log(target) + maxv;
        }
        target = (fo*fo+alpha*alpha*fc*fc)/(beta*epsilon) - target;
      }
      // centric reflection
      else {
        double var = beta*epsilon;
        double arg = fo*alpha*fc/var;
        arg += hl_a*std::cos(pc) + hl_b*std::sin(pc);
        double mabsarg = -std::abs(arg);
        target = (fo*fo + alpha*alpha*fc*fc)/(2*var) + mabsarg
               - std::log((1 + std::exp(2*mabsarg))/2);
      }
      return target;
    }

    inline
    std::complex<double>
    d_target_d_f_calc_one_h(
      double fo,
      double fc,
      double pc,
      double ac,
      double bc,
      double alpha,
      double beta,
      double epsilon,
      bool cf,
      cctbx::hendrickson_lattman<double> const& abcd,
      const af::tiny_plain<double, 4>* cos_sin_table,
      int n_steps,
      double integration_step_size,
      double* workspace)
    {
      const double small = 1.e-9;
      if (fc < small || alpha <= 0 || beta <= 0) {
        return std::complex<double>(0,0);
      }
      double derfc = 0;
      double derpc = 0;
      double cos_pc = std::cos(pc);
      double sin_pc = std::sin(pc);
      double hl_a = abcd.a();
      double hl_b = abcd.b();
      // acentric reflection
      if (!cf) {
        double hl_c = abcd.c();
        double hl_d = abcd.d();
        double arg = 2*alpha*fo/(beta*epsilon);
        double a_prime = arg * fc * cos_pc + hl_a;
        double b_prime = arg * fc * sin_pc + hl_b;
        if (std::abs(hl_c) < small && std::abs(hl_d) < small) {
          double val = std::sqrt(a_prime*a_prime + b_prime*b_prime);
          if(val < small) {
            derfc = 0;
            derpc = 0;
          }
          else {
            double sim = similar(val);
            derfc = sim*arg*(arg*fc + hl_a*cos_pc + hl_b*sin_pc)/val;
            derpc = sim*arg*fc*(hl_a*sin_pc - hl_b*cos_pc)/val;
          }
        }
        // calculate gradients numerically
        else {
          double maxv = 0;
          for(int i=0;i<n_steps;i++) {
            const double* tab = cos_sin_table[i].begin();
            double term = a_prime * tab[0]
                        + b_prime * tab[1]
                        + hl_c    * tab[2]
                        + hl_d    * tab[3];
            if (maxv < term) maxv = term;
            workspace[i] = term;
          }
          double target = 0;
          for(int i=0;i<n_steps;i++) {
            target += std::exp(-maxv+workspace[i]);
          }
          target *= integration_step_size;
          target = -std::log(target) - maxv;
          double deranot = 0;
          double derbnot = 0;
          for(int i=0;i<n_steps;i++) {
            double exp_t_w = std::exp(target+workspace[i]);
            const double* tab = cos_sin_table[i].begin();
            deranot += tab[0] * exp_t_w;
            derbnot += tab[1] * exp_t_w;
          }
          deranot *= integration_step_size;
          derbnot *= integration_step_size;
          derfc = arg*(deranot*cos_pc + derbnot*sin_pc);
          derpc = arg*(deranot*sin_pc - derbnot*cos_pc)*fc;
        }
        derfc = 2*alpha*alpha*fc/(beta*epsilon) - derfc;
      }
      // centric reflection
      else {
        double var = beta*epsilon;
        double arg = hl_a*cos_pc + hl_b*sin_pc + fo*alpha*fc/var;
        double mtwo_arg = -2*arg;
        if(mtwo_arg > 30.) mtwo_arg = 30.0;
        double exp_2_arg = std::exp(mtwo_arg);
        double tmp_tanh = (1-exp_2_arg) / (1+exp_2_arg);
        derfc = alpha*alpha*fc/var - tmp_tanh*fo*alpha/var;
        derpc = 2*tmp_tanh*(hl_a*sin_pc - hl_b*cos_pc);
      }
      return std::complex<double>(
         (derfc*ac - derpc*bc/fc)/fc,
        -(derfc*bc + derpc*ac/fc)/fc);
    }

  //! Maximum-likelihood target function and gradients.
  /*! Incorporates experimental phase information as HL coefficients ABCD.
      As described by Pannu et al, Acta Cryst. (1998). D54, 1285-1294.
      All the equations are reformulated in terms of alpha/beta.
      Pavel Afonine // 14-DEC-2004
   */
  class target_and_gradients : public common_results
  {
    public:
      target_and_gradients(
        af::const_ref<double> const& f_obs,
        af::const_ref<bool> const& r_free_flags,
        af::const_ref<cctbx::hendrickson_lattman<double> > const&
          experimental_phases,
        af::const_ref<std::complex<double> > const& f_calc,
        af::const_ref<double> const& alpha,
        af::const_ref<double> const& beta,
        af::const_ref<double> const& epsilons,
        af::const_ref<bool> const& centric_flags,
        double integration_step_size,
        bool compute_gradients)
      :
        common_results(f_obs.size())
      {
        CCTBX_ASSERT(r_free_flags.size() == 0
                  || r_free_flags.size() == f_obs.size());
        CCTBX_ASSERT(experimental_phases.size() == f_obs.size());
        CCTBX_ASSERT(f_calc.size() == f_obs.size());
        CCTBX_ASSERT(alpha.size() == f_obs.size());
        CCTBX_ASSERT(beta.size() == f_obs.size());
        CCTBX_ASSERT(epsilons.size() == f_obs.size());
        CCTBX_ASSERT(centric_flags.size() == f_obs.size());
        CCTBX_ASSERT(integration_step_size > 0);
        if (f_obs.size() == 0) return;
        detail::r_free_flags_stats rffs(f_obs.size(), r_free_flags.begin());
        CCTBX_ASSERT(rffs.n_work != 0);
        double one_over_n_work = 1./ rffs.n_work;
        if (compute_gradients) {
          gradients_work_.reserve(rffs.n_work);
        }
        hendrickson_lattman<double>::phase_integration_cos_sin_table
          cos_sin_table(scitbx::math::iround(360 / integration_step_size));
        CCTBX_ASSERT(cos_sin_table.n_steps > 0);
        boost::scoped_array<double> workspace(
          new double[cos_sin_table.n_steps]);
        double target_work = 0;
        double target_test = 0;
        for(std::size_t i=0;i<f_obs.size();i++) {
          double fo = f_obs[i];
          double fc = std::abs(f_calc[i]);
          double pc = std::arg(f_calc[i]);
          double ac = std::real(f_calc[i]);
          double bc = std::imag(f_calc[i]);
          double t = target_one_h(
            fo,
            fc,
            pc,
            alpha[i],
            beta[i],
            1.0,
            epsilons[i],
            centric_flags[i],
            experimental_phases[i],
            cos_sin_table.data.get(),
            cos_sin_table.n_steps,
            integration_step_size,
            workspace.get());
          target_per_reflection_[i] = t;
          if (rffs.is_work_refl(i)) {
            target_work += t;
            if (compute_gradients) {
              gradients_work_.push_back(std::conj(
                d_target_d_f_calc_one_h(
                  fo,
                  fc,
                  pc,
                  ac,
                  bc,
                  alpha[i],
                  beta[i],
                  epsilons[i],
                  centric_flags[i],
                  experimental_phases[i],
                  cos_sin_table.data.get(),
                  cos_sin_table.n_steps,
                  integration_step_size,
                  workspace.get())) * one_over_n_work);
            }
          }
          else {
            target_test += t;
          }
        }
        target_work_ = target_work * one_over_n_work;
        if (rffs.n_test != 0) {
          target_test_ = boost::optional<double>(target_test / rffs.n_test);
        }
      }
  };


}}}} // namespace cctbx::xray::targets

#endif // GUARD
