#ifndef CCTBX_RESTRAINTS_ANGLE_H
#define CCTBX_RESTRAINTS_ANGLE_H

#include <cctbx/restraints/utils.h>

namespace cctbx { namespace restraints {

  struct angle_proxy
  {
    angle_proxy() {}

    angle_proxy(
      af::tiny<std::size_t, 3> const& i_seqs_,
      double angle_ideal_,
      double weight_)
    :
      i_seqs(i_seqs_),
      angle_ideal(angle_ideal_),
      weight(weight_)
    {}

    af::tiny<std::size_t, 3> i_seqs;
    double angle_ideal;
    double weight;
  };

  class angle
  {
    public:
      angle() {}

      angle(
        af::tiny<scitbx::vec3<double>, 3> const& sites_,
        double angle_ideal_,
        double weight_)
      :
        sites(sites_),
        angle_ideal(angle_ideal_),
        weight(weight_)
      {
        init_angle_model();
      }

      angle(
        af::const_ref<scitbx::vec3<double> > const& sites_cart,
        angle_proxy const& proxy)
      :
        angle_ideal(proxy.angle_ideal),
        weight(proxy.weight)
      {
        for(int i=0;i<3;i++) {
          std::size_t i_seq = proxy.i_seqs[i];
          CCTBX_ASSERT(i_seq < sites_cart.size());
          sites[i] = sites_cart[i_seq];
        }
        init_angle_model();
      }

      double
      residual() const { return weight * scitbx::fn::pow2(delta); }

      af::tiny<scitbx::vec3<double>, 3>
      gradients(double epsilon=1.e-100) const
      {
        af::tiny<scitbx::vec3<double>, 3> result;
        if (!have_angle_model) {
          result.fill(scitbx::vec3<double>(0,0,0));
        }
        else {
          double
          one_over_grad_acos = std::sqrt(1-scitbx::fn::pow2(cos_angle_model));
          if (one_over_grad_acos < epsilon) {
            result.fill(scitbx::vec3<double>(0,0,0));
          }
          else {
            double grad_factor = -2 * weight * delta/scitbx::constants::pi_180
                               / one_over_grad_acos;
            result[0] = grad_factor / d_01_abs
                      * (d_01_unit * cos_angle_model - d_21_unit);
            result[2] = grad_factor / d_21_abs
                      * (d_21_unit * cos_angle_model - d_01_unit);
            result[1] = -(result[0] + result[2]);
          }
        }
        return result;
      }

      // Not available in Python.
      void
      add_gradients(
        af::ref<scitbx::vec3<double> > const& gradient_array,
        af::tiny<std::size_t, 3> const& i_seqs) const
      {
        af::tiny<scitbx::vec3<double>, 3> grads = gradients();
        for(int i=0;i<3;i++) {
          gradient_array[i_seqs[i]] += grads[i];
        }
      }

#if defined(__MACH__) && defined(__APPLE_CC__) && __APPLE_CC__ <= 1640
      bool dummy_;
#endif
      af::tiny<scitbx::vec3<double>, 3> sites;
      double angle_ideal;
      double weight;
      bool have_angle_model;
    protected:
      double d_01_abs;
      double d_21_abs;
      scitbx::vec3<double> d_01_unit;
      scitbx::vec3<double> d_21_unit;
      double cos_angle_model;
    public:
      double angle_model;
      double delta;

    protected:
      void
      init_angle_model()
      {
        have_angle_model = false;
        d_01_abs = 0;
        d_21_abs = 0;
        d_01_unit.fill(0);
        d_21_unit.fill(0);
        cos_angle_model = -9;
        angle_model = angle_ideal;
        delta = 0;
        scitbx::vec3<double> d_01 = sites[0] - sites[1];
        d_01_abs = d_01.length();
        if (d_01_abs > 0) {
          scitbx::vec3<double> d_21 = sites[2] - sites[1];
          d_21_abs = d_21.length();
          if (d_21_abs > 0) {
            d_01_unit = d_01 / d_01_abs;
            d_21_unit = d_21 / d_21_abs;
            cos_angle_model = std::max(-1.,std::min(1.,d_01_unit*d_21_unit));
            angle_model = std::acos(cos_angle_model)
                        / scitbx::constants::pi_180;
            have_angle_model = true;
            delta = angle_delta_deg(angle_model, angle_ideal);
          }
        }
      }
  };

  inline
  af::shared<double>
  angle_deltas(
    af::const_ref<scitbx::vec3<double> > const& sites_cart,
    af::const_ref<angle_proxy> const& proxies)
  {
    return detail::generic_deltas<angle_proxy, angle>::get(
      sites_cart, proxies);
  }

  inline
  af::shared<double>
  angle_residuals(
    af::const_ref<scitbx::vec3<double> > const& sites_cart,
    af::const_ref<angle_proxy> const& proxies)
  {
    return detail::generic_residuals<angle_proxy, angle>::get(
      sites_cart, proxies);
  }

  inline
  double
  angle_residual_sum(
    af::const_ref<scitbx::vec3<double> > const& sites_cart,
    af::const_ref<angle_proxy> const& proxies,
    af::ref<scitbx::vec3<double> > const& gradient_array)
  {
    return detail::generic_residual_sum<angle_proxy, angle>::get(
      sites_cart, proxies, gradient_array);
  }

}} // namespace cctbx::restraints

#endif // CCTBX_RESTRAINTS_ANGLE_H
