#ifndef CCTBX_RESTRAINTS_CHIRALITY_H
#define CCTBX_RESTRAINTS_CHIRALITY_H

#include <cctbx/restraints/utils.h>

namespace cctbx { namespace restraints {

  struct chirality_proxy
  {
    chirality_proxy() {}

    chirality_proxy(
      af::tiny<unsigned, 4> const& i_seqs_,
      double volume_ideal_,
      bool both_signs_,
      double weight_)
    :
      i_seqs(i_seqs_),
      volume_ideal(volume_ideal_),
      both_signs(both_signs_),
      weight(weight_)
    {}

    af::tiny<unsigned, 4> i_seqs;
    double volume_ideal;
    bool both_signs;
    double weight;
  };

  class chirality
  {
    public:
      chirality() {}

      chirality(
        af::tiny<scitbx::vec3<double>, 4> const& sites_,
        double volume_ideal_,
        bool both_signs_,
        double weight_)
      :
        sites(sites_),
        volume_ideal(volume_ideal_),
        both_signs(both_signs_),
        weight(weight_)
      {
        init_volume_model();
      }

      chirality(
        af::const_ref<scitbx::vec3<double> > const& sites_cart,
        chirality_proxy const& proxy)
      :
        volume_ideal(proxy.volume_ideal),
        both_signs(proxy.both_signs),
        weight(proxy.weight)
      {
        for(int i=0;i<4;i++) {
          std::size_t i_seq = proxy.i_seqs[i];
          CCTBX_ASSERT(i_seq < sites_cart.size());
          sites[i] = sites_cart[i_seq];
        }
        init_volume_model();
      }

      double
      residual() const { return weight * scitbx::fn::pow2(delta); }

      af::tiny<scitbx::vec3<double>, 4>
      gradients() const
      {
        af::tiny<scitbx::vec3<double>, 4> result;
        double f = delta_sign * 2 * weight * delta;
        result[1] = f * d_02_cross_d_03;
        result[2] = f * d_03.cross(d_01);
        result[3] = f * d_01.cross(d_02);
        result[0] = -result[1]-result[2]-result[3];
        return result;
      }

      // Not available in Python.
      void
      add_gradients(
        af::ref<scitbx::vec3<double> > const& gradient_array,
        af::tiny<unsigned, 4> const& i_seqs) const
      {
        af::tiny<scitbx::vec3<double>, 4> grads = gradients();
        for(int i=0;i<4;i++) {
          gradient_array[i_seqs[i]] += grads[i];
        }
      }

#if defined(__MACH__) && defined(__APPLE_CC__) && __APPLE_CC__ <= 1640
      bool dummy_;
#endif
      af::tiny<scitbx::vec3<double>, 4> sites;
      double volume_ideal;
      bool both_signs;
      double weight;
    protected:
      scitbx::vec3<double> d_01;
      scitbx::vec3<double> d_02;
      scitbx::vec3<double> d_03;
      scitbx::vec3<double> d_02_cross_d_03;
    public:
      double volume_model;
      double delta_sign;
      double delta;

    protected:
      void
      init_volume_model()
      {
        d_01 = sites[1] - sites[0];
        d_02 = sites[2] - sites[0];
        d_03 = sites[3] - sites[0];
        d_02_cross_d_03 = d_02.cross(d_03);
        volume_model = d_01 * d_02_cross_d_03;
        delta_sign = -1;
        if (both_signs && volume_model < 0) delta_sign = 1;
        delta = volume_ideal + delta_sign * volume_model;
      }
  };

  inline
  af::shared<double>
  chirality_deltas(
    af::const_ref<scitbx::vec3<double> > const& sites_cart,
    af::const_ref<chirality_proxy> const& proxies)
  {
    return detail::generic_deltas<chirality_proxy, chirality>::get(
      sites_cart, proxies);
  }

  inline
  af::shared<double>
  chirality_residuals(
    af::const_ref<scitbx::vec3<double> > const& sites_cart,
    af::const_ref<chirality_proxy> const& proxies)
  {
    return detail::generic_residuals<chirality_proxy, chirality>::get(
      sites_cart, proxies);
  }

  inline
  double
  chirality_residual_sum(
    af::const_ref<scitbx::vec3<double> > const& sites_cart,
    af::const_ref<chirality_proxy> const& proxies,
    af::ref<scitbx::vec3<double> > const& gradient_array)
  {
    return detail::generic_residual_sum<chirality_proxy, chirality>::get(
      sites_cart, proxies, gradient_array);
  }

}} // namespace cctbx::restraints

#endif // CCTBX_RESTRAINTS_CHIRALITY_H
