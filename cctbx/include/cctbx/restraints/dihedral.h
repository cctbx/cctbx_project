#ifndef CCTBX_RESTRAINTS_DIHEDRAL_H
#define CCTBX_RESTRAINTS_DIHEDRAL_H

#include <cctbx/restraints/utils.h>

namespace cctbx { namespace restraints {

  struct dihedral_proxy
  {
    dihedral_proxy() {}

    dihedral_proxy(
      af::tiny<std::size_t, 4> const& i_seqs_,
      double angle_ideal_,
      double weight_,
      int periodicity_=0)
    :
      i_seqs(i_seqs_),
      angle_ideal(angle_ideal_),
      weight(weight_),
      periodicity(periodicity_)
    {}

    af::tiny<std::size_t, 4> i_seqs;
    double angle_ideal;
    double weight;
    int periodicity;
  };

  /*! CCP4 mon_lib torsion angles and XPLOR/CNS dihedrals
      are identical to MODELLER dihedrals:

      http://salilab.org/modeller/modeller.html

      van Schaik, R. C., Berendsen, H. J., & Torda, A. E. (1993).
      J.Mol.Biol. 234, 751-762.
   */
  class dihedral
  {
    public:
      dihedral() {}

      dihedral(
        af::tiny<scitbx::vec3<double>, 4> const& sites_,
        double angle_ideal_,
        double weight_,
        int periodicity_=0)
      :
        sites(sites_),
        angle_ideal(angle_ideal_),
        weight(weight_),
        periodicity(periodicity_)
      {
        init_angle_model();
      }

      dihedral(
        af::const_ref<scitbx::vec3<double> > const& sites_cart,
        dihedral_proxy const& proxy)
      :
        angle_ideal(proxy.angle_ideal),
        weight(proxy.weight),
        periodicity(proxy.periodicity)
      {
        for(int i=0;i<4;i++) {
          std::size_t i_seq = proxy.i_seqs[i];
          CCTBX_ASSERT(i_seq < sites_cart.size());
          sites[i] = sites_cart[i_seq];
        }
        init_angle_model();
      }

      double
      residual() const { return weight * scitbx::fn::pow2(delta); }

      af::tiny<scitbx::vec3<double>, 4>
      gradients(double epsilon=1.e-100) const
      {
        af::tiny<scitbx::vec3<double>, 4> result;
        double d_21_norm = d_21.length_sq();
        if (   !have_angle_model
            || d_21_norm < epsilon
            || n_0121_norm < epsilon
            || n_2123_norm < epsilon) {
          result.fill(scitbx::vec3<double>(0,0,0));
        }
        else {
          double grad_factor = 2 * weight * delta / scitbx::constants::pi_180
                             * d_21.length();
          result[0] = -grad_factor/n_0121_norm * n_0121;
          result[3] = grad_factor/n_2123_norm * n_2123;
          double d_01_dot_d_21 = d_01 * d_21;
          double d_21_dot_d_23 = d_21 * d_23;
          result[1] = (d_01_dot_d_21/d_21_norm-1) * result[0]
                    - d_21_dot_d_23/d_21_norm * result[3];
          result[2] = (d_21_dot_d_23/d_21_norm-1) * result[3]
                    - d_01_dot_d_21/d_21_norm * result[0];
        }
        return result;
      }

      // Not available in Python.
      void
      add_gradients(
        af::ref<scitbx::vec3<double> > const& gradient_array,
        af::tiny<std::size_t, 4> const& i_seqs) const
      {
        af::tiny<scitbx::vec3<double>, 4> grads = gradients();
        for(int i=0;i<4;i++) {
          gradient_array[i_seqs[i]] += grads[i];
        }
      }

      af::tiny<scitbx::vec3<double>, 4> sites;
      double angle_ideal;
      double weight;
      int periodicity;
      bool have_angle_model;
    protected:
      scitbx::vec3<double> d_01;
      scitbx::vec3<double> d_21;
      scitbx::vec3<double> d_23;
      scitbx::vec3<double> n_0121;
      double n_0121_norm;
      scitbx::vec3<double> n_2123;
      double n_2123_norm;
    public:
      double angle_model;
      double delta;

    protected:
      void
      init_angle_model()
      {
        have_angle_model = false;
        angle_model = angle_ideal;
        delta = 0;
        d_01 = sites[0] - sites[1];
        d_21 = sites[2] - sites[1];
        d_23 = sites[2] - sites[3];
        n_0121 = d_01.cross(d_21);
        n_0121_norm = n_0121.length_sq();
        if (n_0121_norm == 0) return;
        n_2123 = d_21.cross(d_23);
        n_2123_norm = n_2123.length_sq();
        if (n_2123_norm == 0) return;
        double cos_angle_model = std::max(-1.,std::min(1.,
          n_0121 * n_2123 / std::sqrt(n_0121_norm * n_2123_norm)));
        angle_model = std::acos(cos_angle_model)
                    / scitbx::constants::pi_180;
        if (d_21 * (n_0121.cross(n_2123)) < 0) {
          angle_model *= -1;
        }
        have_angle_model = true;
        delta = angle_delta_deg(angle_model, angle_ideal, periodicity);
      }
  };

  inline
  af::shared<double>
  dihedral_deltas(
    af::const_ref<scitbx::vec3<double> > const& sites_cart,
    af::const_ref<dihedral_proxy> const& proxies)
  {
    return detail::generic_deltas<dihedral_proxy, dihedral>::get(
      sites_cart, proxies);
  }

  inline
  af::shared<double>
  dihedral_residuals(
    af::const_ref<scitbx::vec3<double> > const& sites_cart,
    af::const_ref<dihedral_proxy> const& proxies)
  {
    return detail::generic_residuals<dihedral_proxy, dihedral>::get(
      sites_cart, proxies);
  }

  inline
  double
  dihedral_residual_sum(
    af::const_ref<scitbx::vec3<double> > const& sites_cart,
    af::const_ref<dihedral_proxy> const& proxies,
    af::ref<scitbx::vec3<double> > const& gradient_array)
  {
    return detail::generic_residual_sum<dihedral_proxy, dihedral>::get(
      sites_cart, proxies, gradient_array);
  }

}} // namespace cctbx::restraints

#endif // CCTBX_RESTRAINTS_DIHEDRAL_H
