#ifndef CCTBX_GEOMETRY_RESTRAINTS_DIHEDRAL_H
#define CCTBX_GEOMETRY_RESTRAINTS_DIHEDRAL_H

#include <cctbx/geometry_restraints/utils.h>

namespace cctbx { namespace geometry_restraints {

  //! Grouping of indices into array of sites (i_seqs) and parameters.
  struct dihedral_proxy
  {
    //! Default constructor. Some data members are not initialized!
    dihedral_proxy() {}

    //! Constructor.
    dihedral_proxy(
      af::tiny<unsigned, 4> const& i_seqs_,
      double angle_ideal_,
      double weight_,
      int periodicity_=0)
    :
      i_seqs(i_seqs_),
      angle_ideal(angle_ideal_),
      weight(weight_),
      periodicity(periodicity_)
    {}

    //! Sorts i_seqs such that i_seq[0] < i_seq[3] and i_seq[1] < i_seq[2].
    dihedral_proxy
    sort_i_seqs() const
    {
      dihedral_proxy result(*this);
      if (result.i_seqs[0] > result.i_seqs[3]) {
        std::swap(result.i_seqs[0], result.i_seqs[3]);
        result.angle_ideal *= -1;
      }
      if (result.i_seqs[1] > result.i_seqs[2]) {
        std::swap(result.i_seqs[1], result.i_seqs[2]);
        result.angle_ideal *= -1;
      }
      return result;
    }

    //! Indices into array of sites.
    af::tiny<unsigned, 4> i_seqs;
    //! Parameter.
    double angle_ideal;
    //! Parameter.
    double weight;
    //! Parameter.
    int periodicity;
  };

  //! Residual and gradient calculations for dihedral %angle restraint.
  /*! angle_model is defined as the %angle between the plane through
      the three points sites[0],sites[1],sites[2] and the plane through
      the three points sites[1],sites[2],sites[3]. This definition is
      compatible with CCP4 Monomer library torsion %angles, XPLOR/CNS
      dihedrals and MODELLER dihedrals.

      See also:
        http://salilab.org/modeller/manual/manual.html,
        "Features and their derivatives",

        van Schaik, R. C., Berendsen, H. J., & Torda, A. E. (1993).
        J.Mol.Biol. 234, 751-762.
   */
  class dihedral
  {
    public:
      //! Default constructor. Some data members are not initialized!
      dihedral() {}

      //! Constructor.
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

      /*! \brief Coordinates are copied from sites_cart according to
          proxy.i_seqs, parameters are copied from proxy.
       */
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

      //! weight * delta**2.
      /*! See also: Hendrickson, W.A. (1985). Meth. Enzym. 115, 252-270.
       */
      double
      residual() const { return weight * scitbx::fn::pow2(delta); }

      //! Gradients with respect to the four sites.
      /*! The formula for the gradients is singular if certain vectors
          are collinear. However, the gradients converge to zero near
          these singularities. To avoid numerical problems, the
          gradients are set to zero exactly if the norms of certain
          vectors are smaller than epsilon.

          See also:
            http://salilab.org/modeller/manual/manual.html,
            "Features and their derivatives"
       */
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

      //! Support for dihedral_residual_sum.
      /*! Not available in Python.
       */
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

      //! Cartesian coordinates of the sites defining the dihedral %angle.
      af::tiny<scitbx::vec3<double>, 4> sites;
      //! Parameter (usually as passed to the constructor).
      double angle_ideal;
      //! Parameter (usually as passed to the constructor).
      double weight;
      //! Parameter (usually as passed to the constructor).
      int periodicity;
      //! false in singular situations.
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
      //! Value of the dihedral %angle formed by the sites.
      double angle_model;
      /*! \brief Smallest difference between angle_model and angle_ideal
          taking the periodicity into account.
       */
      /*! See also: angle_delta_deg
       */
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

  //! Fast computation of dihedral::delta given an array of dihedral proxies.
  inline
  af::shared<double>
  dihedral_deltas(
    af::const_ref<scitbx::vec3<double> > const& sites_cart,
    af::const_ref<dihedral_proxy> const& proxies)
  {
    return detail::generic_deltas<dihedral_proxy, dihedral>::get(
      sites_cart, proxies);
  }

  /*! Fast computation of dihedral::residual() given an array of
      dihedral proxies.
   */
  inline
  af::shared<double>
  dihedral_residuals(
    af::const_ref<scitbx::vec3<double> > const& sites_cart,
    af::const_ref<dihedral_proxy> const& proxies)
  {
    return detail::generic_residuals<dihedral_proxy, dihedral>::get(
      sites_cart, proxies);
  }

  /*! Fast computation of sum of dihedral::residual() and gradients
      given an array of dihedral proxies.
   */
  /*! The dihedral::gradients() are added to the gradient_array if
      gradient_array.size() == sites_cart.size().
      gradient_array must be initialized before this function
      is called.
      No gradient calculations are performed if gradient_array.size() == 0.
   */
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

}} // namespace cctbx::geometry_restraints

#endif // CCTBX_GEOMETRY_RESTRAINTS_DIHEDRAL_H
