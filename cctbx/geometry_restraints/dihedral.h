#ifndef CCTBX_GEOMETRY_RESTRAINTS_DIHEDRAL_H
#define CCTBX_GEOMETRY_RESTRAINTS_DIHEDRAL_H

#include <cctbx/sgtbx/rt_mx.h>
#include <cctbx/geometry_restraints/utils.h>
#include <cctbx/restraints.h>
#include <scitbx/math/dihedral.h>
#include <boost/format.hpp>

namespace cctbx { namespace geometry_restraints {

  typedef optional_container<af::small<double, 6> > alt_angle_ideals_type;

  //! Grouping of indices into array of sites (i_seqs) and dihedral_params.
  struct dihedral_proxy
  {
    //! Support for shared_proxy_select.
    typedef af::tiny<unsigned, 4> i_seqs_type;

    //! Default constructor. Some data members are not initialized!
    dihedral_proxy() {}

    //! Constructor.
    dihedral_proxy(
      i_seqs_type const& i_seqs_,
      double angle_ideal_,
      double weight_,
      int periodicity_=0,
      alt_angle_ideals_type const& alt_angle_ideals_=alt_angle_ideals_type(),
      double limit_=-1.0,
      bool top_out_=false)
    :
      i_seqs(i_seqs_),
      angle_ideal(angle_ideal_),
      weight(weight_),
      periodicity(periodicity_),
      alt_angle_ideals(alt_angle_ideals_),
      limit(limit_),
      top_out(top_out_)
    {}

    //! Constructor.
    dihedral_proxy(
      i_seqs_type const& i_seqs_,
      optional_container<af::shared<sgtbx::rt_mx> > const& sym_ops_,
      double angle_ideal_,
      double weight_,
      int periodicity_=0,
      alt_angle_ideals_type const& alt_angle_ideals_=alt_angle_ideals_type(),
      double limit_=-1.0,
      bool top_out_=false)
    :
      i_seqs(i_seqs_),
      sym_ops(sym_ops_),
      angle_ideal(angle_ideal_),
      weight(weight_),
      periodicity(periodicity_),
      alt_angle_ideals(alt_angle_ideals_),
      limit(limit_),
      top_out(top_out_)
    {
      if ( sym_ops.get() != 0 ) {
        CCTBX_ASSERT(sym_ops.get()->size() == i_seqs.size());
      }
    }

    //! Support for proxy_select (and similar operations).
    dihedral_proxy(
      i_seqs_type const& i_seqs_,
      dihedral_proxy const& proxy)
    :
      i_seqs(i_seqs_),
      sym_ops(proxy.sym_ops),
      angle_ideal(proxy.angle_ideal),
      weight(proxy.weight),
      periodicity(proxy.periodicity),
      alt_angle_ideals(proxy.alt_angle_ideals),
      limit(proxy.limit),
      top_out(proxy.top_out)
    {
      if ( sym_ops.get() != 0 ) {
        CCTBX_ASSERT(sym_ops.get()->size() == i_seqs.size());
      }
    }

    dihedral_proxy
    scale_weight(
      double factor) const
    {
      return dihedral_proxy(
        i_seqs, sym_ops, angle_ideal, weight*factor,
        periodicity, alt_angle_ideals, limit, top_out);
    }

    //! Sorts i_seqs such that i_seq[0] < i_seq[3] and i_seq[1] < i_seq[2].
    dihedral_proxy
    sort_i_seqs() const
    {
      dihedral_proxy result(*this);
      if (result.i_seqs[0] > result.i_seqs[3]) {
        std::swap(result.i_seqs[0], result.i_seqs[3]);
        if ( sym_ops.get() != 0 ) {
          std::swap(result.sym_ops[0], result.sym_ops[3]);
        }
        result.flip_angle_ideal();
      }
      if (result.i_seqs[1] > result.i_seqs[2]) {
        std::swap(result.i_seqs[1], result.i_seqs[2]);
        if ( sym_ops.get() != 0 ) {
          std::swap(result.sym_ops[1], result.sym_ops[2]);
        }
        result.flip_angle_ideal();
      }
      return result;
    }

    protected:

      //! Helper for sort_i_seqs().
      void
      flip_angle_ideal()
      {
        angle_ideal *= -1;
        if (alt_angle_ideals) {
          alt_angle_ideals_type::value_type& aai = *alt_angle_ideals;
          for(unsigned i=0;i<aai.size();i++) {
            aai[i] *= -1;
          }
        }
      }

      public:

    //! Indices into array of sites.
    i_seqs_type i_seqs;
    //! Optional array of symmetry operations.
    optional_container<af::shared<sgtbx::rt_mx> > sym_ops;
    //! Parameter.
    double angle_ideal;
    //! Parameter.
    double weight;
    //! Parameter.
    int periodicity;
    //! Optional array of alternative angle_ideal.
    alt_angle_ideals_type alt_angle_ideals;
    //! Optional value to set a range for which the angle should be
    //! restrained
    double limit;
    //! Use top-out function or not for residual/gradient
    bool top_out;
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
  class dihedral : protected scitbx::math::dihedral
  {
    public:
      //! Default constructor. Some data members are not initialized!
      dihedral() {}

      //! Constructor.
      dihedral(
        af::tiny<scitbx::vec3<double>, 4> const& sites_,
        double angle_ideal_,
        double weight_,
        int periodicity_=0,
        alt_angle_ideals_type const& alt_angle_ideals_=alt_angle_ideals_type(),
        double limit_=-1.0,
        bool top_out_=false)
      :
        sites(sites_),
        angle_ideal(angle_ideal_),
        weight(weight_),
        periodicity(periodicity_),
        alt_angle_ideals(alt_angle_ideals_),
        limit(limit_),
        top_out(top_out_)
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
        periodicity(proxy.periodicity),
        alt_angle_ideals(proxy.alt_angle_ideals),
        limit(proxy.limit),
        top_out(proxy.top_out)
      {
        for(int i=0;i<4;i++) {
          std::size_t i_seq = proxy.i_seqs[i];
          CCTBX_ASSERT(i_seq < sites_cart.size());
          sites[i] = sites_cart[i_seq];
        }
        init_angle_model();
      }

      /*! \brief Coordinates are obtained from sites_cart according
          to proxy.i_seqs by applying proxy.sym_ops and unit_cell,
          parameters are copied from proxy.
       */
      dihedral(
        uctbx::unit_cell const& unit_cell,
        af::const_ref<scitbx::vec3<double> > const& sites_cart,
        dihedral_proxy const& proxy)
      :
        angle_ideal(proxy.angle_ideal),
        weight(proxy.weight),
        periodicity(proxy.periodicity),
        alt_angle_ideals(proxy.alt_angle_ideals),
        limit(proxy.limit),
        top_out(proxy.top_out)
      {
        for(int i=0;i<4;i++) {
          std::size_t i_seq = proxy.i_seqs[i];
          CCTBX_ASSERT(i_seq < sites_cart.size());
          sites[i] = sites_cart[i_seq];
          if ( proxy.sym_ops.get() != 0 ) {
            sgtbx::rt_mx rt_mx = proxy.sym_ops[i];
            if ( !rt_mx.is_unit_mx() ) {
              sites[i] = unit_cell.orthogonalize(
                rt_mx * unit_cell.fractionalize(sites[i]));
            }
          }
        }
        init_angle_model();
      }

      //! Sinusoidal or harmonic function of delta.
      /*! With periodicity <= 0, the simple harmonic function

            weight * delta**2

          is used (Hendrickson, W.A. (1985). Meth. Enzym. 115, 252-270).
          This function has singularities at angle_ideal+-180/periodicity.

          With periodicity > 0, the sinusoidal function

            weight * 120**2 / (1 - cos(120)) / (periodicity * periodicity)
                   * (1 - cos(periodicity * delta))

          is used, similar to functions used in CHARMM and CNS
          (www.charmm.org, cns-online.org). This function has no
          singularities, is a good approximation of the harmonic
          function around angle_ideal+-120/periodicity, and also
          approximates results from QM calculations reasonably well.

          Run cctbx/geometry_restraints/tst_ext.py to obtain plot files
          for visually comparing the sinusoidal or harmonic functions.
       */
      double
      residual() const
      {
        using scitbx::constants::pi_180;
        double term, delta_local;
        double top;
        delta_local = delta;
        if (limit >= 0 && top_out == false){
          if(delta > limit){
            delta_local = limit;
          }
          else if(delta < (-1.0*limit)){
            delta_local = (-1.0*limit);
          }
        }
        if (periodicity > 0) {
          term = 9600. / (periodicity * periodicity)
               * (1 - std::cos(periodicity * delta_local * pi_180));
        }
        else if (top_out) {
          top = limit * limit;
          //top*(1-exp(-weight*x**2/top))
          term = top * (1.0-std::exp(-weight*delta_local*delta_local/top));
          return term;
        }
        else {
          term = delta_local * delta_local;
        }
        return weight * term;
      }
      //! Gradient of delta with respect to the four sites.
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
      grad_delta(double epsilon=1e-100) const
      {
        af::tiny<scitbx::vec3<double>, 4> result;
        if(limit >= 0 && top_out==false){
          if (std::fabs(delta) > limit){
            result.fill(scitbx::vec3<double>(0,0,0));
            return result;
          }
        }
        double d_21_norm = d_21.length_sq();
        if (   !have_angle_model
            || d_21_norm < epsilon
            || n_0121_norm < epsilon
            || n_2123_norm < epsilon) {
          result.fill(scitbx::vec3<double>(0,0,0));
        }
        else {
          using scitbx::constants::pi_180;
          double grad_factor = d_21.length() / pi_180;
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
      //! Gradient of R = w * delta^2 with respect to the four sites.
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
      gradients(double epsilon=1e-100) const
      {
        using scitbx::constants::pi_180;
        double grad_factor;
        double top;
        if (periodicity > 0) {
            grad_factor = 9600. * weight / periodicity * pi_180
                        * std::sin(periodicity * delta * pi_180);
        }
        else if (top_out) {
            //(2*weight*x/top)*exp(-(weight*x**2)/top)
            top = limit * limit;
            grad_factor = (2.0*weight*delta/top)*std::exp(-(weight*delta*delta)/top);
        }
        else {
            grad_factor = 2 * weight * delta;
        }
        af::tiny<scitbx::vec3<double>, 4> result = grad_delta(epsilon);
        for (std::size_t i=0; i<4; i++) {
          result[i] *= grad_factor;
        }
        return result;
      }

      //! Support for dihedral_residual_sum.
      /*! Not available in Python.
       */
      void
      add_gradients(
        af::ref<scitbx::vec3<double> > const& gradient_array,
        dihedral_proxy::i_seqs_type const& i_seqs) const
      {
        af::tiny<scitbx::vec3<double>, 4> grads = gradients();
        for(int i=0;i<4;i++) {
          gradient_array[i_seqs[i]] += grads[i];
        }
      }

      //! Support for dihedral_residual_sum.
      /*! Not available in Python.

          Inefficient implementation, r_inv_cart is not cached.
          TODO: use asu_mappings to take advantage of caching of r_inv_cart.
       */
      void
      add_gradients(
        uctbx::unit_cell const& unit_cell,
        af::ref<scitbx::vec3<double> > const& gradient_array,
        dihedral_proxy const& proxy) const
      {
        dihedral_proxy::i_seqs_type const& i_seqs = proxy.i_seqs;
        optional_container<af::shared<sgtbx::rt_mx> > const&
          sym_ops = proxy.sym_ops;
        af::tiny<scitbx::vec3<double>, 4> grads = gradients();
        for(int i=0;i<4;i++) {
          if ( sym_ops.get() != 0 && !sym_ops[i].is_unit_mx() ) {
            scitbx::mat3<double>
              r_inv_cart_ = r_inv_cart(unit_cell, sym_ops[i]);
            gradient_array[i_seqs[i]] += grads[i] * r_inv_cart_;
          }
          else { gradient_array[i_seqs[i]] += grads[i]; }
        }
      }

      void
      linearise(
        uctbx::unit_cell const& unit_cell,
        cctbx::restraints::linearised_eqns_of_restraint<double> &linearised_eqns,
        cctbx::xray::parameter_map<cctbx::xray::scatterer<double> > const &parameter_map,
        dihedral_proxy const& proxy) const
      {
        dihedral_proxy::i_seqs_type const& i_seqs = proxy.i_seqs;
        optional_container<af::shared<sgtbx::rt_mx> > const&
          sym_ops = proxy.sym_ops;
        af::tiny<scitbx::vec3<double>, 4> grads = grad_delta();
        std::size_t row_i = linearised_eqns.next_row();
        linearised_eqns.weights[row_i] = proxy.weight;
        linearised_eqns.deltas[row_i] = delta;
        for(int i=0;i<4;i++) {
          grads[i] = unit_cell.fractionalize_gradient(grads[i]);
          if ( sym_ops.get() != 0 && !sym_ops[i].is_unit_mx() ) {
            scitbx::mat3<double> r_inv
              = sym_ops[i].r().inverse().as_double();
            grads[i] = grads[i] * r_inv;
          }
          cctbx::xray::parameter_indices const &ids_i
            = parameter_map[i_seqs[i]];
          if (ids_i.site == -1) continue;
          for (int j=0;j<3;j++) {
            linearised_eqns.design_matrix(row_i, ids_i.site+j) = grads[i][j];
          }
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
      //! Optional array of alternative angle_ideal.
      alt_angle_ideals_type alt_angle_ideals;
      //! Optional limit on deviation for restraint calculation
      double limit;
      //! Use top-out function or notabilities
      bool top_out;
      //! false in singular situations.
      bool have_angle_model;
    public:
      //! Value of the dihedral %angle formed by the sites.
      double angle_model;
      /*! \brief Smallest difference between angle_model and angle_ideal
          taking the periodicity and alt_angle_ideals into account.
       */
      /*! See also: angle_delta_deg
       */
      double delta;

    protected:
      void
      init_angle_model()
      {
        scitbx::math::dihedral::init(sites.begin());
        boost::optional<double> angle_deg = angle(/* deg */ true);
        have_angle_model = bool(angle_deg);
        if (!have_angle_model) return;
        angle_model = *angle_deg;
        delta = angle_delta_deg(angle_model, angle_ideal, periodicity);
        if (alt_angle_ideals) {
          using scitbx::fn::absolute;
          alt_angle_ideals_type::value_type& ai = *alt_angle_ideals;
          delta = angle_delta_deg(angle_model, angle_ideal, periodicity);
          for(unsigned i_ai=0;i_ai<ai.size();i_ai++){
            double a = ai[i_ai];
            if (absolute(a) > 360.){
              throw std::runtime_error((
                boost::format(
                  " dihedral geometry restraint: invalid alt_angle_ideal:"
                  " range = 0-360, alt = %f") % a).str());
            }
            if (a < 0) a += 360;
            double delta_i = angle_delta_deg(
              angle_model, a);
            if (absolute(delta) > absolute(delta_i)){
              delta = delta_i;
            }
          }
        }
      }
  };

  //! Number of proxies with periodicity <= 0.
  inline
  std::size_t
  dihedral_count_harmonic(
    af::const_ref<dihedral_proxy> const& proxies)
  {
    std::size_t result = 0;
    for(std::size_t i=0;i<proxies.size();i++) {
      if (proxies[i].periodicity <= 0) result++;
    }
    return result;
  }

  /*! Fast computation of dihedral::delta given an array of dihedral proxies,
      ignoring proxy.sym_ops.
   */
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
      dihedral proxies, ignoring proxy.sym_ops.
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
      given an array of dihedral proxies, ignoring proxy.sym_ops.
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

  /*! Fast computation of dihedral::delta given an array of dihedral
      proxies, taking into account proxy.sym_ops.
   */
  inline
  af::shared<double>
  dihedral_deltas(
    uctbx::unit_cell const& unit_cell,
    af::const_ref<scitbx::vec3<double> > const& sites_cart,
    af::const_ref<dihedral_proxy> const& proxies)
  {
    return detail::generic_deltas<dihedral_proxy, dihedral>::get(
      unit_cell, sites_cart, proxies);
  }

  /*! Fast computation of dihedral::residual() given an array of
      dihedral proxies, taking into account proxy.sym_ops.
   */
  inline
  af::shared<double>
  dihedral_residuals(
    uctbx::unit_cell const& unit_cell,
    af::const_ref<scitbx::vec3<double> > const& sites_cart,
    af::const_ref<dihedral_proxy> const& proxies)
  {
    return detail::generic_residuals<dihedral_proxy, dihedral>::get(
      unit_cell, sites_cart, proxies);
  }

  /*! Fast computation of sum of dihedral::residual() and gradients
      given an array of dihedral proxies, taking into account
      proxy.sym_ops.
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
    uctbx::unit_cell const& unit_cell,
    af::const_ref<scitbx::vec3<double> > const& sites_cart,
    af::const_ref<dihedral_proxy> const& proxies,
    af::ref<scitbx::vec3<double> > const& gradient_array)
  {
    return detail::generic_residual_sum<dihedral_proxy, dihedral>::get(
      unit_cell, sites_cart, proxies, gradient_array);
  }

}} // namespace cctbx::geometry_restraints

#endif // CCTBX_GEOMETRY_RESTRAINTS_DIHEDRAL_H
