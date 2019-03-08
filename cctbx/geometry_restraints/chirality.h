#ifndef CCTBX_GEOMETRY_RESTRAINTS_CHIRALITY_H
#define CCTBX_GEOMETRY_RESTRAINTS_CHIRALITY_H

#include <cctbx/geometry_restraints/utils.h>
#include <cctbx/restraints.h>
#include <cctbx/geometry/geometry.h>
#include <scitbx/constants.h>

namespace cctbx { namespace geometry_restraints {

  //! Grouping of indices into array of sites (i_seqs) and parameters.
  struct chirality_proxy
  {
    //! Support for shared_proxy_select.
    typedef af::tiny<unsigned, 4> i_seqs_type;

    //! Default constructor. Some data members are not initialized!
    chirality_proxy() {}

    //! Constructor.
    chirality_proxy(
      i_seqs_type const& i_seqs_,
      double volume_ideal_,
      bool both_signs_,
      double weight_,
      unsigned char origin_id_)
    :
      i_seqs(i_seqs_),
      volume_ideal(volume_ideal_),
      both_signs(both_signs_),
      weight(weight_),
      origin_id(origin_id_)
    {}

    //! Constructor.
    chirality_proxy(
      i_seqs_type const& i_seqs_,
      optional_container<af::shared<sgtbx::rt_mx> > const& sym_ops_,
      double volume_ideal_,
      bool both_signs_,
      double weight_,
      unsigned char origin_id_)
    :
      i_seqs(i_seqs_),
      sym_ops(sym_ops_),
      volume_ideal(volume_ideal_),
      both_signs(both_signs_),
      weight(weight_),
      origin_id(origin_id_)
    {
      if ( sym_ops.get() != 0 ) {
        CCTBX_ASSERT(sym_ops.get()->size() == i_seqs.size());
      }
    }

    //! Support for proxy_select (and similar operations).
    chirality_proxy(
      i_seqs_type const& i_seqs_,
      chirality_proxy const& proxy)
    :
      i_seqs(i_seqs_),
      sym_ops(proxy.sym_ops),
      volume_ideal(proxy.volume_ideal),
      both_signs(proxy.both_signs),
      weight(proxy.weight),
      origin_id(proxy.origin_id)
    {
      if ( sym_ops.get() != 0 ) {
        CCTBX_ASSERT(sym_ops.get()->size() == i_seqs.size());
      }
    }

    chirality_proxy
    scale_weight(
      double factor) const
    {
      return chirality_proxy(i_seqs, sym_ops, volume_ideal, both_signs, weight*factor,
          origin_id);
    }

    //! Sorts i_seqs such that i_seq[1] < i_seq[2] < i_seq[3].
    chirality_proxy
    sort_i_seqs() const
    {
      chirality_proxy result(*this);
      for(unsigned i=1;i<3;i++) {
        for(unsigned j=i+1;j<4;j++) {
          if (result.i_seqs[i] > result.i_seqs[j]) {
            std::swap(result.i_seqs[i], result.i_seqs[j]);
            if ( sym_ops.get() != 0 ) {
              std::swap(result.sym_ops[i], result.sym_ops[j]);
            }
            if (!both_signs) result.volume_ideal *= -1;
          }
        }
      }
      return result;
    }

    //! Indices into array of sites.
    i_seqs_type i_seqs;
    //! Optional array of symmetry operations.
    optional_container<af::shared<sgtbx::rt_mx> > sym_ops;
    //! Parameter.
    double volume_ideal;
    //! Parameter.
    bool both_signs;
    //! Parameter.
    double weight;
    unsigned char origin_id;
  };

  //! Residual and gradient calculations for chirality restraint.
  class chirality
  {
    public:
      //! Default constructor. Some data members are not initialized!
      chirality() {}

      //! Constructor.
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

      /*! \brief Coordinates are copied from sites_cart according to
          proxy.i_seqs, parameters are copied from proxy.
       */
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

      chirality(
        uctbx::unit_cell const& unit_cell,
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
          if ( proxy.sym_ops.get() != 0 ) {
            sgtbx::rt_mx rt_mx = proxy.sym_ops[i];
            if ( !rt_mx.is_unit_mx() ) {
              sites[i] = unit_cell.orthogonalize(
                rt_mx * unit_cell.fractionalize(sites[i]));
            }
          }
        }
        init_volume_model();
      }

      //! weight * delta**2.
      /*! See also: Hendrickson, W.A. (1985). Meth. Enzym. 115, 252-270.
       */
      double
      residual() const { return weight * scitbx::fn::pow2(delta); }

      //! Gradients with respect to the four sites.
      /*! See also:
            Spiegel, M. R. & Liu, J. (1998). Mathematical Handbook of
            Formulas and Tables. New York: McGraw-Hill.
       */
      af::tiny<scitbx::vec3<double>, 4>
      gradients() const
      {
        af::tiny<scitbx::vec3<double>, 4> result;
        double f = delta_sign * 2.0 * delta * weight;
        result[1] = f * d_02_cross_d_03;
        result[2] = f * d_03.cross(d_01);
        result[3] = f * d_01.cross(d_02);
        result[0] = -result[1]-result[2]-result[3];
        return result;
      }

      //! Support for chirality_residual_sum.
      /*! Not available in Python.
       */
      void
      add_gradients(
        af::ref<scitbx::vec3<double> > const& gradient_array,
        chirality_proxy::i_seqs_type const& i_seqs) const
      {
        af::tiny<scitbx::vec3<double>, 4> grads = gradients();
        for(int i=0;i<4;i++) {
          gradient_array[i_seqs[i]] += grads[i];
        }
      }

      void
      linearise(
        uctbx::unit_cell const& unit_cell,
        cctbx::restraints::linearised_eqns_of_restraint<double> &linearised_eqns,
        cctbx::xray::parameter_map<cctbx::xray::scatterer<double> > const &parameter_map,
        chirality_proxy const& proxy) const
      {
        chirality_proxy::i_seqs_type const& i_seqs = proxy.i_seqs;
        af::tiny<scitbx::vec3<double>, 4> grads = gradients();
        std::size_t row_i = linearised_eqns.next_row();
        double f = 1.0 /( 2.0 * delta * weight );
        for(int i=0;i<4;i++) {
          grads[i] = unit_cell.fractionalize_gradient(grads[i]);
          if ( sym_ops.get() != 0 && !sym_ops[i].is_unit_mx() ) {
            scitbx::mat3<double> r_inv
              = sym_ops[i].r().inverse().as_double();
            grads[i] = grads[i] * r_inv;
          }
          cctbx::xray::parameter_indices const &ids_i =
            parameter_map[i_seqs[i]];
          if (ids_i.site == -1) continue;
          for (int j=0;j<3;j++) {
            linearised_eqns.design_matrix(row_i, ids_i.site+j) += f * grads[i][j];
          }

          linearised_eqns.weights[row_i] = proxy.weight;
          linearised_eqns.deltas[row_i] = volume_ideal + delta_sign * volume_model;
        }
      }

      //! Cartesian coordinates of the sites defining the chiral center.
      af::tiny<scitbx::vec3<double>, 4> sites;
      //! Optional array of symmetry operations.
      optional_container<af::shared<sgtbx::rt_mx> > sym_ops;
      //! Parameter (usually as passed to the constructor).
      double volume_ideal;
      //! Parameter (usually as passed to the constructor).
      bool both_signs;
      //! Parameter (usually as passed to the constructor).
      double weight;
    protected:
      scitbx::vec3<double> d_01;
      scitbx::vec3<double> d_02;
      scitbx::vec3<double> d_03;
      scitbx::vec3<double> d_02_cross_d_03;
    public:
      //! Value of the chiral volume defined by the sites.
      double volume_model;
      //! See implementation of init_volume_model().
      double delta_sign;
      //! volume_ideal + delta_sign * volume_model
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

  //! Fast computation of chirality::delta given an array of chirality proxies.
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
  chirality_deltas(
    uctbx::unit_cell const& unit_cell,
    af::const_ref<scitbx::vec3<double> > const& sites_cart,
    af::const_ref<chirality_proxy> const& proxies)
  {
    return detail::generic_deltas<chirality_proxy, chirality>::get(
      unit_cell, sites_cart, proxies);
  }

  /*! Fast computation of chirality::residual() given an array of
      chirality proxies.
   */
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
  af::shared<double>
  chirality_residuals(
    uctbx::unit_cell const& unit_cell,
    af::const_ref<scitbx::vec3<double> > const& sites_cart,
    af::const_ref<chirality_proxy> const& proxies)
  {
    return detail::generic_residuals<chirality_proxy, chirality>::get(
      unit_cell, sites_cart, proxies);
  }

  /*! Fast computation of sum of chirality::residual() and gradients
      given an array of chirality proxies.
   */
  /*! The chirality::gradients() are added to the gradient_array if
      gradient_array.size() == sites_cart.size().
      gradient_array must be initialized before this function
      is called.
      No gradient calculations are performed if gradient_array.size() == 0.
   */
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

}} // namespace cctbx::geometry_restraints

#endif // CCTBX_GEOMETRY_RESTRAINTS_CHIRALITY_H
