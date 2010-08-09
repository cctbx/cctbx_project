#ifndef CCTBX_GEOMETRY_RESTRAINTS_BOND_SIMILARITY_H
#define CCTBX_GEOMETRY_RESTRAINTS_BOND_SIMILARITY_H

#include <cctbx/sgtbx/rt_mx.h>
#include <cctbx/geometry_restraints/utils.h>
#include <cctbx/restraints.h>

namespace cctbx { namespace geometry_restraints {

  //! Grouping of indices into array of sites (i_seqs) symmetry
  //! operations and weights.
  struct bond_similarity_proxy
  {
    //! Support for shared_proxy_select.
    typedef af::shared<af::tiny<std::size_t, 2> > i_seqs_type;

    //! Default constructor. Some data members are not initialized!
    bond_similarity_proxy() {}

    //! Constructor.
    bond_similarity_proxy(
      i_seqs_type const& i_seqs_,
      af::shared<double> const& weights_)
    :
      weights(weights_),
      i_seqs(i_seqs_)
    {}

    //! Constructor.
    bond_similarity_proxy(
      i_seqs_type const& i_seqs_,
      af::shared<sgtbx::rt_mx> const& sym_ops_,
      af::shared<double> const& weights_)
    :
      weights(weights_),
      i_seqs(i_seqs_),
      sym_ops(sym_ops_)
    {}

    //! Array of weights.
    af::shared<double> weights;
    //! Array of pairs of indices into array of sites.
    i_seqs_type i_seqs;

    //! Array of symmetry operations to be applied to i_seqs[1]
    //! for i_seqs in i_seqs.
    /*! The bond is between sites_frac[i_seqs[0]] and
        rt_mx_ji * sites_frac[i_seqs[1]].
     */
    optional_container<af::shared<sgtbx::rt_mx> > sym_ops;
  };

  //! Residual and gradient calculations for harmonically restrained bonds.
  class bond_similarity
  {
    public:
      //! Convenience typedef.
      typedef scitbx::vec3<double> vec3;

      //! Default constructor. Some data members are not initialized!
      bond_similarity() {}

      //! Constructor.
      bond_similarity(
        af::shared<af::tiny<vec3, 2> > const& sites_array_,
        af::shared<double> const& weights_)
      :
        sites_array(sites_array_),
        weights(weights_)
        {
          init_deltas();
        }

      /*! \brief Coordinates are copied from sites_cart according to
          proxy.i_seqs, parameters are copied from proxy.
       */
      bond_similarity(
        af::const_ref<vec3> const& sites_cart,
        bond_similarity_proxy const& proxy)
      :
        weights(proxy.weights)
      {
        sites_array.reserve(proxy.i_seqs.size());
        for(int i=0;i<proxy.i_seqs.size();i++) {
          af::tiny<vec3, 2> sites;
          for(int j=0;j<2;j++) {
            std::size_t i_seq = proxy.i_seqs[i][j];
            CCTBX_ASSERT(i_seq < sites_cart.size());
            sites[j] = sites_cart[i_seq];
          }
          sites_array.push_back(sites);
        }
        init_deltas();
      }

      /*! \brief Coordinates are obtained from sites_cart according
          to proxy.i_seqs by applying proxy.sym_ops and unit_cell,
          parameters are copied from proxy.
       */
      bond_similarity(
        uctbx::unit_cell const& unit_cell,
        af::const_ref<vec3> const& sites_cart,
        bond_similarity_proxy const& proxy)
      :
        weights(proxy.weights)
      {
        sites_array.reserve(proxy.i_seqs.size());
        for(int i=0;i<proxy.i_seqs.size();i++) {
          af::tiny<vec3, 2> sites;
          for(int j=0;j<2;j++) {
            std::size_t i_seq = proxy.i_seqs[i][j];
            CCTBX_ASSERT(i_seq < sites_cart.size());
            sites[j] = sites_cart[i_seq];
          }
          if ( proxy.sym_ops.get() != 0 ) {
            sgtbx::rt_mx rt_mx = proxy.sym_ops[i];
            if ( !rt_mx.is_unit_mx() ) {
              sites[1] = unit_cell.orthogonalize(
                rt_mx * unit_cell.fractionalize(sites[1]));
            }
          }
          sites_array.push_back(sites);
        }
        init_deltas();
      }

      //! Array of deviation from mean bond distance for each pair.
      af::shared<double> const&
      deltas() const { return deltas_; }

      //! sqrt(mean_sq(deltas))
      double
      rms_deltas() const {
        return std::sqrt(af::mean_sq(deltas_.const_ref()));
      }

      //! sum(weight_i * delta_i**2) / sum(weights).
      double
      residual() const
      {
        double result = 0;
        af::const_ref<double> weights_ref = weights.const_ref();
        af::const_ref<double> deltas_ref = deltas_.const_ref();
        for(std::size_t i_site=0;i_site<deltas_ref.size();i_site++) {
          result += weights_ref[i_site]
                  * scitbx::fn::pow2(deltas_ref[i_site])
                  / sum_weights_;
        }
        return result;
      }

      //! Gradients with respect to the sites.
      af::shared<af::tiny<vec3, 2> >
      gradients() const
      {
        af::shared<af::tiny<vec3, 2> > result;
        af::tiny<vec3, 2> pair_grads;
        af::const_ref<double> weights_ref = weights.const_ref();
        af::const_ref<double> deltas_ref = deltas_.const_ref();
        af::const_ref<double> distances_ref
          = bond_distances_.const_ref();
        result.reserve(deltas_ref.size());
        for(std::size_t i_pair=0;i_pair<deltas_ref.size();i_pair++) {
          vec3 grad_0 = (2 * weights_ref[i_pair] * deltas_ref[i_pair])
                      / (distances_ref[i_pair] * sum_weights_)
                      * (sites_array[i_pair][0] - sites_array[i_pair][1]);
          pair_grads[0] = grad_0;
          pair_grads[1] = -grad_0;
          result.push_back(pair_grads);
        }
        return result;
      }

      void
      linearise(
        uctbx::unit_cell const& unit_cell,
        cctbx::restraints::linearised_eqns_of_restraint<double> &linearised_eqns,
        cctbx::xray::parameter_map<cctbx::xray::scatterer<double> > const &parameter_map,
        bond_similarity_proxy const& proxy) const
      {
        bond_similarity_proxy::i_seqs_type const& i_seqs = proxy.i_seqs;
        optional_container<af::shared<sgtbx::rt_mx> > const&
          sym_ops = proxy.sym_ops;
        for(std::size_t i_pair=0;i_pair<deltas_.size();i_pair++) {
          std::size_t row_i = linearised_eqns.next_row();
          linearised_eqns.weights[row_i] = weights[i_pair];
          linearised_eqns.deltas[row_i] = deltas_[i_pair];
          vec3 grad_i = (1 - weights[i_pair]/sum_weights_)
            * (sites_array[i_pair][0] - sites_array[i_pair][1])
            / bond_distances_[i_pair];
          grad_i = unit_cell.fractionalize_gradient(grad_i);
          for(int i=0;i<2;i++) {
            if (i == 1) {
              grad_i = -grad_i;
              if (sym_ops.get() != 0 && !sym_ops[i].is_unit_mx() ) {
                scitbx::mat3<double> r_inv
                  = sym_ops[i].r().inverse().as_double();
                grad_i = grad_i * r_inv;
              }
            }
            cctbx::xray::parameter_indices const &ids_i
              = parameter_map[i_seqs[i_pair][i]];
            if (ids_i.site == -1) continue;
            for (int j=0;j<3;j++) {
              linearised_eqns.design_matrix(row_i, ids_i.site+j) = grad_i[j];
            }
          }
        }
      }

      //! Support for bond_similarity_residual_sum.
      /*! Not available in Python.
       */
      void
      add_gradients(
        af::ref<scitbx::vec3<double> > const& gradient_array,
        bond_similarity_proxy::i_seqs_type const& i_seqs) const
      {
        af::const_ref<af::tiny<std::size_t, 2> > i_seqs_ref = i_seqs.const_ref();
        af::shared<af::tiny<vec3, 2> > grads = gradients();
        af::const_ref<af::tiny<vec3, 2> > grads_ref = grads.const_ref();
        for(std::size_t i=0;i<grads_ref.size();i++) {
          gradient_array[i_seqs_ref[i][0]] += grads_ref[i][0];
          gradient_array[i_seqs_ref[i][1]] += grads_ref[i][1];
        }
      }

      //! Support for bond_similarity_residual_sum.
      /*! Not available in Python.

          Inefficient implementation, r_inv_cart is not cached.
          TODO: use asu_mappings to take advantage of caching of r_inv_cart.
       */
      void
      add_gradients(
        uctbx::unit_cell const& unit_cell,
        af::ref<scitbx::vec3<double> > const& gradient_array,
        bond_similarity_proxy const& proxy) const
      {
        af::const_ref<af::tiny<std::size_t, 2> > i_seqs_ref
          = proxy.i_seqs.const_ref();
        optional_container<af::shared<sgtbx::rt_mx> > const&
          sym_ops = proxy.sym_ops;
        af::shared<af::tiny<vec3, 2> > grads = gradients();
        af::const_ref<af::tiny<vec3, 2> > grads_ref = grads.const_ref();
        for(std::size_t i=0;i<grads_ref.size();i++) {
          gradient_array[i_seqs_ref[i][0]] += grads_ref[i][0];
          if ( sym_ops.get() != 0 && !sym_ops[i].is_unit_mx() ) {
            scitbx::mat3<double>
              r_inv_cart_ = r_inv_cart(unit_cell, sym_ops[i]);
            gradient_array[i_seqs_ref[i][1]] += grads_ref[i][1] * r_inv_cart_;
          }
          else gradient_array[i_seqs_ref[i][1]] += grads_ref[i][1];
        }
      }

      //! Weighted mean of the bond lengths
      double
      mean_distance() const { return mean_distance_; }
      //! Cartesian coordinates of bonded sites.
      af::shared<af::tiny<vec3, 2> > sites_array;
      //! Array of weights for each pair.
      af::shared<double> weights;

    protected:
      double mean_distance_;
      double sum_weights_;
      af::shared<double> deltas_;
      af::shared<double> bond_distances_;

      void
      init_deltas()
      {
        af::const_ref<af::tiny<vec3, 2> > sites_ref
          = sites_array.const_ref();
        af::const_ref<double> weights_ref = weights.const_ref();
        bond_distances_.reserve(sites_ref.size());
        mean_distance_ = 0;
        sum_weights_ = 0;
        for (std::size_t i_pair=0;i_pair<sites_array.size();i_pair++) {
          double w = weights_ref[i_pair];
          af::tiny<vec3, 2> sites = sites_array[i_pair];
          bond_distances_.push_back((sites[0] - sites[1]).length());
          mean_distance_ += w * bond_distances_[i_pair];
          sum_weights_ += w;
        }
        CCTBX_ASSERT(sum_weights_ > 0);
        mean_distance_ /= sum_weights_;
        deltas_.reserve(sites_ref.size());
        for (std::size_t i_pair=0;i_pair<sites_array.size();i_pair++) {
          deltas_.push_back(bond_distances_[i_pair] - mean_distance_);
        }
      }
  };

  /*! \brief Fast computation of bond_similarity::rms_deltas() given an
      array of bond_similarity proxies, ignoring proxy.sym_ops.
   */
  inline
  af::shared<double>
  bond_similarity_deltas_rms(
    af::const_ref<scitbx::vec3<double> > const& sites_cart,
    af::const_ref<bond_similarity_proxy> const& proxies)
  {
    af::shared<double> result((af::reserve(proxies.size())));
    for(std::size_t i=0;i<proxies.size();i++) {
      result.push_back(bond_similarity(sites_cart, proxies[i]).rms_deltas());
    }
    return result;
  }

  /*! \brief Fast computation of bond_similarity::residual() given an
      array of bond_similarity proxies, ignoring proxy.sym_ops.
   */
  inline
  af::shared<double>
  bond_similarity_residuals(
    af::const_ref<scitbx::vec3<double> > const& sites_cart,
    af::const_ref<bond_similarity_proxy> const& proxies)
  {
    return detail::generic_residuals<bond_similarity_proxy,
      bond_similarity>::get(sites_cart, proxies);
  }

  /*! Fast computation of sum of bond_similarity::residual() and gradients given
      an array of bond_similarity proxies, ignoring proxy.sym_ops.
   */
  /*! The bond_similarity::gradients() are added to the gradient_array
      if gradient_array.size() == sites_cart.size().
      gradient_array must be initialized before this function
      is called.
      No gradient calculations are performed if gradient_array.size() == 0.
   */
  inline
  double
  bond_similarity_residual_sum(
    af::const_ref<scitbx::vec3<double> > const& sites_cart,
    af::const_ref<bond_similarity_proxy> const& proxies,
    af::ref<scitbx::vec3<double> > const& gradient_array)
  {
    return detail::generic_residual_sum<bond_similarity_proxy,
      bond_similarity>::get(sites_cart, proxies, gradient_array);
  }

  /*! \brief Fast computation of bond_similarity::rms_deltas() given an
      array of bond_similarity proxies, taking account of proxy.sym_ops.
   */
  inline
  af::shared<double>
  bond_similarity_deltas_rms(
    uctbx::unit_cell const& unit_cell,
    af::const_ref<scitbx::vec3<double> > const& sites_cart,
    af::const_ref<bond_similarity_proxy> const& proxies)
  {
    af::shared<double> result((af::reserve(proxies.size())));
    for(std::size_t i=0;i<proxies.size();i++) {
      result.push_back(bond_similarity(
        unit_cell, sites_cart, proxies[i]).rms_deltas());
    }
    return result;
  }

  /*! \brief Fast computation of bond_similarity::residual() given an
      array of bond_similarity proxies, taking account of proxy.sym_ops.
   */
  inline
  af::shared<double>
  bond_similarity_residuals(
    uctbx::unit_cell const& unit_cell,
    af::const_ref<scitbx::vec3<double> > const& sites_cart,
    af::const_ref<bond_similarity_proxy> const& proxies)
  {
    return detail::generic_residuals<bond_similarity_proxy,
      bond_similarity>::get(
      unit_cell, sites_cart, proxies);
  }

  /*! Fast computation of sum of bond_similarity::residual() and gradients given
      an array of bond_similarity proxies, taking account of proxy.sym_ops.
   */
  /*! The bond_similarity::gradients() are added to the gradient_array
      if gradient_array.size() == sites_cart.size().
      gradient_array must be initialized before this function
      is called.
      No gradient calculations are performed if gradient_array.size() == 0.
   */
  inline
  double
  bond_similarity_residual_sum(
    uctbx::unit_cell const& unit_cell,
    af::const_ref<scitbx::vec3<double> > const& sites_cart,
    af::const_ref<bond_similarity_proxy> const& proxies,
    af::ref<scitbx::vec3<double> > const& gradient_array)
  {
    return detail::generic_residual_sum<bond_similarity_proxy,
      bond_similarity>::get(
      unit_cell, sites_cart, proxies, gradient_array);
  }

}} // namespace cctbx::geometry_restraints

#endif // CCTBX_GEOMETRY_RESTRAINTS_BOND_SIMILARITY_H
