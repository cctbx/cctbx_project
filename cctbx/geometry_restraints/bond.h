#ifndef CCTBX_GEOMETRY_RESTRAINTS_BOND_H
#define CCTBX_GEOMETRY_RESTRAINTS_BOND_H

#include <cctbx/geometry_restraints/utils.h>
#include <cctbx/geometry_restraints/asu_cache.h>
#include <cctbx/geometry_restraints/sorted_asu_proxies.h>
#include <cctbx/geometry/geometry.h>
#include <cctbx/restraints.h>

#include <tbxx/optional_copy.hpp>

namespace cctbx { namespace geometry_restraints {

  //! Grouping of bond parameters distance_ideal and weight.
  struct bond_params
  {
    //! Default constructor. Some data members are not initialized!
    bond_params() {}

    //! Constructor.
    bond_params(
      double distance_ideal_,
      double weight_,
      double slack_=0,
      double limit_=-1.0,
      bool top_out_=false)
    :
      distance_ideal(distance_ideal_),
      weight(weight_),
      slack(slack_),
      limit(limit_),
      top_out(top_out_)
    {}

    bond_params
    scale_weight(
      double factor) const
    {
      return bond_params(distance_ideal, weight*factor, slack);
    }

    //! Parameter.
    double distance_ideal;
    //! Parameter.
    double weight;
    //! Parameter.
    double slack;
    //! Parameter.
    double limit;
    //! Parameter.
    bool top_out;
  };

  //! Dictionary of bond parameters for a given index i_seq.
  typedef std::map<unsigned, bond_params> bond_params_dict;
  //! Array of dictionary of bond parameters. The index i_seq is implied.
  typedef af::shared<bond_params_dict> bond_params_table;

  //! Grouping of indices into array of sites (i_seqs) and bond_params.
  struct bond_simple_proxy : bond_params
  {
    //! Default constructor. Some data members are not initialized!
    bond_simple_proxy() {}

    //! Constructor.
    bond_simple_proxy(
      af::tiny<unsigned, 2> const& i_seqs_,
      double distance_ideal_,
      double weight_,
      double slack_=0,
      double limit_=-1.0,
      bool top_out_=false)
    :
      bond_params(distance_ideal_, weight_, slack_, limit_, top_out_),
      i_seqs(i_seqs_)
    {}

    //! Constructor.
    bond_simple_proxy(
      af::tiny<unsigned, 2> const& i_seqs_,
      sgtbx::rt_mx const& rt_mx_ji_,
      double distance_ideal_,
      double weight_,
      double slack_=0,
      double limit_=-1.0,
      bool top_out_=false)
    :
      bond_params(distance_ideal_, weight_, slack_, limit_, top_out_),
      i_seqs(i_seqs_),
      rt_mx_ji(rt_mx_ji_)
    {}

    //! Constructor.
    /*! Not available in Python.
     */
    bond_simple_proxy(
      af::tiny<unsigned, 2> const& i_seqs_,
      bond_params const& params)
    :
      bond_params(params),
      i_seqs(i_seqs_)
    {}

    //! Constructor.
    /*! Not available in Python.
     */
    bond_simple_proxy(
      af::tiny<unsigned, 2> const& i_seqs_,
      sgtbx::rt_mx const& rt_mx_ji_,
      bond_params const& params)
    :
      bond_params(params),
      i_seqs(i_seqs_),
      rt_mx_ji(rt_mx_ji_)
    {}

    //! Sorts i_seqs such that i_seq[0] < i_seq[1].
    bond_simple_proxy
    sort_i_seqs() const
    {
      bond_simple_proxy result(*this);
      if (result.i_seqs[0] > result.i_seqs[1]) {
        std::swap(result.i_seqs[0], result.i_seqs[1]);
      }
      return result;
    }

    //! Indices into array of sites.
    af::tiny<unsigned, 2> i_seqs;

    //! Optional symmetry operation operating on i_seqs[1]
    tbxx::optional_copy<sgtbx::rt_mx> rt_mx_ji;
  };

  //! Grouping of bond_simple_proxy and symmetry operation (rt_mx_ji).
  struct bond_sym_proxy : bond_params
  {
    //! Default constructor. Some data members are not initialized!
    bond_sym_proxy() {}

    //! Constructor.
    bond_sym_proxy(
      af::tiny<unsigned, 2> const& i_seqs_,
      sgtbx::rt_mx const& rt_mx_ji_,
      double distance_ideal_,
      double weight_,
      double slack_=0,
      double limit_=-1.0,
      bool top_out_=false)
    :
      bond_params(distance_ideal_, weight_, slack_, limit_, top_out_),
      i_seqs(i_seqs_),
      rt_mx_ji(rt_mx_ji_)
    {}

    //! Constructor.
    /*! Not available in Python.
     */
    bond_sym_proxy(
      af::tiny<unsigned, 2> const& i_seqs_,
      sgtbx::rt_mx const& rt_mx_ji_,
      bond_params const& params)
    :
      bond_params(params),
      i_seqs(i_seqs_),
      rt_mx_ji(rt_mx_ji_)
    {}

    //! Indices into array of sites.
    af::tiny<unsigned, 2> i_seqs;

    //! Symmetry operation to be applied to i_seqs[1].
    /*! The bond is between sites_frac[i_seqs[0]] and
        rt_mx_ji * sites_frac[i_seqs[1]].
     */
    sgtbx::rt_mx rt_mx_ji;
  };

  //! Grouping of asu_mapping_index_pair and bond_params.
  struct bond_asu_proxy : bond_params, asu_mapping_index_pair
  {
    //! Default constructor. Some data members are not initialized!
    bond_asu_proxy() {}

    //! Constructor.
    bond_asu_proxy(
      asu_mapping_index_pair const& pair_,
      double distance_ideal_,
      double weight_,
      double slack_=0,
      double limit_=-1.0,
      bool top_out_=false)
    :
      bond_params(distance_ideal_, weight_, slack_, limit_, top_out_),
      asu_mapping_index_pair(pair_)
    {}

    //! Constructor.
    bond_asu_proxy(
      asu_mapping_index_pair const& pair_,
      bond_params const& params)
    :
      bond_params(params),
      asu_mapping_index_pair(pair_)
    {}

    //! Conversion to bond_simple_proxy.
    bond_simple_proxy
    as_simple_proxy() const
    {
      return bond_simple_proxy(
        af::tiny<unsigned, 2>(i_seq, j_seq),
        distance_ideal,
        weight,
        slack);
    }
  };

  //! Residual and gradient calculations for harmonically restrained bonds.
  class bond : public cctbx::geometry::distance<double>, public bond_params
  {
    public:
      //! Convenience typedef.
      typedef scitbx::vec3<double> vec3;

      //! Default constructor. Some data members are not initialized!
      bond() {}

      //! Constructor.
      bond(
        af::tiny<scitbx::vec3<double>, 2> const& sites_,
        double distance_ideal_,
        double weight_,
        double slack_=0,
        double limit_=-1.0,
        bool top_out_=false)
      :
        cctbx::geometry::distance<double>(sites_),
        bond_params(distance_ideal_, weight_, slack_, limit_, top_out_)
      {
        init_deltas();
      }

      /*! \brief Coordinates are copied from sites_cart according to
          proxy.i_seqs, parameters are copied from proxy.
       */
      bond(
        af::const_ref<scitbx::vec3<double> > const& sites_cart,
        bond_simple_proxy const& proxy)
      :
        bond_params(proxy.distance_ideal, proxy.weight, proxy.slack,
          proxy.limit, proxy.top_out)
      {
        for(int i=0;i<2;i++) {
          std::size_t i_seq = proxy.i_seqs[i];
          CCTBX_ASSERT(i_seq < sites_cart.size());
          sites[i] = sites_cart[i_seq];
        }
        init_distance_model();
        init_deltas();
      }

      /*! \brief Coordinates are copied from sites_cart according to
          proxy.i_seqs, parameters are copied from proxy.
       */
      bond(
        uctbx::unit_cell const& unit_cell,
        af::const_ref<scitbx::vec3<double> > const& sites_cart,
        bond_simple_proxy const& proxy)
      :
        bond_params(proxy.distance_ideal, proxy.weight, proxy.slack,
          proxy.limit, proxy.top_out)
      {
        for(int i=0;i<2;i++) {
          std::size_t i_seq = proxy.i_seqs[i];
          CCTBX_ASSERT(i_seq < sites_cart.size());
          sites[i] = sites_cart[i_seq];
        }
        if ( proxy.rt_mx_ji ) {
          sites[1] = unit_cell.orthogonalize(
            *proxy.rt_mx_ji * unit_cell.fractionalize(sites[1]));
        }
        init_distance_model();
        init_deltas();
      }

      /*! \brief Coordinates are copied from sites_cart according to
          proxy.i_seqs and proxy.rt_mx_ji, parameters are copied from
          proxy.
       */
      bond(
        uctbx::unit_cell const& unit_cell,
        af::const_ref<scitbx::vec3<double> > const& sites_cart,
        bond_sym_proxy const& proxy)
      :
        bond_params(proxy.distance_ideal, proxy.weight, proxy.slack,
          proxy.limit, proxy.top_out)
      {
        for(int i=0;i<2;i++) {
          std::size_t i_seq = proxy.i_seqs[i];
          CCTBX_ASSERT(i_seq < sites_cart.size());
          sites[i] = sites_cart[i_seq];
        }
        sites[1] = unit_cell.orthogonalize(
          proxy.rt_mx_ji * unit_cell.fractionalize(sites[1]));
        init_distance_model();
        init_deltas();
      }

      /*! \brief Coordinates are copied from sites_cart according to
          proxy.i_seq, proxy.j_seq, parameters are copied from proxy.
       */
      bond(
        af::const_ref<scitbx::vec3<double> > const& sites_cart,
        direct_space_asu::asu_mappings<> const& asu_mappings,
        bond_asu_proxy const& proxy)
      :
        bond_params(proxy.distance_ideal, proxy.weight, proxy.slack,
          proxy.limit, proxy.top_out)
      {
        sites[0] = asu_mappings.map_moved_site_to_asu(
          sites_cart[proxy.i_seq], proxy.i_seq, 0);
        sites[1] = asu_mappings.map_moved_site_to_asu(
          sites_cart[proxy.j_seq], proxy.j_seq, proxy.j_sym);
        init_distance_model();
        init_deltas();
      }

      //! For fast processing. Not available in Python.
      bond(
        asu_cache<> const& cache,
        bond_asu_proxy const& proxy)
      :
        bond_params(proxy.distance_ideal, proxy.weight, proxy.slack,
          proxy.limit, proxy.top_out)
      {
        sites[0] = cache.sites[proxy.i_seq][0];
        sites[1] = cache.sites[proxy.j_seq][proxy.j_sym];
        init_distance_model();
        init_deltas();
      }

      //! weight * delta_slack**2.
      /*! See also: Hendrickson, W.A. (1985). Meth. Enzym. 115, 252-270.
       */
      double
      residual() const {
        // unlike the dihedral angle restraint, the harmonic potential is
        // always used if distance_model < distance_ideal, to compensate for
        // the lack of a nonbonded restraint for the bonded atoms.
        if ((top_out) && (delta_slack < 0)) {
          double top = weight * limit * limit;
          //top*(1-exp(-weight*x**2/top))
          return top * (1.0-std::exp(-weight*delta_slack*delta_slack/top));
        } else {
          return weight * scitbx::fn::pow2(delta_slack);
        }
      }

      //! Gradient of R = w * sum(deltas) with respect to sites[0].
      /*! Not available in Python.
       */
      scitbx::vec3<double>
      gradient_0(double epsilon=1.e-100) const
      {
        if (distance_model < epsilon) return scitbx::vec3<double>(0,0,0);
        if (delta < -slack || delta > slack) {
          if ((top_out) && (delta_slack < 0)) {
            double top = weight * limit * limit;
            return weight * 2 * delta_slack * \
              std::exp(-(weight*delta_slack*delta_slack)/top) * \
              d_distance_d_site_0(epsilon);
          } else {
            return weight * 2 * delta_slack * d_distance_d_site_0(epsilon);
          }
        }
        return scitbx::vec3<double>(0,0,0);
      }

      //! Gradient of R = w * deltas^2 with respect to each site.
      af::tiny<scitbx::vec3<double>, 2>
      gradients() const
      {
        af::tiny<vec3, 2> result;
        result[0] = gradient_0();
        result[1] = -result[0];
        return result;
      }

      //! Support for bond_residual_sum.
      /*! Not available in Python.
       */
      void
      add_gradients(
        af::ref<scitbx::vec3<double> > const& gradient_array,
        af::tiny<unsigned, 2> const& i_seqs) const
      {
        vec3 g0 = gradient_0();
        gradient_array[i_seqs[0]] += g0;
        gradient_array[i_seqs[1]] += -g0;
      }

      //! Support for bond_residual_sum.
      /*! Not available in Python.
       */
      void
      add_gradients(
        af::ref<scitbx::vec3<double> > const& gradient_array,
        direct_space_asu::asu_mappings<> const& asu_mappings,
        asu_mapping_index_pair const& pair) const
      {
        vec3 grad_asu = gradient_0();
        vec3 grad_i_seq = asu_mappings.r_inv_cart(pair.i_seq, 0) * grad_asu;
        gradient_array[pair.i_seq] += grad_i_seq;
        if (pair.j_sym == 0) {
          vec3 grad_j_seq = asu_mappings.r_inv_cart(pair.j_seq, 0) * grad_asu;
          gradient_array[pair.j_seq] -= grad_j_seq;
        }
      }

      //! Support for bond_residual_sum.
      /*! Not available in Python.
       */
      void
      add_gradients(
        asu_cache<>& cache,
        asu_mapping_index_pair const& pair) const
      {
        vec3 grad_asu = gradient_0();
        cache.gradients[pair.i_seq] += grad_asu;
        if (pair.j_sym == 0) {
          cache.gradients[pair.j_seq] -= grad_asu;
        }
      }

      void
      //  The linearised equation of restraint.
      linearise(
        uctbx::unit_cell const& unit_cell,
        cctbx::restraints::linearised_eqns_of_restraint<double> &linearised_eqns,
        cctbx::xray::parameter_map<cctbx::xray::scatterer<double> > const &parameter_map,
        bond_simple_proxy const& proxy) const
      {
        af::tiny<unsigned, 2> const& i_seqs = proxy.i_seqs;
        af::tiny<scitbx::vec3<double>, 2> grads;
        if (delta < -slack || delta > slack) {
          grads[0] = d_distance_d_site_0();
        }
        else { grads[0] = scitbx::vec3<double>(0,0,0); }

        grads[1] = -grads[0];
        std::size_t row_i = linearised_eqns.next_row();
        linearised_eqns.weights[row_i] = proxy.weight;
        linearised_eqns.deltas[row_i] = delta_slack;
        for(int i=0;i<2;i++) {
          grads[i] = unit_cell.fractionalize_gradient(grads[i]);
          if (i == 1 && proxy.rt_mx_ji) {
            scitbx::mat3<double> r_inv
              = proxy.rt_mx_ji->r().inverse().as_double();
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

      //! Difference distance_ideal - distance_model.
      double delta;
      //! sign(delta) * max(0, (abs(delta) - slack)).
      double delta_slack;

    protected:
      void
      init_deltas()
      {
        delta = distance_ideal - distance_model;
        CCTBX_ASSERT(slack >= 0);
        if (delta > slack) {
          delta_slack = delta - slack;
        }
        else if (delta >= -slack) {
          delta_slack = 0;
        }
        else {
          delta_slack = delta + slack;
        }
      }
  };

  //! Extracts bond parameters from array of simple proxies.
  inline
  bond_params_table
  extract_bond_params(
    std::size_t n_seq,
    af::const_ref<bond_simple_proxy> const& bond_simple_proxies)
  {
    bond_params_table tab(n_seq);
    af::ref<bond_params_dict> tab_ref = tab.ref();
    for(std::size_t i_proxy=0;i_proxy<bond_simple_proxies.size();i_proxy++) {
      af::tiny<unsigned, 2> const& i_seqs=bond_simple_proxies[i_proxy].i_seqs;
      CCTBX_ASSERT(i_seqs[0] < tab_ref.size());
      CCTBX_ASSERT(i_seqs[1] < tab_ref.size());
      if (i_seqs[0] < i_seqs[1]) {
        tab_ref[i_seqs[0]][i_seqs[1]] = bond_simple_proxies[i_proxy];
      }
      else {
        tab_ref[i_seqs[1]][i_seqs[0]] = bond_simple_proxies[i_proxy];
      }
    }
    return tab;
  }

  //! Fast computation of bond::distance_model given an array of bond proxies.
  inline
  af::shared<double>
  bond_distances_model(
    af::const_ref<scitbx::vec3<double> > const& sites_cart,
    af::const_ref<bond_simple_proxy> const& proxies)
  {
    af::shared<double> result((af::reserve(proxies.size())));
    for(std::size_t i=0;i<proxies.size();i++) {
      bond restraint(sites_cart, proxies[i]);
      result.push_back(restraint.distance_model);
    }
    return result;
  }

  //! Fast computation of bond::delta given an array of bond proxies.
  inline
  af::shared<double>
  bond_deltas(
    af::const_ref<scitbx::vec3<double> > const& sites_cart,
    af::const_ref<bond_simple_proxy> const& proxies)
  {
    return detail::generic_deltas<bond_simple_proxy, bond>::get(
      sites_cart, proxies);
  }

  //! Fast computation of bond::delta given an array of bond proxies.
  inline
  af::shared<double>
  bond_deltas(
    uctbx::unit_cell const& unit_cell,
    af::const_ref<scitbx::vec3<double> > const& sites_cart,
    af::const_ref<bond_simple_proxy> const& proxies)
  {
    return detail::generic_deltas<bond_simple_proxy, bond>::get(
      unit_cell, sites_cart, proxies);
  }

  //! Fast computation of bond::residual() given an array of bond proxies.
  inline
  af::shared<double>
  bond_residuals(
    af::const_ref<scitbx::vec3<double> > const& sites_cart,
    af::const_ref<bond_simple_proxy> const& proxies)
  {
    return detail::generic_residuals<bond_simple_proxy, bond>::get(
      sites_cart, proxies);
  }

  //! Fast computation of bond::residual() given an array of bond proxies.
  inline
  af::shared<double>
  bond_residuals(
    uctbx::unit_cell const& unit_cell,
    af::const_ref<scitbx::vec3<double> > const& sites_cart,
    af::const_ref<bond_simple_proxy> const& proxies)
  {
    return detail::generic_residuals<bond_simple_proxy, bond>::get(
      unit_cell, sites_cart, proxies);
  }

  /*! Fast computation of sum of bond::residual() and gradients
      given an array of bond proxies.
   */
  /*! The bond::gradients() are added to the gradient_array if
      gradient_array.size() == sites_cart.size().
      gradient_array must be initialized before this function
      is called.
      No gradient calculations are performed if gradient_array.size() == 0.
   */
  inline
  double
  bond_residual_sum(
    af::const_ref<scitbx::vec3<double> > const& sites_cart,
    af::const_ref<bond_simple_proxy> const& proxies,
    af::ref<scitbx::vec3<double> > const& gradient_array)
  {
    return detail::generic_residual_sum<bond_simple_proxy, bond>::get(
      sites_cart, proxies, gradient_array);
  }

  //! Managed group of bond_simple_proxy and bond_asu_proxy arrays.
  typedef sorted_asu_proxies<bond_simple_proxy, bond_asu_proxy>
    bond_sorted_asu_proxies_base;

  //! Fast computation of bond::distance_model given managed proxies.
  inline
  af::shared<double>
  bond_distances_model(
    af::const_ref<scitbx::vec3<double> > const& sites_cart,
    bond_sorted_asu_proxies_base const& sorted_asu_proxies)
  {
    af::shared<double> result = bond_distances_model(
      sites_cart, sorted_asu_proxies.simple.const_ref());
    af::const_ref<bond_asu_proxy> sym = sorted_asu_proxies.asu.const_ref();
    if (sym.size() > 0) {
      result.reserve(sorted_asu_proxies.simple.size() + sym.size());
      direct_space_asu::asu_mappings<> const&
        asu_mappings = *sorted_asu_proxies.asu_mappings();
      for(std::size_t i=0;i<sym.size();i++) {
        bond restraint(sites_cart, asu_mappings, sym[i]);
        result.push_back(restraint.distance_model);
      }
    }
    return result;
  }

  //! Fast computation of bond::delta given managed proxies.
  inline
  af::shared<double>
  bond_deltas(
    af::const_ref<scitbx::vec3<double> > const& sites_cart,
    bond_sorted_asu_proxies_base const& sorted_asu_proxies)
  {
    af::shared<double> result = bond_deltas(
      sites_cart, sorted_asu_proxies.simple.const_ref());
    af::const_ref<bond_asu_proxy> sym = sorted_asu_proxies.asu.const_ref();
    if (sym.size() > 0) {
      result.reserve(sorted_asu_proxies.simple.size() + sym.size());
      direct_space_asu::asu_mappings<> const&
        asu_mappings = *sorted_asu_proxies.asu_mappings();
      for(std::size_t i=0;i<sym.size();i++) {
        bond restraint(sites_cart, asu_mappings, sym[i]);
        result.push_back(restraint.delta);
      }
    }
    return result;
  }

  //! Fast computation of bond::residual() given managed proxies.
  inline
  af::shared<double>
  bond_residuals(
    af::const_ref<scitbx::vec3<double> > const& sites_cart,
    bond_sorted_asu_proxies_base const& sorted_asu_proxies)
  {
    af::shared<double> result = bond_residuals(
      sites_cart, sorted_asu_proxies.simple.const_ref());
    af::const_ref<bond_asu_proxy> sym = sorted_asu_proxies.asu.const_ref();
    if (sym.size() > 0) {
      result.reserve(sorted_asu_proxies.simple.size() + sym.size());
      direct_space_asu::asu_mappings<> const&
        asu_mappings = *sorted_asu_proxies.asu_mappings();
      for(std::size_t i=0;i<sym.size();i++) {
        bond restraint(sites_cart, asu_mappings, sym[i]);
        result.push_back(restraint.residual());
      }
    }
    return result;
  }

  namespace detail {

    inline
    double
    bond_residual_sum(
      af::const_ref<scitbx::vec3<double> > const& sites_cart,
      direct_space_asu::asu_mappings<> const& asu_mappings,
      af::const_ref<bond_asu_proxy> const& proxies,
      std::vector<bool> const& sym_active_flags,
      af::ref<scitbx::vec3<double> > const& gradient_array,
      bool disable_cache=false)
    {
      double result = 0;
      if (!disable_cache) {
        asu_cache<> cache(
          sites_cart,
          asu_mappings,
          sym_active_flags,
          gradient_array.size() != 0);
        for(std::size_t i=0;i<proxies.size();i++) {
          bond restraint(cache, proxies[i]);
          if (proxies[i].j_sym == 0) result += restraint.residual();
          else                       result += restraint.residual()*.5;
          if (gradient_array.size() != 0) {
            restraint.add_gradients(cache, proxies[i]);
          }
        }
        if (gradient_array.size() != 0) {
          cache.add_gradients(gradient_array, asu_mappings);
        }
      }
      else {
        for(std::size_t i=0;i<proxies.size();i++) {
          bond restraint(sites_cart, asu_mappings, proxies[i]);
          if (proxies[i].j_sym == 0) result += restraint.residual();
          else                       result += restraint.residual()*.5;
          if (gradient_array.size() != 0) {
            restraint.add_gradients(gradient_array, asu_mappings, proxies[i]);
          }
        }
      }
      return result;
    }

  } // namespace detail

  /*! Fast computation of sum of bond::residual() and gradients
      given managed proxies.
   */
  /*! The bond::gradients() are added to the gradient_array if
      gradient_array.size() == sites_cart.size().
      gradient_array must be initialized before this function
      is called.
      No gradient calculations are performed if gradient_array.size() == 0.

      Intermediate results are accumulated in an asu cache until
      disable_cache=true. The accumulated results are mapped back
      to the original sites at the end of the calculation. This is
      faster but requires more memory.
   */
  inline
  double
  bond_residual_sum(
    af::const_ref<scitbx::vec3<double> > const& sites_cart,
    bond_sorted_asu_proxies_base const& sorted_asu_proxies,
    af::ref<scitbx::vec3<double> > const& gradient_array,
    bool disable_cache=false)
  {
    double result = bond_residual_sum(
      sites_cart,
      sorted_asu_proxies.simple.const_ref(),
      gradient_array);
    if (sorted_asu_proxies.asu.size() > 0) {
      result += detail::bond_residual_sum(
        sites_cart,
        *sorted_asu_proxies.asu_mappings(),
        sorted_asu_proxies.asu.const_ref(),
        sorted_asu_proxies.asu_active_flags,
        gradient_array,
        disable_cache);
    }
    return result;
  }

}} // namespace cctbx::geometry_restraints

#endif // CCTBX_GEOMETRY_RESTRAINTS_BOND_H
