#ifndef CCTBX_GEOMETRY_RESTRAINTS_NONBONDED_H
#define CCTBX_GEOMETRY_RESTRAINTS_NONBONDED_H

#include <cctbx/geometry_restraints/asu_cache.h>
#include <cctbx/geometry_restraints/sorted_asu_proxies.h>

namespace cctbx { namespace geometry_restraints {

  //! Dictionary of VdW distances (element of nonbonded_distance_table).
  typedef std::map<std::string, double>
    nonbonded_distance_dict;
  //! Table of VdW distances (given two energy types).
  typedef std::map<std::string, std::map<std::string, double> >
    nonbonded_distance_table;

  //! Table of VdW radii (given one energy type).
  typedef std::map<std::string, double>
    nonbonded_radius_table;

  /*! \brief Grouping of parameters for the generation of nonbonded
      pair interactions.
   */
  struct nonbonded_params
  {
    //! Initialization.
    nonbonded_params(
      double factor_1_4_interactions_=2/3.,
      double const_shrink_1_4_interactions_=0,
      double default_distance_=0,
      double minimum_distance_=0)
    :
      factor_1_4_interactions(factor_1_4_interactions_),
      const_shrink_1_4_interactions(const_shrink_1_4_interactions_),
      default_distance(default_distance_),
      minimum_distance(minimum_distance_)
    {}

    //! Table of VdW distances (given two energy types).
    nonbonded_distance_table distance_table;
    //! Table of VdW radii (given one energy type).
    nonbonded_radius_table radius_table;
    //! Multiplicative attenuation factor for 1-4 interactions.
    double factor_1_4_interactions;
    //! Constant reduction of VdW distances for 1-4 interactions.
    double const_shrink_1_4_interactions;
    //! Default VdW distances if table lookup fails.
    /*! An exception is raised if the table lookup fails and
        default_distance == 0.
     */
    double default_distance;
    //! Global minimum VdW distance. May be zero.
    double minimum_distance;
  };

  //! Grouping of indices into array of sites (i_seqs) and vdw_distance.
  struct nonbonded_simple_proxy
  {
    //! Default constructor. Some data members are not initialized!
    nonbonded_simple_proxy() {}

    //! Constructor.
    nonbonded_simple_proxy(
      af::tiny<unsigned, 2> const& i_seqs_,
      double vdw_distance_)
    :
      i_seqs(i_seqs_),
      vdw_distance(vdw_distance_)
    {}

    //! Indices into array of sites.
    af::tiny<unsigned, 2> i_seqs;
    //! VDW distance.
    double vdw_distance;
  };

  //! Grouping of asu_mapping_index_pair and vdw_distance.
  struct nonbonded_asu_proxy : asu_mapping_index_pair
  {
    //! Default constructor. Some data members are not initialized!
    nonbonded_asu_proxy() {}

    //! Constructor.
    nonbonded_asu_proxy(
      asu_mapping_index_pair const& pair_,
      double vdw_distance_)
    :
      asu_mapping_index_pair(pair_),
      vdw_distance(vdw_distance_)
    {}

    //! Constructor.
    /*! Not available in Python.
     */
    nonbonded_simple_proxy
    as_simple_proxy() const
    {
      return nonbonded_simple_proxy(
        af::tiny<unsigned, 2>(i_seq, j_seq),
        vdw_distance);
    }

    //! VDW distance.
    double vdw_distance;
  };

  //! General repulsive function (see PROLSQ and CNS).
  /*! energy(delta) = c_rep*(max(0,(k_rep*vdw_distance)**irexp
                             -delta**irexp))**rexp
   */
  struct repulsion_function
  {
    //! Definition of coefficients.
    repulsion_function(
      double c_rep_=16,
      double k_rep_=1,
      double irexp_=1,
      double rexp_=4)
    :
      c_rep(c_rep_),
      k_rep(k_rep_),
      irexp(irexp_),
      rexp(rexp_)
    {
      CCTBX_ASSERT(rexp > 0);
    }

    //! Support for respulsion class.
    /*! Not available in Python.
     */
    double
    term(double vdw_distance, double delta) const
    {
      if (irexp == 1) return k_rep*vdw_distance - delta;
      return std::pow(k_rep*vdw_distance, irexp) - std::pow(delta, irexp);
    }

    //! Support for respulsion class.
    /*! Not available in Python.
     */
    double
    residual(double term) const
    {
      if (term <= 0) return 0;
      if (rexp == 4) {
        double term_sq = term * term;
        return c_rep * term_sq * term_sq;
      }
      return c_rep * std::pow(term, rexp);
    }

    //! Support for respulsion class.
    /*! Not available in Python.
     */
    double
    gradient_factor(double delta, double term) const
    {
      if (term <= 0 || delta == 0) return 0;
      double d_term_d_r;
      if (irexp == 1) d_term_d_r = -1;
      else            d_term_d_r = -irexp * std::pow(delta, irexp-1);
      if (rexp == 4) {
        return c_rep * rexp * term * term * term * d_term_d_r / delta;
      }
      return c_rep * rexp * std::pow(term, rexp-1) * d_term_d_r / delta;
    }

    double c_rep;
    double k_rep;
    double irexp;
    double rexp;
  };

  //! Residual and gradient calculations for nonbonded restraints.
  class nonbonded
  {
    public:
      //! Convenience typedef.
      typedef scitbx::vec3<double> vec3;

      //! Default constructor. Some data members are not initialized!
      nonbonded() {}

      //! Constructor.
      nonbonded(
        af::tiny<scitbx::vec3<double>, 2> const& sites_,
        double vdw_distance_,
        repulsion_function const& function_=repulsion_function())
      :
        sites(sites_),
        vdw_distance(vdw_distance_),
        function(function_)
      {
        init_term();
      }

      /*! \brief Coordinates are copied from sites_cart according to
          proxy.i_seqs, parameters are copied from proxy.
       */
      nonbonded(
        af::const_ref<scitbx::vec3<double> > const& sites_cart,
        nonbonded_simple_proxy const& proxy,
        repulsion_function const& function_=repulsion_function())
      :
        vdw_distance(proxy.vdw_distance),
        function(function_)
      {
        for(int i=0;i<2;i++) {
          std::size_t i_seq = proxy.i_seqs[i];
          CCTBX_ASSERT(i_seq < sites_cart.size());
          sites[i] = sites_cart[i_seq];
        }
        init_term();
      }

      /*! \brief Coordinates are copied from sites_cart according to
          proxy.i_seq, proxy.j_seq, parameters are copied from proxy.
       */
      nonbonded(
        af::const_ref<scitbx::vec3<double> > const& sites_cart,
        direct_space_asu::asu_mappings<> const& asu_mappings,
        nonbonded_asu_proxy const& proxy,
        repulsion_function const& function_=repulsion_function())
      :
        vdw_distance(proxy.vdw_distance),
        function(function_)
      {
        sites[0] = asu_mappings.map_moved_site_to_asu(
          sites_cart[proxy.i_seq], proxy.i_seq, 0);
        sites[1] = asu_mappings.map_moved_site_to_asu(
          sites_cart[proxy.j_seq], proxy.j_seq, proxy.j_sym);
        init_term();
      }

      //! For fast processing. Not available in Python.
      nonbonded(
        asu_cache<> const& cache,
        nonbonded_asu_proxy const& proxy,
        repulsion_function const& function_=repulsion_function())
      :
        vdw_distance(proxy.vdw_distance),
        function(function_)
      {
        sites[0] = cache.sites[proxy.i_seq][0];
        sites[1] = cache.sites[proxy.j_seq][proxy.j_sym];
        init_term();
      }

      //! Uses function.residual(function.ter(...)).
      double
      residual() const { return function.residual(term_); }

      //! Gradient with respect to sites[0].
      /*! Not available in Python.
       */
      scitbx::vec3<double>
      gradient_0() const
      {
        return diff_vec * function.gradient_factor(delta, term_);
      }

      //! Gradients with respect to both sites.
      af::tiny<scitbx::vec3<double>, 2>
      gradients() const
      {
        af::tiny<scitbx::vec3<double>, 2> result;
        result[0] = gradient_0();
        result[1] = -result[0];
        return result;
      }

      // Not available in Python.
      void
      add_gradients(
        af::ref<scitbx::vec3<double> > const& gradient_array,
        af::tiny<unsigned, 2> const& i_seqs) const
      {
        vec3 g0 = gradient_0();
        gradient_array[i_seqs[0]] += g0;
        gradient_array[i_seqs[1]] += -g0;
      }

      //! Support for nonbonded_residual_sum.
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

      //! Support for nonbonded_residual_sum.
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

      //! Cartesian coordinates of nonbonded sites.
      af::tiny<scitbx::vec3<double>, 2> sites;
      //! Parameter (usually as passed to the constructor).
      double vdw_distance;
      //! Function (usually as passed to the constructor).
      repulsion_function function;
      //! Difference vector sites[0] - sites[1].
      scitbx::vec3<double> diff_vec;
      //! Length of diff_vec.
      double delta;
    protected:
      double term_;

      void
      init_term()
      {
        diff_vec = sites[0] - sites[1];
        delta = diff_vec.length();
        term_ = function.term(vdw_distance, delta);
      }
  };

  //! Fast computation of nonbonded::delta given an array of nonbonded proxies.
  inline
  af::shared<double>
  nonbonded_deltas(
    af::const_ref<scitbx::vec3<double> > const& sites_cart,
    af::const_ref<nonbonded_simple_proxy> const& proxies,
    repulsion_function const& function=repulsion_function())
  {
    af::shared<double> result((af::reserve(proxies.size())));
    for(std::size_t i=0;i<proxies.size();i++) {
      nonbonded restraint(sites_cart, proxies[i], function);
      result.push_back(restraint.delta);
    }
    return result;
  }

  /*! \brief Fast computation of nonbonded::residual() given an array of
      nonbonded proxies.
   */
  inline
  af::shared<double>
  nonbonded_residuals(
    af::const_ref<scitbx::vec3<double> > const& sites_cart,
    af::const_ref<nonbonded_simple_proxy> const& proxies,
    repulsion_function const& function=repulsion_function())
  {
    af::shared<double> result((af::reserve(proxies.size())));
    for(std::size_t i=0;i<proxies.size();i++) {
      nonbonded restraint(sites_cart, proxies[i], function);
      result.push_back(restraint.residual());
    }
    return result;
  }

  /*! Fast computation of sum of nonbonded::residual() and gradients
      given an array of nonbonded proxies.
   */
  /*! The nonbonded::gradients() are added to the gradient_array if
      gradient_array.size() == sites_cart.size().
      gradient_array must be initialized before this function
      is called.
      No gradient calculations are performed if gradient_array.size() == 0.
   */
  inline
  double
  nonbonded_residual_sum(
    af::const_ref<scitbx::vec3<double> > const& sites_cart,
    af::const_ref<nonbonded_simple_proxy> const& proxies,
    af::ref<scitbx::vec3<double> > const& gradient_array,
    repulsion_function const& function=repulsion_function())
  {
    double result = 0;
    for(std::size_t i=0;i<proxies.size();i++) {
      nonbonded restraint(sites_cart, proxies[i], function);
      result += restraint.residual();
      if (gradient_array.size() != 0) {
        restraint.add_gradients(gradient_array, proxies[i].i_seqs);
      }
    }
    return result;
  }

  //! Managed group of nonbonded_simple_proxy and nonbonded_asu_proxy arrays.
  typedef sorted_asu_proxies<nonbonded_simple_proxy, nonbonded_asu_proxy>
    nonbonded_sorted_asu_proxies;

  //! Fast computation of nonbonded::delta given managed proxies.
  inline
  af::shared<double>
  nonbonded_deltas(
    af::const_ref<scitbx::vec3<double> > const& sites_cart,
    nonbonded_sorted_asu_proxies const& sorted_asu_proxies,
    repulsion_function const& function=repulsion_function())
  {
    af::shared<double> result = nonbonded_deltas(
      sites_cart, sorted_asu_proxies.simple.const_ref());
    af::const_ref<nonbonded_asu_proxy> sym = sorted_asu_proxies.asu.const_ref();
    if (sym.size() > 0) {
      result.reserve(sorted_asu_proxies.simple.size() + sym.size());
      direct_space_asu::asu_mappings<> const&
        asu_mappings = *sorted_asu_proxies.asu_mappings();
      for(std::size_t i=0;i<sym.size();i++) {
        nonbonded restraint(sites_cart, asu_mappings, sym[i], function);
        result.push_back(restraint.delta);
      }
    }
    return result;
  }

  //! Fast computation of nonbonded::residual() given managed proxies.
  inline
  af::shared<double>
  nonbonded_residuals(
    af::const_ref<scitbx::vec3<double> > const& sites_cart,
    nonbonded_sorted_asu_proxies const& sorted_asu_proxies,
    repulsion_function const& function=repulsion_function())
  {
    af::shared<double> result = nonbonded_residuals(
      sites_cart, sorted_asu_proxies.simple.const_ref());
    af::const_ref<nonbonded_asu_proxy> sym = sorted_asu_proxies.asu.const_ref();
    if (sym.size() > 0) {
      result.reserve(sorted_asu_proxies.simple.size() + sym.size());
      direct_space_asu::asu_mappings<> const&
        asu_mappings = *sorted_asu_proxies.asu_mappings();
      for(std::size_t i=0;i<sym.size();i++) {
        nonbonded restraint(sites_cart, asu_mappings, sym[i], function);
        result.push_back(restraint.residual());
      }
    }
    return result;
  }

  namespace detail {

    inline
    double
    nonbonded_residual_sum(
      af::const_ref<scitbx::vec3<double> > const& sites_cart,
      direct_space_asu::asu_mappings<> const& asu_mappings,
      af::const_ref<nonbonded_asu_proxy> const& proxies,
      std::vector<bool> const& sym_active_flags,
      af::ref<scitbx::vec3<double> > const& gradient_array,
      repulsion_function const& function=repulsion_function(),
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
          nonbonded restraint(cache, proxies[i], function);
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
          nonbonded restraint(sites_cart, asu_mappings, proxies[i], function);
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

  /*! Fast computation of sum of nonbonded::residual() and gradients
      given managed proxies.
   */
  /*! The nonbonded::gradients() are added to the gradient_array if
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
  nonbonded_residual_sum(
    af::const_ref<scitbx::vec3<double> > const& sites_cart,
    nonbonded_sorted_asu_proxies const& sorted_asu_proxies,
    af::ref<scitbx::vec3<double> > const& gradient_array,
    repulsion_function const& function=repulsion_function(),
    bool disable_cache=false)
  {
    double result = nonbonded_residual_sum(
      sites_cart,
      sorted_asu_proxies.simple.const_ref(),
      gradient_array,
      function);
    if (sorted_asu_proxies.asu.size() > 0) {
      result += detail::nonbonded_residual_sum(
        sites_cart,
        *sorted_asu_proxies.asu_mappings(),
        sorted_asu_proxies.asu.const_ref(),
        sorted_asu_proxies.asu_active_flags,
        gradient_array,
        function,
        disable_cache);
    }
    return result;
  }

}} // namespace cctbx::geometry_restraints

#endif // CCTBX_GEOMETRY_RESTRAINTS_NONBONDED_H
