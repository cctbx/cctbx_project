#ifndef CCTBX_GEOMETRY_RESTRAINTS_NONBONDED_H
#define CCTBX_GEOMETRY_RESTRAINTS_NONBONDED_H

#include <cctbx/geometry_restraints/asu_cache.h>
#include <cctbx/geometry_restraints/sorted_asu_proxies.h>
#include <tbxx/optional_copy.hpp>
#include <set>

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
      double minimum_distance_=0,
      double const_shrink_donor_acceptor_=0)
    :
      factor_1_4_interactions(factor_1_4_interactions_),
      const_shrink_1_4_interactions(const_shrink_1_4_interactions_),
      default_distance(default_distance_),
      minimum_distance(minimum_distance_),
      const_shrink_donor_acceptor(const_shrink_donor_acceptor_)
    {}

    //! Find largest possible VdW distance given nonbonded_types.
    double
    find_max_vdw_distance(
      af::const_ref<std::string> const& nonbonded_types) const
    {
      double result = -1;
      std::set<std::string> unique_types(
        nonbonded_types.begin(), nonbonded_types.end());
      for(std::set<std::string>::const_iterator
            unique_types_i =  unique_types.begin();
            unique_types_i != unique_types.end();
            unique_types_i++) {
        for(std::set<std::string>::const_iterator
              unique_types_j =  unique_types_i;
              unique_types_j != unique_types.end();
              unique_types_j++) {
          double distance = get_nonbonded_distance(
            *unique_types_i,
            *unique_types_j,
            /*donor_acceptor_adjust*/ false);
          if (distance < 0) distance = default_distance;
          if (result < distance) result = distance;
        }
      }
      return std::max(minimum_distance, result);
    }

    //! Determine VdW distance by lookup in distance_table and radius_table.
    /*! Not available in Python.
     */
    double
    get_nonbonded_distance(
      std::string const& type_i,
      std::string const& type_j,
      bool donor_acceptor_adjust,
      int charge_i=0,
      int charge_j=0) const
    {
      nonbonded_distance_table::const_iterator
        distance_dict = distance_table.find(type_i);
      double return_vdw;
      if (distance_dict != distance_table.end()) {
        nonbonded_distance_dict::const_iterator
          dict_entry = distance_dict->second.find(type_j);
        if (dict_entry != distance_dict->second.end()) {
          return dict_entry->second;
        }
      }
      distance_dict = distance_table.find(type_j);
      if (distance_dict != distance_table.end()) {
        nonbonded_distance_dict::const_iterator
          dict_entry = distance_dict->second.find(type_i);
        if (dict_entry != distance_dict->second.end()) {
          return dict_entry->second;
        }
      }
      geometry_restraints::nonbonded_radius_table::const_iterator radius_i;
      if (charge_i != 0) {
        radius_i = ionic_radius_table.find(type_i);
      }
      if ((charge_i == 0) || (radius_i == ionic_radius_table.end())) {
        radius_i = radius_table.find(type_i);
      }
      if (radius_i != radius_table.end()) {
        geometry_restraints::nonbonded_radius_table::const_iterator radius_j;
        if (charge_j != 0) {
          radius_j = ionic_radius_table.find(type_j);
        }
        if ((charge_j == 0) || (radius_j == ionic_radius_table.end())) {
          radius_j = radius_table.find(type_j);
        }
        if (radius_j != radius_table.end()) {
          //return radius_i->second + radius_j->second;
          return_vdw = radius_i->second + radius_j->second;
          // development code for adjusting nonbonded radii
          // for H-bonding atoms - off pending further testing
          if (donor_acceptor_adjust == true) {
            geometry_restraints::nonbonded_radius_table::const_iterator
              h_bond_state_i = donor_acceptor_table.find(type_i);
            if (h_bond_state_i != donor_acceptor_table.end()) {
              geometry_restraints::nonbonded_radius_table::const_iterator
                h_bond_state_j = donor_acceptor_table.find(type_j);
              if (h_bond_state_j != donor_acceptor_table.end()) {
                if ( (h_bond_state_i->second == 1 &&
                      h_bond_state_j->second == 2) || //donor & acceptor
                     (h_bond_state_i->second == 2 &&
                      h_bond_state_j->second == 1) || //acceptor & donor
                     (h_bond_state_i->second == 1 &&
                      h_bond_state_j->second == 3) || //donor & both
                     (h_bond_state_i->second == 3 &&
                      h_bond_state_j->second == 1) || //both & donor
                     (h_bond_state_i->second == 2 &&
                      h_bond_state_j->second == 3) || //acceptor & both
                     (h_bond_state_i->second == 3 &&
                      h_bond_state_j->second == 2) || //both & acceptor
                     (h_bond_state_i->second == 3 &&
                      h_bond_state_j->second == 3) || //both & both
                     (h_bond_state_i->second == 2 &&
                      h_bond_state_j->second == 4) || //acceptor & hydrogen
                     (h_bond_state_i->second == 4 &&
                      h_bond_state_j->second == 2) || //hydrogen & acceptor
                     (h_bond_state_i->second == 3 &&
                      h_bond_state_j->second == 4) || //both & hydrogen
                     (h_bond_state_i->second == 4 &&
                      h_bond_state_j->second == 3) ){ //hydrogen & both
                    //subtract const_shrink_donor_acceptor for H-bond
                   return_vdw -= const_shrink_donor_acceptor;
                 }
              }
            }
          }
          return std::max(minimum_distance, return_vdw);
        }
      }
      return -1;
    }

    /*! \brief Adjusts distance considering factor_1_4_interactions,
        const_shrink_1_4_interactions, and minimum_distance.
     */
    /*! Not available in Python.
     */
    double
    adjust_nonbonded_distance(
      double distance,
      bool is_1_4_interaction) const
    {
      if (is_1_4_interaction) {
        distance *= factor_1_4_interactions;
        distance -= const_shrink_1_4_interactions;
      }
      return std::max(minimum_distance, distance);
    }

    //! Table of VdW distances (given two energy types).
    nonbonded_distance_table distance_table;
    //! Table of VdW radii (given one energy type).
    nonbonded_radius_table radius_table;
    //! Table of VdW radii for ions
    nonbonded_radius_table ionic_radius_table;
    //! Table of Donor/Acceptor status given one energy type
    nonbonded_radius_table donor_acceptor_table; //re-use table type
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
    // subtraction constant for H-bond vdW radii
    double const_shrink_donor_acceptor;
  };

  //! Grouping of indices into array of sites (i_seqs) and vdw_distance.
  struct nonbonded_simple_proxy
  {
    typedef af::tiny<unsigned, 2> i_seqs_type;

    //! Default constructor. Some data members are not initialized!
    nonbonded_simple_proxy() {}

    //! Constructor.
    nonbonded_simple_proxy(
      i_seqs_type const& i_seqs_,
      double vdw_distance_)
    :
      i_seqs(i_seqs_),
      vdw_distance(vdw_distance_)
    {}

    //! Constructor.
    nonbonded_simple_proxy(
      i_seqs_type const& i_seqs_,
      sgtbx::rt_mx const& rt_mx_ji_,
      double vdw_distance_)
    :
      i_seqs(i_seqs_),
      rt_mx_ji(rt_mx_ji_),
      vdw_distance(vdw_distance_)
    {}

    //! Indices into array of sites.
    i_seqs_type i_seqs;

    //! Optional symmetry operation operating on i_seqs[1]
    tbxx::optional_copy<sgtbx::rt_mx> rt_mx_ji;

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
      vdw_distance(vdw_distance_),
      init_pair(pair_)
    {}

    //! Constructor.
    /*! Not available in Python.
     */
    nonbonded_simple_proxy
    as_simple_proxy() const
    {
      return nonbonded_simple_proxy(
        nonbonded_simple_proxy::i_seqs_type(i_seq, j_seq),
        vdw_distance);
    }

    //! VDW distance.
    double vdw_distance;
    asu_mapping_index_pair init_pair; // used for pickling in getinitargs()
  };

  //! General repulsive function (see PROLSQ and CNS).
  /*! energy(delta) = c_rep*(max(0,(k_rep*vdw_distance)**irexp
                             -delta**irexp))**rexp
   */
  struct prolsq_repulsion_function
  {
    //! Definition of coefficients.
    prolsq_repulsion_function(
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

    //! Support for nonbonded class.
    /*! Not available in Python.
     */
    double
    term(double vdw_distance, double delta) const
    {
      if (irexp == 1) return k_rep*vdw_distance - delta;
      return std::pow(k_rep*vdw_distance, irexp) - std::pow(delta, irexp);
    }

    //! Support for nonbonded class.
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

    //! Residual.
    double
    residual(double vdw_distance, double delta) const
    {
      return residual(term(vdw_distance, delta));
    }

    //! Support for nonbonded class.
    /*! Not available in Python.
     */
    double
    d_residual_d_delta_over_delta(
      double delta, double /* vdw_distance */, double term) const
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

  //! energy(delta) = k_rep*vdw_distance/(delta**irexp) repulsive function.
  struct inverse_power_repulsion_function
  {
    //! Definition of coefficients.
    inverse_power_repulsion_function(
      double nonbonded_distance_cutoff_,
      double k_rep_=1,
      double irexp_=1)
    :
      nonbonded_distance_cutoff(nonbonded_distance_cutoff_),
      k_rep(k_rep_),
      irexp(irexp_)
    {}

    //! Support for nonbonded class.
    /*! Not available in Python.
     */
    double
    term(double vdw_distance, double delta) const
    {
      CCTBX_ASSERT(delta != 0);
      if (delta >= nonbonded_distance_cutoff) return 0;
      if (irexp == 1) return k_rep*vdw_distance / delta;
      if (irexp == 2) return k_rep*vdw_distance / delta / delta;
      return k_rep*vdw_distance / std::pow(delta, irexp);
    }

    //! Support for nonbonded class.
    /*! Not available in Python.
     */
    double
    residual(double term) const { return term; }

    //! Residual.
    double
    residual(double vdw_distance, double delta) const
    {
      return residual(term(vdw_distance, delta));
    }

    //! Support for nonbonded class.
    /*! Not available in Python.
     */
    double
    d_residual_d_delta_over_delta(
      double delta, double /* vdw_distance */, double term) const
    {
      if (term == 0) return 0;
      return -irexp * term / delta / delta;
    }

    double nonbonded_distance_cutoff;
    double k_rep;
    double irexp;
  };

  //! energy(delta) = max_residual*((cos(pi*delta)+1)/2)**exponent
  struct cos_repulsion_function
  {
    //! Definition of coefficients.
    cos_repulsion_function(
      double max_residual_,
      double exponent_=1)
    :
      max_residual(max_residual_),
      exponent(exponent_)
    {}

    //! Support for nonbonded class.
    /*! Not available in Python.
     */
    double
    term(double vdw_distance, double delta) const
    {
      if (delta >= vdw_distance) return 0;
      using scitbx::constants::pi;
      double w2 = (std::cos(pi * delta / vdw_distance) + 1) / 2;
      if (exponent == 1) return max_residual * w2;
      if (exponent == 2) return max_residual * w2 * w2;
      return max_residual * std::pow(w2, exponent);
    }

    //! Support for nonbonded class.
    /*! Not available in Python.
     */
    double
    residual(double term) const { return term; }

    //! Residual.
    double
    residual(double vdw_distance, double delta) const
    {
      return residual(term(vdw_distance, delta));
    }

    //! Support for nonbonded class.
    /*! Not available in Python.
     */
    double
    d_residual_d_delta_over_delta(
      double delta, double vdw_distance, double /* term */) const
    {
      if (delta == 0 || delta >= vdw_distance) return 0;
      using scitbx::constants::pi;
      double a = pi * delta / vdw_distance;
      double w = std::cos(a) + 1;
      if (exponent == 1) {
        return -(exponent*max_residual*pi*std::sin(a))
             / (2 * vdw_distance * delta);
      }
      if (exponent == 2) {
        return -(exponent*max_residual*pi*w*std::sin(a))
             / (4 * vdw_distance * delta);
      }
      return -(exponent*max_residual*pi*std::pow(w,exponent-1)*std::sin(a))
           / (std::pow(2.0,exponent) * vdw_distance * delta);
    }

    double max_residual;
    double exponent;
  };

  //! energy(delta) = max_residual*exp(-delta**2/f_sq)
  /*! f_sq = -vdw_distance**2 / log(norm_height_at_vdw_distance)
   */
  struct gaussian_repulsion_function
  {
    //! Definition of coefficients.
    gaussian_repulsion_function(
      double max_residual_,
      double norm_height_at_vdw_distance=0.1)
    :
      max_residual(max_residual_)
    {
      CCTBX_ASSERT(norm_height_at_vdw_distance < 1);
      CCTBX_ASSERT(norm_height_at_vdw_distance > 0);
      log_norm_height_at_vdw_distance = std::log(norm_height_at_vdw_distance);
      CCTBX_ASSERT(log_norm_height_at_vdw_distance < 0);
    }

    double
    norm_height_at_vdw_distance() const
    {
      return std::exp(log_norm_height_at_vdw_distance);
    }

    //! Support for nonbonded class.
    /*! Not available in Python.
     */
    double
    term(double vdw_distance, double delta) const
    {
      double minus_f_sq = (vdw_distance * vdw_distance)
                        / log_norm_height_at_vdw_distance;
      CCTBX_ASSERT(minus_f_sq != 0);
      return max_residual * std::exp(delta*delta / minus_f_sq);
    }

    //! Support for nonbonded class.
    /*! Not available in Python.
     */
    double
    residual(double term) const { return term; }

    //! Residual.
    double
    residual(double vdw_distance, double delta) const
    {
      return residual(term(vdw_distance, delta));
    }

    //! Support for nonbonded class.
    /*! Not available in Python.
     */
    double
    d_residual_d_delta_over_delta(
      double /* delta */, double vdw_distance, double term) const
    {
      double minus_f_sq = (vdw_distance * vdw_distance)
                        / log_norm_height_at_vdw_distance;
      CCTBX_ASSERT(minus_f_sq != 0);
      return 2 * term / minus_f_sq;
    }

    double max_residual;
    double log_norm_height_at_vdw_distance;
  };

  //! Residual and gradient calculations for nonbonded restraints.
  template <typename NonbondedFunction>
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
        NonbondedFunction const& function_)
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
        NonbondedFunction const& function_)
      :
        vdw_distance(proxy.vdw_distance),
        function(function_)
      {
        CCTBX_ASSERT(!proxy.rt_mx_ji);
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
        NonbondedFunction const& function_)
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
        NonbondedFunction const& function_)
      :
        vdw_distance(proxy.vdw_distance),
        function(function_)
      {
        sites[0] = cache.sites[proxy.i_seq][0];
        sites[1] = cache.sites[proxy.j_seq][proxy.j_sym];
        init_term();
      }

      //! Uses function.residual(function.term(...)).
      double
      residual() const { return function.residual(term_); }

      //! Gradient with respect to sites[0].
      /*! Not available in Python.
       */
      scitbx::vec3<double>
      gradient_0() const
      {
        return diff_vec
             * function.d_residual_d_delta_over_delta(
                 delta, vdw_distance, term_);
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
        nonbonded_simple_proxy::i_seqs_type const& i_seqs) const
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
      NonbondedFunction function;
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
    af::const_ref<nonbonded_simple_proxy> const& proxies)
  {
    af::shared<double> result((af::reserve(proxies.size())));
    prolsq_repulsion_function function; // not actually used
    for(std::size_t i=0;i<proxies.size();i++) {
      nonbonded<prolsq_repulsion_function> restraint(
        sites_cart, proxies[i], function);
      result.push_back(restraint.delta);
    }
    return result;
  }

  /*! \brief Fast computation of nonbonded::residual() given an array of
      nonbonded proxies.
   */
  template <typename NonbondedFunction>
  af::shared<double>
  nonbonded_residuals(
    af::const_ref<scitbx::vec3<double> > const& sites_cart,
    af::const_ref<nonbonded_simple_proxy> const& proxies,
    NonbondedFunction const& function)
  {
    af::shared<double> result((af::reserve(proxies.size())));
    for(std::size_t i=0;i<proxies.size();i++) {
      nonbonded<NonbondedFunction> restraint(sites_cart, proxies[i], function);
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
  template <typename NonbondedFunction>
  double
  nonbonded_residual_sum(
    af::const_ref<scitbx::vec3<double> > const& sites_cart,
    af::const_ref<nonbonded_simple_proxy> const& proxies,
    af::ref<scitbx::vec3<double> > const& gradient_array,
    NonbondedFunction const& function)
  {
    double result = 0;
    for(std::size_t i=0;i<proxies.size();i++) {
      nonbonded<NonbondedFunction> restraint(sites_cart, proxies[i], function);
      result += restraint.residual();
      if (gradient_array.size() != 0) {
        restraint.add_gradients(gradient_array, proxies[i].i_seqs);
      }
    }
    return result;
  }

  //! Managed group of nonbonded_simple_proxy and nonbonded_asu_proxy arrays.
  typedef sorted_asu_proxies<nonbonded_simple_proxy, nonbonded_asu_proxy>
    nonbonded_sorted_asu_proxies_base;

  //! Fast computation of nonbonded::delta given managed proxies.
  inline
  af::shared<double>
  nonbonded_deltas(
    af::const_ref<scitbx::vec3<double> > const& sites_cart,
    nonbonded_sorted_asu_proxies_base const& sorted_asu_proxies)
  {
    af::shared<double> result = nonbonded_deltas(
      sites_cart, sorted_asu_proxies.simple.const_ref());
    af::const_ref<nonbonded_asu_proxy>
      sym = sorted_asu_proxies.asu.const_ref();
    if (sym.size() > 0) {
      result.reserve(sorted_asu_proxies.simple.size() + sym.size());
      direct_space_asu::asu_mappings<> const&
        asu_mappings = *sorted_asu_proxies.asu_mappings();
      prolsq_repulsion_function function; // not actually used
      for(std::size_t i=0;i<sym.size();i++) {
        nonbonded<prolsq_repulsion_function> restraint(
          sites_cart, asu_mappings, sym[i], function);
        result.push_back(restraint.delta);
      }
    }
    return result;
  }

  //! Fast computation of nonbonded::residual() given managed proxies.
  template <typename NonbondedFunction>
  af::shared<double>
  nonbonded_residuals(
    af::const_ref<scitbx::vec3<double> > const& sites_cart,
    nonbonded_sorted_asu_proxies_base const& sorted_asu_proxies,
    NonbondedFunction const& function)
  {
    af::shared<double> result = nonbonded_residuals(
      sites_cart, sorted_asu_proxies.simple.const_ref(), function);
    af::const_ref<nonbonded_asu_proxy>
      sym = sorted_asu_proxies.asu.const_ref();
    if (sym.size() > 0) {
      result.reserve(sorted_asu_proxies.simple.size() + sym.size());
      direct_space_asu::asu_mappings<> const&
        asu_mappings = *sorted_asu_proxies.asu_mappings();
      for(std::size_t i=0;i<sym.size();i++) {
        nonbonded<NonbondedFunction> restraint(
          sites_cart, asu_mappings, sym[i], function);
        result.push_back(restraint.residual());
      }
    }
    return result;
  }

  namespace detail {

    template <typename NonbondedFunction>
    double
    nonbonded_residual_sum(
      af::const_ref<scitbx::vec3<double> > const& sites_cart,
      direct_space_asu::asu_mappings<> const& asu_mappings,
      af::const_ref<nonbonded_asu_proxy> const& proxies,
      std::vector<bool> const& sym_active_flags,
      af::ref<scitbx::vec3<double> > const& gradient_array,
      NonbondedFunction const& function,
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
          nonbonded<NonbondedFunction> restraint(cache, proxies[i], function);
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
          nonbonded<NonbondedFunction> restraint(
            sites_cart, asu_mappings, proxies[i], function);
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
  template <typename NonbondedFunction>
  double
  nonbonded_residual_sum(
    af::const_ref<scitbx::vec3<double> > const& sites_cart,
    nonbonded_sorted_asu_proxies_base const& sorted_asu_proxies,
    af::ref<scitbx::vec3<double> > const& gradient_array,
    NonbondedFunction const& function,
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
