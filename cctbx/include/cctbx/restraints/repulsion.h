#ifndef CCTBX_RESTRAINTS_REPULSION_H
#define CCTBX_RESTRAINTS_REPULSION_H

#include <cctbx/restraints/asu_cache.h>
#include <cctbx/restraints/sorted_asu_proxies.h>

namespace cctbx { namespace restraints {

  struct repulsion_simple_proxy
  {
    repulsion_simple_proxy() {}

    repulsion_simple_proxy(
      af::tiny<unsigned, 2> const& i_seqs_,
      double vdw_radius_)
    :
      i_seqs(i_seqs_),
      vdw_radius(vdw_radius_)
    {}

    af::tiny<unsigned, 2> i_seqs;
    double vdw_radius;
  };

  struct repulsion_asu_proxy : asu_mapping_index_pair
  {
    repulsion_asu_proxy() {}

    repulsion_asu_proxy(
      asu_mapping_index_pair const& pair_,
      double vdw_radius_)
    :
      asu_mapping_index_pair(pair_),
      vdw_radius(vdw_radius_)
    {}

    // Not available in Python.
    repulsion_simple_proxy
    as_simple_proxy() const
    {
      return repulsion_simple_proxy(
        af::tiny<unsigned, 2>(i_seq, j_seq),
        vdw_radius);
    }

    double vdw_radius;
  };

  struct repulsion_function
  {
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

    // Not available in Python.
    double
    term(double vdw_radius, double delta) const
    {
      if (irexp == 1) return k_rep*vdw_radius - delta;
      return std::pow(k_rep*vdw_radius, irexp) - std::pow(delta, irexp);
    }

    // Not available in Python.
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

    // Not available in Python.
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

  class repulsion
  {
    public:
      typedef scitbx::vec3<double> vec3;

      repulsion() {}

      repulsion(
        af::tiny<scitbx::vec3<double>, 2> const& sites_,
        double vdw_radius_,
        repulsion_function const& function_=repulsion_function())
      :
        sites(sites_),
        vdw_radius(vdw_radius_),
        function(function_)
      {
        init_term();
      }

      repulsion(
        af::const_ref<scitbx::vec3<double> > const& sites_cart,
        repulsion_simple_proxy const& proxy,
        repulsion_function const& function_=repulsion_function())
      :
        vdw_radius(proxy.vdw_radius),
        function(function_)
      {
        for(int i=0;i<2;i++) {
          std::size_t i_seq = proxy.i_seqs[i];
          CCTBX_ASSERT(i_seq < sites_cart.size());
          sites[i] = sites_cart[i_seq];
        }
        init_term();
      }

      repulsion(
        af::const_ref<scitbx::vec3<double> > const& sites_cart,
        direct_space_asu::asu_mappings<> const& asu_mappings,
        repulsion_asu_proxy const& proxy,
        repulsion_function const& function_=repulsion_function())
      :
        vdw_radius(proxy.vdw_radius),
        function(function_)
      {
        sites[0] = asu_mappings.map_moved_site_to_asu(
          sites_cart[proxy.i_seq], proxy.i_seq, 0);
        sites[1] = asu_mappings.map_moved_site_to_asu(
          sites_cart[proxy.j_seq], proxy.j_seq, proxy.j_sym);
        init_term();
      }

      // Not available in Python.
      repulsion(
        asu_cache<> const& cache,
        repulsion_asu_proxy const& proxy,
        repulsion_function const& function_=repulsion_function())
      :
        vdw_radius(proxy.vdw_radius),
        function(function_)
      {
        sites[0] = cache.sites[proxy.i_seq][0];
        sites[1] = cache.sites[proxy.j_seq][proxy.j_sym];
        init_term();
      }

      double
      residual() const { return function.residual(term_); }

      // Not available in Python.
      scitbx::vec3<double>
      gradient_0() const
      {
        return diff_vec * function.gradient_factor(delta, term_);
      }

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

      // Not available in Python.
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

      // Not available in Python.
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

      af::tiny<scitbx::vec3<double>, 2> sites;
      double vdw_radius;
      repulsion_function function;
      scitbx::vec3<double> diff_vec;
      double delta;
    protected:
      double term_;

      void
      init_term()
      {
        diff_vec = sites[0] - sites[1];
        delta = diff_vec.length();
        term_ = function.term(vdw_radius, delta);
      }
  };

  inline
  af::shared<double>
  repulsion_deltas(
    af::const_ref<scitbx::vec3<double> > const& sites_cart,
    af::const_ref<repulsion_simple_proxy> const& proxies,
    repulsion_function const& function=repulsion_function())
  {
    af::shared<double> result((af::reserve(sites_cart.size())));
    for(std::size_t i=0;i<proxies.size();i++) {
      repulsion restraint(sites_cart, proxies[i], function);
      result.push_back(restraint.delta);
    }
    return result;
  }

  inline
  af::shared<double>
  repulsion_residuals(
    af::const_ref<scitbx::vec3<double> > const& sites_cart,
    af::const_ref<repulsion_simple_proxy> const& proxies,
    repulsion_function const& function=repulsion_function())
  {
    af::shared<double> result((af::reserve(sites_cart.size())));
    for(std::size_t i=0;i<proxies.size();i++) {
      repulsion restraint(sites_cart, proxies[i], function);
      result.push_back(restraint.residual());
    }
    return result;
  }

  inline
  double
  repulsion_residual_sum(
    af::const_ref<scitbx::vec3<double> > const& sites_cart,
    af::const_ref<repulsion_simple_proxy> const& proxies,
    af::ref<scitbx::vec3<double> > const& gradient_array,
    repulsion_function const& function=repulsion_function())
  {
    double result = 0;
    for(std::size_t i=0;i<proxies.size();i++) {
      repulsion restraint(sites_cart, proxies[i], function);
      result += restraint.residual();
      if (gradient_array.size() != 0) {
        restraint.add_gradients(gradient_array, proxies[i].i_seqs);
      }
    }
    return result;
  }

  inline
  af::shared<double>
  repulsion_deltas(
    af::const_ref<scitbx::vec3<double> > const& sites_cart,
    direct_space_asu::asu_mappings<> const& asu_mappings,
    af::const_ref<repulsion_asu_proxy> const& proxies,
    repulsion_function const& function=repulsion_function())
  {
    af::shared<double> result((af::reserve(sites_cart.size())));
    for(std::size_t i=0;i<proxies.size();i++) {
      repulsion restraint(sites_cart, asu_mappings, proxies[i], function);
      result.push_back(restraint.delta);
    }
    return result;
  }

  inline
  af::shared<double>
  repulsion_residuals(
    af::const_ref<scitbx::vec3<double> > const& sites_cart,
    direct_space_asu::asu_mappings<> const& asu_mappings,
    af::const_ref<repulsion_asu_proxy> const& proxies,
    repulsion_function const& function=repulsion_function())
  {
    af::shared<double> result((af::reserve(sites_cart.size())));
    for(std::size_t i=0;i<proxies.size();i++) {
      repulsion restraint(sites_cart, asu_mappings, proxies[i], function);
      result.push_back(restraint.residual());
    }
    return result;
  }

  // Not available in Python.
  inline
  double
  repulsion_residual_sum(
    af::const_ref<scitbx::vec3<double> > const& sites_cart,
    direct_space_asu::asu_mappings<> const& asu_mappings,
    af::const_ref<repulsion_asu_proxy> const& proxies,
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
        repulsion restraint(cache, proxies[i], function);
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
        repulsion restraint(sites_cart, asu_mappings, proxies[i], function);
        if (proxies[i].j_sym == 0) result += restraint.residual();
        else                       result += restraint.residual()*.5;
        if (gradient_array.size() != 0) {
          restraint.add_gradients(gradient_array, asu_mappings, proxies[i]);
        }
      }
    }
    return result;
  }

  inline
  double
  repulsion_residual_sum(
    af::const_ref<scitbx::vec3<double> > const& sites_cart,
    direct_space_asu::asu_mappings<> const& asu_mappings,
    af::const_ref<repulsion_asu_proxy> const& proxies,
    af::ref<scitbx::vec3<double> > const& gradient_array,
    repulsion_function const& function=repulsion_function(),
    bool disable_cache=false)
  {
    std::vector<bool> sym_active_flags;
    if (!disable_cache) {
      sym_active_flags.resize(asu_mappings.mappings_const_ref().size(), true);
    }
    return repulsion_residual_sum(
      sites_cart,
      asu_mappings,
      proxies,
      sym_active_flags,
      gradient_array,
      function,
      disable_cache);
  }

  typedef sorted_asu_proxies<repulsion_simple_proxy, repulsion_asu_proxy>
    repulsion_sorted_asu_proxies;

  inline
  double
  repulsion_residual_sum(
    af::const_ref<scitbx::vec3<double> > const& sites_cart,
    repulsion_sorted_asu_proxies const& sorted_asu_proxies,
    af::ref<scitbx::vec3<double> > const& gradient_array,
    repulsion_function const& function=repulsion_function(),
    bool disable_cache=false)
  {
    double result = repulsion_residual_sum(
      sites_cart,
      sorted_asu_proxies.simple.const_ref(),
      gradient_array,
      function);
    result += repulsion_residual_sum(
      sites_cart,
      *sorted_asu_proxies.asu_mappings(),
      sorted_asu_proxies.sym.const_ref(),
      sorted_asu_proxies.sym_active_flags,
      gradient_array,
      function,
      disable_cache);
    return result;
  }

}} // namespace cctbx::restraints

#endif // CCTBX_RESTRAINTS_REPULSION_H
