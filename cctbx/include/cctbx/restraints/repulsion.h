#ifndef CCTBX_RESTRAINTS_REPULSION_H
#define CCTBX_RESTRAINTS_REPULSION_H

#include <cctbx/restraints/asu_cache.h>

namespace cctbx { namespace restraints {

  struct repulsion_sym_proxy
  {
    repulsion_sym_proxy() {}

    repulsion_sym_proxy(
      direct_space_asu::asu_mapping_index_pair const& pair_,
      double vdw_radius_)
    :
      pair(pair_),
      vdw_radius(vdw_radius_)
    {}

    direct_space_asu::asu_mapping_index_pair pair;
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
        af::tiny<scitbx::vec3<double>, 2> const& sites,
        double vdw_radius,
        repulsion_function const& function_=repulsion_function())
      :
        function(function_)
      {
        diff_vec = sites[0] - sites[1];
        delta = diff_vec.length();
        term_ = function.term(vdw_radius, delta);
      }

      // Not available in Python.
      repulsion(
        af::const_ref<scitbx::vec3<double> > const& sites_cart,
        direct_space_asu::asu_mappings<> const& asu_mappings,
        repulsion_sym_proxy const& proxy,
        repulsion_function const& function_=repulsion_function())
      :
        function(function_)
      {
        vec3 mapped_site_0 = asu_mappings.map_moved_site_to_asu(
          sites_cart[proxy.pair.i_seq], proxy.pair.i_seq, 0);
        vec3 mapped_site_1 = asu_mappings.map_moved_site_to_asu(
          sites_cart[proxy.pair.j_seq], proxy.pair.j_seq, proxy.pair.j_sym);
        diff_vec = mapped_site_0 - mapped_site_1;
        delta = diff_vec.length();
        term_ = function.term(proxy.vdw_radius, delta);
      }

      // Not available in Python.
      repulsion(
        af::const_ref<scitbx::vec3<double> > const& sites_cart,
        asu_cache<> const& cache,
        repulsion_sym_proxy const& proxy,
        repulsion_function const& function_=repulsion_function())
      :
        function(function_)
      {
        vec3 mapped_site_0 = cache.sites[proxy.pair.i_seq][0];
        vec3 mapped_site_1 = cache.sites[proxy.pair.j_seq][proxy.pair.j_sym];
        diff_vec = mapped_site_0 - mapped_site_1;
        delta = diff_vec.length();
        term_ = function.term(proxy.vdw_radius, delta);
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
        direct_space_asu::asu_mappings<> const& asu_mappings,
        direct_space_asu::asu_mapping_index_pair const& pair) const
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
        direct_space_asu::asu_mapping_index_pair const& pair) const
      {
        vec3 grad_asu = gradient_0();
        cache.gradients[pair.i_seq] += grad_asu;
        if (pair.j_sym == 0) {
          cache.gradients[pair.j_seq] -= grad_asu;
        }
      }

      repulsion_function function;
      scitbx::vec3<double> diff_vec;
      double delta;
    protected:
      double term_;
  };

  inline
  af::shared<double>
  repulsion_deltas(
    af::const_ref<scitbx::vec3<double> > const& sites_cart,
    direct_space_asu::asu_mappings<> const& asu_mappings,
    af::const_ref<repulsion_sym_proxy> const& proxies,
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
    af::const_ref<repulsion_sym_proxy> const& proxies,
    repulsion_function const& function=repulsion_function())
  {
    af::shared<double> result((af::reserve(sites_cart.size())));
    for(std::size_t i=0;i<proxies.size();i++) {
      repulsion restraint(sites_cart, asu_mappings, proxies[i], function);
      result.push_back(restraint.residual());
    }
    return result;
  }

  inline
  double
  repulsion_residual_sum(
    af::const_ref<scitbx::vec3<double> > const& sites_cart,
    direct_space_asu::asu_mappings<> const& asu_mappings,
    af::const_ref<repulsion_sym_proxy> const& proxies,
    af::ref<scitbx::vec3<double> > const& gradient_array,
    repulsion_function const& function=repulsion_function(),
    bool disable_cache=false)
  {
    double result = 0;
    if (!disable_cache) {
      asu_cache<> cache(
        sites_cart, asu_mappings, gradient_array.size() != 0);
      for(std::size_t i=0;i<proxies.size();i++) {
        repulsion restraint(sites_cart, cache, proxies[i], function);
        if (proxies[i].pair.j_sym == 0) result += restraint.residual();
        else                            result += restraint.residual()*.5;
        if (gradient_array.size() != 0) {
          restraint.add_gradients(cache, proxies[i].pair);
        }
      }
      if (gradient_array.size() != 0) {
        cache.add_gradients(gradient_array, asu_mappings);
      }
    }
    else {
      for(std::size_t i=0;i<proxies.size();i++) {
        repulsion restraint(sites_cart, asu_mappings, proxies[i], function);
        if (proxies[i].pair.j_sym == 0) result += restraint.residual();
        else                            result += restraint.residual()*.5;
        if (gradient_array.size() != 0) {
          restraint.add_gradients(gradient_array,asu_mappings,proxies[i].pair);
        }
      }
    }
    return result;
  }

}} // namespace cctbx::restraints

#endif // CCTBX_RESTRAINTS_REPULSION_H
