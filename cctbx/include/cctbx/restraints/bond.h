#ifndef CCTBX_RESTRAINTS_BOND_H
#define CCTBX_RESTRAINTS_BOND_H

#include <cctbx/restraints/utils.h>
#include <cctbx/restraints/asu_cache.h>
#include <cctbx/restraints/sorted_asu_proxies.h>
#include <set>

namespace cctbx { namespace restraints {

  struct bond_params
  {
    bond_params() {}

    bond_params(double distance_ideal_, double weight_)
    :
      distance_ideal(distance_ideal_), weight(weight_)
    {}

    double distance_ideal;
    double weight;
  };

  struct bond_simple_proxy : bond_params
  {
    bond_simple_proxy() {}

    bond_simple_proxy(
      af::tiny<std::size_t, 2> const& i_seqs_,
      double distance_ideal_,
      double weight_)
    :
      bond_params(distance_ideal_, weight_),
      i_seqs(i_seqs_)
    {}

    af::tiny<std::size_t, 2> i_seqs;
  };

  struct bond_asu_proxy : bond_params, asu_mapping_index_pair
  {
    bond_asu_proxy() {}

    bond_asu_proxy(
      asu_mapping_index_pair const& pair_,
      double distance_ideal_,
      double weight_)
    :
      bond_params(distance_ideal_, weight_),
      asu_mapping_index_pair(pair_)
    {}

    bond_simple_proxy
    as_simple_proxy() const
    {
      return bond_simple_proxy(
        af::tiny<std::size_t, 2>(i_seq, j_seq),
        distance_ideal,
        weight);
    }
#if defined(__MACH__) && defined(__APPLE_CC__) && __APPLE_CC__ <= 1640
      bool dummy_;
#endif
  };

  class bond : public bond_params
  {
    public:
      typedef scitbx::vec3<double> vec3;

      bond() {}

      bond(
        af::tiny<scitbx::vec3<double>, 2> const& sites_,
        double distance_ideal_,
        double weight_)
      :
        bond_params(distance_ideal_, weight_),
        sites(sites_)
      {
        init_distance_model();
      }

      bond(
        af::const_ref<scitbx::vec3<double> > const& sites_cart,
        bond_simple_proxy const& proxy)
      :
        bond_params(proxy.distance_ideal, proxy.weight)
      {
        for(int i=0;i<2;i++) {
          std::size_t i_seq = proxy.i_seqs[i];
          CCTBX_ASSERT(i_seq < sites_cart.size());
          sites[i] = sites_cart[i_seq];
        }
        init_distance_model();
      }

      bond(
        af::const_ref<scitbx::vec3<double> > const& sites_cart,
        direct_space_asu::asu_mappings<> const& asu_mappings,
        bond_asu_proxy const& proxy)
      :
        bond_params(proxy.distance_ideal, proxy.weight)
      {
        sites[0] = asu_mappings.map_moved_site_to_asu(
          sites_cart[proxy.i_seq], proxy.i_seq, 0);
        sites[1] = asu_mappings.map_moved_site_to_asu(
          sites_cart[proxy.j_seq], proxy.j_seq, proxy.j_sym);
        init_distance_model();
      }

      // Not available in Python.
      bond(
        asu_cache<> const& cache,
        bond_asu_proxy const& proxy)
      :
        bond_params(proxy.distance_ideal, proxy.weight)
      {
        sites[0] = cache.sites[proxy.i_seq][0];
        sites[1] = cache.sites[proxy.j_seq][proxy.j_sym];
        init_distance_model();
      }

      double
      residual() const { return weight * scitbx::fn::pow2(delta); }

      // Not available in Python.
      scitbx::vec3<double>
      gradient_0() const
      {
        return -weight * 2 * delta / distance_model * (sites[0] - sites[1]);
      }

      af::tiny<scitbx::vec3<double>, 2>
      gradients() const
      {
        af::tiny<vec3, 2> result;
        result[0] = gradient_0();
        result[1] = -result[0];
        return result;
      }

      // Not available in Python.
      void
      add_gradients(
        af::ref<scitbx::vec3<double> > const& gradient_array,
        af::tiny<std::size_t, 2> const& i_seqs) const
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
      double distance_model;
      double delta;

    protected:
      void
      init_distance_model()
      {
        distance_model = (sites[0] - sites[1]).length();
        delta = distance_ideal - distance_model;
      }
  };

  inline
  af::shared<double>
  bond_deltas(
    af::const_ref<scitbx::vec3<double> > const& sites_cart,
    af::const_ref<bond_simple_proxy> const& proxies)
  {
    return detail::generic_deltas<bond_simple_proxy, bond>::get(
      sites_cart, proxies);
  }

  inline
  af::shared<double>
  bond_residuals(
    af::const_ref<scitbx::vec3<double> > const& sites_cart,
    af::const_ref<bond_simple_proxy> const& proxies)
  {
    return detail::generic_residuals<bond_simple_proxy, bond>::get(
      sites_cart, proxies);
  }

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

  inline
  af::shared<double>
  bond_deltas(
    af::const_ref<scitbx::vec3<double> > const& sites_cart,
    direct_space_asu::asu_mappings<> const& asu_mappings,
    af::const_ref<bond_asu_proxy> const& proxies)
  {
    af::shared<double> result((af::reserve(sites_cart.size())));
    for(std::size_t i=0;i<proxies.size();i++) {
      bond restraint(sites_cart, asu_mappings,  proxies[i]);
      result.push_back(restraint.delta);
    }
    return result;
  }

  inline
  af::shared<double>
  bond_residuals(
    af::const_ref<scitbx::vec3<double> > const& sites_cart,
    direct_space_asu::asu_mappings<> const& asu_mappings,
    af::const_ref<bond_asu_proxy> const& proxies)
  {
    af::shared<double> result((af::reserve(sites_cart.size())));
    for(std::size_t i=0;i<proxies.size();i++) {
      bond restraint(sites_cart, asu_mappings,  proxies[i]);
      result.push_back(restraint.residual());
    }
    return result;
  }

  //! Not available in Python.
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

  inline
  double
  bond_residual_sum(
    af::const_ref<scitbx::vec3<double> > const& sites_cart,
    direct_space_asu::asu_mappings<> const& asu_mappings,
    af::const_ref<bond_asu_proxy> const& proxies,
    af::ref<scitbx::vec3<double> > const& gradient_array,
    bool disable_cache=false)
  {
    std::vector<bool> sym_active_flags;
    if (!disable_cache) {
      sym_active_flags.resize(asu_mappings.mappings_const_ref().size(), true);
    }
    return bond_residual_sum(
      sites_cart,
      asu_mappings,
      proxies,
      sym_active_flags,
      gradient_array,
      disable_cache);
  }

  typedef sorted_asu_proxies<bond_simple_proxy, bond_asu_proxy>
    bond_sorted_asu_proxies;

  inline
  double
  bond_residual_sum(
    af::const_ref<scitbx::vec3<double> > const& sites_cart,
    bond_sorted_asu_proxies const& sorted_asu_proxies,
    af::ref<scitbx::vec3<double> > const& gradient_array,
    bool disable_cache=false)
  {
    double result = bond_residual_sum(
      sites_cart,
      sorted_asu_proxies.simple.const_ref(),
      gradient_array);
    result += bond_residual_sum(
      sites_cart,
      *sorted_asu_proxies.asu_mappings(),
      sorted_asu_proxies.sym.const_ref(),
      sorted_asu_proxies.sym_active_flags,
      gradient_array,
      disable_cache);
    return result;
  }

  inline
  af::shared<std::set<std::size_t> >
  bond_sets(
    std::size_t n_sites,
    af::const_ref<restraints::bond_simple_proxy> const& proxies)
  {
    af::shared<std::set<std::size_t> > result;
    result.resize(n_sites);
    for(std::size_t i=0;i<proxies.size();i++) {
      restraints::bond_simple_proxy const& proxy = proxies[i];
      CCTBX_ASSERT(proxy.i_seqs[0] < n_sites);
      CCTBX_ASSERT(proxy.i_seqs[1] < n_sites);
      result[proxy.i_seqs[0]].insert(proxy.i_seqs[1]);
      result[proxy.i_seqs[1]].insert(proxy.i_seqs[0]);
    }
    return result;
  }

  inline
  af::shared<std::set<asu_mapping_index> >
  bond_asu_sets(
    std::size_t n_sites,
    af::const_ref<restraints::bond_asu_proxy> const& proxies)
  {
    af::shared<std::set<asu_mapping_index> > result;
    result.resize(n_sites);
    for(std::size_t i=0;i<proxies.size();i++) {
      asu_mapping_index_pair const& pair = proxies[i];
      CCTBX_ASSERT(pair.i_seq < n_sites);
      CCTBX_ASSERT(pair.j_seq < n_sites);
      result[pair.i_seq].insert(asu_mapping_index(pair.j_seq, pair.j_sym));
      if (pair.j_sym == 0) {
        result[pair.j_seq].insert(asu_mapping_index(pair.i_seq, 0));
      }
    }
    return result;
  }

}} // namespace cctbx::restraints

#endif // CCTBX_RESTRAINTS_BOND_H
