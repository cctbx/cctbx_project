#ifndef CCTBX_GEOMETRY_RESTRAINTS_BOND_MISC_H
#define CCTBX_GEOMETRY_RESTRAINTS_BOND_MISC_H

#include <cctbx/geometry_restraints/bond.h>
#include <cctbx/crystal/pair_tables.h>

namespace cctbx { namespace geometry_restraints {

  inline
  double
  home_restraints_summation_skip_special_positions(
    af::const_ref<scitbx::vec3<double> > const& sites_cart,
    af::ref<scitbx::vec3<double> > const& gradients,
    af::const_ref<std::size_t> const& site_symmetry_table_indices,
    af::const_ref<scitbx::vec3<double> > const& home_sites_cart,
    af::const_ref<std::size_t> const& iselection,
    double weight,
    double slack)
  {
    CCTBX_ASSERT(gradients.size() == 0
              || gradients.size() == sites_cart.size());
    CCTBX_ASSERT(site_symmetry_table_indices.size() == 0
              || site_symmetry_table_indices.size() == sites_cart.size());
    CCTBX_ASSERT(home_sites_cart.size() == sites_cart.size());
    double residual_sum = 0;
    for(std::size_t ii=0;ii<iselection.size();ii++) {
      std::size_t i_seq = iselection[ii];
      if (   site_symmetry_table_indices.size() != 0
          && site_symmetry_table_indices[i_seq] != 0) {
        continue;
      }
      af::tiny<scitbx::vec3<double>, 2> sites(
        sites_cart[i_seq], home_sites_cart[i_seq]);
      bond b(sites, /* distance_ideal */ 0, weight, slack);
      residual_sum += b.residual();
      if (gradients.size() != 0) {
        gradients[i_seq] += b.gradient_0();
      }
    }
    return residual_sum;
  }

}} // namespace cctbx::geometry_restraints

#endif // GUARD
