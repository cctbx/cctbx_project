#ifndef MMTBX_SECONDARY_STRUCTURE_H
#define MMTBX_SECONDARY_STRUCTURE_H

#include <mmtbx/error.h>
#include <scitbx/array_family/versa.h>
#include <scitbx/array_family/shared.h>
#include <scitbx/vec3.h>

#include <string>

namespace mmtbx { namespace secondary_structure {

  namespace af = scitbx::af;

  double delta_distance_squared (
    af::const_ref<size_t> base_1,
    af::const_ref<size_t> base_2,
    std::string name_1,
    std::string name_2,
    af::const_ref<std::string> atom_names,
    af::const_ref<scitbx::vec3<double> > sites_cart,
    double distance_ideal)
  {
    for (unsigned i_seq = 0; i_seq < base_1.size(); i_seq++) {
      size_t idx_i = base_1[i_seq];
      if (atom_names[idx_i].compare(name_1) == 0) {
        scitbx::vec3<double> site_i = sites_cart[idx_i];
        for (unsigned j_seq = 0; j_seq < base_2.size(); j_seq++) {
          size_t idx_j = base_2[j_seq];
          if (atom_names[idx_j].compare(name_2) == 0) {
            scitbx::vec3<double> site_j = sites_cart[idx_j];
            double delta_r_ij = (site_i - site_j).length() - distance_ideal;
            if (delta_r_ij < 0) return 0.0;
            return delta_r_ij * delta_r_ij;
          }
        }
        break;
      }
    }
    return -1;
  }

}} // namespace mmtbx::secondary_structure

#endif
