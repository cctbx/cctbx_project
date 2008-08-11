#ifndef CCTBX_ADP_RESTRAINTS_RIGID_BOND_H
#define CCTBX_ADP_RESTRAINTS_RIGID_BOND_H

#include <cctbx/error.h>
#include <cctbx/adptbx.h>

#include <scitbx/array_family/versa.h>
#include <vector>
#include <scitbx/array_family/shared.h>
#include <cctbx/import_scitbx_af.h>

namespace cctbx { namespace adp_restraints {

using scitbx::vec3;
using scitbx::mat3;
using scitbx::sym_mat3;

class rigid_bond_pair {
public:
    rigid_bond_pair(vec3<double> const& site1,
                    vec3<double> const& site2,
                    sym_mat3<double> const& ustar1,
                    sym_mat3<double> const& ustar2,
                    cctbx::uctbx::unit_cell const& uc);
    double z_12() { return z_12_; }
    double z_21() { return z_21_; }
    double delta_z() { return delta_z_; }
private:
    double z_12_, z_21_, delta_z_;
};

}} // namespace cctbx::apd_restraints

#endif // CCTBX_ADP_RESTRAINTS_RIGID_BOND_H
