#include <assert.h>
#include <math.h>
#include <iostream>
#include <cctbx/adptbx.h>
#include <cctbx/error.h>
#include <cctbx/adp_restraints/rigid_bond.h>

using namespace std;
namespace cctbx { namespace adp_restraints {

rigid_bond_pair::rigid_bond_pair(vec3<double> const& site1,
                                 vec3<double> const& site2,
                                 sym_mat3<double> const& ustar1,
                                 sym_mat3<double> const& ustar2,
                                 cctbx::uctbx::unit_cell const& uc)
{
    sym_mat3<double> g = uc.metrical_matrix();
    vec3<double> l_12 = site1 - site2;
    vec3<double> l_21 = site2 - site1;
    double bond_length_sq = l_12 * g * l_12;
    z_12_ = (g * l_12) * ustar1 * (g * l_12);
    z_21_ = (g * l_21) * ustar2 * (g * l_21);
    delta_z_ = (z_12_ - z_21_) / bond_length_sq;
}

}} // namespace cctbx::apd_restraints
