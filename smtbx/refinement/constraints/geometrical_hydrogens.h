#ifndef SMTBX_REFINEMENT_CONSTRAINTS_GEOMETRICAL_HYDROGENS_H
#define SMTBX_REFINEMENT_CONSTRAINTS_GEOMETRICAL_HYDROGENS_H

#include <smtbx/refinement/constraints/reparametrisation.h>
#include <scitbx/constants.h>
#include <scitbx/math/orthonormal_basis.h>
#include <boost/optional.hpp>
#include <cmath>


namespace smtbx { namespace refinement { namespace constraints {

namespace constants {
  using namespace scitbx::constants;
  static double const cos_tetrahedral_angle = -1./3.;
  static double const tetrahedral_angle = std::acos(cos_tetrahedral_angle);
  static double const sin_tetrahedral_angle
    = std::sqrt(1. - cos_tetrahedral_angle*cos_tetrahedral_angle);
  static double const half_sqrt_3 = std::sqrt(3.)/2;
}


/// Model of Y-XHn with tetrahedral angles
/**
  X is referred to as the "pivot" and Y as the "pivot neighbour".

  All angles Hi-X-Hj and Hi-X-Y are tetrahedral.
  All distances X-Hi are equal. That unique distance may be a variable
  parameter if stretching is allowed.
  A free rotation around the bond Y-X is allowed, characterized by an
  azimuthal angle.

  The Hydrogen sites ride on the pivot site.
*/
class terminal_tetrahedral_xhn_sites : public crystallographic_parameter
{
public:
  /// Construct sites depending on the given parameter groups
  /** e_zero_azimuth is the vector defining azimuth = 0. It shall be such that
      it can never become nearly colinear to the bond between the pivot
      and its neighbour.
   */
  terminal_tetrahedral_xhn_sites(site_parameter *pivot,
                                 site_parameter *pivot_neighbour,
                                 independent_scalar_parameter *azimuth,
                                 independent_scalar_parameter *length,
                                 cart_t const &e_zero_azimuth,
                                 af::small<scatterer_pointer, 3> const &hydrogen)
    : crystallographic_parameter(4),
      e_zero_azimuth(e_zero_azimuth),
      hydrogen(hydrogen),
      x_h(hydrogen.size())
  {
    set_arguments(pivot, pivot_neighbour, azimuth, length);
  }

  virtual std::size_t size() const;

  virtual void linearise(uctbx::unit_cell const &unit_cell,
                         sparse_matrix_type *jacobian_transpose);

  virtual void store(uctbx::unit_cell const &unit_cell) const;

private:
  af::small<scatterer_pointer, 3> hydrogen;
  cart_t e_zero_azimuth;
  af::small<cart_t, 3> x_h;
};


}}}

#endif // GUARD
