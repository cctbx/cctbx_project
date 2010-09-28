#ifndef SMTBX_REFINEMENT_CONSTRAINTS_GEOMETRICAL_HYDROGENS_H
#define SMTBX_REFINEMENT_CONSTRAINTS_GEOMETRICAL_HYDROGENS_H

#include <smtbx/refinement/constraints/reparametrisation.h>
#include <scitbx/constants.h>
#include <scitbx/math/orthonormal_basis.h>
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


/// Base class for geometrical hydrogen
template <int n_hydrogens>
class geometrical_hydrogen_sites
  : public crystallographic_parameter
{
public:
  geometrical_hydrogen_sites(int n_parameters,
                             scatterer_type *h)
  : crystallographic_parameter(n_parameters), hydrogen(h)
  {}

  geometrical_hydrogen_sites(int n_parameters,
                             scatterer_type *h0, scatterer_type *h1)
  : crystallographic_parameter(n_parameters), hydrogen(h0, h1)
  {}

  geometrical_hydrogen_sites(int n_parameters,
                             af::tiny<scatterer_type *, n_hydrogens> h)
  : crystallographic_parameter(n_parameters), hydrogen(h)
  {}

  virtual scatterer_sequence_type scatterers() const {
    return hydrogen.const_ref();
  }

  virtual std::size_t size() const {
    return 3*n_hydrogens;
  }

  virtual void store(uctbx::unit_cell const &unit_cell) const {
    for (int i=0; i<hydrogen.size(); ++i) {
      hydrogen[i]->site = unit_cell.fractionalize(x_h[i]);
    }
  }

protected:
  af::tiny<scatterer_type *, n_hydrogens> hydrogen;
  af::tiny<cart_t, n_hydrogens> x_h;
};


/// Model of Y-XHn with tetrahedral angles
/**
  X is referred to as the "pivot" and Y as the "pivot neighbour".

  All angles Hi-X-Hj and Hi-X-Y are tetrahedral.
  All distances X-Hi are equal. That unique distance may be a variable
  parameter if stretching is allowed.
  The Hydrogen sites ride on the pivot site.
*/
template <int n_hydrogens, bool staggered>
class terminal_tetrahedral_xhn_sites
  : public geometrical_hydrogen_sites<n_hydrogens>
{
public:
  /// Construct Hydrogens freely rotating about the bond X-Y
  /** e_zero_azimuth is the vector such that the plane (e_zero_azimuth, XY)
      defines azimuth = 0. Thus it shall be such that it can never become nearly
      colinear to the bond between the pivot and its neighbour.
   */
  terminal_tetrahedral_xhn_sites(site_parameter *pivot,
                                 site_parameter *pivot_neighbour,
                                 independent_scalar_parameter *azimuth,
                                 independent_scalar_parameter *length,
                                 cart_t const &e_zero_azimuth,
                                 af::tiny<crystallographic_parameter::scatterer_type *,
                                          n_hydrogens> const &hydrogen)
    : geometrical_hydrogen_sites<n_hydrogens>(4, hydrogen),
      e_zero_azimuth(e_zero_azimuth)
  {
    SMTBX_ASSERT(!staggered);
    this->set_arguments(pivot, pivot_neighbour, azimuth, length);
  }

  /// Construct Hydrogens staggered on the specified site
  /** stagger shall be a neighbour of the pivot neighbour onto which to
      stagger the Hydrogen's.
   */
  terminal_tetrahedral_xhn_sites(site_parameter *pivot,
                                 site_parameter *pivot_neighbour,
                                 site_parameter *stagger,
                                 independent_scalar_parameter *length,
                                 af::tiny<crystallographic_parameter::scatterer_type *,
                                          n_hydrogens> const &hydrogen)
  : geometrical_hydrogen_sites<n_hydrogens>(4, hydrogen)
  {
    SMTBX_ASSERT(staggered);
    this->set_arguments(pivot, pivot_neighbour, stagger, length);
  }

  virtual void linearise(uctbx::unit_cell const &unit_cell,
                         sparse_matrix_type *jacobian_transpose);

private:
  cart_t e_zero_azimuth;
};


/// A parameter representing an angle in radians starting as tetrahedral
class angle_starting_tetrahedral : public independent_scalar_parameter
{
public:
  angle_starting_tetrahedral(bool variable=true)
    : independent_scalar_parameter(constants::tetrahedral_angle, variable)
  {}
};


/// Model of X-CH2-Y
/**
  C is referred to as the "pivot" and X and Y as pivot's neighbour 1 and 2.

  All angles Hi-C-X and Hi-C-Y are equal.
  The angle H-C-H is refinable (flapping).
*/
class secondary_ch2_sites
  : public geometrical_hydrogen_sites<2>
{
public:
  secondary_ch2_sites(site_parameter *pivot,
                      site_parameter *pivot_neighbour_0,
                      site_parameter *pivot_neighbour_1,
                      independent_scalar_parameter *length,
                      angle_starting_tetrahedral *h_c_h,
                      scatterer_type *hydrogen_0,
                      scatterer_type *hydrogen_1)
  : geometrical_hydrogen_sites<2>(5, hydrogen_0, hydrogen_1)
  {
    set_arguments(pivot, pivot_neighbour_0, pivot_neighbour_1, length, h_c_h);
  }

  virtual void linearise(uctbx::unit_cell const &unit_cell,
                         sparse_matrix_type *jacobian_transpose);
};


/// Model of tertiary CH
/** All angles Hi-C-X are equal.
 */
class tertiary_ch_site
  : public geometrical_hydrogen_sites<1>
{
public:
  tertiary_ch_site(site_parameter *pivot,
                   site_parameter *pivot_neighbour_0,
                   site_parameter *pivot_neighbour_1,
                   site_parameter *pivot_neighbour_2,
                   independent_scalar_parameter *length,
                   scatterer_type *hydrogen)
  : geometrical_hydrogen_sites<1>(5, hydrogen)
  {
    set_arguments(pivot,
                  pivot_neighbour_0, pivot_neighbour_1, pivot_neighbour_2,
                  length);
  }

  virtual void linearise(uctbx::unit_cell const &unit_cell,
                         sparse_matrix_type *jacobian_transpose);
};


/// Model of aromatic C-H or amide N-H
/** Denoting C or N as X and the two neighbours of X as Y and Z,
    Z-X-Y is bisected by X-H.
*/
class secondary_planar_xh_site:
  public geometrical_hydrogen_sites<1>
{
public:
  secondary_planar_xh_site(site_parameter *pivot,
                           site_parameter *pivot_neighbour_0,
                           site_parameter *pivot_neighbour_1,
                           independent_scalar_parameter *length,
                           scatterer_type *hydrogen)
  : geometrical_hydrogen_sites<1>(4, hydrogen)
  {
    set_arguments(pivot, pivot_neighbour_0, pivot_neighbour_1, length);
  }

  virtual void linearise(uctbx::unit_cell const &unit_cell,
                         sparse_matrix_type *jacobian_transpose);
};


/// Model of terminal Z-Y=XH2 (ethylenic CH2 or amide NH2)
/**
    X is referred to as the "pivot" whereas Y is the pivot's neighbour
    and Z is the pivot's neighbour's substituent.

    The two hydrogen atoms are in the plane ZYX,
    and XY bissects the 120-degree angle H1-X-H2.
*/
class terminal_planar_xh2_sites
  : public geometrical_hydrogen_sites<2>
{
public:
  terminal_planar_xh2_sites(site_parameter *pivot,
                            site_parameter *pivot_neighbour,
                            site_parameter *pivot_neighbour_substituent,
                            independent_scalar_parameter *length,
                            scatterer_type *hydrogen_0,
                            scatterer_type *hydrogen_1)
    : geometrical_hydrogen_sites<2>(4, hydrogen_0, hydrogen_1)
  {
    set_arguments(pivot, pivot_neighbour, pivot_neighbour_substituent, length);
  }

  virtual void linearise(uctbx::unit_cell const &unit_cell,
                         sparse_matrix_type *jacobian_transpose);
};


/// Model of acetylenic X-CH
/**
    X-C-H is linear
*/
class terminal_linear_ch_site
  : public geometrical_hydrogen_sites<1>
{
public:
  terminal_linear_ch_site(site_parameter *pivot,
                          site_parameter *pivot_neighbour,
                          independent_scalar_parameter *length,
                          scatterer_type *hydrogen)
    : geometrical_hydrogen_sites<1>(3, hydrogen)
  {
    set_arguments(pivot, pivot_neighbour, length);
  }

  virtual void linearise(uctbx::unit_cell const &unit_cell,
                         sparse_matrix_type *jacobian_transpose);
};

}}}

#endif // GUARD
