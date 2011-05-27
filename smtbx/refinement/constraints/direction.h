#ifndef SMTBX_REFINEMENT_CONSTRAINTS_DIRECTION_H
#define SMTBX_REFINEMENT_CONSTRAINTS_DIRECTION_H
#include <scitbx/math/least_squares_plane.h>
#include <scitbx/matrix/eigensystem.h>
#include <smtbx/refinement/constraints/reparametrisation.h>

namespace smtbx { namespace refinement { namespace constraints {


/** direction interface
 */
class direction_base {
public:
  virtual ~direction_base() {}
  // returns normalised vector
  virtual cart_t direction(uctbx::unit_cell const &unit_cell) const = 0;
};

/** static direction
 */
class static_direction : public direction_base {
  cart_t direction_;
protected:
  static af::shared<cart_t> sites_to_points(uctbx::unit_cell const &unit_cell,
    af::shared<site_parameter *> const &sites)
  {
    af::shared<cart_t> points(sites.size());
    for (int i=0; i < sites.size(); i++)
      points[i] = unit_cell.orthogonalize(sites[i]->value);
    return points;
  }
public:
  static_direction(cart_t const &direction)
  : direction_(direction.normalize())
  {}

  virtual cart_t direction(uctbx::unit_cell const &) const {
    return direction_;
  }

  static static_direction best_line(af::shared<cart_t> const &points) {
    if (points.size() == 2)
      return static_direction((points[1]-points[0]).normalize());
    return static_direction(calc_plane_normal(points, 2));
  }
  static cart_t calc_best_line(uctbx::unit_cell const &unit_cell,
    af::shared<site_parameter *> const &sites)
  {
    return best_line(sites_to_points(unit_cell, sites)).direction_;
  }

  static static_direction best_plane_normal(af::shared<cart_t> const &points) {
    return static_direction(calc_plane_normal(points, 0));
  }
  static cart_t calc_best_plane_normal(uctbx::unit_cell const &unit_cell,
    af::shared<site_parameter *> const &sites)
  {
    return best_plane_normal(sites_to_points(unit_cell, sites)).direction_;
  }

  static cart_t calc_plane_normal(uctbx::unit_cell const &unit_cell,
    af::shared<site_parameter *> const &sites,
    int eigen_value_index)
  {
    return calc_plane_normal(
      sites_to_points(unit_cell, sites), eigen_value_index);
  }

  static cart_t calc_plane_normal(af::shared<cart_t> const &points,
    int eigen_value_index)
  {
    SMTBX_ASSERT(!(eigen_value_index < 0 || eigen_value_index > 2));
    // plane normal
    if (points.size() == 3 && eigen_value_index == 0) {
      return (points[0]-points[1]).cross(points[2]-points[1]).normalize();
    }
    cart_t cnt(0,0,0);
    for (int i=0; i < points.size(); i++)
      cnt += points[i];
    cnt /= double(points.size());
    scitbx::sym_mat3<double> smv(0,0,0,0,0,0);
    for (int i=0; i < points.size(); i++) {
      cart_t p = points[i]-cnt;
      for (int j=0; j < 3; j++)
        for (int k=j; k < 3; k++)
          smv(j,k) += p[j]*p[k];
    }
    scitbx::matrix::eigensystem::real_symmetric<double> es(smv);
    /* sort the eigenvalues and vectors in acending order:
    this leaves index 0 for the best plane and index 2 - for the worst
    plane - i.e. best line */
    int indices[3] = {0,1,2};
    bool swaps = true;
    while (swaps) {
      swaps = false;
      for (int i=0; i < 2; i++ )  {
        if (es.values()[indices[i]] > es.values()[indices[i+1]] ) {
          std::swap(indices[i], indices[i+1]);
          swaps = true;
        }
      }
    }
    af::versa<double, af::c_grid<2> > ev = es.vectors();
    int off = indices[eigen_value_index]*3;
    return cart_t(ev[off], ev[off+1], ev[off+2]);
  }
};


/** vector direction, best line direction
 */
class vector_direction : public direction_base {
  af::shared<site_parameter *> sites;
public:
  vector_direction(site_parameter *from, site_parameter *to)
  : sites(2)
  {
    sites[0] = from;
    sites[1] = to;
  }

  vector_direction(af::shared<site_parameter *> const &sites_)
  : sites(sites_)
  {
    SMTBX_ASSERT(!(sites.size() < 2));
  }

  virtual cart_t direction(uctbx::unit_cell const &unit_cell) const {
    return static_direction::calc_best_line(unit_cell, sites);
  }
};

/** cross-product/best plane direction
*/
class normal_direction : public direction_base {
  af::shared<site_parameter *> sites;
public:
  normal_direction(af::shared<site_parameter *> const &sites_)
  : sites(sites_)
  {
    SMTBX_ASSERT(!(sites.size() < 3));
  }

  virtual cart_t direction(uctbx::unit_cell const &unit_cell) const {
    return static_direction::calc_best_plane_normal(unit_cell, sites);
  }
};

}}}

#endif // GUARD
