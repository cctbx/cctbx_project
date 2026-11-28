#ifndef MMTBX_ROTAMER_H
#define MMTBX_ROTAMER_H

#include <boost/python/list.hpp>
#include <boost/python/extract.hpp>
#include <cctbx/maptbx/target_and_gradients.h>
#include <scitbx/math/rotate_around_axis.h>
#include <scitbx/vec3.h>
#include <vector>
#include <algorithm>

#include <cctbx/maptbx/interpolation.h>

namespace mmtbx { namespace rotamer {
namespace af=scitbx::af;
namespace rs_simple=cctbx::maptbx::target_and_gradients::simple;

template <typename FloatType=double>
class moving
{
  public:
    af::shared<scitbx::vec3<FloatType> > sites_cart;
    af::shared<scitbx::vec3<FloatType> > sites_cart_start;
    af::shared<FloatType> radii;
    af::shared<FloatType> weights;
    af::shared<af::tiny<std::size_t, 2> > bonded_pairs;
    af::shared<FloatType> ref_map_max;
    af::shared<FloatType> ref_map_min;

    moving() { // Can I remove this?
     sites_cart.fill(0.0);
     sites_cart_start.fill(0.0);
     radii.fill(0.0);
     weights.fill(0.0);
    }

    moving(af::shared<scitbx::vec3<FloatType> > const& sites_cart_,
           af::shared<scitbx::vec3<FloatType> > const& sites_cart_start_,
           af::shared<FloatType> const& radii_,
           af::shared<FloatType> const& weights_,
           boost::python::list const& bonded_pairs_,
           af::shared<FloatType> const& ref_map_max_,
           af::shared<FloatType> const& ref_map_min_
           )
    :
      sites_cart(sites_cart_),
      sites_cart_start(sites_cart_start_),
      radii(radii_),
      weights(weights_),
      ref_map_max(ref_map_max_),
      ref_map_min(ref_map_min_)
    {
      for(std::size_t i=0;i<boost::python::len(bonded_pairs_);i++) {
         af::shared<size_t> p =
           boost::python::extract<af::shared<size_t> >(bonded_pairs_[i])();
         bonded_pairs.push_back( af::tiny<std::size_t, 2>(p[0],p[1]) );
      }
    }
};

template <typename FloatType=double>
class fixed
{
  public:
    af::shared<scitbx::vec3<FloatType> > sites_cart;
    af::shared<FloatType> radii;
    fixed() {
     sites_cart.fill(0.0);
     radii.fill(0.0);
    }

    fixed(af::shared<scitbx::vec3<FloatType> > const& sites_cart_,
          af::shared<FloatType> const& radii_)
    :
      sites_cart(sites_cart_), radii(radii_)
    {}

};

template <typename FloatType=double>
class fit {
public:
  af::shared<af::shared<std::size_t> > axes;
  af::shared<af::shared<std::size_t> > rotatable_points_indices;
  // Sets of angles (one angle per axis in each) one after another. Tables of
  // rotamers have up to 5.e5 sets: given as a Python list of lists they took
  // longer to convert (for every residue) than to use; one array of numbers
  // is used as it is.
  af::shared<FloatType> angles_owned;
  af::const_ref<FloatType> angles_flat;
  std::size_t n_sets;
  af::shared<scitbx::vec3<FloatType> > all_points_result;
  FloatType score_;
  FloatType score_start_;

  fit() {}

  void set_axes(
    boost::python::list const& axes_,
    boost::python::list const& rotatable_points_indices_)
  {
    SCITBX_ASSERT(boost::python::len(axes_)==
                  boost::python::len(rotatable_points_indices_));
    for(std::size_t i=0;i<boost::python::len(axes_);i++) {
      axes.push_back(
        boost::python::extract<af::shared<std::size_t > >(axes_[i])());
      rotatable_points_indices.push_back(
        boost::python::extract<af::shared<std::size_t> >(
          rotatable_points_indices_[i])());
    }
  }

  // list of lists
  void set_angles(boost::python::list const& angles_array_)
  {
    n_sets = boost::python::len(angles_array_);
    for(std::size_t j=0;j<n_sets;j++) {
      af::shared<FloatType> angles =
        boost::python::extract<af::shared<FloatType> >(angles_array_[j])();
      SCITBX_ASSERT(angles.size() <= axes.size());
      for(std::size_t i=0;i<axes.size();i++) {
        // sets shorter than the number of axes: the rest is not rotated
        angles_owned.push_back(i<angles.size() ? angles[i] : 0);
      }
    }
    angles_flat = angles_owned.const_ref();
  }

  // one array: n_sets * (number of axes) values
  void set_angles(af::const_ref<FloatType> const& angles_flat_)
  {
    angles_flat = angles_flat_;
    n_sets = 0;
    if(axes.size()>0) {
      SCITBX_ASSERT(angles_flat.size() % axes.size() == 0);
      n_sets = angles_flat.size() / axes.size();
    }
  }

  // Search state. Sets of angles that share leading angles share the
  // rotations about the first axes: sites after rotation about axes 0..i-1
  // are kept (levels[i]; levels[0] are the sites given) and only what follows
  // the first angle that differs from the previous set is redone. Tables
  // sorted by angles (rotamer_manager.get_chi_angles_flat) make consecutive
  // sets differ mostly in the last angle. Any order gives the same sites as
  // rotating each set from scratch.
  std::vector<af::shared<scitbx::vec3<FloatType> > > levels;
  std::size_t n_valid; // levels[1..n_valid] are those of the previous set

  void start(af::const_ref<scitbx::vec3<FloatType> > const& sites)
  {
    levels.clear();
    for(std::size_t i=0;i<=axes.size();i++) {
      levels.push_back(af::shared<scitbx::vec3<FloatType> >(
        sites.begin(), sites.end()));
    }
    n_valid = 0;
  }

  // number of levels that set j shares with the set before it
  std::size_t shared_levels(std::size_t j)
  {
    if(j==0) return 0;
    std::size_t n_axes = axes.size();
    const FloatType* a = &angles_flat[j*n_axes];
    const FloatType* b = a - n_axes;
    std::size_t p = 0;
    while(p<n_valid && a[p]==b[p]) p++;
    return p;
  }

  // levels[i+1]: levels[i] rotated about axis i by the angle of set j
  void rotate_level(
    std::size_t j,
    std::size_t i,
    af::const_ref<FloatType> const& sin_table,
    af::const_ref<FloatType> const& cos_table,
    FloatType const& step,
    int const& n)
  {
    std::copy(levels[i].begin(), levels[i].end(), levels[i+1].begin());
    scitbx::math::rotate_points_around_axis(
      axes[i][0],
      axes[i][1],
      levels[i+1].ref(),
      rotatable_points_indices[i].const_ref(),
      angles_flat[j*axes.size()+i],
      sin_table,
      cos_table,
      step,
      n);
  }

  // sites rotated by set j of angles (sets are to be asked for in order)
  af::shared<scitbx::vec3<FloatType> > const& rotated(
    std::size_t j,
    af::const_ref<FloatType> const& sin_table,
    af::const_ref<FloatType> const& cos_table,
    FloatType const& step,
    int const& n)
  {
    for(std::size_t i=shared_levels(j);i<axes.size();i++) {
      rotate_level(j, i, sin_table, cos_table, step, n);
    }
    n_valid = axes.size();
    return levels[axes.size()];
  }

  // AnglesType: boost::python::list (of lists) or af::const_ref<FloatType>
  template <typename AnglesType>
  fit(
    FloatType target_value,
    boost::python::list const& axes_,
    boost::python::list const& rotatable_points_indices_,
    AnglesType const& angles_array_,
    af::const_ref<FloatType, af::c_grid_padded<3> > const& density_map,
    af::shared<scitbx::vec3<FloatType> > all_points,
    cctbx::uctbx::unit_cell const& unit_cell,
    af::const_ref<std::size_t> const& selection,
    af::const_ref<FloatType> const& sin_table,
    af::const_ref<FloatType> const& cos_table,
    FloatType const& step,
    int const& n)
  :
  score_(target_value), score_start_(target_value)
  {
    set_axes(axes_, rotatable_points_indices_);
    set_angles(angles_array_);
    start(all_points.const_ref());
    for(std::size_t j=0;j<n_sets;j++) {
      af::shared<scitbx::vec3<FloatType> > const& sites = rotated(
        j, sin_table, cos_table, step, n);
      FloatType mv = rs_simple::target<FloatType, FloatType>(
        unit_cell,
        density_map,
        sites.const_ref(),
        selection);
      if(mv>score_) {
        all_points_result = sites.deep_copy();
        score_ = mv;
      }
    }
  }

  // A set of angles is a candidate if no atom of selection_clash is closer
  // to a bumper than the sum of radii and map values at all atoms of
  // selection_rsr are between their min and max; the best candidate has the
  // largest sum of weighted map values, if larger than that of
  // sites_cart_start.
  // An atom is in place once the last axis that moves it is done, so it is
  // checked there, once for all sets that share the angles up to that axis:
  // if it fails none of them is a candidate and they are skipped. Same
  // result as checking all atoms for every set.
  template <typename AnglesType>
  fit(
    fixed<double> const& xyzrad_bumpers,
    boost::python::list const& axes_,
    boost::python::list const& rotatable_points_indices_,
    AnglesType const& angles_array_,
    af::const_ref<FloatType, af::c_grid_padded<3> > const& density_map,
    moving<double> const& all_points,
    cctbx::uctbx::unit_cell const& unit_cell,
    af::const_ref<std::size_t> const& selection_clash,
    af::const_ref<std::size_t> const& selection_rsr,
    af::const_ref<FloatType> const& sin_table,
    af::const_ref<FloatType> const& cos_table,
    FloatType const& step,
    int const& n)
  {
    set_axes(axes_, rotatable_points_indices_);
    set_angles(angles_array_);
    // initial score
    af::tiny<FloatType, 2> r = rs_simple::score<FloatType, FloatType>(
      unit_cell,
      density_map,
      all_points.sites_cart_start.const_ref(),
      selection_rsr,
      all_points.bonded_pairs,
      all_points.weights.const_ref(),
      all_points.ref_map_max.const_ref(),
      all_points.ref_map_min.const_ref());
    score_ = r[1];
    score_start_ = score_;
    std::size_t n_axes = axes.size();
    std::size_t n_atoms = all_points.sites_cart.size();
    // level at which an atom is in place: 0 if no axis moves it
    std::vector<std::size_t> level_of(n_atoms, 0);
    for(std::size_t i=0;i<n_axes;i++) {
      for(std::size_t k=0;k<rotatable_points_indices[i].size();k++) {
        std::size_t a = rotatable_points_indices[i][k];
        SCITBX_ASSERT(a < n_atoms);
        level_of[a] = i+1;
      }
    }
    std::vector<std::vector<std::size_t> > clash_at(n_axes+1);
    std::vector<std::vector<std::size_t> > rsr_at(n_axes+1);
    for(std::size_t k=0;k<selection_clash.size();k++) {
      clash_at[level_of[selection_clash[k]]].push_back(selection_clash[k]);
    }
    for(std::size_t k=0;k<selection_rsr.size();k++) {
      rsr_at[level_of[selection_rsr[k]]].push_back(selection_rsr[k]);
    }
    af::const_ref<scitbx::vec3<FloatType> >
      bumpers = xyzrad_bumpers.sites_cart.const_ref();
    af::const_ref<FloatType> bumper_radii = xyzrad_bumpers.radii.const_ref();
    af::const_ref<FloatType> weights = all_points.weights.const_ref();
    af::const_ref<FloatType> map_max = all_points.ref_map_max.const_ref();
    af::const_ref<FloatType> map_min = all_points.ref_map_min.const_ref();
    std::vector<FloatType> value(n_atoms, 0); // weighted map value at atom
    // Bumpers that atoms checked at a level can reach. Rotation about the
    // axis of the level keeps the distance of an atom to a point of that
    // axis, and bumpers farther than sqrt(14.45) from an atom do not count:
    // only bumpers within (largest such distance + that) of the axis point
    // can matter. They are listed again when the axis point moves, that is
    // when the level before is redone.
    const FloatType reach = std::sqrt(14.45) + 1.e-3;
    std::vector<std::vector<std::size_t> > near(n_axes+1);
    std::vector<bool> near_valid(n_axes+1, false);
    for(std::size_t m=0;m<bumpers.size();m++) near[0].push_back(m);
    near_valid[0] = true;
    start(all_points.sites_cart.const_ref());
    // npos: none; else the first level with an atom that fails
    const std::size_t npos = n_axes+2;
    std::size_t failed = npos;
    for(std::size_t j=0;j<=n_sets;j++) { // j=0: level 0, then set j-1
      std::size_t first, last;
      if(j==0) { first = 0; last = 0; }
      else {
        std::size_t p = shared_levels(j-1);
        n_valid = p;
        if(failed <= p) continue;
        failed = npos;
        first = p+1;
        last = n_axes;
      }
      for(std::size_t level=first;level<=last;level++) {
        if(level>0) {
          rotate_level(j-1, level-1, sin_table, cos_table, step, n);
          n_valid = level;
          if(level<n_axes) near_valid[level+1] = false;
        }
        af::const_ref<scitbx::vec3<FloatType> >
          sites = levels[level].const_ref();
        bool ok = true;
        std::vector<std::size_t> const& ca = clash_at[level];
        std::vector<std::size_t>& nb = near[level];
        if(!near_valid[level]) {
          scitbx::vec3<FloatType> const& b = sites[axes[level-1][1]];
          FloatType r_max = 0;
          for(std::size_t k=0;k<ca.size();k++) {
            r_max = std::max(r_max, (sites[ca[k]]-b).length());
          }
          r_max = (r_max + reach)*(r_max + reach);
          nb.clear();
          for(std::size_t m=0;m<bumpers.size();m++) {
            if((bumpers[m]-b).length_sq() <= r_max) nb.push_back(m);
          }
          near_valid[level] = true;
        }
        for(std::size_t k=0;k<ca.size() && ok;k++) {
          scitbx::vec3<FloatType> const& s1 = sites[ca[k]];
          FloatType rad = all_points.radii[ca[k]];
          for(std::size_t i_m=0;i_m<nb.size();i_m++) {
            std::size_t m = nb[i_m];
            scitbx::vec3<FloatType> const& s2 = bumpers[m];
            FloatType d = (s1[0]-s2[0])*(s1[0]-s2[0]) +
                          (s1[1]-s2[1])*(s1[1]-s2[1]) +
                          (s1[2]-s2[2])*(s1[2]-s2[2]);
            if(d>14.45) continue;
            if(std::sqrt(d)<bumper_radii[m]+rad) {
              ok = false;
              break;
            }
          }
        }
        std::vector<std::size_t> const& ra = rsr_at[level];
        for(std::size_t k=0;k<ra.size() && ok;k++) {
          std::size_t a = ra[k];
          FloatType mv = cctbx::maptbx::eight_point_interpolation(
            density_map, unit_cell.fractionalize(sites[a]));
          if(mv < map_min[a] || mv > map_max[a]) ok = false;
          value[a] = mv * weights[a];
        }
        if(!ok) {
          failed = level;
          break;
        }
      }
      if(j==0 && failed != npos) break; // an atom that does not move fails
      if(j==0 || failed != npos) continue;
      FloatType mv = 0;
      for(std::size_t k=0;k<selection_rsr.size();k++) {
        mv += value[selection_rsr[k]];
      }
      if(mv>score_) {
        all_points_result = levels[n_axes].deep_copy();
        score_ = mv;
      }
    }
  }

  template <typename AnglesType>
  fit(
    af::shared<scitbx::vec3<double> > const& sites_cart_start,
    boost::python::list const& axes_,
    boost::python::list const& rotatable_points_indices_,
    AnglesType const& angles_array_,
    af::shared<scitbx::vec3<double> > all_points,
    af::const_ref<double> const& sin_table,
    af::const_ref<double> const& cos_table,
    double const& step,
    int const& n)
  {
    set_axes(axes_, rotatable_points_indices_);
    set_angles(angles_array_);
    double mv_best = 1.e+9;//target_value;
    start(all_points.const_ref());
    for(std::size_t j=0;j<n_sets;j++) {
      af::shared<scitbx::vec3<FloatType> > const& sites = rotated(
        j, sin_table, cos_table, step, n);
      double mv = 0;
      for(std::size_t i=0;i<sites.size();i++) {
        scitbx::vec3<double> delta = sites[i]-sites_cart_start[i];
        mv += std::sqrt(delta[0]*delta[0]+delta[1]*delta[1]+delta[2]*delta[2]);
      }
      if(mv<mv_best) {
        all_points_result = sites.deep_copy();
        mv_best = mv;
      }
    }
  }

  af::shared<scitbx::vec3<double> > result() { return all_points_result; }

  FloatType score() { return score_;}
  FloatType score_start() { return score_start_;}

};

}} // namespace mmtbx::rotamer

#endif // MMTBX_ROTAMER_H
