#ifndef MMTBX_ROTAMER_H
#define MMTBX_ROTAMER_H

#include <boost/python/list.hpp>
#include <boost/python/extract.hpp>
#include <cctbx/maptbx/target_and_gradients.h>
#include <scitbx/math/rotate_around_axis.h>
#include <scitbx/vec3.h>

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
  af::shared<af::shared<FloatType> > angles_array;
  af::shared<scitbx::vec3<FloatType> > all_points_result;
  FloatType score_;
  FloatType score_start_;

  fit() {}

  fit(
    FloatType target_value,
    boost::python::list const& axes_,
    boost::python::list const& rotatable_points_indices_,
    boost::python::list const& angles_array_,
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
    SCITBX_ASSERT(boost::python::len(axes_)==
                  boost::python::len(rotatable_points_indices_));
    for(std::size_t i=0;i<boost::python::len(axes_);i++) {
      axes.push_back(
        boost::python::extract<af::shared<std::size_t > >(axes_[i])());
      rotatable_points_indices.push_back(
        boost::python::extract<af::shared<std::size_t> >(
          rotatable_points_indices_[i])());
    }
    for(std::size_t i=0;i<boost::python::len(angles_array_);i++) {
      angles_array.push_back(
        boost::python::extract<af::shared<FloatType> >(angles_array_[i])());
    }
    for(std::size_t j=0;j<angles_array.size();j++) {
      af::shared<FloatType> angles = angles_array[j];
      af::shared<scitbx::vec3<FloatType> > all_points_cp = all_points.deep_copy();
      for(std::size_t i=0;i<angles.size();i++) {
       scitbx::math::rotate_points_around_axis(
          axes[i][0],
          axes[i][1],
          all_points_cp.ref(),
          rotatable_points_indices[i].const_ref(),
          angles[i],
          sin_table,
          cos_table,
          step,
          n);
      }
      FloatType mv = rs_simple::target<FloatType, FloatType>(
        unit_cell,
        density_map,
        all_points_cp.ref(),
        selection);
      if(mv>score_) {
        all_points_result = all_points_cp.deep_copy();
        score_ = mv;
      }
    }
  }

  fit(
    fixed<double> const& xyzrad_bumpers,
    boost::python::list const& axes_,
    boost::python::list const& rotatable_points_indices_,
    boost::python::list const& angles_array_,
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
    SCITBX_ASSERT(boost::python::len(axes_)==
                  boost::python::len(rotatable_points_indices_));
    for(std::size_t i=0;i<boost::python::len(axes_);i++) {
      axes.push_back(
        boost::python::extract<af::shared<std::size_t > >(axes_[i])());
      rotatable_points_indices.push_back(
        boost::python::extract<af::shared<std::size_t> >(
          rotatable_points_indices_[i])());
    }
    for(std::size_t i=0;i<boost::python::len(angles_array_);i++) {
      angles_array.push_back(
        boost::python::extract<af::shared<FloatType> >(angles_array_[i])());
    }
    af::tiny<int, 3> a = density_map.accessor().all();
    af::flex_grid<> const& g = af::flex_grid<>(a[0],a[1],a[2]);
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
    for(std::size_t j=0;j<angles_array.size();j++) {
      af::shared<FloatType> angles = angles_array[j];
      af::shared<scitbx::vec3<FloatType> >
        all_points_cp = all_points.sites_cart.deep_copy();
      for(std::size_t i=0;i<angles.size();i++) {
        scitbx::math::rotate_points_around_axis(
          axes[i][0],
          axes[i][1],
          all_points_cp.ref(),
          rotatable_points_indices[i].const_ref(),
          angles[i],
          sin_table,
          cos_table,
          step,
          n);
      }
      bool ignore = false;
      for(std::size_t k=0;k<selection_clash.size();k++) {
        if(ignore) break;
        cctbx::cartesian<FloatType> s1 = all_points_cp[ selection_clash[k] ];
        FloatType rad = all_points.radii[ selection_clash[k] ];
        for(std::size_t m=0;m<xyzrad_bumpers.sites_cart.size();m++) {
          cctbx::cartesian<FloatType> s2 = xyzrad_bumpers.sites_cart[m];
          FloatType d = (s1[0]-s2[0])*(s1[0]-s2[0]) +
                        (s1[1]-s2[1])*(s1[1]-s2[1]) +
                        (s1[2]-s2[2])*(s1[2]-s2[2]);
          if(d>14.45) continue;
          if(std::sqrt(d)<xyzrad_bumpers.radii[m]+rad) {
            ignore = true;
            break;
          }
        }
      }
      if(!ignore) {
        af::tiny<FloatType, 2> r = rs_simple::score<FloatType, FloatType>(
          unit_cell,
          density_map,
          all_points_cp.ref(),
          selection_rsr,
          all_points.bonded_pairs,
          all_points.weights.const_ref(),
          all_points.ref_map_max.const_ref(),
          all_points.ref_map_min.const_ref());
        if(r[0]>0 && r[1]>score_) {
          all_points_result = all_points_cp;
          score_ = r[1];
        }
      }
    }
  }

  fit(
    af::shared<scitbx::vec3<double> > const& sites_cart_start,
    boost::python::list const& axes_,
    boost::python::list const& rotatable_points_indices_,
    boost::python::list const& angles_array_,
    af::shared<scitbx::vec3<double> > all_points,
    af::const_ref<double> const& sin_table,
    af::const_ref<double> const& cos_table,
    double const& step,
    int const& n)
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
    for(std::size_t i=0;i<boost::python::len(angles_array_);i++) {
      angles_array.push_back(
        boost::python::extract<af::shared<double> >(angles_array_[i])());
    }
    double mv_best = 1.e+9;//target_value;
    for(std::size_t j=0;j<angles_array.size();j++) {
      af::shared<double> angles = angles_array[j];
      af::shared<scitbx::vec3<double> > all_points_cp = all_points.deep_copy();
      for(std::size_t i=0;i<angles.size();i++) {
       scitbx::math::rotate_points_around_axis(
          axes[i][0],
          axes[i][1],
          all_points_cp.ref(),
          rotatable_points_indices[i].const_ref(),
          angles[i],
          sin_table,
          cos_table,
          step,
          n);
      }
      double mv = 0;
      for(std::size_t i=0;i<all_points_cp.size();i++) {
        scitbx::vec3<double> delta = all_points_cp[i]-sites_cart_start[i];
        mv += std::sqrt(delta[0]*delta[0]+delta[1]*delta[1]+delta[2]*delta[2]);
      }
      if(mv<mv_best) {
        all_points_result = all_points_cp.deep_copy();
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
