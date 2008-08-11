#ifndef SCITBX_MATH_MINIMUM_COVERING_SPHERE
#define SCITBX_MATH_MINIMUM_COVERING_SPHERE

#include <scitbx/vec3.h>
#include <scitbx/error.h>
#include <vector>

namespace scitbx { namespace math {

  template <typename FloatType=double>
  class sphere_3d
  {
    public:
      //! Default constructor. Some data members are not initialized!
      sphere_3d() {}

      sphere_3d(
        vec3<FloatType> const& center,
        FloatType const& radius)
      :
        center_(center),
        radius_(radius)
      {}

      //! Center of the minimum covering sphere.
      vec3<FloatType>
      center() const { return center_; }

      //! Radius of the minimum covering sphere.
      FloatType
      radius() const { return radius_; }

      //! New sphere with expanded radius.
      sphere_3d
      expand(FloatType const& additional_radius) const
      {
        return sphere_3d(center_, radius_+additional_radius);
      }

      //! New sphere with radius*(1+additional_relative_radius).
      sphere_3d
      expand_relative(FloatType const& additional_relative_radius) const
      {
        return sphere_3d(center_, radius_*(1+additional_relative_radius));
      }

      //! True if point is inside sphere.
      bool
      is_inside(vec3<FloatType> const& point)
      {
        if ((point - center_).length_sq() <= radius_ * radius_) {
          return true;
        }
        return false;
      }

      //! Coordinates of lower-left corner of covering box.
      vec3<FloatType>
      box_min() const { return center_ - radius_; }

      //! Coordinates of upper-right corner of covering box.
      vec3<FloatType>
      box_max() const { return center_ + radius_; }

    protected:
      vec3<FloatType> center_;
      FloatType radius_;
  };

  //! Minimum covering sphere of a set of 3-dimensional points.
  /*! Algorithm due to

      C.L. Lawson, SIAM Review, Vol. 7, No. 3 (Jul., 1965), 415-416.

      This algorithm is compact and numerically stable but converges
      very slowly when many of the points are near the surface of the
      sphere. For a more involved but faster algorithm see:

      http://www.mel.nist.gov/msidlibrary/doc/hopp95.pdf

      Theodore H. Hopp, Charles P. Reeve;
      An Algorithm for Computing the Minimum Covering Sphere in Any Dimension;
      National Institute of Standards and Technology,
      Gaithersburg, MD 20899-0001;
      (NIST IR 5831).

      Note that Lawson's algorithm also generalizes to sets of N-dimensional
      points, but this implementation is restricted to the 3-dimensional
      case.
   */
  template <typename FloatType=double>
  class minimum_covering_sphere_3d : public sphere_3d<FloatType>
  {
    public:
      static const vec3<FloatType> default_center_if_no_points;

      //! Default constructor. Some data members are not initialized!
      minimum_covering_sphere_3d() {}

      //! Execution of Lawson's algorithm.
      /*! The iterative algorithm is terminated if

          tau - sigma < tau * epsilon,

          using Lawson's notation.

          if radius_if_one_or_no_points is less than 0 an exception
          is thrown if points.size() == 0.
       */
      minimum_covering_sphere_3d(
        af::const_ref<vec3<FloatType> > const& points,
        FloatType const& epsilon=1.e-6,
        FloatType const& radius_if_one_or_no_points=1,
        vec3<FloatType> const& center_if_no_points=default_center_if_no_points)
      :
        n_iterations_(0)
      {
        SCITBX_ASSERT(points.size() > 0 || radius_if_one_or_no_points >= 0);
        SCITBX_ASSERT(epsilon > 0);
        if (points.size() == 0) {
          this->center_ = center_if_no_points;
          this->radius_ = radius_if_one_or_no_points;
          return;
        }
        if (points.size() == 1) {
          this->center_ = points[0];
          this->radius_ = radius_if_one_or_no_points;
          return;
        }
        typedef FloatType f_t;
        std::vector<f_t> weights(points.size(), 1./points.size());
        while (true) {
          this->center_.fill(0);
          for(std::size_t i=0;i<points.size();i++) {
            this->center_ += weights[i] * points[i];
          }
          f_t sum_w_r_sq = 0;
          f_t sum_w_r = 0;
          this->radius_ = 0;
          for(std::size_t i=0;i<points.size();i++) {
            f_t r_sq = (this->center_ - points[i]).length_sq();
            sum_w_r_sq += weights[i] * r_sq;
            f_t r = std::sqrt(r_sq);
            if (this->radius_ < r) this->radius_ = r;
            weights[i] *= r;
            sum_w_r += weights[i];
          }
          if (this->radius_ - std::sqrt(sum_w_r_sq) < this->radius_*epsilon) {
            break;
          }
          SCITBX_ASSERT(sum_w_r != 0);
          for(std::size_t i=0;i<points.size();i++) {
            weights[i] /= sum_w_r;
          }
          n_iterations_++;
        }
      }

      //! Number of iterations.
      std::size_t
      n_iterations() const { return n_iterations_; }

    protected:
      std::size_t n_iterations_;
  };

  template <typename FloatType>
  const
  vec3<FloatType>
  minimum_covering_sphere_3d<FloatType>::default_center_if_no_points(0,0,0);

}} // namespace scitbx::math

#endif // SCITBX_MATH_MINIMUM_COVERING_SPHERE
