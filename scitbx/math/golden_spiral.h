#ifndef SCITBX_MATH_GOLDEN_SPIRAL_H
#define SCITBX_MATH_GOLDEN_SPIRAL_H

#include <boost/math/constants/constants.hpp>

#include <vector>
#include <cmath>

namespace scitbx { namespace math {

  //! Unit vectors spread over the sphere on a golden-angle spiral.
  /*! The heights along y step evenly by 2 / count from next to the -y pole,
      half a step in from each pole, and the azimuth about y advances by the
      golden angle, pi * (3 - sqrt(5)). Each point stands for 4 pi / count of
      the sphere (unit_area).
   */
  template <typename Vector>
  class golden_spiral
  {
    public:
      typedef Vector vector_type;
      typedef typename Vector::value_type value_type;
      //! A type of its own, so that its Python class does not claim
      //! std::vector<Vector>, which other modules convert too.
      struct storage_type : std::vector<vector_type> {};

    private:
      size_t count_;
      value_type unit_area_;
      storage_type points_;

      static const value_type pitch;

    public:
      golden_spiral(size_t count)
      :
        count_(count),
        unit_area_(4.0 / count * boost::math::constants::pi<value_type>())
      {
        points_.reserve(count);
        const value_type step = 2.0 / count;
        const value_type offset = step / 2.0;
        for (size_t index = 0; index < count; ++index) {
          const value_type y = step * index - 1.0 + offset;
          const value_type r = std::sqrt(1.0 - y * y);
          const value_type phi = pitch * index;
          points_.push_back(
            vector_type(std::cos(phi) * r, y, std::sin(phi) * r));
        }
      }

      const value_type&
      unit_area() const { return unit_area_; }

      const size_t&
      count() const { return count_; }

      const storage_type&
      points() const { return points_; }
  };

  template <typename Vector>
  const typename golden_spiral<Vector>::value_type
  golden_spiral<Vector>::pitch =
    boost::math::constants::pi<typename golden_spiral<Vector>::value_type>()
    * (3.0 - std::sqrt(5.0));

}} // namespace scitbx::math

#endif // SCITBX_MATH_GOLDEN_SPIRAL_H
