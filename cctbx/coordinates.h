#ifndef CCTBX_COORDINATES_H
#define CCTBX_COORDINATES_H

#include <scitbx/vec3.h>
#include <cctbx/import_scitbx_af.h>

namespace cctbx {

  template <typename FloatType>
  class fractional;

  template <typename IntType>
  class grid_point;

  //! Class for cartesian (orthogonal, real) coordinates.
  /*! The template parameter FloatType should be a floating point type
      (e.g. float or double).
      <p>
      See also: class fractional
   */
  template <typename FloatType = double>
  class cartesian : public scitbx::vec3<FloatType>
  {
    public:
      typedef scitbx::vec3<FloatType> base_type;

      //! Default constructor: elements are not initialized!
      cartesian() {}

      //! The elements of the coordinate vector are copied from v.
      template <typename OtherFloatType>
      cartesian(af::tiny_plain<OtherFloatType, 3> const& v)
      {
        for(std::size_t i=0;i<3;i++) this->elems[i] = v[i];
      }

      //! The elements of the coordinate vector are copied from xyz.
      explicit
      cartesian(const FloatType* xyz)
      {
        for(std::size_t i=0;i<3;i++) this->elems[i] = xyz[i];
      }

      //! The elements of the coordinate vector are initialized with x,y,z.
      cartesian(FloatType const& x, FloatType const& y, FloatType const& z)
      {
        this->elems[0] = x; this->elems[1] = y; this->elems[2] = z;
      }

    private:
      // disables construction from fractional
      template <typename OtherFloatType>
      cartesian(fractional<OtherFloatType> const&);
      // disables construction from grid_point
      template <typename OtherFloatType>
      cartesian(grid_point<OtherFloatType> const&);
  };

  //! Class for fractional coordinates.
  /*! The template parameter FloatType should be a floating point type
      (e.g. float or double).
      <p>
      See also: class cartesian
   */
  template <typename FloatType = double>
  class fractional : public scitbx::vec3<FloatType>
  {
    public:
      typedef scitbx::vec3<FloatType> base_type;

      //! Default constructor: elements are not initialized!
      fractional() {}

      //! The elements of the coordinate vector are copied from v.
      template <typename OtherFloatType>
      fractional(af::tiny_plain<OtherFloatType, 3> const& v)
      {
        for(std::size_t i=0;i<3;i++) this->elems[i] = v[i];
      }

      //! The elements of the coordinate vector are copied from xyz.
      template <typename OtherFloatType>
      explicit
      fractional(const OtherFloatType* xyz)
      {
        for(std::size_t i=0;i<3;i++) this->elems[i] = xyz[i];
      }

      //! The elements of the coordinate vector are initialized with x,y,z.
      fractional(FloatType const& x, FloatType const& y, FloatType const& z)
      {
        this->elems[0] = x; this->elems[1] = y; this->elems[2] = z;
      }

      /*! \brief Apply modulus operation such that 0.0 <= x < 1.0
          for all elements of the coordinate vector.
       */
      fractional mod_positive() const
      {
        fractional result;
        for(std::size_t i=0;i<3;i++) {
          result[i] = std::fmod(this->elems[i], 1.);
          while (result[i] <  0.) result[i] += 1.;
          while (result[i] >= 1.) result[i] -= 1.;
        }
        return result;
      }

      /*! \brief Apply modulus operation such that -0.5 < x <= 0.5
          for all elements of the coordinate vector.
       */
      fractional mod_short() const
      {
        fractional result;
        for(std::size_t i=0;i<3;i++) {
          result[i] = std::fmod(this->elems[i], 1.);
          if      (result[i] <= -.5) result[i] += 1.;
          else if (result[i] >   .5) result[i] -= 1.;
        }
        return result;
      }

      scitbx::vec3<int>
      unit_shifts() const
      {
        scitbx::vec3<int> result;
        for(std::size_t i=0;i<3;i++) {
          if (this->elems[i] >= 0.) result[i] = int(this->elems[i] + 0.5);
          else                      result[i] = int(this->elems[i] - 0.5);
        }
        return result;
      }

    private:
      // disables construction from cartesian
      template <typename OtherFloatType>
      fractional(cartesian<OtherFloatType> const&);
      // disables construction from grid_point
      template <typename OtherFloatType>
      fractional(grid_point<OtherFloatType> const&);
  };

  //! Class for grid_point coordinates.
  /*! The template parameter IntType should be an integral type
      (e.g. signed int or signed long).
      <p>
      See also: class fractional, class cartesian
   */
  template <typename IntType = signed long>
  class grid_point : public scitbx::vec3<IntType>
  {
    public:
      typedef scitbx::vec3<IntType> base_type;

      //! Default constructor: elements are not initialized!
      grid_point() {}

      //! The elements of the coordinate vector are copied from v.
      template <typename OtherIntType>
      grid_point(af::tiny_plain<OtherIntType, 3> const& v)
      {
        for(std::size_t i=0;i<3;i++) this->elems[i] = v[i];
      }

      //! The elements of the coordinate vector are copied from xyz.
      explicit
      grid_point(const IntType* xyz)
      {
        for(std::size_t i=0;i<3;i++) this->elems[i] = xyz[i];
      }

      //! The elements of the coordinate vector are initialized with x,y,z.
      grid_point(IntType const& x, IntType const& y, IntType const& z)
      {
        this->elems[0] = x; this->elems[1] = y; this->elems[2] = z;
      }

    private:
      // disables construction from fractional
      template <typename OtherIntType>
      grid_point(fractional<OtherIntType> const&);
      // disables construction from cartesian
      template <typename OtherIntType>
      grid_point(cartesian<OtherIntType> const&);
  };

} // namespace cctbx

#endif // CCTBX_COORDINATES_H
