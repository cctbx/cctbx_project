#ifndef SCITBX_MATH_SUPERPOSE_H
#define SCITBX_MATH_SUPERPOSE_H

#include <scitbx/matrix/eigensystem.h>
#include <scitbx/math/r3_rotation.h>
#include <scitbx/array_family/shared_reductions.h>
#include <scitbx/array_family/shared_algebra.h>
#include <scitbx/array_family/ref_algebra.h>

#include <boost/lambda/lambda.hpp>
#include <boost/lambda/bind.hpp>

namespace scitbx {  namespace math {  namespace superpose {

  template <typename FloatType=double>
  class superposition
  {
    private:
      mat3< FloatType > rotation_;
      vec3< FloatType > translation_;

    public:
      superposition(
          const af::const_ref< vec3< FloatType > >& reference,
          const af::const_ref< vec3< FloatType > >& moving,
          const mat3< FloatType > (*rotation_calculation)(
              const af::const_ref< vec3< FloatType > >&,
              const af::const_ref< vec3< FloatType > >& ) =
            superposition::kearsley_rotation );

      // Accessors
      const mat3< FloatType >& get_rotation() const;
      const vec3< FloatType >& get_translation() const;

      // Rotation calculation alternatives
      static const mat3< FloatType > kearsley_rotation(
          const af::const_ref< vec3< FloatType > >& reference,
          const af::const_ref< vec3< FloatType > >& moving );

      // Helper functions
      static const af::tiny< af::shared< FloatType >, 3 >
      decompose_array_of_vec3(
          const af::const_ref< vec3< FloatType > >& array );

      static const mat3< FloatType > quaternion_to_rotmat(
          FloatType const& q1,
          FloatType const& q2,
          FloatType const& q3,
          FloatType const& q4 );
  };


  template <typename FloatType=double>
  class least_squares_fit
  {
    private:
      superposition<FloatType> superposition_;
      af::shared< vec3< FloatType > > reference_sites_;
      af::shared< vec3< FloatType > > other_sites_best_fit_;

    public:
      least_squares_fit(
          const af::shared< vec3< FloatType > >& reference_sites,
          const af::shared< vec3< FloatType > >& moving_sites );

      // Accessors
      const mat3< FloatType >& get_rotation() const;
      const vec3< FloatType >& get_translation() const;
      const af::shared< vec3< FloatType > >& other_sites_best_fit() const;

      // Methods
      FloatType other_sites_rmsd() const;
  };


//
// class superposition
//


// Constructor
template <typename FloatType>
superposition<FloatType>::superposition(
    const af::const_ref< vec3< FloatType > >& reference_coords,
    const af::const_ref< vec3< FloatType > >& moving_coords,
    const mat3< FloatType > (*rotation_calculation)(
        const af::const_ref< vec3< FloatType > >&,
        const af::const_ref< vec3< FloatType > >& ) )
{
  vec3< FloatType > reference_mean = af::mean( reference_coords );
  vec3< FloatType > moving_mean = af::mean( moving_coords );

  rotation_ = rotation_calculation(
      ( reference_coords - reference_mean ).const_ref(),
      ( moving_coords - moving_mean ).const_ref() );

  translation_ = reference_mean - rotation_ * moving_mean;
}


// Accessors
template <typename FloatType>
const mat3< FloatType >&
superposition<FloatType>::get_rotation() const
{
  return rotation_;
}


template <typename FloatType>
const vec3< FloatType >&
superposition<FloatType>::get_translation() const
{
  return translation_;
}


// Rotation calculation altenatives
template <typename FloatType>
const mat3< FloatType >
superposition<FloatType>::kearsley_rotation(
    const af::const_ref< vec3< FloatType > >& reference,
    const af::const_ref< vec3< FloatType > >& moving )
{
  assert ( reference.size() == moving.size() );

  af::tiny< af::shared< FloatType >, 3 > diff(
      decompose_array_of_vec3( ( reference - moving ).const_ref() ) );
  af::tiny< af::shared< FloatType >, 3 > summ(
      decompose_array_of_vec3( ( reference + moving ).const_ref() ) );

  // Fill up 4*4 matrix
  FloatType matrix_memory[16];

  af::ref< FloatType, af::c_grid< 2 > > matrix(
      matrix_memory,
      af::c_grid< 2 >( 4, 4 ) );

  // Diagonal elements
  matrix( 0, 0 ) = af::sum(
      diff[0] * diff[0] + diff[1] * diff[1] + diff[2] * diff[2] );
  matrix( 1, 1 ) = af::sum(
      diff[0] * diff[0] + summ[1] * summ[1] + summ[2] * summ[2] );
  matrix( 2, 2 ) = af::sum(
      summ[0] * summ[0] + diff[1] * diff[1] + summ[2] * summ[2] );
  matrix( 3, 3 ) = af::sum(
      summ[0] * summ[0] + summ[1] * summ[1] + diff[2] * diff[2] );

  // Off-diagonal elements (NB, matrix is symmetric)
  matrix( 0, 1 ) = af::sum( summ[1] * diff[2] - diff[1] * summ[2] );
  matrix( 1, 0 ) = matrix( 0, 1 );
  matrix( 0, 2 ) = af::sum( diff[0] * summ[2] - summ[0] * diff[2] );
  matrix( 2, 0 ) = matrix( 0, 2 );
  matrix( 0, 3 ) = af::sum( summ[0] * diff[1] - diff[0] * summ[1] );
  matrix( 3, 0 ) = matrix( 0, 3 );

  matrix( 1, 2 ) = af::sum( diff[0] * diff[1] - summ[0] * summ[1] );
  matrix( 2, 1 ) = matrix( 1, 2 );
  matrix( 1, 3 ) = af::sum( diff[0] * diff[2] - summ[0] * summ[2] );
  matrix( 3, 1 ) = matrix( 1, 3 );

  matrix( 2, 3 ) = af::sum( diff[1] * diff[2] - summ[1] * summ[2] );
  matrix( 3, 2 ) = matrix( 2, 3 );

  // Solve characteristic equation
  matrix::eigensystem::real_symmetric< FloatType > eigensystem( matrix );

  // The eigenvectors corresponding to the lowest four eigenvalues form
  // a unit quaternion
  return math::r3_rotation::unit_quaternion_as_matrix(
    -eigensystem.vectors()[12], // rotation in opposite direction
     eigensystem.vectors()[13],
     eigensystem.vectors()[14],
     eigensystem.vectors()[15]);
}


// Helper methods
template <typename FloatType>
const af::tiny< af::shared< FloatType >, 3 >
superposition<FloatType>::decompose_array_of_vec3(
    const af::const_ref< vec3< FloatType > >& array_of_vec3 )
{
  af::tiny< af::shared< FloatType >, 3 > array_of_coordinates;
  size_t array_length( array_of_vec3.size() );

  for ( int coord = 0; coord < 3; ++coord )
  {
    array_of_coordinates[ coord ].reserve( array_length );
    std::transform(
        array_of_vec3.begin(),
        array_of_vec3.end(),
        std::back_inserter( array_of_coordinates[ coord ] ),
        boost::lambda::ret< FloatType >( boost::lambda::_1[ coord ] ) );
  }

  return array_of_coordinates;
}


//
// class least_squares_fit
//


// Constructor
template <typename FloatType>
least_squares_fit<FloatType>::least_squares_fit(
    const af::shared< vec3< FloatType > >& reference_sites,
    const af::shared< vec3< FloatType > >& moving_sites )
    : superposition_(
          reference_sites.const_ref(),
          moving_sites.const_ref() ),
      reference_sites_( reference_sites )
{
  other_sites_best_fit_.reserve( reference_sites.size() );

  // Double return-type specification needed otherwise BLL compilation fails
  std::transform(
      moving_sites.begin(),
      moving_sites.end(),
      std::back_inserter( other_sites_best_fit_ ),
      boost::lambda::ret< vec3< FloatType > >(
          boost::lambda::ret< vec3< FloatType > >(
              superposition_.get_rotation() * boost::lambda::_1 )
              + superposition_.get_translation() ) );
}


// Accessors
template <typename FloatType>
const mat3< FloatType >&
least_squares_fit<FloatType>::get_rotation() const
{
  return superposition_.get_rotation();
}


template <typename FloatType>
const vec3< FloatType >&
least_squares_fit<FloatType>::get_translation() const
{
  return superposition_.get_translation();
}


template <typename FloatType>
const af::shared< vec3< FloatType > >&
least_squares_fit<FloatType>::other_sites_best_fit() const
{
  return other_sites_best_fit_;
}


// Methods
template <typename FloatType>
FloatType
least_squares_fit<FloatType>::other_sites_rmsd() const
{
  af::shared< FloatType > length_difference_sq;
  length_difference_sq.reserve( reference_sites_.size() );

  std::transform(
      reference_sites_.begin(),
      reference_sites_.end(),
      other_sites_best_fit_.begin(),
      std::back_inserter( length_difference_sq ),
      boost::lambda::bind< FloatType >(
          &vec3< FloatType >::length_sq,
          boost::lambda::ret< vec3< FloatType > >(
              boost::lambda::_1 - boost::lambda::_2 ) ) );

  return sqrt( af::sum( length_difference_sq ) / length_difference_sq.size() );
}

}  }  } // namespace scitbx::math::superpose

#endif // SCITBX_MATH_SUPERPOSE_H
