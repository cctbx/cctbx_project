#ifndef SCITBX_MATH_SUPERPOSE_H
#define SCITBX_MATH_SUPERPOSE_H

#include <scitbx/vec3.h>
#include <scitbx/mat3.h>
#include <scitbx/array_family/shared.h>
#include <scitbx/array_family/ref.h>
#include <scitbx/array_family/tiny.h>
#include <scitbx/array_family/shared_algebra.h>
#include <scitbx/array_family/ref_algebra.h>
#include <scitbx/array_family/shared_reductions.h>
#include <scitbx/array_family/ref_reductions.h>
#include <scitbx/array_family/accessors/c_grid.h>
#include <scitbx/math/eigensystem.h>

#include <boost/lambda/lambda.hpp>
#include <boost/lambda/bind.hpp>

#include <algorithm>

using namespace scitbx;

namespace scitbx {  namespace math {  namespace superpose {

  class superposition
  {
    private:
      mat3< double > rotation_;
      vec3< double > translation_;

    public:
      superposition(
          const af::const_ref< vec3< double > >& reference,
          const af::const_ref< vec3< double > >& moving,
          const mat3< double > (*rotation_calculation)(
              const af::const_ref< vec3< double > >&,
              const af::const_ref< vec3< double > >& ) =
            superposition::kearsley_rotation );
      ~superposition();

      // Accessors
      const mat3< double >& get_rotation() const;
      const vec3< double >& get_translation() const;

      // Rotation calculation alternatives
      static const mat3< double > kearsley_rotation(
          const af::const_ref< vec3< double > >& reference,
          const af::const_ref< vec3< double > >& moving );

      // Helper functions
      static const af::tiny< af::shared< double >, 3 >
      decompose_array_of_vec3(
          const af::const_ref< vec3< double > >& array );

      static const mat3< double > quaternion_to_rotmat(
          double q1,
          double q2,
          double q3,
          double q4 );
  };


  class least_squares_fit
  {
    private:
      superposition superposition_;
      af::shared< vec3< double > > reference_sites_;
      af::shared< vec3< double > > other_sites_best_fit_;

    public:
      least_squares_fit(
          const af::shared< vec3< double > >& reference_sites,
          const af::shared< vec3< double > >& moving_sites );
      ~least_squares_fit();

      // Accessors
      const mat3< double >& get_rotation() const;
      const vec3< double >& get_translation() const;
      const af::shared< vec3< double > >& other_sites_best_fit() const;

      // Methods
      double other_sites_rmsd() const;
  };


//
// class superposition
//


// Constructor
superposition::superposition(
    const af::const_ref< vec3< double > >& reference_coords,
    const af::const_ref< vec3< double > >& moving_coords,
    const mat3< double > (*rotation_calculation)(
        const af::const_ref< vec3< double > >&,
        const af::const_ref< vec3< double > >& ) )
{
  vec3< double > reference_mean = af::mean( reference_coords );
  vec3< double > moving_mean = af::mean( moving_coords );

  rotation_ = rotation_calculation(
      ( reference_coords - reference_mean ).const_ref(),
      ( moving_coords - moving_mean ).const_ref() );

  translation_ = reference_mean - rotation_ * moving_mean;
}


// Destructor
superposition::~superposition()
{
}


// Accessors
const mat3< double >&
superposition::get_rotation() const
{
  return rotation_;
}


const vec3< double >&
superposition::get_translation() const
{
  return translation_;
}


// Rotation calculation altenatives
const mat3< double >
superposition::kearsley_rotation(
    const af::const_ref< vec3< double > >& reference,
    const af::const_ref< vec3< double > >& moving )
{
  assert ( reference.size() == moving.size() );

  af::tiny< af::shared< double >, 3 > diff(
      decompose_array_of_vec3( ( reference - moving ).const_ref() ) );
  af::tiny< af::shared< double >, 3 > summ(
      decompose_array_of_vec3( ( reference + moving ).const_ref() ) );

  // Fill up 4*4 matrix
  double matrix_memory[16];

  af::ref< double, af::c_grid< 2 > > matrix(
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
  math::eigensystem::real_symmetric< double > eigenvalues( matrix );

  // The eigenvectors corresponding to the lowest four eignevalues form
  // a unit quaternion
  return quaternion_to_rotmat(
      eigenvalues.vectors()[12],
      eigenvalues.vectors()[13],
      eigenvalues.vectors()[14],
      eigenvalues.vectors()[15] );
}


// Helper methods
const af::tiny< af::shared< double >, 3 >
superposition::decompose_array_of_vec3(
    const af::const_ref< vec3< double > >& array_of_vec3 )
{
  af::tiny< af::shared< double >, 3 > array_of_coordinates;
  size_t array_length( array_of_vec3.size() );

  for ( int coord = 0; coord < 3; ++coord )
  {
    array_of_coordinates[ coord ].reserve( array_length );
    std::transform(
        array_of_vec3.begin(),
        array_of_vec3.end(),
        std::back_inserter( array_of_coordinates[ coord ] ),
        boost::lambda::ret< double >( boost::lambda::_1[ coord ] ) );
  }

  return array_of_coordinates;
}


const mat3< double >
superposition::quaternion_to_rotmat( double q1, double q2, double q3, double q4 )
{
  return mat3< double >(
      q1*q1 + q2*q2 - q3*q3 - q4*q4,      // Element ( 1, 1 )
      2.0 * ( q2*q3 + q1*q4 ),            // Element ( 1, 2 )
      2.0 * ( q2*q4 - q1*q3 ),            // Element ( 1, 3 )
      2.0 * ( q2*q3 - q1*q4 ),            // Element ( 2, 1 )
      q1 * q1 + q3*q3 - q2*q2 - q4*q4,    // Element ( 2, 2 )
      2.0 * ( q3*q4 + q1*q2 ),            // Element ( 2, 3 )
      2.0 * ( q2*q4 + q1*q3 ),            // Element ( 3, 1 )
      2.0 * ( q3*q4 - q1*q2 ),            // Element ( 3, 2 )
      q1*q1 + q4*q4 - q2*q2 - q3*q3 );    // Element ( 3, 3 )
}


//
// class least_squares_fit
//


// Constructor
least_squares_fit::least_squares_fit(
    const af::shared< vec3< double > >& reference_sites,
    const af::shared< vec3< double > >& moving_sites )
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
      boost::lambda::ret< vec3< double > >(
          boost::lambda::ret< vec3< double > >(
              superposition_.get_rotation() * boost::lambda::_1 )
              + superposition_.get_translation() ) );
}


// Destructor
least_squares_fit::~least_squares_fit()
{
}


// Accessors
const mat3< double >&
least_squares_fit::get_rotation() const
{
  return superposition_.get_rotation();
}


const vec3< double >&
least_squares_fit::get_translation() const
{
  return superposition_.get_translation();
}


const af::shared< vec3< double > >&
least_squares_fit::other_sites_best_fit() const
{
  return other_sites_best_fit_;
}


// Methods
double
least_squares_fit::other_sites_rmsd() const
{
  af::shared< double > length_difference_sq;
  length_difference_sq.reserve( reference_sites_.size() );

  std::transform(
      reference_sites_.begin(),
      reference_sites_.end(),
      other_sites_best_fit_.begin(),
      std::back_inserter( length_difference_sq ),
      boost::lambda::bind< double >(
          &vec3< double >::length_sq,
          boost::lambda::ret< vec3< double > >(
              boost::lambda::_1 - boost::lambda::_2 ) ) );

  return sqrt( af::sum( length_difference_sq ) / length_difference_sq.size() );
}

}  }  } // namespace scitbx::math::superpose

#endif // SCITBX_MATH_SUPERPOSE_H
