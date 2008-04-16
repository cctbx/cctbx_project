#ifndef SCITBX_MATH_EULER_ANGLES_H
#define SCITBX_MATH_EULER_ANGLES_H

#include <scitbx/constants.h>
#include <scitbx/vec3.h>
#include <scitbx/mat3.h>

namespace scitbx {  namespace math {  namespace euler_angles {
    
template< typename FloatType >
const scitbx::mat3< FloatType >
xyz_matrix( const FloatType x_deg, const FloatType y_deg, const FloatType z_deg )
{
    FloatType x_rad( scitbx::deg_as_rad( x_deg ) );
    FloatType y_rad( scitbx::deg_as_rad( y_deg ) );
    FloatType z_rad( scitbx::deg_as_rad( z_deg ) );
    
    FloatType sin_x( std::sin( x_rad ) );
    FloatType sin_y( std::sin( y_rad ) );
    FloatType sin_z( std::sin( z_rad ) );
    
    FloatType cos_x( std::cos( x_rad ) );
    FloatType cos_y( std::cos( y_rad ) );
    FloatType cos_z( std::cos( z_rad ) );
    
    return scitbx::mat3< FloatType >(
             cos_y * cos_z,
            -cos_y * sin_z,
             sin_y,
            
             cos_x * sin_z  +  sin_x * sin_y * cos_z,
             cos_x * cos_z  -  sin_x * sin_y * sin_z,
            -sin_x * cos_y,
            
             sin_x * sin_z  -  cos_x * sin_y * cos_z,
             sin_x * cos_z  +  cos_x * sin_y * sin_z,
             cos_x * cos_y );
}


template< typename FloatType >
const scitbx::vec3< FloatType >
xyz_angles( const scitbx::mat3< FloatType >& m, const FloatType eps = 1.0e-12 )
{
    if ( m( 0, 2 ) > ( 1 - eps ) )
    {
        return scitbx::vec3< FloatType >(
                scitbx::rad_as_deg( std::atan2( m(2, 1), m(1, 1) ) ),
                scitbx::rad_as_deg( scitbx::constants::pi_2 ),
                scitbx::rad_as_deg( 0.0 ) );
    }
    
    else if ( m( 0, 2 ) < ( -1 + eps ) )
    {
        return scitbx::vec3< FloatType >(
                scitbx::rad_as_deg( std::atan2( m(2, 1), m(1, 1) ) ),
                scitbx::rad_as_deg( -scitbx::constants::pi_2 ),
                scitbx::rad_as_deg( 0.0 ) );
    }
    
    else
    {
        return scitbx::vec3< FloatType >(
                scitbx::rad_as_deg( std::atan2( -m(1, 2), m(2, 2) ) ),
                scitbx::rad_as_deg( std::asin( m(0, 2) ) ),
                scitbx::rad_as_deg( std::atan2( -m(0, 1), m(0, 0) ) ) );
    }
}


template< typename FloatType >
const scitbx::mat3< FloatType >
yzx_matrix( const FloatType y_deg, const FloatType z_deg, const FloatType x_deg )
{
    FloatType y_rad( scitbx::deg_as_rad( y_deg ) );
    FloatType z_rad( scitbx::deg_as_rad( z_deg ) );
    FloatType x_rad( scitbx::deg_as_rad( x_deg ) );
    
    FloatType sin_y( std::sin( y_rad ) );
    FloatType sin_z( std::sin( z_rad ) );
    FloatType sin_x( std::sin( x_rad ) );
    
    FloatType cos_y( std::cos( y_rad ) );
    FloatType cos_z( std::cos( z_rad ) );
    FloatType cos_x( std::cos( x_rad ) );
    
    return scitbx::mat3< FloatType >(
             cos_y * cos_z,
             sin_x * sin_y  -  cos_x * cos_y * sin_z,
             cos_x * sin_y  +  cos_y * sin_x * sin_z,
            
             sin_z,
             cos_x * cos_z,
            -cos_z * sin_x,
            
            -cos_z * sin_y,
             cos_y * sin_x  +  cos_x * sin_y * sin_z,
             cos_x * cos_y  -  sin_x * sin_y * sin_z );
}


template< typename FloatType >
const scitbx::vec3< FloatType >
yzx_angles( const scitbx::mat3< FloatType >& m, const FloatType eps = 1.0e-12 )
{
    if ( m( 1, 0 ) > ( 1 - eps ) )
    {
        return scitbx::vec3< FloatType >(
                scitbx::rad_as_deg( std::atan2( m(0, 2), m(2, 2) ) ),
                scitbx::rad_as_deg( scitbx::constants::pi_2 ),
                scitbx::rad_as_deg( 0.0 ) );
    }
    
    else if ( m( 1, 0 ) < ( -1 + eps ) )
    {
        return scitbx::vec3< FloatType >(
                scitbx::rad_as_deg( std::atan2( m(0, 2), m(2, 2) ) ),
                scitbx::rad_as_deg( -scitbx::constants::pi_2 ),
                scitbx::rad_as_deg( 0.0 ) );
    }
    
    else
    {
        return scitbx::vec3< FloatType >(
                scitbx::rad_as_deg( std::atan2( -m(2, 0), m(0, 0) ) ),
                scitbx::rad_as_deg( std::asin( m(1, 0) ) ),
                scitbx::rad_as_deg( std::atan2( -m(1, 2), m(1, 1) ) ) );
    }
}


template< typename FloatType >
const scitbx::mat3< FloatType >
zyz_matrix( const FloatType z1_deg, const FloatType y2_deg, const FloatType z3_deg )
{
    FloatType z1_rad( scitbx::deg_as_rad( z1_deg ) );
    FloatType y2_rad( scitbx::deg_as_rad( y2_deg ) );
    FloatType z3_rad( scitbx::deg_as_rad( z3_deg ) );
    
    FloatType sin_z1( std::sin( z1_rad ) );
    FloatType sin_y2( std::sin( y2_rad ) );
    FloatType sin_z3( std::sin( z3_rad ) );
    
    FloatType cos_z1( std::cos( z1_rad ) );
    FloatType cos_y2( std::cos( y2_rad ) );
    FloatType cos_z3( std::cos( z3_rad ) );
    
    return scitbx::mat3< FloatType >(
             cos_y2 * cos_z1 * cos_z3  -  sin_z1 * sin_z3,
            -cos_z3 * sin_z1  -  cos_y2 * cos_z1 * sin_z3,
             cos_z1 * sin_y2,
            
             cos_y2 * cos_z3 * sin_z1  +  cos_z1 * sin_z3,
             cos_z1 * cos_z3  -  cos_y2 * sin_z1 * sin_z3,
             sin_y2 * sin_z1,
            
            -cos_z3 * sin_y2,
             sin_y2 * sin_z3,
             cos_y2 );
}


template< typename FloatType >
const scitbx::vec3< FloatType >
zyz_angles( const scitbx::mat3< FloatType >& m, const FloatType eps = 1.0e-12 )
{
    if ( ( 1 - eps ) < m( 2, 2 ) )
    {
        return scitbx::vec3< FloatType >(
                scitbx::rad_as_deg( std::atan2( -m(0, 1), m(1, 1) ) ),
                scitbx::rad_as_deg( 0.0 ),
                scitbx::rad_as_deg( 0.0 ) );
    }
    
    else if ( m( 2, 2 ) < ( -1 + eps ) )
    {
        return scitbx::vec3< FloatType >(
                scitbx::rad_as_deg( std::atan2( -m(0, 1), m(1, 1) ) ),
                scitbx::rad_as_deg( scitbx::constants::pi ),
                scitbx::rad_as_deg( 0.0 ) );
    }
    
    else
    {
        return scitbx::vec3< FloatType >(
                scitbx::rad_as_deg( std::atan2( m(1, 2), m(0, 2) ) ),
                scitbx::rad_as_deg( std::acos( m(2, 2) ) ),
                scitbx::rad_as_deg( std::atan2( m(2, 1), -m(2, 0) ) ) );
    }
}
    
}  }  } //namespace scitbx::math::euler_angles

#endif /*SCITBX_MATH_EULER_ANGLES_H*/
