#ifndef RSTBX_A_G_CONVERSION_H
#define RSTBX_A_G_CONVERSION_H

#include <scitbx/array_family/flex_types.h>
#include <scitbx/array_family/shared.h>
#include <scitbx/vec3.h>
#include <scitbx/mat3.h>
#include <cctbx/crystal_orientation.h>
#include <scitbx/sym_mat3.h>
#include <cctbx/uctbx.h>

namespace rstbx {
typedef scitbx::af::shared<double> farray;
typedef scitbx::af::shared<int>    iarray;
typedef scitbx::af::shared<scitbx::vec3<double> > vec3array;
typedef scitbx::vec3<double>                      vec3;
typedef scitbx::mat3<double>                      mat3;
typedef scitbx::af::shared<cctbx::miller::index<> > marray;
//! Shorthand for default sym_mat3 type in unit cell toolbox.
typedef scitbx::sym_mat3<double> uc_sym_mat3;

namespace symmetry {

struct AG { // convert orientation matrix A to metrical matrix g & reverse
  /*The general orientation matrix A is re-expressed in terms of the
    upper-triangular fractionalization matrix F by means of the following
    transformation:
                           F = (D * C * B * A).transpose()
    where D,C,B are three rotation matrices.
  */

  cctbx::crystal_orientation orientation;
  double phi,psi,theta; //in radians
  mat3 B,C,D,F;
  uc_sym_mat3 G;

  void forward(cctbx::crystal_orientation const& ori){
    orientation = ori;
    mat3 A( ori.reciprocal_matrix() );

    phi = std::atan2(A[2],-A[8]);
    B = mat3(std::cos(phi),0.,std::sin(phi),
                    0.,           1.,0.,
             -std::sin(phi),0.,std::cos(phi));
    mat3 BA (B * A);

    psi = std::atan2(-BA[5],BA[8]);
    C = mat3(1.,0.,0.,
             0., std::cos(psi),std::sin(psi),
             0.,-std::sin(psi),std::cos(psi));
    mat3 CBA (C * BA);

    theta = std::atan2(-CBA[1],CBA[4]);
    D = mat3(  std::cos(theta),std::sin(theta),0.,
              -std::sin(theta),std::cos(theta),0.,
               0.,         0.,        1.);
    F = (D * CBA).transpose();

    mat3 G9 (A.transpose()*A); /*  3x3 form of metrical matrix
      Note 3/23/2009 we are using the reciprocal space metrical matrix
       (a*.a*, a*.b*, a*.c*, a*.b*, b*.b*, b*.c*, a*.c*, b*.c*, c*.c*)
    */

    G = uc_sym_mat3(G9[0],G9[4],G9[8],G9[1],G9[2],G9[5]);
    // (a*.a*, b*.b*, c*.c*, a*.b*, a*.c*, b*.c*)
  }

  void setAngles(double const& phi, double const& psi, double const& theta){
    B = mat3(std::cos(phi),0.,std::sin(phi),
                  0.,           1.,0.,
            -std::sin(phi),0.,std::cos(phi));
    C = mat3(1.,0.,0.,
             0.,std::cos(psi),std::sin(psi),
             0.,-std::sin(psi),std::cos(psi));
    D = mat3( std::cos(theta),std::sin(theta),0.,
             -std::sin(theta),std::cos(theta),0.,
                  0.,         0.,        1.);

  }

  void validate_and_setG(uc_sym_mat3 const&g){
    double g0=g[0],g1=g[1],g2=g[2],g3=g[3],g4=g[4],g5=g[5];

    //Note:g0 = a*.a*  g1 = b*.b*  g2 = c*.c*  g3 = a*.b*  g4 = a*.c*  g5 = b*.c*

    if (g2 <= 0.){
      throw scitbx::error("g2 <= 0."); g2 = 1.E-7; }
    double cstrz = std::sqrt(g2);
    double bstrz = g5/cstrz;
    double astrz = g4/cstrz;
    if (g1-bstrz*bstrz <= 0.){
      throw scitbx::error("g1-bstrz*bstrz <= 0."); g1 = bstrz*bstrz + 1.E-7; }
    double bstry = std::sqrt(g1-bstrz*bstrz);
    double astry = (g3-astrz*bstrz)/bstry;
    if (g0 - astry*astry -astrz*astrz <= 0.){
      throw scitbx::error("g0 - astry*astry -astrz*astrz <= 0."); g0 = astry*astry + astrz*astrz + 1.E-7; }
    //double astrx = std::sqrt(g0 - astry*astry -astrz*astrz);

    G = uc_sym_mat3(g0,g1,g2,g3,g4,g5);
    cctbx::uctbx::unit_cell ersatz_uc ( cctbx::uctbx::unit_cell(G).reciprocal() );

    if ( ersatz_uc.volume() <= 40.){ throw SCITBX_ERROR(
      "Metrical matrix g is expected to be in the reciprocal setting;this appears to be direct space");
    }

  };

  mat3 back()const{

    cctbx::uctbx::unit_cell ersatz_uc ( cctbx::uctbx::unit_cell(G).reciprocal() );

    // ersatz F is the fract. matrix, PDB Convention, compatible with CCTBX
    mat3 ersatzF (ersatz_uc.fractionalization_matrix());

    // Fback is the lower-triangular matrix compatible with the Rsymop
    // paper, equation (3) = {{a*_x,0,0},{a*_y,b*_y,0},{a*_z,b*_z,c*_z}}
    mat3 Fback (ersatzF.transpose());

    return (B.inverse() * C.inverse() * D.inverse() * Fback);

  }

  cctbx::crystal_orientation
  back_as_orientation() const{
    return cctbx::crystal_orientation( back(), true );
  }
};
}
} //namespace rstbx
#endif// RSTBX_A_G_CONVERSION_H
