#ifndef SCITBX_RIGID_BODY_SPATIAL_LIB_H
#define SCITBX_RIGID_BODY_SPATIAL_LIB_H

#include <scitbx/rotr3.h>
#include <scitbx/array_family/versa_matrix.h>
#include <scitbx/rigid_body/matrix_helpers.h>

namespace scitbx { namespace rigid_body {

//! See essence/spatial_lib.py
namespace spatial_lib {

  template <typename ElementType>
  af::tiny<ElementType, 6>
  as_tiny_6(
    vec3<ElementType> const& a,
    vec3<ElementType> const& l)
  {
    af::tiny<ElementType, 6> result;
    for(unsigned i=0;i<3;i++) result[i] = a[i];
    for(unsigned i=0;i<3;i++) result[i+3] = l[i];
    return result;
  }

  //! RBDA Tab. 2.2, p. 23.
  /*! Spatial coordinate transform (rotation around origin).
      Calculates the coordinate transform matrix from A to B coordinates
      for spatial motion vectors, in which frame B is rotated relative to
      frame A.
   */
  template <typename FloatType>
  af::versa<FloatType, af::mat_grid>
  xrot(
    mat3<FloatType> const& e)
  {
    FloatType const coeffs[] = {
      e[0], e[1], e[2], 0, 0, 0,
      e[3], e[4], e[5], 0, 0, 0,
      e[6], e[7], e[8], 0, 0, 0,
      0, 0, 0,          e[0], e[1], e[2],
      0, 0, 0,          e[3], e[4], e[5],
      0, 0, 0,          e[6], e[7], e[8]};
    return af::versa_mat_grid(coeffs, 6, 6);
  }

  //! RBDA Tab. 2.2, p. 23.
  /*! Spatial coordinate transform (translation of origin).
      Calculates the coordinate transform matrix from A to B coordinates
      for spatial motion vectors, in which frame B is translated by an
      amount r (3D vector) relative to frame A.
   */
  template <typename FloatType>
  af::versa<FloatType, af::mat_grid>
  xtrans(
    vec3<FloatType> const& r)
  {
    FloatType const coeffs[] = {
          1,     0,     0, 0, 0, 0,
          0,     1,     0, 0, 0, 0,
          0,     0,     1, 0, 0, 0,
          0,  r[2], -r[1], 1, 0, 0,
      -r[2],     0,  r[0], 0, 1, 0,
       r[1], -r[0],     0, 0, 0, 1};
    return af::versa_mat_grid(coeffs, 6, 6);
  }

  //! RBDA Eq. 2.28, p. 22.
  template <typename FloatType>
  af::versa<FloatType, af::mat_grid>
  cb_as_spatial_transform(
    rotr3<FloatType> const& cb)
  {
    return af::matrix_multiply(
      xrot(cb.r).const_ref(),
      xtrans(-cb.t * cb.r).const_ref());
  }

  //! RBDA Eq. 2.31, p. 25.
  /*! Spatial cross-product operator (motion).
      Calculates the 6x6 matrix such that the expression crm(v)*m is the
      cross product of the spatial motion vectors v and m.
   */
  template <typename FloatType>
  af::versa<FloatType, af::mat_grid>
  crm(
    af::tiny<FloatType, 6> const& v)
  {
    FloatType const coeffs[] = {
          0, -v[2],  v[1],     0,     0,     0,
       v[2],     0, -v[0],     0,     0,     0,
      -v[1],  v[0],     0,     0,     0,     0,
          0, -v[5],  v[4],     0, -v[2],  v[1],
       v[5],     0, -v[3],  v[2],     0, -v[0],
      -v[4],  v[3],     0, -v[1],  v[0],     0};
    return af::versa_mat_grid(coeffs, 6, 6);
  }

  //! RBDA Eq. 2.32, p. 25.
  /*! Spatial cross-product operator (force).
      Calculates the 6x6 matrix such that the expression crf(v)*f is the
      cross product of the spatial motion vector v with the spatial force
      vector f.
   */
  template <typename FloatType>
  af::versa<FloatType, af::mat_grid>
  crf(
    af::tiny<FloatType, 6> const& v)
  {
    return -af::matrix_transpose(crm(v).const_ref());
  }

  //! RBDA Eq. 2.63, p. 33.
  /*! Spatial rigid-body inertia from mass, CoM and rotational inertia.
      Calculates the spatial inertia matrix of a rigid body from its
      mass, centre of mass (3D vector) and rotational inertia (3x3 matrix)
      about its centre of mass.
   */
  template <typename FloatType>
  af::versa<FloatType, af::mat_grid>
  mci(
    FloatType const& m,
    vec3<FloatType> const& c,
    sym_mat3<FloatType> const& i)
  {
    typedef FloatType ft;
    mat3<ft> cx = mat3<ft>::cross_product_matrix(c);
    sym_mat3<ft> imcc = i + m * cx.self_times_self_transpose();
    ft const coeffs[] = {
      imcc[0], imcc[3], imcc[4],  m*cx[0], m*cx[1], m*cx[2],
      imcc[3], imcc[1], imcc[5],  m*cx[3], m*cx[4], m*cx[5],
      imcc[4], imcc[5], imcc[2],  m*cx[6], m*cx[7], m*cx[8],
      m*cx[0], m*cx[3], m*cx[6],  m, 0, 0,
      m*cx[1], m*cx[4], m*cx[7],  0, m, 0,
      m*cx[2], m*cx[5], m*cx[8],  0, 0, m};
    return af::versa_mat_grid(coeffs, 6, 6);
  }

  //! RBDA Eq. 2.67, p. 35.
  template <typename FloatType>
  FloatType
  kinetic_energy(
    af::const_ref<FloatType, af::mat_grid> const& i_spatial,
    af::tiny<FloatType, 6> const& v_spatial)
  {
    af::tiny<FloatType, 6> iv;
    matrix_mul(iv, i_spatial, v_spatial.const_ref());
    return 0.5 * dot_product(v_spatial, iv);
  }

}}} // namespace scitbx::rigid_body::spatial_lib

#endif // GUARD
