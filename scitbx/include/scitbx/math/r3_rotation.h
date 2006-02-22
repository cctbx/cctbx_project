#ifndef SCITBX_MATH_R3_ROTATION_H
#define SCITBX_MATH_R3_ROTATION_H

#include <scitbx/matrix/row_echelon.h>
#include <scitbx/mat3.h>
#include <scitbx/constants.h>

namespace scitbx { namespace math { namespace r3_rotation {

  //! Conversion of axis and angle to a rotation matrix.
  /*! http://skal.planet-d.net/demo/matrixfaq.htm
   */
  template <typename FloatType>
  scitbx::mat3<FloatType>
  axis_and_angle_as_matrix(
    scitbx::vec3<FloatType> const& axis,
    FloatType angle,
    bool deg=false,
    FloatType const& min_axis_length=1.e-15)
  {
    SCITBX_ASSERT(min_axis_length > 0);
    FloatType u = axis[0];
    FloatType v = axis[1];
    FloatType w = axis[2];
    FloatType l = std::sqrt(u*u+v*v+w*w);
    if (l < min_axis_length) {
      throw std::runtime_error(
        "Very short rotation axis vector may lead to"
        " numerical instabilities.");
    }
    u /= l;
    v /= l;
    w /= l;
    if (deg) angle *= scitbx::constants::pi_180;
    FloatType c = std::cos(angle);
    FloatType s = std::sin(angle);
    FloatType oc = 1-c;
    FloatType uoc = u*oc;
    FloatType voc = v*oc;
    FloatType woc = w*oc;
    FloatType us = u*s;
    FloatType vs = v*s;
    FloatType ws = w*s;
    return scitbx::mat3<FloatType>(
       c + u*uoc,
     -ws + u*voc,
      vs + u*woc,
      ws + v*uoc,
       c + v*voc,
     -us + v*woc,
     -vs + w*uoc,
      us + w*voc,
       c + w*woc);
  }

  //! Numerically robust computation of rotation axis and angle.
  /*! The rotation axis is determined by solving the system
      R-I=0 using row echelon reduction with full pivoting
      (matrix::row_echelon::full_pivoting_small<>) and
      back-substitution. R is the input rotation matrix,
      I the identity matrix.

      The rotation angle is determined in three steps. First,
      the cross product of the normalized axis vector and
      each basis vector is determined. Each cross product
      yields a vector perpendicular to the rotation axis,
      or the null vector if the rotation axis is parallel
      to one of the basis vectors. The cross product with
      the largest length is selected for the second step.
      Let this vector be "perp". It is multiplied with the
      rotation matrix to yield a second vector "r_perp"
      which is also perpendicular to the rotation axis. The
      angle between perp and r_perp is the rotation angle.
      In the third step, the magnitude of the angle is determined
      via the scalar product. The sign relative to the rotation
      axis vector is determined via axis.dot(perp.cross(r_perp)).

      The algorithm works reliably for all rotation matrices
      even if the rotation angle is very small (the corresponding
      unit test is tst_r3_rotation.py). The accuracy of the results
      is determined purely by the precision of the floating-point
      type. In contrast to quaternion-based alternative algorithms
      (http://skal.planet-d.net/demo/matrixfaq.htm),
      arbitrary cutoff tolerances are not needed at any stage.
      At the same time, the algorithm is very fast.

      If the input matrix is not a rotation matrix the results
      are meaningless. To check for this condition, use the
      as_matrix() member function to reconstruct a rotation
      matrix and compare with the input matrix element-by-element.
      Large discrepancies are an indication of an improper input
      matrix.
   */
  template <typename FloatType=double>
  struct axis_and_angle_from_matrix
  {
    //! Normalized rotation axis.
    scitbx::vec3<FloatType> axis;
    //! Rotation angle in radians.
    FloatType angle_rad;

    //! Default constructor. Some data members are not initialized!
    axis_and_angle_from_matrix() {}

    //! Computation of rotation axis and angle.
    axis_and_angle_from_matrix(scitbx::mat3<FloatType> const& r)
    {
      // obtain axis by solving the system r-i=0
      scitbx::mat3<FloatType> m_work(
        r[0]-1, r[1],   r[2],
        r[3],   r[4]-1, r[5],
        r[6],   r[7],   r[8]-1);
      matrix::row_echelon::full_pivoting_small<double, 3, 3>
        row_echelon_form(
          af::ref<FloatType, af::c_grid<2> >(
            m_work.begin(), af::c_grid<2>(3,3)),
          /*min_abs_pivot*/ 0,
          /*max_rank*/ 2);
      axis = row_echelon_form.back_substitution(
        af::small<double, 3>(row_echelon_form.free_cols.size(), 1));
      FloatType& u = axis[0];
      FloatType& v = axis[1];
      FloatType& w = axis[2];
      // normalize axis
      FloatType uu = u*u;
      FloatType vv = v*v;
      FloatType ww = w*w;
      FloatType l = std::sqrt(uu+vv+ww);
      axis /= l;
      // determine basis vector b leading to maximum length of axis x b,
      // with b=(1,0,0), b=(0,1,0), b=(0,0,1). This leads to a vector
      // perpendicular to axis. This vector is normalized.
      scitbx::vec3<FloatType> perp;
      if (vv >= uu) {
        if (ww >= uu) {
          l = std::sqrt(v*v+w*w);
          perp[0] = 0; perp[1] = w/l; perp[2] = -v/l;
        }
        else {
          l = std::sqrt(u*u+v*v);
          perp[0] = v/l; perp[1] = -u/l; perp[2] = 0;
        }
      }
      else if (ww >= vv) {
        l = std::sqrt(u*u+w*w);
        perp[0] = -w/l; perp[1] = 0; perp[2] = u/l;
      }
      else {
        l = std::sqrt(u*u+v*v);
        perp[0] = v/l; perp[1] = -u/l; perp[2] = 0;
      }
      // apply rotation to the perpendicular vector
      scitbx::vec3<FloatType> r_perp = r * perp;
      // determine the angle
      FloatType cos_angle = perp * r_perp;
      if      (cos_angle <= -1) angle_rad = scitbx::constants::pi;
      else if (cos_angle >=  1) angle_rad = 0;
      else {
        angle_rad = std::acos(cos_angle);
        if (axis * perp.cross(r_perp) < 0) {
          angle_rad *= -1;
        }
      }
    }

    //! Rotation angle in radians or degrees.
    FloatType
    angle(bool deg=false) const
    {
      if (!deg) return angle_rad;
      return angle_rad / scitbx::constants::pi_180;
    }

    //! Reconstructs the rotation matrix using axis_and_angle_as_matrix().
    scitbx::mat3<FloatType>
    as_matrix() const
    {
      return axis_and_angle_as_matrix(axis, angle_rad);
    }
  };

}}} // namespace scitbx::math::r3_rotation

#endif // SCITBX_MATH_R3_ROTATION_H
