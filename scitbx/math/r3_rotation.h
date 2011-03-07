#ifndef SCITBX_MATH_R3_ROTATION_H
#define SCITBX_MATH_R3_ROTATION_H

#include <scitbx/matrix/row_echelon_full_pivoting_small.h>
#include <scitbx/mat3.h>
#include <scitbx/constants.h>

namespace scitbx { namespace math {

  //! Algorithms for R3 (i.e. 3-dimensional space) rotation matrices.
  namespace r3_rotation {

  namespace detail {

    // hiding static const char* in a function to avoid gcc 3.2 -Wall warning
    inline
    const char*
    very_short_axis_message()
    {
      static const char* result =
        "Very short rotation axis vector may lead to numerical instabilities.";
      return result;
    }

  } // namespace detail

  //! Conversion of axis and angle to a rotation matrix.
  /*! http://skal.planet-d.net/demo/matrixfaq.htm
   */
  template <typename FloatType>
  mat3<FloatType>
  axis_and_angle_as_matrix(
    vec3<FloatType> const& axis,
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
      throw std::runtime_error(detail::very_short_axis_message());
    }
    u /= l;
    v /= l;
    w /= l;
    if (deg) angle *= constants::pi_180;
    FloatType c = std::cos(angle);
    FloatType s = std::sin(angle);
    FloatType oc = 1-c;
    FloatType uoc = u*oc;
    FloatType voc = v*oc;
    FloatType woc = w*oc;
    FloatType us = u*s;
    FloatType vs = v*s;
    FloatType ws = w*s;
    return mat3<FloatType>(
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

  //! Conversion without validation of inputs.
  template <typename FloatType>
  af::tiny<FloatType, 4>
  normalized_axis_and_angle_rad_as_unit_quaternion(
    const FloatType* axis,
    FloatType const& angle)
  {
    FloatType h = angle * 0.5;
    FloatType ca = std::cos(h);
    FloatType sa = std::sin(h);
    return af::tiny<FloatType, 4>(ca, axis[0]*sa, axis[1]*sa, axis[2]*sa);
  }

  //! Conversion with validation of inputs.
  template <typename FloatType>
  af::tiny<FloatType, 4>
  axis_and_angle_as_unit_quaternion(
    vec3<FloatType> const& axis,
    FloatType angle,
    bool deg=false,
    FloatType const& min_axis_length=1.e-15)
  {
    SCITBX_ASSERT(min_axis_length > 0);
    FloatType l = axis.length();
    if (l < min_axis_length) {
      throw std::runtime_error(detail::very_short_axis_message());
    }
    if (deg) angle *= constants::pi_180;
    return normalized_axis_and_angle_rad_as_unit_quaternion(
      (axis / l).begin(), angle);
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
    vec3<FloatType> axis;
    //! Rotation angle in radians.
    FloatType angle_rad;

    //! Default constructor. Some data members are not initialized!
    axis_and_angle_from_matrix() {}

    //! Computation of rotation axis and angle.
    axis_and_angle_from_matrix(mat3<FloatType> const& r)
    {
      // obtain axis by solving the system r-i=0
      mat3<FloatType> a_work(
        r[0]-1, r[1],   r[2],
        r[3],   r[4]-1, r[5],
        r[6],   r[7],   r[8]-1);
      matrix::row_echelon::full_pivoting_small<double, 3, 3>
        row_echelon_form(
          af::ref<FloatType, af::c_grid<2> >(
            a_work.begin(), af::c_grid<2>(3,3)),
          /*min_abs_pivot*/ 0,
          /*max_rank*/ 2);
      axis = row_echelon_form.back_substitution(
        a_work.begin(),
        af::small<double, 3>(row_echelon_form.nullity, 1));
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
      vec3<FloatType> perp;
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
      vec3<FloatType> r_perp = r * perp;
      // determine the angle
      FloatType cos_angle = perp * r_perp;
      if      (cos_angle <= -1) angle_rad = constants::pi;
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
      return angle_rad / constants::pi_180;
    }

    //! Reconstructs the rotation matrix using axis_and_angle_as_matrix().
    mat3<FloatType>
    as_matrix() const
    {
      return axis_and_angle_as_matrix(axis, angle_rad);
    }

    /*! \brief Returns (q0,q1,q2,q3), where q0 is the scalar part
        of the unit quaternion.
     */
    af::tiny<FloatType, 4>
    as_unit_quaternion() const
    {
      return normalized_axis_and_angle_rad_as_unit_quaternion(
        axis.begin(), angle_rad);
    }
  };

  //! Returns rotation matrix mapping given_unit_vector onto target_unit_vector.
  /*! It is not checked if the input vectors are unit vectors.
      The result is meaningless if this is not true.

      See also: axis_and_angle_as_matrix()
   */
  template <typename FloatType>
  mat3<FloatType>
  vector_to_vector(
    vec3<FloatType> const& given_unit_vector,
    vec3<FloatType> const& target_unit_vector,
    FloatType const& sin_angle_is_zero_threshold=1.e-10)
  {
    typedef FloatType ft;
    typedef mat3<ft> m3;
    vec3<ft> perp = given_unit_vector.cross(target_unit_vector);
    ft c = given_unit_vector * target_unit_vector;
    ft s = perp.length();
    if (s < sin_angle_is_zero_threshold) {
      if (c > 0) {
        return m3(1,0,0,0,1,0,0,0,1);
      }
      perp = target_unit_vector.ortho(/* normalize */ true);
      // specialization of general code below, with s=0 and c=-1
      ft u = perp[0];
      ft v = perp[1];
      ft w = perp[2];
      ft tu = 2*u;
      ft tv = 2*v;
      return m3(tu*u-1, tu*v, tu*w, tu*v, tv*v-1, tv*w, tu*w, tv*w, 2*w*w-1);
    }
    ft us = perp[0];
    ft vs = perp[1];
    ft ws = perp[2];
    ft u = us / s;
    ft v = vs / s;
    ft w = ws / s;
    ft oc = 1-c;
    ft uoc = u*oc;
    ft voc = v*oc;
    ft woc = w*oc;
    ft uvoc = u*voc;
    ft uwoc = u*woc;
    ft vwoc = v*woc;
    return m3(
       c + u*uoc,
     -ws + uvoc,
      vs + uwoc,
      ws + uvoc,
       c + v*voc,
     -us + vwoc,
     -vs + uwoc,
      us + vwoc,
       c + w*woc);
  }

  //! Returns rotation matrix mapping given_unit_vector onto (0,0,1).
  /*! It is not checked if the input vector is a unit vector.
      The result is meaningless if this is not true.

      The implementation is a simplification of vector_to_vector().
   */
  template <typename FloatType>
  mat3<FloatType>
  vector_to_001(
    vec3<FloatType> const& given_unit_vector,
    FloatType const& sin_angle_is_zero_threshold=1.e-10)
  {
    typedef FloatType ft;
    typedef mat3<ft> m3;
    ft x = given_unit_vector[0];
    ft y = given_unit_vector[1];
    ft c = given_unit_vector[2];
    ft s = std::sqrt(x*x + y*y);
    if (s < sin_angle_is_zero_threshold) {
      if (c > 0) {
        return m3(1,0,0,0,1,0,0,0,1);
      }
      return m3(1,0,0,0,-1,0,0,0,-1);
    }
    ft us = y;
    ft vs = -x;
    ft u = us / s;
    ft v = vs / s;
    ft oc = 1-c;
    ft uoc = u*oc;
    ft voc = v*oc;
    ft uvoc = u*voc;
    return m3(c + u*uoc, uvoc, vs, uvoc, c + v*voc, -us, -vs, us, c);
  }

  //! Returns rotation matrix mapping given_unit_vector onto (0,1,0).
  /*! It is not checked if the input vector is a unit vector.
      The result is meaningless if this is not true.

      The implementation is a simplification of vector_to_vector().
   */
  template <typename FloatType>
  mat3<FloatType>
  vector_to_010(
    vec3<FloatType> const& given_unit_vector,
    FloatType const& sin_angle_is_zero_threshold=1.e-10)
  {
    typedef FloatType ft;
    typedef mat3<ft> m3;
    ft x = given_unit_vector[0];
    ft c = given_unit_vector[1];
    ft z = given_unit_vector[2];
    ft s = std::sqrt(x*x + z*z);
    if (s < sin_angle_is_zero_threshold) {
      if (c > 0) {
        return m3(1,0,0,0,1,0,0,0,1);
      }
      return m3(1,0,0,0,-1,0,0,0,-1);
    }
    ft us = -z;
    ft ws = x;
    ft u = us / s;
    ft w = ws / s;
    ft oc = 1-c;
    ft uoc = u*oc;
    ft woc = w*oc;
    ft uwoc = u*woc;
    return m3(c + u*uoc, -ws, uwoc, ws, c, -us, uwoc, us, c + w*woc);
  }

  //! Returns rotation matrix mapping given_unit_vector onto (1,0,0).
  /*! It is not checked if the input vector is a unit vector.
      The result is meaningless if this is not true.

      The implementation is a simplification of vector_to_vector().
   */
  template <typename FloatType>
  mat3<FloatType>
  vector_to_100(
    vec3<FloatType> const& given_unit_vector,
    FloatType const& sin_angle_is_zero_threshold=1.e-10)
  {
    typedef FloatType ft;
    typedef mat3<ft> m3;
    ft c = given_unit_vector[0];
    ft y = given_unit_vector[1];
    ft z = given_unit_vector[2];
    ft s = std::sqrt(y*y + z*z);
    if (s < sin_angle_is_zero_threshold) {
      if (c > 0) {
        return m3(1,0,0,0,1,0,0,0,1);
      }
      return m3(-1,0,0,0,1,0,0,0,-1);
    }
    ft vs = z;
    ft ws = -y;
    ft v = vs / s;
    ft w = ws / s;
    ft oc = 1-c;
    ft voc = v*oc;
    ft woc = w*oc;
    ft vwoc = v*woc;
    return m3(c, -ws, vs, ws, c + v*voc, vwoc, -vs, vwoc, c + w*woc);
  }

  //! Unit quaternion (a.k.a. Euler parameters) as matrix.
  /*! The unit quaternion elements must satisfy the normalization condition
        q0**2+q1**2+q2**2+q3**3 = 1
      but this is not checked.

      Also implemented in Python: scitbx/matrix/__init__.py
   */
  template <typename FloatType>
  mat3<FloatType>
  unit_quaternion_as_matrix(
    FloatType const& q0,
    FloatType const& q1,
    FloatType const& q2,
    FloatType const& q3)
  {
    FloatType q0_q0 = q0*q0;
    FloatType q0_q1 = q0*q1;
    FloatType q0_q2 = q0*q2;
    FloatType q0_q3 = q0*q3;
    FloatType q1_q2 = q1*q2;
    FloatType q1_q3 = q1*q3;
    FloatType q2_q3 = q2*q3;
    return mat3<FloatType>(
      2*(q0_q0+q1*q1)-1, 2*(q1_q2-q0_q3),   2*(q1_q3+q0_q2),
      2*(q1_q2+q0_q3),   2*(q0_q0+q2*q2)-1, 2*(q2_q3-q0_q1),
      2*(q1_q3-q0_q2),   2*(q2_q3+q0_q1),   2*(q0_q0+q3*q3)-1);
  }

  template <typename FloatType>
  mat3<FloatType>
  unit_quaternion_as_matrix(
    af::tiny<FloatType, 4> const& q)
  {
    return unit_quaternion_as_matrix(q[0], q[1], q[2], q[3]);
  }

  //! Matrix as unit quaternion (a.k.a. Euler parameters).
  /*! The matrix elements must satisfy the orthogonality condition
        r.transpose()*r = identity
      but this is not thoroughly checked.

      Based on work by:
        Shepperd (1978), J. Guidance and Control, 1, 223-224.
        Sam Buss, http://math.ucsd.edu/~sbuss/MathCG
        Robert Hanson, jmol/Jmol/src/org/jmol/util/Quaternion.java

      Also implemented in Python: scitbx/matrix/__init__.py
   */
  template <typename FloatType>
  af::tiny<FloatType, 4>
  matrix_as_unit_quaternion(
    mat3<FloatType> const& r)
  {
    typedef FloatType ft;
    ft w, x, y, z;
    ft trace = r[0] + r[4] + r[8];
    if (trace >= 0.5) {
      w = std::sqrt(1 + trace);
      ft d = w + w;
      w *= 0.5;
      x = (r[7] - r[5]) / d;
      y = (r[2] - r[6]) / d;
      z = (r[3] - r[1]) / d;
    }
    else {
      unsigned mx = 2;
      if (r[0] > r[4]) {
        if (r[0] > r[8]) mx = 0;
      }
      else if (r[4] > r[8]) mx = 1;
      ft invalid_cutoff = 0.8; // not critical; true value is closer to 0.83
      const char* invalid_message = "Not a r3_rotation matrix.";
      if (mx == 0) {
        ft x_sq = 1 + r[0] - r[4] - r[8];
        if (x_sq < invalid_cutoff) throw std::runtime_error(invalid_message);
        x = std::sqrt(x_sq);
        ft d = x + x;
        x *= 0.5;
        w = (r[7] - r[5]) / d;
        y = (r[3] + r[1]) / d;
        z = (r[6] + r[2]) / d;
      }
      else if (mx == 1) {
        ft y_sq = 1 + r[4] - r[0] - r[8];
        if (y_sq < invalid_cutoff) throw std::runtime_error(invalid_message);
        y = std::sqrt(y_sq);
        ft d = y + y;
        y *= 0.5;
        w = (r[2] - r[6]) / d;
        x = (r[3] + r[1]) / d;
        z = (r[7] + r[5]) / d;
      }
      else {
        ft z_sq = 1 + r[8] - r[0] - r[4];
        if (z_sq < invalid_cutoff) throw std::runtime_error(invalid_message);
        z = std::sqrt(z_sq);
        ft d = z + z;
        z *= 0.5;
        w = (r[3] - r[1]) / d;
        x = (r[6] + r[2]) / d;
        y = (r[7] + r[5]) / d;
      }
    }
    return af::tiny<FloatType, 4>(w, x, y, z);
  }

}}} // namespace scitbx::math::r3_rotation

#endif // SCITBX_MATH_R3_ROTATION_H
