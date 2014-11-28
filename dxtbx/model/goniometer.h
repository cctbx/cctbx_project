/*
 * goniometer.h
 *
 *  Copyright (C) 2013 Diamond Light Source
 *
 *  Author: James Parkhurst
 *
 *  This code is distributed under the BSD license, a copy of which is
 *  included in the root directory of this package.
 */
#ifndef DXTBX_MODEL_GONIOMETER_H
#define DXTBX_MODEL_GONIOMETER_H

#include <iostream>
#include <scitbx/vec3.h>
#include <scitbx/mat3.h>
#include <scitbx/math/r3_rotation.h>
#include <scitbx/constants.h>
#include <scitbx/array_family/simple_io.h>
#include <scitbx/array_family/simple_tiny_io.h>
#include <dxtbx/error.h>
#include "model_helpers.h"

namespace dxtbx { namespace model {

  using scitbx::vec3;
  using scitbx::mat3;
  using scitbx::constants::pi;
  using scitbx::math::r3_rotation::axis_and_angle_as_matrix;

  /** A goniometer base class */
  class GoniometerBase {};

  /**
   * A class to represent the rotation axis for a standard rotation
   * geometry diffraction data set.
   *
   * The rotation axis assumed to have it's origin at the origin of the
   * laboratory coordinate system. The rotation axis vector is normalized
   * when set using either the constructor or the rotation axis 'setter'
   *
   * The fixed rotation matrix, F, represents the additional fixed rotation of
   * the sample attached to the rotating axis - for example the effects of
   * kappa and phi for an omega scan on a kappa goniometer. A second rotation
   * matrix, S, may be applied which is used to align the rotation axis,
   * for example the omega & kappa axes for a phi scan. These rotation matrices
   * should be applied as:
   *
   *  A = [S][R][F][U][B]
   *
   * in a standard orientation matrix.
   *
   * In the document detailing the conventions used:
   *    rotation_axis -> m2
   *    fixed_rotation_matrix -> *unspecified*
   */
  class Goniometer : public GoniometerBase {

  public:

    /**
     * Initialise the goniometer. The fixed rotation matrix is set to the
     * identity matrix and the rotation axis is set to the x axis
     * @param rotation_axis The goniometer rotation axis
     */
    Goniometer()
      : rotation_axis_(1.0, 0.0, 0.0),
        fixed_rotation_(
          1.0, 0.0, 0.0,
          0.0, 1.0, 0.0,
          0.0, 0.0, 1.0),
        setting_rotation_(
          1.0, 0.0, 0.0,
          0.0, 1.0, 0.0,
          0.0, 0.0, 1.0) {}

    /**
     * Initialise the goniometer. The fixed rotation matrix is set to the
     * identity matrix.
     * @param rotation_axis The goniometer rotation axis
     */
    Goniometer(vec3 <double> rotation_axis)
      : fixed_rotation_(
          1.0, 0.0, 0.0,
          0.0, 1.0, 0.0,
          0.0, 0.0, 1.0),
        setting_rotation_(
          1.0, 0.0, 0.0,
          0.0, 1.0, 0.0,
          0.0, 0.0, 1.0) {
      DXTBX_ASSERT(rotation_axis.length() > 0);
      rotation_axis_ = rotation_axis.normalize();
    }

    /**
     * Initialise the goniometer.
     * @param rotation_axis The goniometer rotation axis
     * @param fixed_rotation The additional fixed rotation of the sample
     */
    Goniometer(vec3 <double> rotation_axis,
               mat3 <double> fixed_rotation)
      : fixed_rotation_(fixed_rotation),
        setting_rotation_(
          1.0, 0.0, 0.0,
          0.0, 1.0, 0.0,
          0.0, 0.0, 1.0) {
      DXTBX_ASSERT(rotation_axis.length() > 0);
      rotation_axis_ = rotation_axis.normalize();
    }

    /**
     * Initialise the goniometer.
     * @param rotation_axis The goniometer rotation axis
     * @param fixed_rotation The additional fixed rotation of the sample
     * @param setting_rotation The additional setting rotation of the axis
     */
    Goniometer(vec3 <double> rotation_axis,
               mat3 <double> fixed_rotation,
               mat3 <double> setting_rotation)
      : fixed_rotation_(fixed_rotation),
        setting_rotation_(setting_rotation) {
      DXTBX_ASSERT(rotation_axis.length() > 0);
      rotation_axis_ = rotation_axis.normalize();
    }

    /** Virtual destructor */
    virtual ~Goniometer() {}

    /** Get the rotation axis */
    vec3 <double> get_rotation_axis() const {
      return rotation_axis_;
    }

    /** Get the fixed rotation matrix */
    mat3 <double> get_fixed_rotation() const {
      return fixed_rotation_;
    }

    /** Get the setting matrix */
    mat3 <double> get_setting_rotation() const {
      return setting_rotation_;
    }

    /** Set the rotation axis */
    void set_rotation_axis(vec3 <double> rotation_axis) {
      DXTBX_ASSERT(rotation_axis.length() > 0);
      rotation_axis_ = rotation_axis.normalize();
    }

    /** Set the fixed rotation matrix */
    void set_fixed_rotation(mat3 <double> fixed_rotation) {
      fixed_rotation_ = fixed_rotation;
    }

    /** Set the setting rotation matrix */
    void set_setting_rotation(mat3 <double> setting_rotation) {
      setting_rotation_ = setting_rotation;
    }

    /** Check rotation axes are (almost) the same */
    bool operator==(const Goniometer &b) {
      double eps = 1.0e-6;

      return std::abs(angle_safe(rotation_axis_, b.rotation_axis_)) <= eps
      && fixed_rotation_.const_ref().all_approx_equal(b.fixed_rotation_.const_ref(), eps)
      && setting_rotation_.const_ref().all_approx_equal(b.setting_rotation_.const_ref(), eps);
    }

    /** Check rotation axes are not (almost) the same */
    bool operator!=(const Goniometer &goniometer) {
      return !(*this == goniometer);
    }

    friend std::ostream& operator<<(std::ostream& os, const Goniometer &gonio);

  protected:

    vec3 <double> rotation_axis_;
    mat3 <double> fixed_rotation_;
    mat3 <double> setting_rotation_;
  };

  /** Print goniometer data */
  inline
  std::ostream& operator<<(std::ostream& os, const Goniometer &g) {
    os << "Goniometer:\n";
    os << "    Rotation axis:   " << g.get_rotation_axis().const_ref() << "\n";
    os << "    Fixed rotation:  " << g.get_fixed_rotation().const_ref() << "\n";
    os << "    Setting rotation:" << g.get_setting_rotation().const_ref() << "\n";
    return os;
  }

}} // namespace dxtbx::model

#endif // DXTBX_MODEL_GONIOMETER_H
