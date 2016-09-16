/*
 * multi_axis_goniometer.h
 *
 *  Copyright (C) 2016 Diamond Light Source
 *
 *  Author: Richard Gildea
 *
 *  This code is distributed under the BSD license, a copy of which is
 *  included in the root directory of this package.
 */
#ifndef DXTBX_MODEL_MULTI_AXIS_GONIOMETER_H
#define DXTBX_MODEL_MULTI_AXIS_GONIOMETER_H

#include <iostream>
#include <scitbx/array_family/shared.h>
#include <scitbx/vec3.h>
#include <scitbx/mat3.h>
#include <scitbx/math/r3_rotation.h>
#include <scitbx/constants.h>
#include <dxtbx/error.h>
#include "goniometer.h"

namespace dxtbx { namespace model {

  using scitbx::vec3;
  using scitbx::mat3;
  using scitbx::constants::pi;
  using scitbx::math::r3_rotation::axis_and_angle_as_matrix;

  /*
   * A class representing a generic multi-axis goniometer. The axes should be
   * given in the order as viewed from the crystal to the base of the goniometer.
   * I.e. for a kappa goniometer the axes (and angles) would be provided in the
   * order phi -> kappa -> omega. The scan axis is identified as the index into
   * the list of axes provided.
   */
  class MultiAxisGoniometer : public Goniometer {
  public:

    /* Default constructor */
    MultiAxisGoniometer()
      : axes_(scitbx::af::shared<vec3<double> >(1, vec3<double>(1.0,0.0,0.0))),
        angles_(scitbx::af::shared<double>((0.0))),
        scan_axis_(0) {}

    /*
     * Initialise the goniometer.
     * @param axes The axes viewed from the crystal to the goniometer base
     * @param angles The corresponding axis angles
     * @param scan_axis The index into the list of axes identifying the scan axis
     */
    MultiAxisGoniometer(const scitbx::af::const_ref<vec3<double> > &axes,
                        const scitbx::af::const_ref<double> &angles,
                        const scitbx::af::const_ref<std::string> &names,
                        std::size_t scan_axis)
      : axes_(axes.begin(), axes.end()),
        angles_(angles.begin(), angles.end()),
        names_(names.begin(), names.end()),
        scan_axis_(scan_axis)
    {
      DXTBX_ASSERT(axes.size() >= 1);
      DXTBX_ASSERT(scan_axis < axes.size());
      init();
    }

    /* Virtual destructor */
    virtual ~MultiAxisGoniometer() {}

    /* Get the goniometer axes */
    scitbx::af::shared<vec3<double> > get_axes() const {
      return scitbx::af::shared<vec3<double> >(axes_.begin(), axes_.end());
    }

    /* Set the goniometer axes */
    void set_axes(scitbx::af::const_ref<vec3<double> > axes) {
      DXTBX_ASSERT(axes.size() == axes_.size());
      for (std::size_t i = 0; i < axes.size(); i++) {
        axes_[i] = axes[i];
      }
      init();
    }

    /* Get the goniometer angles */
    scitbx::af::shared<double> get_angles() const {
      return scitbx::af::shared<double>(angles_.begin(), angles_.end());
    }

    /* Set the goniometer angles */
    void set_angles(scitbx::af::shared<double> angles) {
      DXTBX_ASSERT(angles.size() == angles_.size());
      for (std::size_t i = 0; i < angles.size(); i++) {
        angles_[i] = angles[i];
      }
      init();
    }

    /* Get the axis names */
    scitbx::af::shared<std::string> get_names() const {
      return scitbx::af::shared<std::string>(names_.begin(), names_.end());
    }

    /* Get the scan axis */
    std::size_t get_scan_axis() const {
      return scan_axis_;
    }

    friend std::ostream& operator<<(std::ostream& os, const MultiAxisGoniometer &g);

  protected:

    void init() {
      setting_rotation_ = calculate_setting_rotation();
      fixed_rotation_ = calculate_fixed_rotation();
      // if we are using the setting rotation then this should be
      // applied before this axis => do not include the setting here...
      rotation_axis_ = axes_[scan_axis_];
    }

    /* Calculate the fixed rotation */
    mat3 <double> calculate_fixed_rotation() {
      mat3<double> fixed_rotation(
        1.0, 0.0, 0.0,
        0.0, 1.0, 0.0,
        0.0, 0.0, 1.0);
      for (std::size_t i = 0; i < scan_axis_; i++) {
        mat3 <double> R = axis_and_angle_as_matrix(axes_[i], angles_[i], true);
        fixed_rotation = R * fixed_rotation;
      }
      return fixed_rotation;
    }

    /* Calculate the setting rotation */
    mat3 <double> calculate_setting_rotation() {
      mat3<double> setting_rotation(
        1.0, 0.0, 0.0,
        0.0, 1.0, 0.0,
        0.0, 0.0, 1.0);
      for (std::size_t i = scan_axis_+1; i < axes_.size(); i++) {
        mat3 <double> R = axis_and_angle_as_matrix(axes_[i], angles_[i], true);
        setting_rotation = R * setting_rotation;
      }
      return setting_rotation;
    }

    scitbx::af::shared<vec3<double> > axes_;
    scitbx::af::shared<double> angles_;
    scitbx::af::shared<std::string> names_;
    std::size_t scan_axis_;
  };

  /* Print goniometer info */
  inline
  std::ostream& operator<<(std::ostream& os, const MultiAxisGoniometer &g) {
    os << "Goniometer:\n";
    os << "    Rotation axis:   " << g.get_rotation_axis().const_ref() << "\n";
    os << "    Fixed rotation:  " << g.get_fixed_rotation().const_ref() << "\n";
    os << "    Setting rotation:" << g.get_setting_rotation().const_ref() << "\n";
    for (std::size_t i=0; i < g.get_axes().size(); i++) {
      os << "    Axis #" << i << " (" << g.get_names()[i] << ")" << ":  "
         << g.get_axes()[i].const_ref() << "\n";
    }
    os << "    scan axis: #" << g.get_scan_axis()
       << " (" << g.get_names()[g.get_scan_axis()] << ")" << "\n";
    return os;
  }

}} // namespace dxtbx::model

#endif // DXTBX_MODEL_MULTI_AXIS_GONIOMETER_H
