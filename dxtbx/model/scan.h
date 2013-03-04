/*
 * scan.h
 *
 *  Copyright (C) 2013 Diamond Light Source
 *
 *  Author: James Parkhurst
 *
 *  This code is distributed under the BSD license, a copy of which is
 *  included in the root directory of this package.
 */
#ifndef DXTBX_MODEL_SCAN_H
#define DXTBX_MODEL_SCAN_H

#include <iostream>
#include <scitbx/vec2.h>
#include <scitbx/array_family/flex_types.h>
#include <scitbx/array_family/simple_io.h>
#include <scitbx/array_family/simple_tiny_io.h>
#include <dxtbx/error.h>
#include "scan_helpers.h"

namespace dxtbx { namespace model {

  using scitbx::vec2;
  using scitbx::af::flex_double;

  /** A scan base class */
  class ScanBase {};

  /** A class to represent a scan */
  class ScanData : public ScanBase {
  public:

    /** The default constructor */
    ScanData() {}

    /**
     * Initialise the class
     * @param image_range The range of images covered by the scan
     * @param starting_angle The starting rotation angle
     * @param oscillation_range The oscillation range of each frame
     */
    ScanData(vec2 <int> image_range,
             vec2 <double> oscillation,
             double exposure_time)
      : image_range_(image_range),
        oscillation_(oscillation),
        exposure_time_(exposure_time),
        num_images_(1 + image_range_[1] - image_range_[0]) {
      DXTBX_ASSERT(num_images_ >= 0);
      DXTBX_ASSERT(exposure_time_ >= 0.0);
    }

    /**
     * Initialise the class
     * @param image_range The range of images covered by the scan
     * @param starting_angle The starting rotation angle
     * @param oscillation_range The oscillation range of each frame
     */
    ScanData(vec2 <int> image_range,
             vec2 <double> oscillation,
             double exposure_time,
             const flex_double &epochs)
      : image_range_(image_range),
        oscillation_(oscillation),
        exposure_time_(exposure_time),
        epochs_(epochs),
        num_images_(1 + image_range_[1] - image_range_[0]) {
      DXTBX_ASSERT(num_images_ >= 0);
      DXTBX_ASSERT(exposure_time_ >= 0.0);
      DXTBX_ASSERT(epochs_.size() == num_images_);
    }

    /** Virtual destructor */
    virtual ~ScanData() {}

    /** Get the image range */
    vec2 <int> get_image_range() const {
      return image_range_;
    }

    /** Get the oscillation */
    vec2 <double> get_oscillation() const {
      return oscillation_;
    }

    /** Get the number of images */
    int get_num_images() const {
      return num_images_;
    }

    /** Get the exposure time */
    double get_exposure_time() const {
      return exposure_time_;
    }

    /** Get the image epochs */
    flex_double get_epochs() {
      return epochs_;
    }

    /** Set the image range */
    void set_image_range(vec2 <int> image_range) {
      image_range_ = image_range;
      num_images_ = image_range_[1] - image_range_[0];
      DXTBX_ASSERT(num_images_ >= 0);
    }

    /** Set the oscillation */
    void set_oscillation(vec2 <double> oscillation) {
      oscillation_ = oscillation;
    }

    /** Set the exposure time */
    void set_exposure_time(double exposure_time) {
      exposure_time_ = exposure_time;
    }

    /** Set the image epochs */
    void set_epochs(const flex_double &epochs) {
      epochs_ = epochs;
    }

    /** Get the total oscillation range of the scan */
    vec2 <double> get_oscillation_range() const {
      return vec2 <double> (
        oscillation_[0],
        oscillation_[0] + num_images_ * oscillation_[1]);
    }

    /** Get the image angle and oscillation width as a tuple */
    vec2 <double> get_image_oscillation(int index) const {
      return vec2 <double> (
        oscillation_[0] + (index - image_range_[0]) * oscillation_[1],
        oscillation_[1]);
    }

    /** Get the image epoch */
    double get_image_epoch(int index) const {
      DXTBX_ASSERT(0 <= index && index < epochs_.size());
      return epochs_[index];
    }

    /** Check the scans are the same */
    bool operator==(const ScanData &scan) const {
      double eps = 1.0e-6;
      double d_angle = std::abs(oscillation_[0] - scan.oscillation_[0]);
      double d_range = std::abs(oscillation_[1] - scan.oscillation_[1]);
      double d_expos = std::abs(exposure_time_ - scan.exposure_time_);
      return image_range_ == scan.image_range_ &&
             d_angle <= eps && d_range <= eps && d_expos <= eps;
    }

    /** Check the scans are not the same */
    bool operator!=(const ScanData &scan) const {
      return !(*this == scan);
    }

    /**
     * Check if the angle is the range of angles coverd by the scan.
     */
    bool is_angle_valid(double angle) const {
      return is_angle_in_range(get_oscillation_range(), angle);
    }

    /** Check if the frame is valid */
    bool is_frame_valid(double frame) const {
      return (image_range_[0] <= frame && frame < image_range_[1]);
    }

    /**
     * Calculate the angle corresponding to the given frame
     * @param frame The frame number
     * @returns The angle at the given frame
     */
    double get_angle_from_frame(double frame) const {
      return oscillation_[0] + (frame - image_range_[0]) * oscillation_[1];
    }

    /**
     * Calculate the frame corresponding to the given angle
     * @param angle The angle
     * @returns The frame at the given angle
     */
    double get_frame_from_angle(double angle) const {
      return image_range_[0] + (angle - oscillation_[0]) / oscillation_[1];
    }

    /**
     * A function to calculate all the frames in the scan at which an observation
     * with a given angle will be observed. I.e. for a given angle, find all the
     * equivalent angles (i.e. mod 2pi) within the scan range and calculate the
     * frame number for each angle.
     * Calculate and return an array of frame numbers at which a reflection
     * with a given rotation angle will be observed.
     * @param angle The rotation angle of the reflection
     * @returns The array of frame numbers
     */
    flex_double get_frames_with_angle(double angle) const {
      flex_double result = get_mod2pi_angles_in_range(
        get_oscillation_range(), angle);
      for (std::size_t i = 0; i < result.size(); ++i) {
        result[i] = get_frame_from_angle(result[i]);
      }
      return result;
    }

    friend std::ostream& operator<<(std::ostream &os, const ScanData &s);

  private:

    vec2 <int> image_range_;
    vec2 <double> oscillation_;
    double exposure_time_;
    flex_double epochs_;
    int num_images_;
  };

  /** Print ScanData information */
  inline
  std::ostream& operator<<(std::ostream &os, const ScanData &s) {
    os << "Scan:\n";
    os << "    image range:   " << s.get_image_range().const_ref() << "\n";
    os << "    oscillation:   " << s.get_oscillation().const_ref() << "\n";
    os << "    exposure time: " << s.get_exposure_time() << "\n";
    return os;
  }

}} // namespace dxtbx::model

#endif // DXTBX_MODEL_SCAN_H
