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

#include <scitbx/vec2.h>
#include <scitbx/array_family/flex_types.h>
#include <dxtbx/error.h>

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
             double exposure_time,
             const flex_double &epochs)
      : image_range_(image_range),
        oscillation_(oscillation),
        num_images_(1 + image_range_[1] - image_range_[0]),
        exposure_time_(exposure_time),
        epochs_(epochs) {
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

    /** Get the total oscillation range of the scan */
    vec2 <double> get_oscillation_range() const {
      return vec2 <double> (
        oscillation_[0],
        oscillation_[0] + num_images_ * oscillation_[1]);
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
    bool operator==(const ScanData &scan) {
      double eps = 1.0e-6;
      double d_angle = std::abs(oscillation_[0] - scan.oscillation_[0]);
      double d_range = std::abs(oscillation_[1] - scan.oscillation_[1]);
      double d_expos = std::abs(exposure_time_ - scan.exposure_time_);
      return image_range_ == scan.image_range_ &&
             d_angle <= eps && d_range <= eps && d_expos <= eps;
    }

    /** Check the scans are not the same */
    bool operator!=(const ScanData &scan) {
      return !(*this == scan);
    }

  private:

    vec2 <int> image_range_;
    vec2 <double> oscillation_;
    int num_images_;
    double exposure_time_;
    flex_double epochs_;
  };

}} // namespace dxtbx::model

#endif // DXTBX_MODEL_SCAN_H
