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

#include <cmath>
#include <iostream>
#include <scitbx/vec2.h>
#include <scitbx/array_family/flex_types.h>
#include <scitbx/array_family/simple_io.h>
#include <scitbx/array_family/simple_tiny_io.h>
#include <dxtbx/error.h>
#include "scan_helpers.h"

namespace dxtbx { namespace model {

  using std::fabs;
  using scitbx::vec2;
  using scitbx::af::flex_double;
  using scitbx::rad_as_deg;

  /** A scan base class */
  class ScanBase {};

  /** A class to represent a scan */
  class Scan : public ScanBase {
  public:

    /** The default constructor */
    Scan() :
      image_range_(0, 0),
      oscillation_(0.0, 0.0),
      exposure_time_(0.0),
      num_images_(0) {}

    /**
     * Initialise the class
     * @param image_range The range of images covered by the scan
     * @param starting_angle The starting rotation angle
     * @param oscillation_range The oscillation range of each frame
     */
    Scan(vec2 <int> image_range,
         vec2 <double> oscillation,
         double exposure_time)
      : image_range_(image_range),
        oscillation_(oscillation),
        exposure_time_(exposure_time),
        num_images_(1 + image_range_[1] - image_range_[0]),
        epochs_(num_images_) {
      DXTBX_ASSERT(num_images_ >= 0);
      DXTBX_ASSERT(exposure_time_ >= 0.0);
    }

    /**
     * Initialise the class
     * @param image_range The range of images covered by the scan
     * @param starting_angle The starting rotation angle
     * @param oscillation_range The oscillation range of each frame
     */
    Scan(vec2 <int> image_range,
         vec2 <double> oscillation,
         double exposure_time,
         const flex_double &epochs)
      : image_range_(image_range),
        oscillation_(oscillation),
        exposure_time_(exposure_time),
        num_images_(1 + image_range_[1] - image_range_[0]),
        epochs_(epochs) {
      DXTBX_ASSERT(num_images_ >= 0);
      DXTBX_ASSERT(exposure_time_ >= 0.0);
      DXTBX_ASSERT(epochs_.size() == num_images_);
    }

    /** Virtual destructor */
    virtual ~Scan() {}

    /** Get the image range */
    vec2 <int> get_image_range() const {
      return image_range_;
    }

    /** Get the array range (zero based) */
    vec2<int> get_array_range() const {
      return vec2<int>(0, image_range_[1]);
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
    flex_double get_epochs() const {
      return epochs_;
    }

    /** Set the image range */
    void set_image_range(vec2 <int> image_range) {
      image_range_ = image_range;
      num_images_ = 1 + image_range_[1] - image_range_[0];
      epochs_.resize(num_images_);
      DXTBX_ASSERT(num_images_ > 0);
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
      DXTBX_ASSERT(epochs.size() == num_images_);
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
      DXTBX_ASSERT(image_range_[0] <= index && index <= image_range_[1]);
      return epochs_[index - image_range_[0]];
    }

    /** Check the scans are the same */
    bool operator==(const Scan &scan) const {
      double eps = 1.0e-6;
      double d_angle = std::abs(oscillation_[0] - scan.oscillation_[0]);
      double d_range = std::abs(oscillation_[1] - scan.oscillation_[1]);
      double d_expos = std::abs(exposure_time_ - scan.exposure_time_);
      return image_range_ == scan.image_range_ &&
             d_angle <= eps && d_range <= eps && d_expos <= eps;
    }

    /** Check the scans are not the same */
    bool operator!=(const Scan &scan) const {
      return !(*this == scan);
    }

    /** Comparison operator */
    bool operator<(const Scan &scan) const {
      return image_range_[0] < scan.image_range_[0];
    }

    /** Comparison operator */
    bool operator<=(const Scan &scan) const {
      return image_range_[0] <= scan.image_range_[0];
    }

    /** Comparison operator */
    bool operator>(const Scan &scan) const {
      return image_range_[0] > scan.image_range_[0];
    }

    /** Comparison operator */
    bool operator>=(const Scan &scan) const {
      return image_range_[0] >= scan.image_range_[0];
    }

    /**
     * Return a new sweep which cosists of the contents of this sweep and
     * the contents of the other sweep, provided that they are consistent -
     * if they are not consistent then an AssertionError will result.
     */
    Scan operator+(const Scan &scan) const {

      // Ensurw scans are consistent
      double eps = 1e-6;
      DXTBX_ASSERT(exposure_time_ == scan.exposure_time_);
      DXTBX_ASSERT(image_range_[1] + 1 == scan.image_range_[0]);
      DXTBX_ASSERT(fabs(oscillation_[1] - scan.oscillation_[1]) < eps);
      DXTBX_ASSERT(fabs(get_oscillation_range()[1] -
                        scan.get_oscillation_range()[0]) < eps);

      // Create the new image range
      vec2<int> new_image_range(image_range_[0], scan.image_range_[1]);

      // Create the new array of epochs
      flex_double new_epochs(epochs_.size() + scan.epochs_.size());
      for (std::size_t i = 0; i < epochs_.size(); ++i) {
        new_epochs[i] = epochs_[i];
      }
      for (std::size_t i = 0; i < scan.epochs_.size(); ++i) {
        new_epochs[i + epochs_.size()] = scan.epochs_[i];
      }

      // Return a new scan
      return Scan(new_image_range, oscillation_, exposure_time_, new_epochs);
    }

    /**
     * Check if the angle is the range of angles coverd by the scan.
     */
    bool is_angle_valid(double angle) const {
      return is_angle_in_range(get_oscillation_range(), angle);
    }

    /** Check if the index is valid */
    bool is_image_index_valid(double index) const {
      return (image_range_[0] <= index && index <= image_range_[1]);
    }

    /** Check if the array index is valid */
    bool is_array_index_valid(double index) const {
      return is_image_index_valid(index + 1);
    }

    /**
     * Calculate the angle corresponding to the given frame
     * @param index The frame number
     * @returns The angle at the given frame
     */
    double get_angle_from_image_index(double index) const {
      return oscillation_[0] + (index - image_range_[0]) * oscillation_[1];
    }

    /**
     * Calculate the angle corresponding to the given zero based frame
     * @param index The frame number
     * @returns The angle at the given frame
     */
    double get_angle_from_array_index(double index) const {
      return get_angle_from_image_index(index + 1);
    }

    /**
     * Calculate the frame corresponding to the given angle
     * @param angle The angle
     * @returns The frame at the given angle
     */
    double get_image_index_from_angle(double angle) const {
      return image_range_[0] + (angle - oscillation_[0]) / oscillation_[1];
    }

    /**
     * Calculate the zero based frame corresponding to the given angle
     * @param angle The angle
     * @returns The frame at the given angle
     */
    double get_array_index_from_angle(double angle) const {
      return get_image_index_from_angle(angle) - 1;
    }

    /**
     * A function to calculate all the frames in the scan at which an
     * observation with a given angle will be observed. I.e. for a given angle,
     * find all the equivalent angles (i.e. mod 2pi) within the scan range and#
     * calculate the frame number for each angle.
     * Calculate and return an array of frame numbers at which a reflection
     * with a given rotation angle will be observed.
     * @param angle The rotation angle of the reflection
     * @returns The array of frame numbers
     */
    flex_double get_image_indices_with_angle(double angle) const {
      flex_double result = get_mod2pi_angles_in_range(
        get_oscillation_range(), angle);
      for (std::size_t i = 0; i < result.size(); ++i) {
        result[i] = get_image_index_from_angle(result[i]);
      }
      return result;
    }

    /**
     * Calculate and return an array of zero based frame numbers at which a
     * reflection with a given rotation angle will be observed.
     * @param angle The rotation angle of the reflection
     * @returns The array of frame numbers
     */
    flex_double get_array_indices_with_angle(double angle) const {
      flex_double result = get_mod2pi_angles_in_range(
        get_oscillation_range(), angle);
      for (std::size_t i = 0; i < result.size(); ++i) {
        result[i] = get_array_index_from_angle(result[i]);
      }
      return result;
    }

    friend std::ostream& operator<<(std::ostream &os, const Scan &s);

  private:

    vec2 <int> image_range_;
    vec2 <double> oscillation_;
    double exposure_time_;
    int num_images_;
    flex_double epochs_;
  };

  /** Print Scan information */
  inline
  std::ostream& operator<<(std::ostream &os, const Scan &s) {
    // Print oscillation as degrees!
    vec2<double> oscillation = s.get_oscillation();
    oscillation[0] = rad_as_deg(oscillation[0]);
    oscillation[1] = rad_as_deg(oscillation[1]);
    os << "Scan:\n";
    os << "    image range:   " << s.get_image_range().const_ref() << "\n";
    os << "    oscillation:   " << oscillation.const_ref() << "\n";
    os << "    exposure time: " << s.get_exposure_time() << "\n";
    return os;
  }

}} // namespace dxtbx::model

#endif // DXTBX_MODEL_SCAN_H
