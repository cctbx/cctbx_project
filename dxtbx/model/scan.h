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
#include <scitbx/array_family/shared.h>
#include <scitbx/array_family/simple_io.h>
#include <scitbx/array_family/simple_tiny_io.h>
#include <dxtbx/error.h>
#include "scan_helpers.h"

namespace dxtbx { namespace model {

  using scitbx::vec2;
  using scitbx::rad_as_deg;
  using scitbx::constants::pi;

  /** A scan base class */
  class ScanBase {};

  /** A class to represent a scan */
  class Scan : public ScanBase {
  public:

    /** The default constructor */
    Scan() :
      image_range_(0, 0),
      oscillation_(0.0, 0.0),
      num_images_(0),
      batch_offset_(0) {}

    /**
     * Initialise the class
     * @param image_range The range of images covered by the scan
     * @param starting_angle The starting rotation angle
     * @param oscillation_range A tuple containing the start angle of the first
     *                          image and the oscillation range of each frame
     */
    Scan(vec2 <int> image_range,
         vec2 <double> oscillation,
         int batch_offset=0)
      : image_range_(image_range),
        oscillation_(oscillation),
        num_images_(1 + image_range_[1] - image_range_[0]),
        batch_offset_(batch_offset),
        exposure_times_(num_images_, 0.0),
        epochs_(num_images_, 0.0) {
      DXTBX_ASSERT(num_images_ >= 0);
    }

    /**
     * Initialise the class
     * @param image_range The range of images covered by the scan
     * @param batch_offset The batch offset for the scan
     * @param starting_angle The starting rotation angle
     * @param oscillation_range The oscillation range of each frame
     */
    Scan(vec2 <int> image_range,
         vec2 <double> oscillation,
         const scitbx::af::shared<double> &exposure_times,
         const scitbx::af::shared<double> &epochs,
         int batch_offset=0)
      : image_range_(image_range),
        oscillation_(oscillation),
        num_images_(1 + image_range_[1] - image_range_[0]),
        batch_offset_(batch_offset),
        exposure_times_(exposure_times),
        epochs_(epochs) {
      DXTBX_ASSERT(num_images_ >= 0);
      if (exposure_times_.size() == 1 && num_images_ > 1) {
        // assume same exposure time for all images - there is
        // probably a better way of coding this...
        scitbx::af::shared<double> expanded_exposure_times;
        expanded_exposure_times.reserve(num_images_);
        for (int j = 0; j < num_images_; j++) {
          expanded_exposure_times.push_back(exposure_times[0]);
          exposure_times_ = expanded_exposure_times;
        }
      }
      DXTBX_ASSERT(exposure_times_.size() == num_images_);
      DXTBX_ASSERT(epochs_.size() == num_images_);
    }

    /** Copy */
    Scan(const Scan &rhs)
      : image_range_(rhs.image_range_),
        oscillation_(rhs.oscillation_),
        num_images_(rhs.num_images_),
        batch_offset_(rhs.batch_offset_),
        exposure_times_(scitbx::af::reserve(rhs.exposure_times_.size())),
        epochs_(scitbx::af::reserve(rhs.epochs_.size())) {
      std::copy(rhs.epochs_.begin(), rhs.epochs_.end(),
        std::back_inserter(epochs_));
      std::copy(rhs.exposure_times_.begin(), rhs.exposure_times_.end(),
        std::back_inserter(exposure_times_));
    }

    /** Virtual destructor */
    virtual ~Scan() {}

    /** Get the image range */
    vec2 <int> get_image_range() const {
      return image_range_;
    }

    /** Get the batch offset */
    int get_batch_offset() const {
      return batch_offset_;
    }

    /** Get the batch number for a given image index */
    int get_batch_for_image_index(int index) const {
      return index + batch_offset_;
    }

    /** Get the batch number for a given array index */
    int get_batch_for_array_index(int index) const {
      return index + batch_offset_ + 1;
    }

    /** Get the batch range */
    vec2 <int> get_batch_range() const {
      return vec2<int>(
        image_range_[0] + batch_offset_, image_range_[1] + batch_offset_);
    }

    /** Get the array range (zero based) */
    vec2<int> get_array_range() const {
      return vec2<int>(image_range_[0]-1, image_range_[1]);
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
    scitbx::af::shared<double> get_exposure_times() const {
      return exposure_times_;
    }

    /** Get the image epochs */
    scitbx::af::shared<double> get_epochs() const {
      return epochs_;
    }

    /** Set the image range */
    void set_image_range(vec2 <int> image_range) {
      image_range_ = image_range;
      num_images_ = 1 + image_range_[1] - image_range_[0];
      epochs_.resize(num_images_);
      exposure_times_.resize(num_images_);
      DXTBX_ASSERT(num_images_ > 0);
    }

    /** Set the batch_offset */
    void set_batch_offset(int batch_offset) {
      batch_offset_ = batch_offset;
    }

    /** Set the oscillation */
    void set_oscillation(vec2 <double> oscillation) {
      oscillation_ = oscillation;
    }

    /** Set the exposure time */
    void set_exposure_times(scitbx::af::shared<double> exposure_times) {
      DXTBX_ASSERT(exposure_times.size() == num_images_);
      exposure_times_ = exposure_times;
    }

    /** Set the image epochs */
    void set_epochs(const scitbx::af::shared<double> &epochs) {
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

    double get_image_exposure_time(int index) const {
      DXTBX_ASSERT(image_range_[0] <= index && index <= image_range_[1]);
      return exposure_times_[index - image_range_[0]];
    }

    /** Check the scans are the same */
    bool operator==(const Scan &rhs) const {
      double eps = 1e-7;
      return image_range_ == rhs.image_range_
          && batch_offset_ == rhs.batch_offset_
          && std::abs(oscillation_[0] - rhs.oscillation_[0]) < eps
          && std::abs(oscillation_[1] - rhs.oscillation_[1]) < eps
          && exposure_times_.const_ref().all_approx_equal(rhs.exposure_times_.const_ref(), eps)
          && epochs_.const_ref().all_approx_equal(rhs.epochs_.const_ref(), eps);
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
     * Append the rhs scan onto the current scan
     */
    void append(const Scan &rhs, double scan_tolerance) {
      double eps = scan_tolerance * std::abs(oscillation_[1]);
      DXTBX_ASSERT(eps > 0);
      DXTBX_ASSERT(std::abs(oscillation_[1]) > 0.0);
      DXTBX_ASSERT(image_range_[1] + 1 == rhs.image_range_[0]);
      DXTBX_ASSERT(std::abs(oscillation_[1] - rhs.oscillation_[1]) < eps);
      DXTBX_ASSERT(batch_offset_ == rhs.batch_offset_);
      // sometimes ticking through 0 the first difference is not helpful
      double diff_2pi = std::abs(mod_2pi(get_oscillation_range()[1]) -
                                 mod_2pi(rhs.get_oscillation_range()[0]));
      double diff_abs = std::abs(get_oscillation_range()[1] -
                                 rhs.get_oscillation_range()[0]);
      DXTBX_ASSERT(std::min(diff_2pi, diff_abs) < eps * get_num_images());
      image_range_[1] = rhs.image_range_[1];
      num_images_ = 1 + image_range_[1] - image_range_[0];
      exposure_times_.reserve(exposure_times_.size() + exposure_times_.size());
      epochs_.reserve(epochs_.size() + epochs_.size());
      std::copy(rhs.exposure_times_.begin(), rhs.exposure_times_.end(),
        std::back_inserter(exposure_times_));
      std::copy(rhs.epochs_.begin(), rhs.epochs_.end(),
        std::back_inserter(epochs_));
    }

    /**
     * Append the rhs scan onto the current scan
     */
    Scan& operator+=(const Scan &rhs) {
      // Set the epsilon to 1% of oscillation range
      append(rhs, 0.01);
      return *this;
    }

    /**
     * Return a new sweep which consists of the contents of this sweep and
     * the contents of the other sweep, provided that they are consistent -
     * if they are not consistent then an AssertionError will result.
     */
    Scan operator+(const Scan &rhs) const {
      Scan lhs(*this);
      lhs += rhs;
      return lhs;
    }

    /**
     * Check if the angle is in the range of angles covered by the scan.
     */
    bool is_angle_valid(double angle) const {
      return is_angle_in_range(get_oscillation_range(), angle);
    }

    /** Check if the index is valid */
    bool is_image_index_valid(double index) const {
      return (image_range_[0] <= index && index <= image_range_[1]);
    }

    /** Check if a given batch is valid */
    bool is_batch_valid(int batch) const {
      vec2<int> batch_range = get_batch_range();
      return (batch_range[0] <= batch && batch <= batch_range[1]);
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
    scitbx::af::shared< vec2<double> > get_image_indices_with_angle(
        double angle) const {
      scitbx::af::shared<double> angles = get_mod2pi_angles_in_range(
        get_oscillation_range(), angle);
      scitbx::af::shared< vec2<double> > result(angles.size());
      for (std::size_t i = 0; i < result.size(); ++i) {
        result[i][0] = angles[i];
        result[i][1] = get_image_index_from_angle(angles[i]);
      }
      return result;
    }

    /**
     * Calculate and return an array of zero based frame numbers at which a
     * reflection with a given rotation angle will be observed.
     * @param angle The rotation angle of the reflection
     * @returns The array of frame numbers
     */
    scitbx::af::shared< vec2<double> > get_array_indices_with_angle(
        double angle, double padding=0, bool deg=false) const {
      DXTBX_ASSERT(padding >= 0);
      if (deg == true) {
        padding = padding * pi / 180.0;
      }
      vec2<double> range = get_oscillation_range();
      range[0] -= padding;
      range[1] += padding;
      scitbx::af::shared<double> angles = get_mod2pi_angles_in_range(
        range, angle);
      scitbx::af::shared< vec2<double> > result(angles.size());
      for (std::size_t i = 0; i < result.size(); ++i) {
        result[i][0] = angles[i];
        result[i][1] = get_array_index_from_angle(angles[i]);
      }
      return result;
    }

    Scan operator[](int index) const {
      // Check index
      DXTBX_ASSERT((index >= 0) && (index < get_num_images()));
      int image_index = get_image_range()[0] + index;

      // Create the new epoch array
      scitbx::af::shared<double> new_epochs(1);
      new_epochs[0] = get_image_epoch(image_index );
      scitbx::af::shared<double> new_exposure_times(1);
      new_exposure_times[0] = get_image_exposure_time(image_index );

      // Return scan
      return Scan(vec2<int>(image_index, image_index),
        get_image_oscillation(image_index ),
        new_exposure_times, new_epochs,
        get_batch_offset());
    }

    friend std::ostream& operator<<(std::ostream &os, const Scan &s);

  private:

    vec2 <int> image_range_;
    vec2 <double> oscillation_;
    int num_images_;
    int batch_offset_;
    scitbx::af::shared<double> exposure_times_;
    scitbx::af::shared<double> epochs_;
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
    if (s.num_images_ > 0) {
      os << "    exposure time: " << s.exposure_times_.const_ref()[0] << "\n";
    }
    return os;
  }

}} // namespace dxtbx::model

#endif // DXTBX_MODEL_SCAN_H
