/*
* panel_data.h
*
*  Copyright (C) 2013 Diamond Light Source
*
*  Author: James Parkhurst
*
*  This code is distributed under the BSD license, a copy of which is
*  included in the root directory of this package.
*/
#ifndef DXTBX_MODEL_PANEL_DATA_H
#define DXTBX_MODEL_PANEL_DATA_H

#include <string>
#include <scitbx/array_family/tiny_types.h>
#include <dxtbx/model/model_helpers.h>
#include <dxtbx/model/virtual_panel.h>
#include <dxtbx/error.h>

namespace dxtbx { namespace model {

  using scitbx::af::int4;
  using scitbx::af::tiny;

  /**
   * A panel class.
   */
  class PanelData : public VirtualPanel {
  public:

    /** Construct the panel with the simple px->mm strategy */
    PanelData()
      : pixel_size_(0.0, 0.0),
        image_size_(0, 0),
        trusted_range_(0.0, 0.0),
        thickness_(0.0) {}

    /** Construct with data */
    PanelData(std::string type,
          std::string name,
          tiny<double,3> fast_axis,
          tiny<double,3> slow_axis,
          tiny<double,3> origin,
          tiny<double,2> pixel_size,
          tiny<std::size_t,2> image_size,
          tiny<double,2> trusted_range,
          double thickness,
          std::string material)
      : pixel_size_(pixel_size),
        image_size_(image_size),
        trusted_range_(trusted_range),
        thickness_(thickness),
        material_(material) {
      set_type(type);
      set_name(name);
      set_local_frame(fast_axis, slow_axis, origin);
    }

    virtual ~PanelData() {}

    /** Get the pixel size */
    tiny<double,2> get_pixel_size() const {
      return pixel_size_;
    }

    /** Set the pixel size */
    void set_pixel_size(tiny<double,2> pixel_size) {
      pixel_size_ = pixel_size;
    }

    /** Get the image size */
    tiny<std::size_t,2> get_image_size() const {
      return image_size_;
    }

    /** Set the image size */
    void set_image_size(tiny<std::size_t,2> image_size) {
      image_size_ = image_size;
    }

    /** Get the trusted range */
    tiny<double,2> get_trusted_range() const {
      return trusted_range_;
    }

    /** Set the trusted range */
    void set_trusted_range(tiny<double,2> trusted_range) {
      trusted_range_ = trusted_range;
    }

    /** Get the thickness */
    double get_thickness() const {
      return thickness_;
    }

    /** Set the thickness */
    void set_thickness(double thickness) {
      thickness_ = thickness;
    }

    /** Get the material */
    std::string get_material() const {
      return material_;
    }

    /** Set the material */
    void set_material(const std::string &material) {
      material_ = material;
    }

    /** Get the mask array */
    scitbx::af::shared<int4> get_mask() const {
      scitbx::af::shared<int4> result((scitbx::af::reserve(mask_.size())));
      std::copy(mask_.begin(), mask_.end(), std::back_inserter(result));
      return result;
    }

    /** Set the mask */
    void set_mask(const scitbx::af::const_ref<int4> &mask) {
      mask_.clear();
      std::copy(mask.begin(), mask.end(), std::back_inserter(mask_));
    }

    /** Add an element to the mask */
    void add_mask(int f0, int s0, int f1, int s1) {
      mask_.push_back(int4(f0, f1, s0, s1));
    }

    /** @returns True/False this is the same as the other */
    bool operator==(const PanelData &rhs) const {
      return VirtualPanel::operator==(rhs)
          && image_size_.const_ref().all_eq(rhs.image_size_.const_ref())
          && pixel_size_.const_ref().all_approx_equal(rhs.pixel_size_.const_ref(), 1e-6)
          && trusted_range_.const_ref().all_approx_equal(rhs.trusted_range_.const_ref(), 1e-6);
    }

    /** @returns True/False this is not the same as the other */
    bool operator!=(const PanelData &rhs) const {
      return !(*this == rhs);
    }

    /** @returns True/False the panels are similar */
    bool is_similar_to(const PanelData &rhs) const {
      return image_size_.const_ref().all_eq(
              rhs.image_size_.const_ref())
          && pixel_size_.const_ref().all_approx_equal(
              rhs.pixel_size_.const_ref(), 1e-7);
//          && trusted_range_.const_ref().all_approx_equal(
//              rhs.trusted_range_.const_ref(), 1e-7);
    }

  protected:
    tiny<double,2> pixel_size_;
    tiny<std::size_t,2> image_size_;
    tiny<double,2> trusted_range_;
    double thickness_;
    std::string material_;
    scitbx::af::shared<int4> mask_;
  };

}} // namespace dxtbx::model

#endif // DXTBX_MODEL_PANEL_H

