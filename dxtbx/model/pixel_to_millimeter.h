/*
* pixel_to_millimeter.h
*
*  Copyright (C) 2013 Diamond Light Source
*
*  Author: James Parkhurst
*
*  This code is distributed under the BSD license, a copy of which is
*  included in the root directory of this package.
*/
#ifndef DXTBX_MODEL_PIXEL_TO_MILLIMETER_H
#define DXTBX_MODEL_PIXEL_TO_MILLIMETER_H

#include <scitbx/vec2.h>
#include <dxtbx/model/parallax_correction.h>
#include <dxtbx/model/panel_data.h>

namespace dxtbx { namespace model {

  using scitbx::vec2;

  /**
   * Base class for the pixel to millimeter strategy
   */
  class PxMmStrategy {
  public:

    /** Virtual desctructor */
    virtual ~PxMmStrategy() {}

    /** @returns the name */
    virtual std::string name() const = 0;

    /**
     * Convert a pixel coordinate to a millimeter coordinate
     * @param panel The panel structure
     * @param xy The (x, y) pixel coordinate
     * @return The (x, y) millimeter coordinate
     */
    virtual vec2<double> to_millimeter(const PanelData &panel,
      vec2<double> xy) const = 0;

    /**
     * Convert a millimeter coordinate to a pixel coordinate
     * @param panel The panel structure
     * @param xy The (x, y) millimeter coordinate
     * @return The (x, y) pixel coordinate
     */
    virtual vec2<double> to_pixel(const PanelData &panel,
      vec2<double> xy) const = 0;
  };

  /**
   * The simple pixel to millimeter strategy. Multiply the pixel coordinate
   * by the pixel size and vice versa.
   */
  class SimplePxMmStrategy : public PxMmStrategy {
  public:

    /** Virtual desctructor */
    virtual ~SimplePxMmStrategy() {}

    /** @returns the name */
    virtual std::string name() const {
      return "SimplePxMmStrategy";
    }

    /**
     * Convert a pixel coordinate to a millimeter coordinate
     * @param panel The panel structure
     * @param xy The (x, y) pixel coordinate
     * @return The (x, y) millimeter coordinate
     */
    vec2<double> to_millimeter(const PanelData &panel,
        vec2<double> xy) const {
      vec2<double> pixel_size = panel.get_pixel_size();
      return vec2<double> (xy[0] * pixel_size[0], xy[1] * pixel_size[1]);
    }

    /**
     * Convert a millimeter coordinate to a pixel coordinate
     * @param panel The panel structure
     * @param xy The (x, y) millimeter coordinate
     * @return The (x, y) pixel coordinate
     */
    vec2<double> to_pixel(const PanelData &panel,
        vec2<double> xy) const {
      vec2<double> pixel_size = panel.get_pixel_size();
      return vec2<double> (xy[0] / pixel_size[0], xy[1] / pixel_size[1]);
    }
  };

  /**
   * The parallax corrected strategy. From the simple conversion, then
   * perform a parallax correction.
   */
  class ParallaxCorrectedPxMmStrategy : public SimplePxMmStrategy {
  public:
    ParallaxCorrectedPxMmStrategy(double la) : la_(la) {}

    /** Virtual desctructor */
    virtual ~ParallaxCorrectedPxMmStrategy() {}

    /** @returns the name */
    virtual std::string name() const {
      return "ParallaxCorrectedPxMmStrategy";
    }

    /** @returns the attenuation length */
    double attenuation_length() const {
      return la_;
    }

    /**
     * Convert a pixel coordinate to a millimeter coordinate
     * @param panel The panel structure
     * @param xy The (x, y) pixel coordinate
     * @return The (x, y) millimeter coordinate
     */
    vec2<double> to_millimeter(const PanelData &panel,
        vec2<double> xy) const {
      return parallax_correction_inv(
        panel.get_distance(), la_, panel.get_normal_origin(),
        SimplePxMmStrategy::to_millimeter(panel, xy));
    }

    /**
     * Convert a millimeter coordinate to a pixel coordinate
     * @param panel The panel structure
     * @param xy The (x, y) millimeter coordinate
     * @return The (x, y) pixel coordinate
     */
    vec2<double> to_pixel(const PanelData &panel,
        vec2<double> xy) const {
      return SimplePxMmStrategy::to_pixel(panel,
        parallax_correction(panel.get_distance(), la_,
          panel.get_normal_origin(), xy));
    }

  protected:
    double la_;
  };

}} // namespace dxtbx::model

#endif /* DXTBX_MODEL_PIXEL_TO_MILLIMETER_H */
