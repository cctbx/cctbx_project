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

namespace dxtbx { namespace model {

  using scitbx::vec2;

  /** Pre-declare the panel class */
  class Panel;

  /**
   * Base class for the pixel to millimeter strategy
   */
  class PxMmStrategy {
  public:

    /**
     * Convert a pixel coordinate to a millimeter coordinate
     * @param panel The panel structure
     * @param xy The (x, y) pixel coordinate
     * @return The (x, y) millimeter coordinate
     */
    virtual vec2<double> to_millimeter(const Panel &panel,
      vec2<double> xy) const = 0;

    /**
     * Convert a millimeter coordinate to a pixel coordinate
     * @param panel The panel structure
     * @param xy The (x, y) millimeter coordinate
     * @return The (x, y) pixel coordinate
     */
    virtual vec2<double> to_pixel(const Panel &panel,
      vec2<double> xy) const= 0;
  };

  /**
   * The simple pixel to millimeter strategy. Multiply the pixel coordinate
   * by the pixel size and vice versa.
   */
  class SimplePxMmStrategy : public PxMmStrategy {
  public:

    virtual vec2<double> to_millimeter(const Panel &panel,
        vec2<double> xy) const;

    virtual vec2<double> to_pixel(const Panel &panel, vec2<double> xy) const;
  };

  /**
   * The parallax corrected strategy. From the simple conversion, then
   * perform a parallax correction.
   */
  class ParallaxCorrectedPxMmStrategy : public SimplePxMmStrategy {
  public:
    ParallaxCorrectedPxMmStrategy(double la) : la_(la) {}

    virtual vec2<double> to_millimeter(const Panel &panel,
        vec2<double> xy) const;

    virtual vec2<double> to_pixel(const Panel &panel, vec2<double> xy) const;

  protected:
    double la_;
  };

}} // namespace dxtbx::model

#endif /* DXTBX_MODEL_PIXEL_TO_MILLIMETER_H */
