/*
* pixel_to_millimeteir_interface.h
*
*  Copyright (C) 2014 Diamond Light Source
*
*  Author: James Parkhurst
*
*  This code is distributed under the BSD license, a copy of which is
*  included in the root directory of this package.
*/
#ifndef DXTBX_MODEL_PIXEL_TO_MILLIMETER_INTERFACE_H
#define DXTBX_MODEL_PIXEL_TO_MILLIMETER_INTERFACE_H

#include <scitbx/vec2.h>

namespace dxtbx { namespace model {

  using scitbx::vec2;

  /** Pre-declare the panel class */
  class Panel;

  /**
   * Base class for the pixel to millimeter strategy
   */
  class PxMmStrategy {
  public:

    /** Virtual desctructor */
    virtual ~PxMmStrategy() {}

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
      vec2<double> xy) const = 0;
  };

}} // namespace dxtbx::model

#endif /* DXTBX_MODEL_PIXEL_TO_MILLIMETER_INTERFACE_H */
