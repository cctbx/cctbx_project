/*
* pixel_to_millimeter.cc
*
*  Copyright (C) 2013 Diamond Light Source
*
*  Author: James Parkhurst
*
*  This code is distributed under the BSD license, a copy of which is
*  included in the root directory of this package.
*/

#include <dxtbx/model/panel.h>
#include <dxtbx/model/pixel_to_millimeter.h>

namespace dxtbx { namespace model {

  /**
   * Convert a pixel coordinate to a millimeter coordinate
   * @param panel The panel structure
   * @param xy The (x, y) pixel coordinate
   * @return The (x, y) millimeter coordinate
   */  
  vec2<double> SimplePxMmStrategy::to_millimeter(const Panel &panel, 
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
  vec2<double> SimplePxMmStrategy::to_pixel(const Panel &panel, 
      vec2<double> xy) const {
    vec2<double> pixel_size = panel.get_pixel_size();
    return vec2<double> (xy[0] / pixel_size[0], xy[1] / pixel_size[1]);
  }
  
  /**
   * Convert a pixel coordinate to a millimeter coordinate
   * @param panel The panel structure
   * @param xy The (x, y) pixel coordinate
   * @return The (x, y) millimeter coordinate
   */  
  vec2<double> ParallaxCorrectedPxMmStrategy::to_millimeter(const Panel &panel, 
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
  vec2<double> ParallaxCorrectedPxMmStrategy::to_pixel(const Panel &panel, 
      vec2<double> xy) const {
    return SimplePxMmStrategy::to_pixel(panel, 
      parallax_correction(panel.get_distance(), la_, 
        panel.get_normal_origin(), xy));
  }

}} // namespace dxtbx::model
