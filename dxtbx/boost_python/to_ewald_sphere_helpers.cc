/*
 * to_ewald_sphere_helpers.cc
 *
 *  Copyright (C) 2013 Diamond Light Source
 *
 *  Author: James Parkhurst
 *
 *  This code is distributed under the BSD license, a copy of which is
 *  included in the root directory of this package.
 */
#include <boost/python.hpp>
#include <boost/python/def.hpp>
#include <scitbx/vec3.h>
#include <scitbx/array_family/flex_types.h>
#include <dxtbx/model/beam.h>
#include <dxtbx/model/detector.h>
#include <dxtbx/model/goniometer.h>
#include <dxtbx/model/scan.h>

namespace dxtbx { namespace boost_python {

  using namespace boost::python;
  
  using scitbx::vec2;
  using scitbx::vec3;
  using scitbx::af::flex_grid;
  using dxtbx::model::Beam;
  using dxtbx::model::Detector;
  using dxtbx::model::Goniometer;
  using dxtbx::model::Scan;

  typedef scitbx::af::flex<vec3<double> >::type flex_vec3_double;

  class ImageToEwaldSphere {
  public:
    ImageToEwaldSphere(const Beam &beam, const Detector &detector,
                       const Goniometer &gonio, const Scan &scan)
      : beam_(beam),
        detector_(detector),
        gonio_(gonio),
        scan_(scan) {}
    
    flex_vec3_double operator()(int frame, std::size_t panel) {

      // Check panel
      DXTBX_ASSERT(panel < detector_.size());

      // Get size and create array
      std::size_t slow_size = detector_[0].get_image_size()[1];
      std::size_t fast_size = detector_[0].get_image_size()[0];
      flex_vec3_double x(flex_grid<>(slow_size, fast_size));
      
      // Get rotation angle
      double phi = scan_.get_angle_from_array_index(frame - 0.5);    
      
      // Get coordinate for each pixel
      for (std::size_t j = 0; j < slow_size; ++j) {
        for (std::size_t i = 0; i < fast_size; ++i) {
          vec3<double> s1 = detector_[panel].get_pixel_lab_coord(vec2<double>(i, j));
          x(j, i) = s1.normalize().unit_rotate_around_origin(
            gonio_.get_rotation_axis(), phi) / beam_.get_wavelength();
        }
      }
      
      // Return array
      return x;
    }
    
    Beam beam_;
    Detector detector_;
    Goniometer gonio_;
    Scan scan_;
  };


  void export_to_ewald_sphere_helpers()
  {
    class_<ImageToEwaldSphere>("ImageToEwaldSphere", no_init)
      .def(init<const Beam&,
                const Detector&,
                const Goniometer&,
                const Scan&>((
          arg("beam"), 
          arg("detector"), 
          arg("goniometer"), 
          arg("scan"))))
      .def("__call__", &ImageToEwaldSphere::operator());  
  }

}} // namespace dxtbx::boost_python
