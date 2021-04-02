#include <simtbx/nanoBragg/nanoBragg.h>

//Contributed by Billy Poon, LBNL.

// function declaration from nanoBraggCUDA.cu
using simtbx::nanoBragg::shapetype;
extern "C"
void nanoBraggSpotsCUDA(int deviceId, int spixels, int fpixels, int roi_xmin, int roi_xmax,
                        int roi_ymin, int roi_ymax, int oversample,
                        int point_pixel, double pixel_size, double subpixel_size,
                        int steps, double detector_thickstep,
                        int detector_thicksteps, double detector_thick,
                        double detector_mu, double sdet_vector[4],
                        double fdet_vector[4], double odet_vector[4],
                        double pix0_vector[4], int curved_detector,
                        double distance, double close_distance,
                        double beam_vector[4], double Xbeam, double Ybeam,
                        double dmin, double phi0, double phistep, int phisteps,
                        double spindle_vector[4], int sources, double *source_X,
                        double *source_Y, double * source_Z, double * source_I,
                        double * source_lambda, double a0[4], double b0[4],
                        double c0[4], shapetype xtal_shape, double mosaic_spread,
                        int mosaic_domains, double * mosaic_umats, double Na,
                        double Nb, double Nc, double V_cell, double water_size,
                        double water_F, double water_MW, double r_e_sqr,
                        double fluence, double Avogadro, int integral_form,
                        double default_F, int interpolate, double *** Fhkl,
                        int h_min, int h_max, int h_range, int k_min, int k_max,
                        int k_range, int l_min, int l_max, int l_range, int hkls,
                        int nopolar, double polar_vector[4], double polarization,
                        double fudge, int unsigned short * maskimage,
                        float * floatimage /*out*/, double * omega_sum/*out*/,
                        int * sumn /*out*/, double * sum /*out*/,
                        double * sumsqr /*out*/, double * max_I/*out*/,
                        double * max_I_x/*out*/, double * max_I_y /*out*/, double spot_scale);

namespace simtbx {
namespace nanoBragg {

/* add spots from nanocrystal simulation */
void
nanoBragg::add_nanoBragg_spots_cuda()
{
  /*
     water_size not defined in class, CLI argument, defaults to 0
  */
  double water_size = 0.0;

  /*
     missing constants
  */
  double water_F = 2.57;
  double water_MW = 18.0;

  /* make sure we are normalizing with the right number of sub-steps */
  steps = phisteps*mosaic_domains*oversample*oversample;
  subpixel_size = pixel_size/oversample;

  /* declare a float version of floatimage for output */
  float* float_floatimage = new float[raw_pixels.size()];
#ifdef HAVE_NANOBRAGG_SPOTS_CUDA

  nanoBraggSpotsCUDA(device_Id, spixels, fpixels, roi_xmin, roi_xmax,
                     roi_ymin, roi_ymax, oversample,
                     point_pixel /* bool */, pixel_size, subpixel_size,
                     steps, detector_thickstep,
                     detector_thicksteps, detector_thick,
                     detector_attnlen /* detector_mu */, sdet_vector,
                     fdet_vector, odet_vector,
                     pix0_vector, curved_detector /* bool */,
                     distance, close_distance,
                     beam_vector, Xbeam, Ybeam,
                     dmin, phi0, phistep, phisteps,
                     spindle_vector, sources, source_X,
                     source_Y, source_Z, source_I,
                     source_lambda, a0, b0,
                     c0, xtal_shape, mosaic_spread,
                     mosaic_domains, mosaic_umats, Na,
                     Nb, Nc, V_cell, water_size,
                     water_F, water_MW, r_e_sqr,
                     fluence, Avogadro, integral_form,
                     default_F, interpolate, Fhkl,
                     h_min, h_max, h_range, k_min, k_max,
                     k_range, l_min, l_max, l_range, hkls,
                     nopolar /* bool */, polar_vector, polarization,
                     fudge, maskimage,
                     /* output */
                     float_floatimage, &omega_sum,
                     &sumn, &sum,
                     &sumsqr, &max_I,
                     &max_I_x, &max_I_y, spot_scale);
#else
  throw SCITBX_ERROR("no CUDA implementation of nanoBragg_add_spots");
#endif
  /* convert float_floatimage to double */
  for (int i=0; i<raw_pixels.size(); i++) {
    raw_pixels[i] = double(float_floatimage[i]);
  }
  delete[] float_floatimage;

  floatimage = raw_pixels.begin();
}
// end of add_nanoBragg_spots_cuda()

}}// namespace simtbx::nanoBragg
