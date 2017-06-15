#include <cctbx/boost_python/flex_fwd.h>

#include <boost/python/module.hpp>
#include <boost/python/scope.hpp>
#include <boost/python/class.hpp>
#include <boost/python/def.hpp>
#include <boost/python/tuple.hpp>
#include <boost/python/enum.hpp>
#include <boost/python/extract.hpp>

#include <simtbx/nanoBragg/nanoBragg.h>



using namespace boost::python;

namespace simtbx {
namespace nanoBragg {

namespace boost_python { namespace {

  boost::python::tuple
  testuple()
  {
    double test = NAN;
    return boost::python::make_tuple(isnan(test),2,3,4);
  }

  /* getter/setter functions for convenient python-side members
     here we apply and scaling required to initialize parameters in convenient units
     internal storage is always in meters and radians
     after most initializations we need to update inter-related parameters sensibly,
     such as pixel size and detector size when the number of detector pixels has already been set */

  /* square pixel size, exposed to Python in mm */
  static double get_pixel_size_mm(nanoBragg const& nanoBragg) {return nanoBragg.pixel_size*1000.;}
  static void   set_pixel_size_mm(nanoBragg& nanoBragg, double const& value) {
      nanoBragg.pixel_size = value/1000.;
      /* re-defines detector size */
      nanoBragg.init_detector();
  }

  /* number of pixels along each edge of detector */
  static scitbx::vec2<int> get_detpixels_fastslow(nanoBragg const& nanoBragg) {
      scitbx::vec2<int> value;
      value[0]=nanoBragg.fpixels;
      value[1]=nanoBragg.spixels;
      return value;
  }
  static void   set_detpixels_fastslow(nanoBragg& nanoBragg, scitbx::vec2<int> const& value) {
      nanoBragg.fpixels = value[0];
      nanoBragg.spixels = value[1];
      /* re-defines detector size */
      nanoBragg.detsize_f=nanoBragg.detsize_s=0.0;
      nanoBragg.init_detector();
  }

  /* size of each edge of detector in mm */
  static vec2 get_detsize_fastslow_mm(nanoBragg const& nanoBragg) {
      vec2 value;
      value[0]=nanoBragg.detsize_f*1000.0;
      value[1]=nanoBragg.detsize_s*1000.0;
      return value;
  }
  static void   set_detsize_fastslow_mm(nanoBragg& nanoBragg, vec2 const& value) {
      nanoBragg.detsize_f = value[0]/1000.0;
      nanoBragg.detsize_s = value[1]/1000.0;
      /* update number of pixels */
      nanoBragg.fpixels=nanoBragg.spixels=0;
      /* free the array? */
      nanoBragg.init_detector();
  }

  /* number of sub-pixels in each direction for oversampling */
  static int  get_oversample(nanoBragg const& nanoBragg) {return nanoBragg.oversample;}
  static void set_oversample(nanoBragg& nanoBragg, int const& value) {
      nanoBragg.oversample = value;
      /* don't keep re-setting this! */
      nanoBragg.user_oversample=true;
      nanoBragg.update_oversample();
      nanoBragg.update_steps();
}


  /* region of interest in pixels */
  static boost::python::tuple get_roi(nanoBragg const& nanoBragg) {
      boost::python::tuple frange;
      boost::python::tuple srange;
      frange = boost::python::make_tuple(nanoBragg.roi_xmin,nanoBragg.roi_xmax);
      srange = boost::python::make_tuple(nanoBragg.roi_ymin,nanoBragg.roi_ymax);

      return boost::python::make_tuple(frange,srange);
  }
  static void   set_roi(nanoBragg& nanoBragg, boost::python::tuple const& tupleoftuples) {
      boost::python::tuple frange;
      boost::python::tuple srange;
      frange = extract<boost::python::tuple>(tupleoftuples[0]);
      srange = extract<boost::python::tuple>(tupleoftuples[1]);
      nanoBragg.roi_xmin=extract<int>(frange[0]);
      nanoBragg.roi_xmax=extract<int>(frange[1]);
      nanoBragg.roi_ymin=extract<int>(srange[0]);
      nanoBragg.roi_ymax=extract<int>(srange[1]);
       /* fix this up so its not out-of-range */
      nanoBragg.init_beamcenter();
  }


  /* debugging option: select pixel to print out to screen */
  static scitbx::vec2<int> get_printout_pixel_fastslow(nanoBragg const& nanoBragg) {
      scitbx::vec2<int> value;
      value[0]=nanoBragg.printout_fpixel;
      value[1]=nanoBragg.printout_spixel;
      return value;
  }
  static void   set_printout_pixel_fastslow(nanoBragg& nanoBragg, scitbx::vec2<int> const& value) {
      nanoBragg.printout_fpixel = value[0];
      nanoBragg.printout_spixel = value[1];
  }



  /* beam center in some arbitrary convention (conventions defined below) */
  static vec2 get_beam_center_mm(nanoBragg const& nanoBragg) {
      vec2 value;
      value[0]=nanoBragg.Xbeam*1000.;
      value[1]=nanoBragg.Ybeam*1000.;
      return value;
  }
  static void   set_beam_center_mm(nanoBragg& nanoBragg, vec2 const& value) {
      nanoBragg.Xbeam = value[0]/1000.;
      nanoBragg.Ybeam = value[1]/1000.;
      nanoBragg.Fbeam = nanoBragg.Sbeam = NAN;
      nanoBragg.Xclose = nanoBragg.Yclose = NAN;
      nanoBragg.Fclose = nanoBragg.Sclose = NAN;
      nanoBragg.ORGX = nanoBragg.ORGY = NAN;
      nanoBragg.reconcile_parameters();
  }
  static vec2 get_XDS_ORGXY(nanoBragg const& nanoBragg) {
      vec2 value;
      value[0]=nanoBragg.ORGX;
      value[1]=nanoBragg.ORGY;
      return value;
  }
  static void   set_XDS_ORGXY(nanoBragg& nanoBragg, vec2 const& value) {
      nanoBragg.ORGX = value[0];
      nanoBragg.ORGY = value[1];
      nanoBragg.Fbeam = nanoBragg.Sbeam = NAN;
      nanoBragg.Xclose = nanoBragg.Yclose = NAN;
      nanoBragg.Fclose = nanoBragg.Sclose = NAN;
      nanoBragg.Xbeam = nanoBragg.Ybeam = NAN;
      nanoBragg.reconcile_parameters();
  }



  /* number of unit cells along edge in each cell axis direction */
  static scitbx::vec3<int> get_Nabc(nanoBragg const& nanoBragg) {
      scitbx::vec3<int> value;
      value[0]=static_cast<int>(nanoBragg.Na);
      value[1]=static_cast<int>(nanoBragg.Nb);
      value[2]=static_cast<int>(nanoBragg.Nc);
      return value;
  }
  static void   set_Nabc(nanoBragg& nanoBragg, scitbx::vec3<int> const& value) {
      nanoBragg.Na = value[0];
      nanoBragg.Nb = value[1];
      nanoBragg.Nc = value[2];
//      init_interpolator();
      nanoBragg.update_oversample();
      nanoBragg.update_steps();
  }


  /* specify overall crystal size instead of number of cells */
  static vec3 get_xtalsize_mm(nanoBragg const& nanoBragg) {
      vec3 value;
      value[0]=nanoBragg.sample_x*1000.;
      value[1]=nanoBragg.sample_y*1000.;
      value[2]=nanoBragg.sample_z*1000.;
      return value;
  }
  static void   set_xtalsize_mm(nanoBragg& nanoBragg, vec3 const& value) {
      nanoBragg.sample_x = value[0]/1000.;
      nanoBragg.sample_y = value[1]/1000.;
      nanoBragg.sample_z = value[2]/1000.;
      /* need to update oversampling */
      nanoBragg.update_oversample();
      nanoBragg.init_interpolator();
      nanoBragg.update_steps();
  }

  /* amorphous material thickness along beam (mm)  */
  static double get_amorphous_thick_mm(nanoBragg const& nanoBragg) {return nanoBragg.amorphous_thick*1000.;}
  static void   set_amorphous_thick_mm(nanoBragg& nanoBragg, double const& value) {
      nanoBragg.amorphous_thick = value/1000.;
  }
  /* amorphous material molecular weight (Da) */
  static double get_amorphous_MW_Da(nanoBragg const& nanoBragg) {return nanoBragg.amorphous_MW;}
  static void   set_amorphous_MW_Da(nanoBragg& nanoBragg, double const& value) {
      nanoBragg.amorphous_MW = value;
  }
  /* amorphous material density (g/cm^3) */
  static double get_amorphous_density_gcm3(nanoBragg const& nanoBragg) {return nanoBragg.amorphous_density*1e6;}
  static void   set_amorphous_density_gcm3(nanoBragg& nanoBragg, double const& value) {
      nanoBragg.amorphous_density = value*1e6;
  }
  /* now we can compute the scale factor of intensity from amorphous material structure factor */

  /* unit cell dimensions */
  static cctbx::uctbx::unit_cell get_cell_Adeg(nanoBragg const& nanoBragg) {
      scitbx::af::double6 six;
      six[0]=nanoBragg.a[0]*1e10; // convert from internal storage of meters and radians
      six[1]=nanoBragg.b[0]*1e10;
      six[2]=nanoBragg.c[0]*1e10;
      six[3]=nanoBragg.alpha*RTD;
      six[4]=nanoBragg.beta*RTD;
      six[5]=nanoBragg.gamma*RTD;
      cctbx::uctbx::unit_cell cell(six);
      return cell;
  }
  static void   set_cell_Adeg(nanoBragg& nanoBragg, cctbx::uctbx::unit_cell const& value) {
      nanoBragg.a_A[0] = value.parameters()[0];  // initialize Angstrom version, update at the end
      nanoBragg.b_A[0] = value.parameters()[1];
      nanoBragg.c_A[0] = value.parameters()[2];
      nanoBragg.alpha = value.parameters()[3]/RTD;
      nanoBragg.beta  = value.parameters()[4]/RTD;
      nanoBragg.gamma = value.parameters()[5]/RTD;
      nanoBragg.user_cell = 1;
      /* re-generate cell structure and apply any missetting angles to it */
      nanoBragg.init_cell();
//      reconcile_parameters();
  }


  /* unit cell dimensions as a tuple */
  static boost::python::tuple get_cell_astuple(nanoBragg const& nanoBragg) {
      double a,b,c,al,be,ga;
      a=nanoBragg.a[0]*1e10; // convert from internal storage of meters and radians
      b=nanoBragg.b[0]*1e10;
      c=nanoBragg.c[0]*1e10;
      al=nanoBragg.alpha*RTD;
      be=nanoBragg.beta*RTD;
      ga=nanoBragg.gamma*RTD;
      return boost::python::make_tuple(a,b,c,al,be,ga);
  }
  static void   set_cell_astuple(nanoBragg& nanoBragg, boost::python::tuple const& value) {
      nanoBragg.a_A[0] = extract<double>(value[0]);  // initialize Angstrom version, update at the end
      nanoBragg.b_A[0] = extract<double>(value[1]);
      nanoBragg.c_A[0] = extract<double>(value[2]);;
      nanoBragg.alpha  = extract<double>(value[3])/RTD;
      nanoBragg.beta   = extract<double>(value[4])/RTD;
      nanoBragg.gamma  = extract<double>(value[5])/RTD;
      nanoBragg.user_cell = 1;
      /* re-generate cell structure and apply any missetting angles to it */
      nanoBragg.init_cell();
//      reconcile_parameters();
  }


  /* use a MOSFLM-style A=UB matrix */
  static mat3 get_Amatrix(nanoBragg const& nanoBragg) {
      mat3 Amatrix;
      Amatrix(0,0)=nanoBragg.a_star[1]; // internal storage is reciprocal Angstrom
      Amatrix(0,1)=nanoBragg.a_star[2];
      Amatrix(0,2)=nanoBragg.a_star[3];
      Amatrix(1,0)=nanoBragg.b_star[1];
      Amatrix(1,1)=nanoBragg.b_star[2];
      Amatrix(1,2)=nanoBragg.b_star[3];
      Amatrix(2,0)=nanoBragg.c_star[1];
      Amatrix(2,1)=nanoBragg.c_star[2];
      Amatrix(2,2)=nanoBragg.c_star[3];
      return Amatrix;
  }
  static void   set_Amatrix(nanoBragg& nanoBragg, mat3 const& value) {
      nanoBragg.a_star[1] = value(0,0);  // initialize Angstrom version, update at the end
      nanoBragg.a_star[2] = value(0,1);
      nanoBragg.a_star[3] = value(0,2);
      nanoBragg.b_star[1] = value(1,0);
      nanoBragg.b_star[2] = value(1,1);
      nanoBragg.b_star[3] = value(1,2);
      nanoBragg.c_star[1] = value(2,0);
      nanoBragg.c_star[2] = value(2,1);
      nanoBragg.c_star[3] = value(2,2);
      nanoBragg.user_cell = 0;
      /* re-generate cell structure and apply any missetting angles to it */
      nanoBragg.init_cell();
//      reconcile_parameters();
  }



  /* crystal misseting angles, will be added after any matrix */
  static vec3 get_misset_deg(nanoBragg const& nanoBragg) {
      vec3 value;
      value[0]=nanoBragg.misset[1]*RTD;
      value[1]=nanoBragg.misset[2]*RTD;
      value[2]=nanoBragg.misset[3]*RTD;
      return value;
  }
  static void   set_misset_deg(nanoBragg& nanoBragg, vec3 const& value) {
      nanoBragg.misset[0] = 1;
      nanoBragg.misset[1] = value[0]/RTD;
      nanoBragg.misset[2] = value[1]/RTD;
      nanoBragg.misset[3] = value[2]/RTD;
      /* need to update A matrix */
      nanoBragg.init_cell();
  }


  /* detector misseting angles, will be applied */
  static vec3 get_detector_rot_deg(nanoBragg const& nanoBragg) {
      vec3 value;
      value[0]=nanoBragg.detector_rotx*RTD;
      value[1]=nanoBragg.detector_roty*RTD;
      value[2]=nanoBragg.detector_rotz*RTD;
      return value;
  }
  static void   set_detector_rot_deg(nanoBragg& nanoBragg, vec3 const& value) {
      nanoBragg.detector_rotx = value[0]/RTD;
      nanoBragg.detector_roty = value[1]/RTD;
      nanoBragg.detector_rotz = value[2]/RTD;
      /* need to update detector stuff */
      nanoBragg.init_detector();
  }

  /* twotheta is a special extra missetting angle, applied in addition to above, but pivoting about sample instead of beam spot */
  static double get_detector_twotheta_deg(nanoBragg const& nanoBragg) {
      return nanoBragg.detector_twotheta*RTD;;
  }
  static void   set_detector_twotheta_deg(nanoBragg& nanoBragg, double const& value) {
      nanoBragg.detector_twotheta = value/RTD;
      nanoBragg.detector_pivot = BEAM;
      /* need to update detector stuff */
      nanoBragg.init_detector();
  }

  /* sample-to-direct-beam-spot distance, exposed to Python in mm */
  static double get_distance_mm(nanoBragg const& nanoBragg) {return nanoBragg.distance*1000.;}
  static void   set_distance_mm(nanoBragg& nanoBragg, double const& value) {
      nanoBragg.distance = value/1000.;
      nanoBragg.close_distance = NAN;
      nanoBragg.reconcile_parameters();
  }
  /* sample-to-nearest-point-on-detector distance, exposed to Python in mm */
  static double get_close_distance_mm(nanoBragg const& nanoBragg) {return nanoBragg.close_distance*1000.;}
  static void   set_close_distance_mm(nanoBragg& nanoBragg, double const& value) {
      nanoBragg.close_distance = value/1000.;
      nanoBragg.distance = NAN;
      nanoBragg.reconcile_parameters();
  }



  /* unit vectors defining camera coordinate system */

  /* first, the X-ray beam, default is (1,0,0) */
  static vec3 get_beam_vector(nanoBragg const& nanoBragg) {
      vec3 value;
      value[0]=nanoBragg.beam_vector[1];
      value[1]=nanoBragg.beam_vector[2];
      value[2]=nanoBragg.beam_vector[3];
      return value;
  }
  static void   set_beam_vector(nanoBragg& nanoBragg, vec3 const& value) {
      nanoBragg.beam_vector[0] = 0;
      nanoBragg.beam_vector[1] = value[0];
      nanoBragg.beam_vector[2] = value[1];
      nanoBragg.beam_vector[3] = value[2];
      /* need to update everything? */
      nanoBragg.init_detector();
  }

  /* the fast direction of the detector pixel array, default is (0,0,1) */
  static vec3 get_fdet_vector(nanoBragg const& nanoBragg) {
      vec3 value;
      value[0]=nanoBragg.fdet_vector[1];
      value[1]=nanoBragg.fdet_vector[2];
      value[2]=nanoBragg.fdet_vector[3];
      return value;
  }
  static void   set_fdet_vector(nanoBragg& nanoBragg, vec3 const& value) {
      nanoBragg.fdet_vector[0] = 0;
      nanoBragg.fdet_vector[1] = value[0];
      nanoBragg.fdet_vector[2] = value[1];
      nanoBragg.fdet_vector[3] = value[2];
      /* need to update everything? */
      nanoBragg.init_detector();
  }

  /* the slow direction of the detector pixel array, default is (0,-1,0) */
  static vec3 get_sdet_vector(nanoBragg const& nanoBragg) {
      vec3 value;
      value[0]=nanoBragg.sdet_vector[1];
      value[1]=nanoBragg.sdet_vector[2];
      value[2]=nanoBragg.sdet_vector[3];
      return value;
  }
  static void   set_sdet_vector(nanoBragg& nanoBragg, vec3 const& value) {
      nanoBragg.sdet_vector[0] = 0;
      nanoBragg.sdet_vector[1] = value[0];
      nanoBragg.sdet_vector[2] = value[1];
      nanoBragg.sdet_vector[3] = value[2];
      /* need to update everything? */
      nanoBragg.init_detector();
  }

  /* the orthogonal vector to the detector pixel array, default is (1,0,0), same as the beam */
  /* note that this is updated every time using sdet_vector and fdet_vector bu init_detector() */
  static vec3 get_odet_vector(nanoBragg const& nanoBragg) {
      vec3 value;
      value[0]=nanoBragg.odet_vector[1];
      value[1]=nanoBragg.odet_vector[2];
      value[2]=nanoBragg.odet_vector[3];
      return value;
  }
  static void   set_odet_vector(nanoBragg& nanoBragg, vec3 const& value) {
      nanoBragg.odet_vector[0] = 0;
      nanoBragg.odet_vector[1] = value[0];
      nanoBragg.odet_vector[2] = value[1];
      nanoBragg.odet_vector[3] = value[2];
      /* need to update everything? */
      nanoBragg.init_detector();
  }

  /* the direction vector of detector_twotheta rotation axis (about the sample), default is (0,1,0) */
  static vec3 get_twotheta_axis(nanoBragg const& nanoBragg) {
      vec3 value;
      value[0]=nanoBragg.twotheta_axis[1];
      value[1]=nanoBragg.twotheta_axis[2];
      value[2]=nanoBragg.twotheta_axis[3];
      return value;
  }
  static void   set_twotheta_axis(nanoBragg& nanoBragg, vec3 const& value) {
      nanoBragg.twotheta_axis[0] = 0;
      nanoBragg.twotheta_axis[1] = value[0];
      nanoBragg.twotheta_axis[2] = value[1];
      nanoBragg.twotheta_axis[3] = value[2];
      /* need to update everything? */
      nanoBragg.init_detector();
  }

  /* direction of the E vector of the X-ray beam. will be made normal to beam direction, default: (0,0,1) */
  static vec3 get_polar_vector(nanoBragg const& nanoBragg) {
      vec3 value;
      value[0]=nanoBragg.polar_vector[1];
      value[1]=nanoBragg.polar_vector[2];
      value[2]=nanoBragg.polar_vector[3];
      return value;
  }
  static void   set_polar_vector(nanoBragg& nanoBragg, vec3 const& value) {
      nanoBragg.polar_vector[0] = 0;
      nanoBragg.polar_vector[1] = value[0];
      nanoBragg.polar_vector[2] = value[1];
      nanoBragg.polar_vector[3] = value[2];
      /* need to update everything? */
      nanoBragg.init_detector();
  }

  /* the direction vector of the phi spindle rotation axis, default is (0,0,1) */
  static vec3 get_spindle_vector(nanoBragg const& nanoBragg) {
      vec3 value;
      value[0]=nanoBragg.spindle_vector[1];
      value[1]=nanoBragg.spindle_vector[2];
      value[2]=nanoBragg.spindle_vector[3];
      return value;
  }
  static void   set_spindle_vector(nanoBragg& nanoBragg, vec3 const& value) {
      nanoBragg.spindle_vector[0] = 0;
      nanoBragg.spindle_vector[1] = value[0];
      nanoBragg.spindle_vector[2] = value[1];
      nanoBragg.spindle_vector[3] = value[2];
      /* need to update everything? */
      nanoBragg.init_detector();
  }




  /* hkl, without F */
  static indices get_indices(nanoBragg nanoBragg) {
      int i,h,k,l;
      /* initialize python from built-in ***Fhkl array */
      nanoBragg.pythony_indices = indices();
      int hkls = nanoBragg.h_range*nanoBragg.k_range*nanoBragg.l_range;
      nanoBragg.pythony_indices.resize(hkls);
      i=0;
      for(h=nanoBragg.h_min;h<=nanoBragg.h_max;++h){
          for(k=nanoBragg.k_min;k<=nanoBragg.k_max;++k){
              for(l=nanoBragg.l_min;l<=nanoBragg.l_max;++l){

                  if ( (h<=nanoBragg.h_max) && (h>=nanoBragg.h_min) && (k<=nanoBragg.k_max)
                    && (k>=nanoBragg.k_min) && (l<=nanoBragg.l_max) && (l>=nanoBragg.l_min)  ) {
                      /* populate only valid variables */
                      miller_t hkl (h,k,l);
                      nanoBragg.pythony_indices[i] = hkl;
                      ++i;
                  }
              }
          }
      }
      return nanoBragg.pythony_indices;
  }
  static void   set_indices(nanoBragg& nanoBragg, indices const& value) {
      if(nanoBragg.verbose>3) printf(" about to initialize indices\n");
      nanoBragg.pythony_indices = value;
      if(nanoBragg.verbose>3) printf(" done\n");

      /* re-initialize */
      nanoBragg.init_Fhkl();
  }


  /* F without hkl : translate F values, hope number is the same as miller indices */
  static af::shared<double> get_amplitudes(nanoBragg nanoBragg) {
      int i,h,k,l;
      nanoBragg.pythony_amplitudes = af::shared<double>();
      int hkls = nanoBragg.h_range*nanoBragg.k_range*nanoBragg.l_range;
      nanoBragg.pythony_amplitudes.resize(hkls);
      i=0;
      for(h=nanoBragg.h_min;h<=nanoBragg.h_max;++h){
          for(k=nanoBragg.k_min;k<=nanoBragg.k_max;++k){
              for(l=nanoBragg.l_min;l<=nanoBragg.l_max;++l){

                  if ( (h<=nanoBragg.h_max) && (h>=nanoBragg.h_min) && (k<=nanoBragg.k_max)
                    && (k>=nanoBragg.k_min) && (l<=nanoBragg.l_max) && (l>=nanoBragg.l_min)  ) {
                      /* populate only valid values */
                      nanoBragg.pythony_amplitudes[i] = nanoBragg.Fhkl[h-nanoBragg.h_min][k-nanoBragg.k_min][l-nanoBragg.l_min];
                      ++i;
                  }
              }
          }
      }
      return nanoBragg.pythony_amplitudes;
  }
  static void   set_amplitudes(nanoBragg& nanoBragg, af::shared<double> const& value) {
      if(nanoBragg.verbose>3) printf(" about to initialize amplitudes\n");
      nanoBragg.pythony_amplitudes = value;
      if(nanoBragg.verbose>3) printf(" done\n");

      /* re-initialize ***Fhkl array */
      nanoBragg.init_Fhkl();
  }


  /* F and hkl as a tuple: translate F values, hope number is the same as miller indices */
  static boost::python::tuple get_Fhkl_tuple(nanoBragg nanoBragg) {
      int h,k,l;
      double temp;
      int hkls = nanoBragg.h_range*nanoBragg.k_range*nanoBragg.l_range;
      if(nanoBragg.verbose>3) printf(" expecting %d hkls\n",hkls);
//      nanoBragg.pythony_indices = indices();
//      nanoBragg.pythony_amplitudes = af::shared<double>();
      nanoBragg.pythony_indices = indices(hkls,af::init_functor_null<scitbx::vec3<int> >());
      nanoBragg.pythony_amplitudes = af::shared<double>(hkls,af::init_functor_null<double>());
//      nanoBragg.pythony_indices = indices(hkls);
//      nanoBragg.pythony_amplitudes = af::shared<double>(hkls);
//      nanoBragg.pythony_indices.resize(hkls);
//      nanoBragg.pythony_amplitudes.resize(hkls);
      int i=0;

      for(h=nanoBragg.h_min;h<=nanoBragg.h_max;++h){
          for(k=nanoBragg.k_min;k<=nanoBragg.k_max;++k){
              for(l=nanoBragg.l_min;l<=nanoBragg.l_max;++l){

                  if ( (h<=nanoBragg.h_max) && (h>=nanoBragg.h_min) && (k<=nanoBragg.k_max)
                    && (k>=nanoBragg.k_min) && (l<=nanoBragg.l_max) && (l>=nanoBragg.l_min)  ) {
                      /* populate all stored values */
                      miller_t hkl (h,k,l);
                      if(nanoBragg.verbose>3) printf(" about to access indices[%d] at %p\n",i,&nanoBragg.pythony_indices[i]);
                      nanoBragg.pythony_indices[i] = hkl;
//                      nanoBragg.pythony_indices.push_back(hkl);
                      if(nanoBragg.verbose>3) printf(" about to access (%d,%d,%d) Fhkl[%d][%d][%d]\n",
                                           h,k,l,h-nanoBragg.h_min,k-nanoBragg.k_min,l-nanoBragg.l_min);
                      temp = nanoBragg.Fhkl[h-nanoBragg.h_min][k-nanoBragg.k_min][l-nanoBragg.l_min];
                      if(nanoBragg.verbose>3) printf(" about to access amplitudes[%d] at %p\n",
                                            i,&nanoBragg.pythony_amplitudes[i]);
                      nanoBragg.pythony_amplitudes[i] = temp;
//                      nanoBragg.pythony_amplitudes.push_back(temp);
                      if(nanoBragg.verbose>3) printf(" done\n");
                      ++i;
                  }
              }
          }
      }
      /* return a tuple so there is no confusion about order of initialization */
      if(nanoBragg.verbose>3) SCITBX_EXAMINE(nanoBragg.pythony_indices.size());
      return boost::python::make_tuple(nanoBragg.pythony_indices,nanoBragg.pythony_amplitudes);
  }
  static void   set_Fhkl_tuple(nanoBragg& nanoBragg, boost::python::tuple const& value) {
//      nanoBragg.pythony_indices = boost::python::tuple::get<0>(value);
//      nanoBragg.pythony_amplitudes = boost::python::tuple::get<1>(value);
//      boost::python::extract( nanoBragg.pythony_indices , nanoBragg.pythony_amplitudes ) = value;
      if(nanoBragg.verbose>3) printf(" about to initialize indices\n");
//printf("GOTHERE extract<indices>(value[0]).size()=%d\n",extract<indices >(value[0]).size());
      nanoBragg.pythony_indices = extract<indices >(value[0]);
      if(nanoBragg.verbose>3) printf(" about to initialize amplitudes\n");
      nanoBragg.pythony_amplitudes = extract<af::shared<double> >(value[1]);

      /* re-initialize ***Fhkl array */
      nanoBragg.init_Fhkl();
  }

  /* allow default value for missing Fs */
  static double get_default_F(nanoBragg const& nanoBragg) {return nanoBragg.default_F;}
  static void   set_default_F(nanoBragg& nanoBragg, double const& value) {
      nanoBragg.default_F = value;
      /* make sure this gets applied! */
      nanoBragg.init_Fhkl();
  }



  /* X-ray wavelength */
  static double get_lambda_A(nanoBragg const& nanoBragg) {return nanoBragg.lambda0*1e10;}
  static void   set_lambda_A(nanoBragg& nanoBragg, double const& value) {nanoBragg.lambda0 = value/1e10;}
  /* X-ray wavelength from energy */
  static double get_energy_eV(nanoBragg const& nanoBragg) {return 12398.42/(nanoBragg.lambda0*1e10);}
  static void   set_energy_eV(nanoBragg& nanoBragg, double const& value) {nanoBragg.lambda0 = (12398.42/value)/1e10;}

  /* X-ray fluence in photons/meter^2 */
  static double get_fluence(nanoBragg const& nanoBragg) {
      return nanoBragg.fluence;
  }
  static void   set_fluence(nanoBragg& nanoBragg, double const& value) {
      nanoBragg.fluence = value;
      /* force update of the flux... */
      nanoBragg.flux = 0.0;
      nanoBragg.init_beam();
  }

  /* X-ray flux in photons/s */
  static double get_flux(nanoBragg const& nanoBragg) {
      return nanoBragg.flux;
  }
  static void   set_flux(nanoBragg& nanoBragg, double const& value) {
      nanoBragg.flux = value;
      /* re-defines fluence */
      nanoBragg.init_beam();
  }

  /* X-ray exposure time, mainly defines fluence */
  static double get_exposure_s(nanoBragg const& nanoBragg) {
      return nanoBragg.exposure;
  }
  static void   set_exposure_s(nanoBragg& nanoBragg, double const& value) {
      nanoBragg.exposure = value;
      /* re-defines fluence */
      nanoBragg.init_beam();
  }

  /* X-ray beam size, mainly defines fluence */
  static double get_beamsize_mm(nanoBragg const& nanoBragg) {
      return nanoBragg.beamsize*1000.0;
  }
  static void   set_beamsize_mm(nanoBragg& nanoBragg, double const& value) {
      nanoBragg.beamsize = value;
      /* re-defines fluence */
      nanoBragg.init_beam();
  }


  /* x-ray beam spectral dispersion in % */
  static double get_dispersion_pct(nanoBragg const& nanoBragg) {
      return nanoBragg.dispersion*100.0;
  }
  static void   set_dispersion_pct(nanoBragg& nanoBragg, double const& value) {
      nanoBragg.dispersion = value/100.0;
      /* need to re-create source table */
      nanoBragg.init_steps();
      nanoBragg.init_sources();
      nanoBragg.update_steps();
  }

  /* number of discrete wavelengths in X-ray source, default: 1 or 2 ; you will want more */
  static int get_dispsteps(nanoBragg const& nanoBragg) {
      return nanoBragg.dispsteps;
  }
  static void   set_dispsteps(nanoBragg& nanoBragg, int const& value) {
      nanoBragg.dispsteps = value;
//      nanoBragg.dispersion=-1.0;
      /* need to re-create source table */
      nanoBragg.init_steps();
      nanoBragg.init_sources();
      nanoBragg.update_steps();
  }




  /* x-ray beam horizontal,vertical divergence in mrad */
  static vec2 get_divergence_hv_mrad(nanoBragg const& nanoBragg) {
      vec2 value;
      value[0]=nanoBragg.hdivrange*1000.0;
      value[1]=nanoBragg.vdivrange*1000.0;
      return value;
  }
  static void   set_divergence_hv_mrad(nanoBragg& nanoBragg, vec2 const& value) {
      nanoBragg.hdivrange = value[0]/1000.0;
      nanoBragg.vdivrange = value[1]/1000.0;
      /* re-determine stepsize */
      nanoBragg.hdivstep=nanoBragg.vdivstep=-1.0;
      /* need to re-create source table */
      nanoBragg.init_steps();
      nanoBragg.init_sources();
      nanoBragg.update_steps();
  }

  /* number of discrete incidence angles in X-ray source, default: 1 or 2 ; you will want more */
  static scitbx::vec2<int> get_divsteps_hv(nanoBragg const& nanoBragg) {
      scitbx::vec2<int> value;
      value[0]=nanoBragg.hdivsteps;
      value[1]=nanoBragg.vdivsteps;
      return value;
  }
  static void   set_divsteps_hv(nanoBragg& nanoBragg, scitbx::vec2<int> const& value) {
      nanoBragg.hdivsteps = value[0];
      nanoBragg.vdivsteps = value[1];
      /* need to re-create source table */
      nanoBragg.vdivstep = nanoBragg.hdivstep = -1.0;
      nanoBragg.init_steps();
      nanoBragg.init_sources();
      nanoBragg.update_steps();
  }

  /* alternative to steps: x-ray beam divergence step size, horizontal,vertical in mrad */
  static vec2 get_divstep_hv_mrad(nanoBragg const& nanoBragg) {
      vec2 value;
      value[0]=nanoBragg.hdivstep*1000.0;
      value[1]=nanoBragg.vdivstep*1000.0;
      return value;
  }
  static void   set_divstep_hv_mrad(nanoBragg& nanoBragg, vec2 const& value) {
      nanoBragg.hdivstep = value[0]/1000.0;
      nanoBragg.vdivstep = value[1]/1000.0;
      /* need to re-create source table */
      nanoBragg.vdivsteps = nanoBragg.hdivsteps = -1;
      nanoBragg.init_steps();
      nanoBragg.init_sources();
      nanoBragg.update_steps();
  }





  /* crystal mosaic spread, in deg */
  static double get_mosaic_spread_deg(nanoBragg const& nanoBragg) {
      return nanoBragg.mosaic_spread*RTD;;
  }
  static void   set_mosaic_spread_deg(nanoBragg& nanoBragg, double const& value) {
      nanoBragg.mosaic_spread = value/RTD;
      /* need to re-create mosaic domain table */
      nanoBragg.init_mosaicity();
      nanoBragg.init_steps();
      nanoBragg.update_steps();
      nanoBragg.show_mosaic_blocks();
  }
  /* number of discrete mosaic domains, default: 1 or 2 ; you will want more */
  static int get_mosaic_domains(nanoBragg const& nanoBragg) {
      return nanoBragg.mosaic_domains;
  }
  static void   set_mosaic_domains(nanoBragg& nanoBragg, int const& value) {
      nanoBragg.mosaic_domains = value;
      /* need to re-create mosaic domain table */
      nanoBragg.init_mosaicity();
      nanoBragg.init_steps();
      nanoBragg.update_steps();
      nanoBragg.show_mosaic_blocks();
  }




  /* spindle starting angle phi, in deg */
  static double get_phi_deg(nanoBragg const& nanoBragg) {
      return nanoBragg.phi0*RTD;;
  }
  static void   set_phi_deg(nanoBragg& nanoBragg, double const& value) {
      nanoBragg.phi0 = value/RTD;
      /* need to re-create phi step table */
      nanoBragg.init_steps();
      nanoBragg.update_steps();
      nanoBragg.show_phisteps();
  }

  /* spindle rotation ange, in deg */
  static double get_osc_deg(nanoBragg const& nanoBragg) {
      return nanoBragg.osc*RTD;;
  }
  static void   set_osc_deg(nanoBragg& nanoBragg, double const& value) {
      nanoBragg.osc = value/RTD;
      /* need to re-create phi step table */
      nanoBragg.phistep=-1.0;
      nanoBragg.init_steps();
      nanoBragg.update_steps();
      nanoBragg.show_phisteps();
  }

  /* number of discrete phi rotations, default: 1 or 2 ; you will want more */
  static int get_phisteps(nanoBragg const& nanoBragg) {
      return nanoBragg.phisteps;
  }
  static void   set_phisteps(nanoBragg& nanoBragg, int const& value) {
      nanoBragg.phisteps = value;
      /* need to re-create phi step table */
      nanoBragg.phistep=-1.0;
      nanoBragg.init_steps();
      nanoBragg.update_steps();
      nanoBragg.show_phisteps();
  }

  /* spindle angle phi step, in deg */
  static double get_phistep_deg(nanoBragg const& nanoBragg) {
      return nanoBragg.phi0*RTD;;
  }
  static void   set_phistep_deg(nanoBragg& nanoBragg, double const& value) {
      nanoBragg.phistep = value/RTD;
      /* need to re-create phi step table */
      nanoBragg.init_steps();
      nanoBragg.update_steps();
      nanoBragg.show_phisteps();
  }


  /* hijack detctor_pivot as an int */
  static int get_hijack_detector_pivot(nanoBragg const& nanoBragg) {
      if(nanoBragg.detector_pivot == SAMPLE) return 1;
      return 0;
  }
  static void   set_hijack_detector_pivot(nanoBragg& nanoBragg, int const& value) {
      if(value == 0) nanoBragg.detector_pivot = BEAM;
      if(value == 1) nanoBragg.detector_pivot = SAMPLE;
  }


  void
  init_module() {
    using namespace boost::python;

    typedef return_value_policy<return_by_value> rbv;
    typedef default_call_policies dcp;

    def("testuple", &testuple,
        "test that basic python stuff works, as well as NAN implementation");

    class_<nanoBragg>("nanoBragg",no_init)
      /* constructor that takes a dials detector model */
      .def(init<const dxtbx::model::Detector&>(
        (arg_("detector")),
        "nanoBragg simulation initialized from a dxtbx detector object"))

       /* constructor that takes any and all parameters with sensible defaults */
      .def(init<scitbx::vec2<int>,
                scitbx::vec3<int>,
                cctbx::uctbx::unit_cell,
                vec3,
                vec2,
                double,
                double,
                double,
                double,
                double,
                double,
                int,
                int>(
        (arg_("detpixels_slowfast")=scitbx::vec2<int>(1024,1024),
         arg_("Ncells_abc")=scitbx::vec3<int>(1,1,1),
         arg_("unit_cell_Adeg")=cctbx::uctbx::unit_cell(scitbx::af::double6(78.0,78.0,38.0,90.0,90.0,90.0)),
         arg_("misset_deg")=vec3(0.0,0.0,0.0),
         arg_("beam_center_mm")=vec2(NAN,NAN),
         arg_("distance_mm")=100,
         arg_("pixel_size_mm")=0.1,
         arg_("wavelength_A")=1,
         arg_("divergence_mrad")=0,
         arg_("dispersion_pct")=0,
         arg_("mosaic_spread_deg")=0,
         arg_("oversample")=0,
         arg_("verbose")=1),
         "nanoBragg simulation initialized with most common parameters. All parameters have sensible defaults."))

      /* toggle printing stuff to screen - does not work completely yet */
      .add_property("verbose",
                     make_getter(&nanoBragg::verbose,rbv()),
                     make_setter(&nanoBragg::verbose,dcp()),
                     "non-zero turns on debugging messages of increasing verbosity.")
      /* toggle detailed print-out of each and every pixel value and parameters */
      .add_property("printout",
                     make_getter(&nanoBragg::printout,rbv()),
                     make_setter(&nanoBragg::printout,dcp()),
                     "debugging option: print to screen all values associated with a pixel")
      /* select a single pixel to be verbose about */
      .add_property("printout_pixel_fastslow",
                     make_function(&get_printout_pixel_fastslow,rbv()),
                     make_function(&set_printout_pixel_fastslow,dcp()),
                     "fast/slow pixel coordinate of particular pixel to print to screen.  Default is all if printout=True.")
      /* toggle screen printing of progrees in percent */
      .add_property("progress_meter",
                     make_getter(&nanoBragg::progress_meter,rbv()),
                     make_setter(&nanoBragg::progress_meter,dcp()),
                     "toggle screen printing of simulation progress in percent.")



      /* number of unit cells in each direction */
      .add_property("Ncells_abc",
                     make_function(&get_Nabc,rbv()),
                     make_function(&set_Nabc,dcp()),
                     "number of unit cells along each unit cell axis")
      /* alternatively, specify crystal size in mm instead of number of cells */
      .add_property("xtalsize_mm",
                     make_function(&get_xtalsize_mm,rbv()),
                     make_function(&set_xtalsize_mm,dcp()),
                     "alternative: specify crystal size (mm) instead of Ncells_abc")



      /* specify properties of amorphous material contribution to background */
      .add_property("amorphous_thick_mm",
                     make_function(&get_amorphous_thick_mm,rbv()),
                     make_function(&set_amorphous_thick_mm,dcp()),
                     "thickness of amorphous material (mm), most important along beam path")
      .add_property("amorphous_density_gcm3",
                     make_function(&get_amorphous_density_gcm3,rbv()),
                     make_function(&set_amorphous_density_gcm3,dcp()),
                     "density of amorphous material (g/cm^3), needed to calculate scale factor")
      .add_property("amorphous_MW_Da",
                     make_function(&get_amorphous_MW_Da,rbv()),
                     make_function(&set_amorphous_MW_Da,dcp()),
                     "molecular weight of amorphous material (Da, or g/mol), needed to calculate scale factor")
//      .add_property("amorphous_F_vs_stol",
//                     make_function(&get_amorphous_F_vs_stol,rbv()),
//                     make_function(&set_amorphous_F_vs_stol,dcp()),
//                     "structure factor of amorphous material (electrons/molecule) vs sin(theta)/lambda (Angstrom), same as tabulated for individual atom form factors in International Tables")



      /* unit cell dimensions */
      .add_property("unit_cell_Adeg",
                     make_function(&get_cell_Adeg,rbv()),
                     make_function(&set_cell_Adeg,dcp()),
                     "uctbx unit cell")
      .add_property("unit_cell_tuple",
                     make_function(&get_cell_astuple,rbv()),
                     make_function(&set_cell_astuple,dcp()),
                     "unit cell as a 6-tuple")

      /* crystal orientation as missetting angles (deg) */
      .add_property("missets_deg",
                     make_function(&get_misset_deg,rbv()),
                     make_function(&set_misset_deg,dcp()),
                     "crystal mis-orientation about x,y,z axes as missetting angles (deg)")
      /* crystal orientation as an A=UB matrix */
      .add_property("Amatrix",
                     make_function(&get_Amatrix,rbv()),
                     make_function(&set_Amatrix,dcp()),
                     "unit cell, wavelength, and crystal orientation as Arndt & Wonacott A=UB matrix")

      /* beam center, in convention specified below */
      .add_property("beam_center_mm",
                     make_function(&get_beam_center_mm,rbv()),
                     make_function(&set_beam_center_mm,dcp()),
                     "x-ray beam center (mm) in specified convention (dials, adxv, mosflm, denzo, XDS)")
      .add_property("XDS_ORGXY",
                     make_function(&get_XDS_ORGXY,rbv()),
                     make_function(&set_XDS_ORGXY,dcp()),
                     "XDS-convention beam center is nearest pixel in detector plane to sample position (fast,slow)")

      /* specify the detector pivot point, note: this is an enum */
      .add_property("detector_pivot",
                     make_getter(&nanoBragg::detector_pivot,rbv()),
                     make_setter(&nanoBragg::detector_pivot,dcp()),
                     "specify the center of detector rotations, either the sample, or the direct-beam spot")

      /* specify the detector pivot point, note: this is an enum */
      .add_property("hijack_detector_pivot",
                     make_function(&get_hijack_detector_pivot,rbv()),
                     make_function(&set_hijack_detector_pivot,dcp()),
                     "hijack setting detector pivot from a string becaues python enums suck ")



      /* specify convention? MOSFLM, XDS, Denzo, ADXV, DIALS */

      /* unit vectors specifying coordinate system */
      .add_property("fdet_vector",
                     make_function(&get_fdet_vector,rbv()),
                     make_function(&set_fdet_vector,dcp()),
                     "unit vector representing 3-space direction of increasing fast-pixel-coordinate values")
      .add_property("sdet_vector",
                     make_function(&get_sdet_vector,rbv()),
                     make_function(&set_sdet_vector,dcp()),
                     "unit vector representing 3-space direction of increasing slow-pixel-coordinate values")
      .add_property("odet_vector",
                     make_function(&get_odet_vector,rbv()),
                     make_function(&set_odet_vector,dcp()),
                     "unit vector representing 3-space direction of increasing detector distance")
      .add_property("beam_vector",
                     make_function(&get_beam_vector,rbv()),
                     make_function(&set_beam_vector,dcp()),
                     "unit vector representing 3-space direction of x-ray beam propagation (the k-vector)")
      .add_property("polar_vector",
                     make_function(&get_polar_vector,rbv()),
                     make_function(&set_polar_vector,dcp()),
                     "unit vector representing 3-space direction of x-ray beam polarization (the E-vector)")
      .add_property("spindle_axis",
                     make_function(&get_spindle_vector,rbv()),
                     make_function(&set_spindle_vector,dcp()),
                     "unit vector representing 3-space direction of crystal phi spindle rotation axis")
      .add_property("twotheta_axis",
                     make_function(&get_twotheta_axis,rbv()),
                     make_function(&set_twotheta_axis,dcp()),
                     "unit vector representing 3-space direction of detector twotheta rotation axis (about sample)")

      /* sample-to-direct-beam-spot distance */
      .add_property("distance_meters",
                     make_getter(&nanoBragg::distance,rbv()),
                     make_setter(&nanoBragg::distance,dcp()),
                     "sample-to-detector distance (meters)")
      .add_property("distance_mm",
                     make_function(&get_distance_mm,rbv()),
                     make_function(&set_distance_mm,dcp()),
                     "sample-to-detector distance (mm)")
      /* sample-to-nearest-point-on-detector distance */
      .add_property("close_distance_mm",
                     make_function(&get_close_distance_mm,rbv()),
                     make_function(&set_close_distance_mm,dcp()),
                     "distance between sample and nearest point in the detector pixel plane (mm)")

      /* rotation of detector about sample */
      .add_property("detector_twotheta_deg",
                     make_function(&get_detector_twotheta_deg,rbv()),
                     make_function(&set_detector_twotheta_deg,dcp()),
                     "rotation of detector about twotheta_axis at sample position")

      /* detector face size in mm */
      .add_property("detsize_fastslow_mm",
                     make_function(&get_detsize_fastslow_mm,rbv()),
                     make_function(&set_detsize_fastslow_mm,dcp()),
                     "overall size of the detector fast and slow pixel directions (mm)")
      /* detector face size in pixels */
      .add_property("detpixels_fastslow",
                     make_function(&get_detpixels_fastslow,rbv()),
                     make_function(&set_detpixels_fastslow,dcp()),
                     "overall size of the detector fast and slow pixel directions (pixels)")

      /* rotation of detector about direct-beam-spot or sample, depending on pivot */
      .add_property("detector_rot_deg",
                     make_function(&get_detector_rot_deg,rbv()),
                     make_function(&set_detector_rot_deg,dcp()),
                     "rotation of detector (deg) about direct-beam-spot (or sample position, depending on pivot)")

      /* True forces all pixels to be same distance from sample */
      .add_property("curved_detector",
                     make_getter(&nanoBragg::curved_detector,rbv()),
                     make_setter(&nanoBragg::curved_detector,dcp()),
                     "True forces all pixels to be same distance from sample")

      /* detector pixel size in mm */
      .add_property("pixel_size_mm",
                     make_function(&get_pixel_size_mm,rbv()),
                     make_function(&set_pixel_size_mm,dcp()),
                     "square detector pixel size (mm)")

      /* True turns solid angle effect off */
      .add_property("point_pixel",
                     make_getter(&nanoBragg::point_pixel,rbv()),
                     make_setter(&nanoBragg::point_pixel,dcp()),
                     "True turns pixel solid angle effect off, pixels behave like infinitely sharp points")

      /* set Kahn polarization factor between -1 and 1, default: 0  - note polarization=0 is DIFFERENT from nopolar */
      .add_property("polarization",
                     make_getter(&nanoBragg::polarization,rbv()),
                     make_setter(&nanoBragg::polarization,dcp()),
                     "set Kahn polarization factor between -1 and 1, default: 0  - note polarization=0 is DIFFERENT from nopolar")
      /* True turns polarization effect off */
      .add_property("nopolar",
                     make_getter(&nanoBragg::nopolar,rbv()),
                     make_setter(&nanoBragg::nopolar,dcp()),
                     "True turns polarization effect off, F^2 scattering in all directions")

     /* override oversampling, default: reciprocal sub-pixel fits xtalsize */
      .add_property("oversample",
                     make_function(&get_oversample,rbv()),
                     make_function(&set_oversample,dcp()),
                     "override pixel oversampling, speed is inversely proportional to square of this quantity, but critical if reciprocal crystal size is smaller than a pixel. default: set oversample so that reciprocal sub-pixel fits xtalsize")

     /* select a region-of-interest: xmin xmax ymin ymax */
      .add_property("region-of-interest",
                     make_function(&get_roi,rbv()),
                     make_function(&set_roi,dcp()),
                     "region of interst on detector: fast_min fast_max slow_min slow_max")

      /* wavelength in Angstrom */
      .add_property("wavelength_A",
                     make_function(&get_lambda_A,rbv()),
                     make_function(&set_lambda_A,dcp()),
                     "centroid x-ray wavelength (Angstrom)")
      /* or photon energy */
      .add_property("energy_eV",
                     make_function(&get_energy_eV,rbv()),
                     make_function(&set_energy_eV,dcp()),
                     "alternately specify x-ray wavelength as a photon energy (eV)")

      /* override value of x-ray fluence in photons/m^2 : default is sufficient for intensity=F^2 */
      .add_property("fluence",
                     make_function(&get_fluence,rbv()),
                     make_function(&set_fluence,dcp()),
                     "x-ray fluence upon sample (photons/meter^2), default makes photons/steradian=F^2 exactly")
      /* specify x-ray flux in photons/s, updates fluence */
      .add_property("flux",
                     make_function(&get_flux,rbv()),
                     make_function(&set_flux,dcp()),
                     "x-ray flux (photons/s), mainly to update fluence")
      /* exposure time, updates fluence */
      .add_property("exposure_s",
                     make_function(&get_exposure_s,rbv()),
                     make_function(&set_exposure_s,dcp()),
                     "x-ray exposure time (seconds), mainly to update fluence")
      /* x-ray beam size in mm, updates fluence */
      .add_property("beamsize_mm",
                     make_function(&get_beamsize_mm,rbv()),
                     make_function(&set_beamsize_mm,dcp()),
                     "x-ray beam size in both directions (mm), mainly to update fluence")

      /* x-ray beam spectral dispersion in percent */
      .add_property("dispersion_pct",
                     make_function(&get_dispersion_pct,rbv()),
                     make_function(&set_dispersion_pct,dcp()),
                     "x-ray beam spectral dispersion (%)")
      /* spectral dispersion steps */
      .add_property("dispsteps",
                     make_function(&get_dispsteps,rbv()),
                     make_function(&set_dispsteps,dcp()),
                     "number of micro-steps across wavelength range, speed is inversely proportional to this")

      /* x-ray beam divergence in mrad */
      .add_property("divergence_hv_mrad",
                     make_function(&get_divergence_hv_mrad,rbv()),
                     make_function(&set_divergence_hv_mrad,dcp()),
                     "x-ray beam divergence in horizontal and vertical directions (mrad)")
      /* beam divergence micro-steps in wavelength sweep */
      .add_property("divsteps_hv",
                     make_function(&get_divsteps_hv,rbv()),
                     make_function(&set_divsteps_hv,dcp()),
                     "number of x-ray beam divergence micro-steps across source in horizontal and vertical directions, speed is inversely proportional to this")
      /* beam divergence internal micro-step size in mrad */
      .add_property("divstep_hv_mrad",
                     make_function(&get_divstep_hv_mrad,rbv()),
                     make_function(&set_divstep_hv_mrad,dcp()),
                     "beam divergence internal micro-step size across source in horizontal and vertical directions (mrad)")
      /* Trim automatic range of sources for divergence to lie within an ellipse instead of a rectangle */
      .add_property("round_div",
                     make_getter(&nanoBragg::round_div,rbv()),
                     make_setter(&nanoBragg::round_div,dcp()),
                     "Trim automatic range of sources for divergence to lie within an ellipse instead of a rectangle")

      /* spindle phi rotation start in deg */
      .add_property("phi_deg",
                     make_function(&get_phi_deg,rbv()),
                     make_function(&set_phi_deg,dcp()),
                     "spindle phi rotation start angle (deg)")
      /* spindle phi rotation range in deg */
      .add_property("osc_deg",
                     make_function(&get_osc_deg,rbv()),
                     make_function(&set_osc_deg,dcp()),
                     "spindle phi rotation range (deg)")
      /* number of internal phi micro-steps in phi sweep */
      .add_property("phisteps",
                     make_function(&get_phisteps,rbv()),
                     make_function(&set_phisteps,dcp()),
                     "number of internal phi micro-steps in phi sweep, speed is inversely proportional to this")
      /* size of the internal phi micro-step, in deg */
      .add_property("phistep_deg",
                     make_function(&get_phistep_deg,rbv()),
                     make_function(&set_phistep_deg,dcp()),
                     "size of the internal phi micro-step (deg)")


      /* crystal mosaic spread spherical cap range in deg */
      .add_property("mosaic_spread_deg",
                     make_function(&get_mosaic_spread_deg,rbv()),
                     make_function(&set_mosaic_spread_deg,dcp()),
                     "crystal mosaic spread spherical cap range, not radius (deg)")
      /* number of discrete mosaic domains to generate, speed is inversely proportional to this */
      .add_property("mosaic_domains",
                     make_function(&get_mosaic_domains,rbv()),
                     make_function(&set_mosaic_domains,dcp()),
                     "number of discrete mosaic domains to generate, speed is inversely proportional to this")

      /* hkl and F */
      .add_property("indices",
                     make_function(&get_indices,rbv()),
                     make_function(&set_indices,dcp()),
                     "miller flex array, use miller array Fhkl instead")
      .add_property("amplitudes",
                     make_function(&get_amplitudes,rbv()),
                     make_function(&set_amplitudes,dcp()),
                     "structure factor amplitude flex array, use miller array Fhkl instead")
      /* hkl and F as a 2-tuple of indices and flex-double */
      .add_property("Fhkl_tuple",
                     make_function(&get_Fhkl_tuple,rbv()),
                     make_function(&set_Fhkl_tuple,dcp()),
                     "hkl and F as a 2-tuple of indices and flex-double, use miller array Fhkl instead")
      /* override value of missing Fs, default 0 */
      .add_property("default_F",
                     make_function(&get_default_F,rbv()),
                     make_function(&set_default_F,dcp()),
                     "override value of missing Fs, default 0 (useful if you just want spots fast)")


      /* toggle interpolation between structure factors */
      .add_property("interpolate",
                     make_getter(&nanoBragg::interpolate,rbv()),
                     make_setter(&nanoBragg::interpolate,dcp()),
                     "toggle interpolation between structure factors for inter-Bragg peaks")
      /* True turns spots into flat-tops same width as sinc function */
      .add_property("xtal_shape",
                     make_getter(&nanoBragg::xtal_shape,rbv()),
                     make_setter(&nanoBragg::xtal_shape,dcp()),
                     "select crystal shape transform: square is exact, round is approximate ellipsoid, gauss for gaussian spots, or tophat for top-hat-spots")
      /* experimental: use integral form instead of oversamling */
      .add_property("integral_form",
                     make_getter(&nanoBragg::integral_form,rbv()),
                     make_setter(&nanoBragg::integral_form,dcp()),
                     "experimental: use integral form instead of oversamling")

      /* random number seed for mosaic domain generation */
      .add_property("mosaic_seed",
                     make_getter(&nanoBragg::mosaic_seed,rbv()),
                     make_setter(&nanoBragg::mosaic_seed,dcp()),
                     "random number seed for mosaic domain generation, keep the same for a given crystal")
      /* random number seed for noise generation (use a different value for different images) */
      .add_property("seed",
                     make_getter(&nanoBragg::seed,rbv()),
                     make_setter(&nanoBragg::seed,dcp()),
                     "random number seed for noise generation, use a different value for different images")
      /* random number seed for detector calibration error (usually same for all frames) */
      .add_property("calib_seed",
                     make_getter(&nanoBragg::calib_seed,rbv()),
                     make_setter(&nanoBragg::calib_seed,dcp()),
                     "random number seed for detector calibration error, keep the same for a given detector")

      /* 2D flex array representing pixel values in expected photons, not neccesarily integers */
      .add_property("raw",
                     make_getter(&nanoBragg::raw,rbv()),
                     make_setter(&nanoBragg::raw,dcp()),
                     "2D flex array representing floating-point pixel values in expected photons, not neccesarily integers")

      /* print to screen a summary of all initialized parameters */
      .def("show_params",&nanoBragg::show_params,
       "print out all simulation parameters, just like the standalone program")

      /* actual run of the spot simulation */
      .def("add_nanoBragg_spots",&nanoBragg::add_nanoBragg_spots,
       "actually run the spot simulation, going pixel-by-pixel over the region-of-interst")

      /* actual run of the background simulation */
      .def("add_background",&nanoBragg::add_background,
       "run the non-Bragg simulation, adding background from speficied amorphous materials")

      /* blur the image with specified point-spread function */
      .def("apply_psf",&nanoBragg::apply_psf,
        (arg_("psf_type")=0,arg_("fwhm_pixels")=0,arg_("psf_radius")=0),
       "blur the image with specified point-spread function, may be done before or after adding noise")

      /* actual run of the noise simulation */
      .def("add_noise",&nanoBragg::add_noise,
       "apply specified Poisson, calibration, flicker and read-out noise to the pixels")

      .def("to_smv_format",&nanoBragg::to_smv_format,
        (arg_("fileout"),arg_("intfile_scale")=0,arg_("adc_offset")=40),
        "interally produce an SMV-format image file on disk from the raw pixel array\nintfile_scale is applied before rounding off to integral pixel values")
    ;
    // end of nanoBragg class definition

    // Expose enums to Python
    enum_<pivot>("pivot",
     "description of detector pivot point for detector rotation angles")
      .value("Beam",BEAM)
      .value("Sample",SAMPLE)
      .export_values();
    enum_<shapetype>("shapetype",
       "select shape of the crystal, spot, or point-spread function" )
      .value("SQUARE",SQUARE)
      .value("ROUND",ROUND)
      .value("GAUSS",GAUSS)
      .value("TOPHAT",TOPHAT)
      .value("FIBER",FIBER)
      .value("UNKNOWN",UNKNOWN)
      .export_values();
    enum_<convention>("convention",
       "beam center convention to use when interpreting Xbeam Ybeam")
      .value("CUSTOM",CUSTOM)
      .value("ADXV",ADXV)
      .value("MOSFLM",MOSFLM)
      .value("XDS",XDS)
      .value("DIALS",DIALS)
      .value("DENZO",DENZO)
      .export_values();

  } // end of init_module

} // end of namespace
} // end of namespace boost_python
} // end of namespace nanoBragg
} // end of namespace simtbx


BOOST_PYTHON_MODULE(simtbx_nanoBragg_ext)
{
  simtbx::nanoBragg::boost_python::init_module();
}


