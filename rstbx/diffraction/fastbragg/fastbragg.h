#ifndef RSTBX_DIFF_FASTBRAGG_H
#define RSTBX_DIFF_FASTBRAGG_H
//34567890        20        30        40        50        60        70        80        90
#include <map>
#include <string>
#include <ctime>
#include <limits>
#include <scitbx/vec3.h>
#include <scitbx/mat3.h>
#include <scitbx/array_family/shared.h>
#include <scitbx/array_family/flex_types.h>
#include <scitbx/array_family/accessors/c_grid.h>
#include <cctbx/crystal_orientation.h>
#include <cctbx/miller.h>

namespace rstbx {
namespace diffraction {
namespace fastbragg {

namespace af = scitbx::af;
typedef scitbx::vec3<double> vec3;
typedef scitbx::mat3<double> mat3;
typedef cctbx::miller::index<> miller_t;
typedef af::shared<miller_t > indices;

/* Fourier transform of a sphere */
double sinc3(double const& x);

//Plan to replace ran1, gammln & poidev with Boost functions in future

/* returns a uniform random deviate between 0 and 1 */
/* Random number generator ran1 from Computers in Physics */
/* Volume 6 No. 5, 1992, 522-524, Press and Teukolsky */
/* To generate real random numbers 0.0-1.0 */
/* Should be seeded with a negative integer */
float ran1(long *idum);

/* ln of the gamma function */
float gammln(float xx);

float poidev(float xm, long *idum);

//! Class detector defines the image dimensions and pixel size
/*! All attributes are set through initializer functions:
      The constructor defines slow x fast pixels and pixel size
      set_region_of_interest() defines a subarea, if desired
      set_oversampling() allows > 1 samplings per pixel, if desires
*/
struct detector {
  int ypixels,zpixels;                     // number of pixels in slow & fast dimensions
  int roi_zmin,roi_zmax,roi_ymin,roi_ymax; // region of interest definitions
  double pixel_sz;                         // pixel size in meters
  int oversample;                          // number of sub-pixels per pixel
  double subpixel_sz;
  double max_I;
  af::versa<double, af::c_grid<2> > raw;

  inline
  detector(int const& slowpixels, int const& fastpixels, double const& pixel_size):
    ypixels(slowpixels),zpixels(fastpixels),
    roi_zmin(0),roi_zmax(fastpixels),
    roi_ymin(0),roi_ymax(slowpixels),
    pixel_sz(pixel_size),oversample(1),
    subpixel_sz(pixel_size),
    raw(af::c_grid<2> (slowpixels,fastpixels))
    {  }

  inline
  void set_region_of_interest(int const& slow_min, int const& slow_max,
                              int const& fast_min, int const& fast_max) {
    roi_ymin = slow_min;
    roi_ymax = slow_max;
    roi_zmin = fast_min;
    roi_zmax = fast_max;
  }

  void set_oversampling(int const& oversampling){
    oversample = oversampling;
    subpixel_sz = pixel_sz/oversample;
  }

  inline
  af::versa<int, af::c_grid<2> >
  intdata(){
    af::versa<int, af::c_grid<2> > result(af::c_grid<2> (ypixels,zpixels),
                                             af::init_functor_null<int>());
    double * dptr = raw.begin();
    int * iptr = result.begin();
    for (int count = 0; count < raw.size(); ++count){
      *iptr++ = (int)std::ceil( (*dptr++)-0.5 );
    }
    return result;
  }
};

//! Class camera defines the experimental geometry
/*! All attributes are addressed directly
      X : direction from crystal to detector
      Y : slow pixel direction
      Z : fast pixel direction
*/
struct camera {
  double distance;    //distance from origin to detector center in meters
  double Ybeam,Zbeam; //slow & fast coordinates of direct beam in meters
  double lambda0;     //central wavelength in meters
  double dispersion;  //fractional full-width spectral dispersion (delta lambda/lambda)
  int dispsteps;      //number of wavelengths in above range
  double hdivrange;   //horizontal angular spread of source points in radians
  double vdivrange;   //vertical angular spread of source points in radians
  double hdivstep;    //number of source points in the horizontal
  double vdivstep;    //number of source points in the vertical
  double source_distance; //distance from origin to source in meters
  double fluence;     //incident x-ray fluence in photons/meter^2

  inline
  af::shared<double>
  get_wavelengths() const{
    af::shared<double> result;
    /* count dispersion steps with sweep over spectral dispersion */
    double dispstep,effective_dispersion=dispersion;
    if(dispsteps > 1){
        dispstep = lambda0*dispersion/(dispsteps-1)-1e-23;
    } else {
        dispstep = 0;
    }
    if(dispstep == 0.0) { dispstep = 1e99; effective_dispersion = 0.0; };

    for(double lambda = lambda0*(1.0 - effective_dispersion/2);
               lambda <= lambda0*(1.0 + effective_dispersion/2);
               lambda += dispstep){
        result.push_back(lambda);
    }
    SCITBX_ASSERT(result.size()==dispsteps);
    return result;
  }

  inline
  int
  get_divsteps() const {

    /* count divsteps sweep over solid angle of beam divergence */
    int divsteps = 0;
    for(double hdiv=-hdivrange/2;hdiv<=hdivrange/2+1e-11;hdiv+=hdivstep){
    for(double vdiv=-vdivrange/2;vdiv<=vdivrange/2+1e-11;vdiv+=vdivstep){
    //   printf("%d divergence steps %g %g\n",divsteps,hdiv,vdiv);
      /* force an elliptical divergence */
      if( hdivrange != 0 && vdivrange != 0) {
        SCITBX_ASSERT( hdivrange != 0); // assertion works around an
        SCITBX_ASSERT( vdivrange != 0); // apparent optimizer bug
        if( (hdiv*hdiv/hdivrange/hdivrange +
             vdiv*vdiv/vdivrange/vdivrange)*4.0 > 1.0 ) continue;
      }
      ++divsteps;
    }
    }
    return divsteps;
  }

};

struct crystal {
  cctbx::crystal_orientation orientation;
  indices miller;
  af::shared<double> amplitudes;
  int Na,Nb,Nc;

  inline
  std::map<miller_t,double>
  get_amplitude_mapping() const {
    std::map<miller_t,double> result;
    for (int i = 0; i < miller.size(); ++i){
      result[miller[i]]=amplitudes[i];
    }
    return result;
  }
};

struct integer_miller_index_policy { // with IndexType == cctbx::miller::index<>

  template <typename IndexType>
  static void push_back_index(af::shared<IndexType>& container, vec3 const& hklvec){

    /* round off to nearest whole index */
    int h0 = static_cast<int>(std::ceil(hklvec[0]-0.5));
    int k0 = static_cast<int>(std::ceil(hklvec[1]-0.5));
    int l0 = static_cast<int>(std::ceil(hklvec[2]-0.5));
    container.push_back(IndexType(h0,k0,l0));
  }
};

struct vec3_double_miller_index_policy { // with IndexType == scitbx::vec3<double>

  template <typename IndexType>
  static void push_back_index(af::shared<IndexType>& container, vec3 const& hklvec){
    container.push_back(hklvec);
  }
};

//! Simulation of nanocrystallography "still".  Contributed by James Holton, LBNL.
/*! Algorithm is described in Kirian, RA, Wang, X, Weierstall, U, Schmidt, KE,
      Spence, JCH, Hunter, M, Fromme, P, White, T, Chapman, HN & Holton, J (2010).
      "Femtosecond protein nanocrystallography: data analysis methods",
      Optics Express 18, 5713-5723.
    This is a perfect-lattice diffraction simulator (no mosaicity).
    Function sweep_over_detector() essentially implements eq. (1) from that reference.
    Still needed:  Fix the aliasing errors introduced by sampling the shape transform
    at the center of each pixel, instead of integrating over the body of the pixel.
    This can rapidly lead to "missing" spots as the crystal becomes large (shape transform
    becomes sharp compared to inter-pixel distance).  nearBragg has the same problem.
    With a symmetric Gaussian spot shape, the over-pixel integral should be do-able
    analytically, but has not been implemented yet.
*/
class fast_bragg_simulation {
 private:
   detector D;
   camera C;
   crystal X;

 public:
  bool printout;
  fast_bragg_simulation():
    D(detector(0,0,0.)),printout(false){}

  inline
  void set_detector(detector const& indata){
    D = indata;
  }

  inline
  void set_camera(camera const& incamera){
    C = incamera;
  }

  inline
  void set_crystal(crystal const& incrystal){
    X = incrystal;
  }

  inline
  void sweep_over_detector(bool const& verbose){
    D.max_I = 0.0;
    int j_image_ptr = 0;
    int progress_pixel = 0;
    af::shared<double> lambdas=C.get_wavelengths();
    int divsteps = C.get_divsteps();
    int steps = divsteps*C.dispsteps*D.oversample*D.oversample;
    mat3 Amat = X.orientation.direct_matrix();
    std::map<miller_t,double> Fhkl = X.get_amplitude_mapping();
    double PI = scitbx::constants::pi;
    double* floatimage(D.raw.begin());
    int progress_pixels = (D.roi_zmax-D.roi_zmin+1)*(D.roi_ymax-D.roi_ymin+1);

    /* Thomson cross section */
    double r_e_sqr = 7.94079248018965e-30;

    /* water background stuff */
    double water_size = 0.0;
    double water_scale = 1.75;
    // 2.57^2*r_e_sqr*1e6*6.022e23/18.0
    for(int ypixel=0;ypixel<D.ypixels;++ypixel){
     for(int zpixel=0;zpixel<D.zpixels;++zpixel){

      if(zpixel < D.roi_zmin || zpixel > D.roi_zmax ||
         ypixel < D.roi_ymin || ypixel > D.roi_ymax) {
           ++j_image_ptr; continue;
      }

      /* reset photon count */
      double I_photon_count = 0;
      double omega_pixel = 0;
      double polar = 0;

      for(int suby=0;suby<D.oversample;++suby){
       for(int subz=0;subz<D.oversample;++subz){

        double Zdet = D.subpixel_sz*(zpixel*D.oversample + subz);
        double Ydet = D.subpixel_sz*(ypixel*D.oversample + suby);
//      Zdet = pixel*zpixel;
//      Ydet = pixel*ypixel;

        /* construct detector pixel position */
        vec3 pixel_xyz( C.distance, Ydet-C.Ybeam, Zdet-C.Zbeam );

        /* construct the unit vector to this pixel */
        double air_path = pixel_xyz.length();
             vec3 S_xyz = pixel_xyz/air_path;

        if (omega_pixel==0.0){
          /* solid angle subtended by a pixel: (pix/air_path)^2*cos(2theta) */
          omega_pixel = D.pixel_sz*D.pixel_sz*C.distance/(air_path*air_path*air_path);


          /* polarization factor for this pixel */
          double costwotheta =
              std::sqrt(pixel_xyz[1]*pixel_xyz[1]+pixel_xyz[2]*pixel_xyz[2])/air_path;
          polar = 0.5*(1.0+costwotheta*costwotheta);
        }

        /* sweep over wavelengths */
        for(int i_lambda=0; i_lambda < lambdas.size(); ++i_lambda){

            /* sweep over solid angle of beam divergence */
            for(double hdiv=-C.hdivrange/2;hdiv<=C.hdivrange/2+1e-11;hdiv+=C.hdivstep){
            for(double vdiv=-C.vdivrange/2;vdiv<=C.vdivrange/2+1e-11;vdiv+=C.vdivstep){

                /* force an elliptical divergence */
                if( C.hdivrange != 0 && C.vdivrange != 0) {
                  if((hdiv*hdiv/C.hdivrange/C.hdivrange +
                      vdiv*vdiv/C.vdivrange/C.vdivrange)*4 > 1.0)
                    {continue;}
                }

                /* construct source position (flat, coherently-emitting plane) */
                vec3 source_xyz( -C.source_distance,
                                 std::atan(hdiv)*C.source_distance,
                                 std::atan(vdiv)*C.source_distance);

                /* construct the incident beam unit vector */
                double source_path = source_xyz.length();
                vec3 S0_xyz = -source_xyz/source_path;

                /* construct the scattering vector for this pixel */
                vec3 s_xyz = (S_xyz-S0_xyz)/lambdas[i_lambda];

                /* construct fractional Miller indices */
                vec3 hklvec = (1.e-10 * Amat) * s_xyz; //Convert Amat to meters first

                /* round off to nearest whole index */
                int h0 = static_cast<int>(std::ceil(hklvec[0]-0.5));
                int k0 = static_cast<int>(std::ceil(hklvec[1]-0.5));
                int l0 = static_cast<int>(std::ceil(hklvec[2]-0.5));

                /* structure factor of the unit cell */
                miller_t hkl (h0, k0, l0);
                double F_cell = Fhkl[hkl];

                /* structure factor of the lattice (paralelpiped crystal)
              F_latt=sin(PI*Na*h)*sin(PI*Nb*k)*sin(PI*Nc*l)/sin(PI*h)/sin(PI*k)/sin(PI*l);
                */
                double F_latt = 1.0;
                double denom;
                denom  = std::sin(PI*hklvec[0]);
                if(denom==0.0){
                    F_latt *= X.Na;
                }
                else
                {
                    F_latt *= std::sin(PI*X.Na*hklvec[0])/denom;
                }
                denom  = std::sin(PI*hklvec[1]);
                if(denom==0.0){
                    F_latt *= X.Nb;
                }
                else
                {
                    F_latt *= std::sin(PI*X.Nb*hklvec[1])/denom;
                }
                denom  = std::sin(PI*hklvec[2]);
                if(denom==0.0){
                    F_latt *= X.Nc;
                }
                else
                {
                    F_latt *= std::sin(PI*X.Nc*hklvec[2])/denom;
                }
               //F_latt=Na*Nb*Nc*sinc(PI*Na*(h-h0))*sinc(PI*Nb*(k-k0))*sinc(PI*Nc*(l-l0));
                F_latt = X.Na * X.Nb * X.Nc * sinc3(PI*X.Na*(hklvec[0]-h0)) *
                         sinc3(PI*X.Nb*(hklvec[1]-k0)) * sinc3(PI*X.Nc*(hklvec[2]-l0));

                /* convert amplitudes into intensity (photons per steradian) */
                I_photon_count += F_cell*F_cell*F_latt*F_latt;
            }
            }
          }
         }
        }

        /* add background from something? */
        double I_water =
          water_scale*C.fluence*polar*water_size*water_size*water_size*omega_pixel;

        floatimage[j_image_ptr]=
              r_e_sqr*C.fluence*polar*I_photon_count/steps*omega_pixel + I_water;
        if(floatimage[j_image_ptr] > D.max_I) D.max_I = floatimage[j_image_ptr];
        if( printout ) {
            printf("%4d %4d   %15.10g\n", zpixel,ypixel,floatimage[j_image_ptr]);
        }else if (verbose){
          if(progress_pixel % ( progress_pixels/20 ) == 0 ||
            ((10*progress_pixel<progress_pixels || 10*progress_pixel>9*progress_pixels) &&
            (progress_pixel % (progress_pixels/100) == 0))) {
              printf("%d%% done\n",progress_pixel*100/progress_pixels);
          }
//        printf("%d %d\n",j_image_ptr,j_image_ptr % ( pixels/100 ));
//     if(j_image_ptr % ( pixels/100 ) == 0) printf("%d%% done\n",j_image_ptr*100/pixels);
        }
//      printf("%d ",intimage[j_image_ptr]);
//      printf("%.10g ",floatimage[j_image_ptr]);
//      printf("%.10g ",max_I);
//      printf("\n");
//;fflush(stdout);
        ++j_image_ptr;
        ++progress_pixel;
      }
    }

  }

  template <typename IndexType, typename MillerConversionPolicy>
  af::shared<IndexType>
  sweep_over_detector_get_indices(bool const& verbose){
    D.max_I = 0.0;
    int j_image_ptr = 0;
    int progress_pixel = 0;
    af::shared<double> lambdas=C.get_wavelengths();
    mat3 Amat = X.orientation.direct_matrix();
    af::shared<IndexType> miller_set;

    for(int ypixel=0;ypixel<D.ypixels;++ypixel){
     for(int zpixel=0;zpixel<D.zpixels;++zpixel){

      if(zpixel < D.roi_zmin || zpixel > D.roi_zmax ||
         ypixel < D.roi_ymin || ypixel > D.roi_ymax) {
           ++j_image_ptr; continue;
      }

      /* reset photon count */
      double omega_pixel = 0;
      SCITBX_ASSERT(D.oversample==1);
      for(int suby=0;suby<D.oversample;++suby){
       for(int subz=0;subz<D.oversample;++subz){

        double Zdet = D.subpixel_sz*(zpixel*D.oversample + subz);
        double Ydet = D.subpixel_sz*(ypixel*D.oversample + suby);
//      Zdet = pixel*zpixel;
//      Ydet = pixel*ypixel;

        /* construct detector pixel position */
        vec3 pixel_xyz( C.distance, Ydet-C.Ybeam, Zdet-C.Zbeam );

        /* construct the unit vector to this pixel */
        double air_path = pixel_xyz.length();
             vec3 S_xyz = pixel_xyz/air_path;

        if (omega_pixel==0.0){
          /* solid angle subtended by a pixel: (pix/air_path)^2*cos(2theta) */
          omega_pixel = D.pixel_sz*D.pixel_sz*C.distance/(air_path*air_path*air_path);
        }

        /* sweep over wavelengths */
        SCITBX_ASSERT(lambdas.size()==1);
        for(int i_lambda=0; i_lambda < lambdas.size(); ++i_lambda){

            /* sweep over solid angle of beam divergence */
            SCITBX_ASSERT(C.hdivrange==0);
            SCITBX_ASSERT(C.vdivrange==0);
            for(double hdiv=-C.hdivrange/2;hdiv<=C.hdivrange/2+1e-11;hdiv+=C.hdivstep){
            for(double vdiv=-C.vdivrange/2;vdiv<=C.vdivrange/2+1e-11;vdiv+=C.vdivstep){

                /* force an elliptical divergence */
                if( C.hdivrange != 0 && C.vdivrange != 0) {
                  if((hdiv*hdiv/C.hdivrange/C.hdivrange +
                      vdiv*vdiv/C.vdivrange/C.vdivrange)*4 > 1.0)
                    {continue;}
                }

                /* construct source position (flat, coherently-emitting plane) */
                vec3 source_xyz( -C.source_distance,
                                 std::atan(hdiv)*C.source_distance,
                                 std::atan(vdiv)*C.source_distance);

                /* construct the incident beam unit vector */
                double source_path = source_xyz.length();
                vec3 S0_xyz = -source_xyz/source_path;

                /* construct the scattering vector for this pixel */
                vec3 s_xyz = (S_xyz-S0_xyz)/lambdas[i_lambda];

                /* construct fractional Miller indices */
                vec3 hklvec = (1.e-10 * Amat) * s_xyz; //Convert Amat to meters first

                MillerConversionPolicy::push_back_index(miller_set,hklvec);
            }
            }
          }
         }
        }
        ++j_image_ptr;
        ++progress_pixel;
      }
    }
    return miller_set;
  }

  inline
  void to_smv_format(std::string const& fileout, double const& intfile_scale_r,
                     int const& saturation, bool const& verbose){

    int pixels = D.ypixels * D.zpixels;
    double * floatimage = D.raw.begin();
    FILE* outfile;

    /* output as doubles */
    /*
    std::string floatfilename("floatimage.bin");
    printf("writing %s as %d %d-byte floats\n",
                   floatfilename.c_str(),pixels,(int)sizeof(double));
    outfile = fopen(floatfilename.c_str(),"w");
    fwrite(floatimage,sizeof(double),pixels,outfile);
    fclose(outfile);
    */

    /* output as ints */
    long seed = -std::time((time_t *)0);
    double OFFSET = 0.0;

    af::versa<unsigned short int, af::c_grid<2> > intimage_v(
      af::c_grid<2> (D.ypixels,D.zpixels));
    unsigned short int * intimage = intimage_v.begin();

    double intfile_scale = intfile_scale_r;
    int j = 0;
    if(intfile_scale_r <= 0.0){
        printf("providing default scaling: max_I = %g\n",D.max_I);
        intfile_scale = 55000.0/D.max_I;
        printf("providing default scaling: intfile_scale = %f\n",intfile_scale);
    }

    double max_value = (double)std::numeric_limits<unsigned short int>::max();
    for(int ypixel=0;ypixel<D.ypixels;++ypixel){
      for(int zpixel=0;zpixel<D.zpixels;++zpixel){
        if(zpixel < D.roi_zmin || zpixel > D.roi_zmax ||
           ypixel < D.roi_ymin || ypixel > D.roi_ymax) {
           ++j; continue;
        }
        intimage[j] = (unsigned short int) (std::min(max_value,
                      floatimage[j] *intfile_scale+0.0*ran1(&seed)+OFFSET) );
        //printf("%d %d = %d\n",zpixel,ypixel,intimage[j]);
        ++j;
      }
    }
    if (verbose){
      printf("writing %s as %d-byte integers\n",fileout.c_str(),
            (int)sizeof(unsigned short int));}
    outfile = fopen(fileout.c_str(),"w");
    fprintf(outfile,"{\nHEADER_BYTES=512;\nDIM=2;\nBYTE_ORDER=little_endian;");
    fprintf(outfile,"\nTYPE=unsigned_short;\nSIZE1=%d;\nSIZE2=%d;",D.zpixels,D.ypixels);
    fprintf(outfile,"\nPIXEL_SIZE=%g;\nDISTANCE=%g;\n",D.pixel_sz*1000.,C.distance*1000.);
    fprintf(outfile,"WAVELENGTH=%g;\nBEAM_CENTER_X=%g;\nBEAM_CENTER_Y=%g;\n",
                     C.lambda0*1e10,C.Zbeam*1000.0,(D.ypixels*D.pixel_sz-C.Ybeam)*1000);
    fprintf(outfile,"PHI=0;\nOSC_START=0;\nOSC_RANGE=0;\n");
    fprintf(outfile,"DETECTOR_SN=777;\nCCD_IMAGE_SATURATION=%d;\n",saturation);
    fprintf(outfile,"}\f");
    while ( ftell(outfile) < 512 ){ fprintf(outfile," "); };
    fwrite(intimage,sizeof(unsigned short int),pixels,outfile);
    fclose(outfile);

    /* simulate Poisson noise */
    /*
    j = 0;
    for(int ypixel=0;ypixel<D.ypixels;++ypixel){
      for(int zpixel=0;zpixel<D.zpixels;++zpixel){
        if(zpixel < D.roi_zmin || zpixel > D.roi_zmax ||
           ypixel < D.roi_ymin || ypixel > D.roi_ymax) {
           ++j; continue;
        }
        intimage[j] = (unsigned short int) ( poidev( floatimage[j], &seed ) +OFFSET );
        ++j;
      }
    }
    std::string noisefilename("noiseimage.img");

    printf("writing %s as %d-byte integers\n",noisefilename.c_str(),
            (int)sizeof(unsigned short int));
    outfile = fopen(noisefilename.c_str(),"w");
    fprintf(outfile,"{\nHEADER_BYTES=512;\nDIM=2;\nBYTE_ORDER=little_endian;");
    fprintf(outfile,"\nTYPE=unsigned_short;\nSIZE1=%d;\nSIZE2=%d;",D.zpixels,D.ypixels);
    fprintf(outfile,"\nPIXEL_SIZE=%g;\nDISTANCE=%g;\n",D.pixel_sz*1000.,C.distance*1000.);
    fprintf(outfile,"WAVELENGTH=%g;\nBEAM_CENTER_X=%g;\nBEAM_CENTER_Y=%g;\n",
                     C.lambda0*1e10,C.Zbeam*1000.0,(D.ypixels*D.pixel_sz-C.Ybeam)*1000);
    fprintf(outfile,"PHI=0;\nOSC_START=0;\nOSC_RANGE=0;\n");
    fprintf(outfile,"DETECTOR_SN=777;\nCCD_IMAGE_SATURATION=%d;\n",saturation);
    fprintf(outfile,"}\f");
    while ( ftell(outfile) < 512 ){ fprintf(outfile," "); };
    fwrite(intimage,sizeof(unsigned short int),pixels,outfile);
    fclose(outfile);
    */

  }

};


}}}// namespace rstbx::diffraction::fastbragg
#endif //RSTBX_DIFF_FASTBRAGG_H
