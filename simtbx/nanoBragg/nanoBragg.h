#ifndef SIMTBX_NANOBRAGG_H
#define SIMTBX_NANOBRAGG_H
//34567890        20        30        40        50        60        70        80        90
#include <map>
#include <string>
#include <ctime>
#include <limits>
#include <scitbx/vec3.h>
#include <scitbx/mat3.h>
#include <scitbx/constants.h>
#include <scitbx/array_family/shared.h>
#include <scitbx/array_family/flex_types.h>
#include <scitbx/array_family/accessors/c_grid.h>
#include <scitbx/array_family/versa.h>
#include <cctbx/uctbx.h>
#include <cctbx/crystal_orientation.h>
#include <cctbx/miller.h>
#include <dxtbx/model/detector.h>
#include <dxtbx/model/beam.h>
#include <boost/math/special_functions/erf.hpp>
#include <boost/math/special_functions/fpclassify.hpp>
#include <boost_adaptbx/python_streambuf.h>
#include <omptbx/omp_or_stubs.h>

using boost::math::erf;
using boost::math::isnan;
#define isnan(X) boost::math::isnan(X)


/* need this on macs */
#define _USE_MATH_DEFINES
/* need NAN and isnan() to define uninitialized values */
#include <cmath>
#include <cfloat>
/* seem to need these on windoze and MacOS */
#ifndef NAN
//#define NAN strtod("NAN",NULL)
/* this works on Windoze */
#define NAN sqrt((long double)-1)
#endif
#ifndef DBL_MIN
//#define DBL_MIN 1e-99
#endif
#ifndef M_PI
#define M_PI scitbx::constants::pi
#endif
/* need this for read_text_file() */
#include <stdarg.h>

namespace simtbx {
namespace nanoBragg {

namespace af = scitbx::af;
typedef scitbx::vec2<double> vec2;
typedef scitbx::vec3<double> vec3;
typedef scitbx::mat3<double> mat3;
typedef cctbx::miller::index<> miller_t;
typedef af::shared<miller_t > indices;


/* fundamental constants of physics that should not be changeable */
/* Avogadro's number */
static const double Avogadro = 6.02214179e23;
/* convert from radians to degrees: 180/pi */
static const double RTD = 180.0/M_PI;
/* Thomson cross section ((e^2)/(4*PI*epsilon0*m*c^2))^2 */
static const double r_e_sqr = 7.94079248018965e-30;


/* Holton's general utility to read in text file into double arrays at provided addresses */
size_t read_text_file(char *filename, size_t nargs, ... );

/* cubic spline interpolation functions */
void polint(double *xa, double *ya, double x, double *y);
void polin2(double *x1a, double *x2a, double **ya, double x1,double x2, double *y);
void polin3(double *x1a, double *x2a, double *x3a, double ***ya, double x1,double x2, double x3, double *y);


/* James Holton's personal linear algebra functions */
/* rotate a 3-vector in space applied in order phix,phiy,phiz*/
double *rotate(double *v, double *newv, double phix, double phiy, double phiz);
/* rotate a 3-vector about a unit vector axis */
double *rotate_axis(double *v, double *newv, double *axis, double phi);
/* rotate a 3-vector using a 9-element unitary matrix */
double *rotate_umat(double *v, double *newv, double *umat);

/* vector cross product where vector magnitude is 0th element */
double *cross_product(double *x, double *y, double *z);
/* vector inner product where vector magnitude is 0th element */
double dot_product(double *x, double *y);
/* compute difference between two vectors */
double vector_diff(double *vector, double *origin_vector, double *new_vector);
/* measure magnitude of vector and put it in 0th element */
double magnitude(double *vector);
/* scale the magnitude of a vector */
double vector_scale(double *vector, double *new_vector, double scale);
/* force the magnitude of vector to given value */
double vector_rescale(double *vector, double *new_vector, double magnitude);
/* make a unit vector pointing in same direction and report magnitude (both args can be same vector) */
double unitize(double *vector, double *new_unit_vector);


/* polarization factor from vectors */
double polarization_factor(double kahn_factor, double *incident, double *diffracted, double *axis);


/* generate random unitary rotation matrix within a spherical cap */
double *mosaic_rotation_umat(double mosaicity, double umat[9], long *idum);
/* convert unitary matrix into missetting angles */
double *umat2misset(double umat[9],double *missets);


/* returns a uniform random deviate between 0 and 1 */
/* Random number generator ran1 from Computers in Physics */
/* Volume 6 No. 5, 1992, 522-524, Press and Teukolsky */
/* To generate real random numbers 0.0-1.0 */
/* Should be seeded with a negative integer */
double ran1(long *idum);

/* ln of the gamma function */
double gammln(double xx);

/* return Gaussian deviate with rms=1 and FWHM = 2/sqrt(log(2)) */
double gaussdev(long *idum);
/* return Poissonian deviate given expectation value */
double poidev(double xm, long *idum);


/* Fourier transform of a truncated lattice */
double sincg(double x, double N);
/* Fourier transform of a sphere */
double sinc3(double x);
/* Fourier transform of a spherically-truncated lattice */
double sinc_conv_sinc3(double x);

/* typedefs to help remember options */
typedef enum { SAMPLE, BEAM } pivot;
typedef enum { UNKNOWN, SQUARE, ROUND, GAUSS, TOPHAT, FIBER } shapetype;
typedef enum { CUSTOM, ADXV, MOSFLM, XDS, DIALS, DENZO } convention;

/* math functions for point spread */
/* 2D Gaussian integral=1 */
double ngauss2D(double x, double y, double fwhm);
/* integral of Gaussian fwhm=1 integral=1 */
double ngauss2D_integ(double x, double y);
/* unit volume integrated over a pixel, fwhm = 1 */
double ngauss2D_pixel(double x,double y,double pix);
double integrate_gauss_over_pixel(double x, double y, double fwhm, double pix);

/* fiber-coupled CCD PSF: g/(2*pi)*(g**2+x**2+y**2)**(-3/2) see Holton JSR (2012) */
double fiber2D(double x,double y,double g);
double fiber2D_integ(double x,double y,double g);
double fiber2D_pixel(double x,double y,double g,double pix);
double integrate_fiber_over_pixel(double x, double y, double g, double pix);

/* median filter tools */
double fmedian(unsigned int n, double arr[]);
double fmedian_with_rejection(unsigned int n, double arr[],double sigma_cutoff, double *final_mad, int *final_n);
double fmedian_absolute_deviation(unsigned int n, double arr[], double median_value);
double fmean_with_rejection(unsigned int starting_points, double arr[], double sigma_cutoff, double *final_rmsd, int *final_n);


//! Simulation of nanocrystal diffraction.  Contributed by James Holton, LBNL.
/*! */
class nanoBragg {
  public:
    /* optional file stuff, to be removed eventually? */
    char *matfilename;  // = NULL;
    char *hklfilename;  // = NULL;
//    char *dumpfilename;  // = "Fdump.bin\0";
    char *stolfilename;  // = NULL;
    char *imginfilename;  // = NULL;
    char *maskfilename;  // = NULL;
    char *stoloutfilename;  // = "output.stol\0";
    char *sourcefilename;  // = NULL;
    char *floatfilename;  // = "floatimage.bin\0";
    //char *sinfilename = "sinimage.bin\0";
    //char *cosfilename = "cosimage.bin\0";
    char *intfilename;  // = "intimage.img\0";
    char *pgmfilename;  // = "image.pgm\0";
    char *noisefilename;  // = "noiseimage.img\0";
    FILE *infile;  // = NULL;
//    FILE *Fdumpfile;  // = NULL;
//    FILE *outfile;  // = NULL;
//    FILE *stoloutfile;  // = NULL;

    /* progress meter stuff */
    long progress_pixel,progress_pixels;
    bool progress_meter;
    bool babble;
    bool printout;
    int printout_spixel,printout_fpixel;
    int verbose;

    /* x-ray beam properties */
    double beam_vector[4];
    bool coherent;
    bool far_source;
    bool round_div;
    double lambda,*lambda_of;
    double mosaic_spread,*mosaic_umats,mosaic_missets[4];
    double dispersion,dispstep,lambda0;
    double hdiv,hdivstep,hdivrange;
    double vdiv,vdivstep,vdivrange;
    double source_path,source_distance;
    int divsteps,hdivsteps,vdivsteps,dispsteps;
    int hdiv_tic,vdiv_tic,disp_tic,mos_tic;
    int mosaic_domains;
    double weight;
    int source,sources,allocated_sources;
    double *source_X,*source_Y,*source_Z,*source_I,*source_lambda;

    /* version of source list to pass back to Python */
    af::shared<vec3> pythony_source_XYZ;
    af::shared<double> pythony_source_intensity;
    af::shared<double> pythony_source_lambda;
    scitbx::af::versa<dxtbx::model::Beam, scitbx::af::flex_grid<> > pythony_beams;

    /* incident x-ray fluence in photons/m^2   default equivalent to unity
        that is, one electron will scatter 1 ph/SR after a fluence of 1.26e29 ph/m^2
        this places the input file on a photons/pixel scale */
    double fluence;
    double flux,exposure,beamsize;

    /* sample size stuff */
    int    N;
    double Na,Nb,Nc;
    double xtalsize_max,xtalsize_a,xtalsize_b,xtalsize_c;
    double reciprocal_pixel_size;

    /* macroscopic xtal properties, cut with beam? */
    shapetype xtal_shape;
    double hrad_sqr,fudge;
    double xtal_size_x;            /* m */
    double xtal_size_y;            /* m */
    double xtal_size_z;            /* m */
    double xtal_density;             /* g/m^3 */
    double xtal_molecular_weight;    /* g/mol */
    double xtal_volume,xtal_molecules;

    /* amorphous material properties */
    double amorphous_sample_x, amorphous_sample_y, amorphous_sample_z;
    double amorphous_volume;
    double amorphous_molecular_weight;
    double amorphous_density;
    double amorphous_molecules;
    /* scale factor = Fbg^2*r_e_sqr*fluence*Avogadro*volume*density/molecular_weight
                             m^2     ph/m^2  /mol      m^3   g/m^3    g/mol   */
    /* water Fbg = 2.57 in forward direction */


    /* detector stuff */
    double pixel_size; // = 0.1e-3;
    double pixel_pos[4];
    int fpixel,spixel,fpixels,spixels,pixels;
    int allocated_pixels,allocated_stols; // used to decide if we need to erase and re-calloc
    double distance; // = 100.0e-3;
    double detsize_f; // = 102.4e-3;
    double detsize_s; // = 102.4e-3;
    double detector_attnlen; //= 234 um;
    double detector_thick; // =0.0;
    double detector_thickstep,parallax,capture_fraction;
    int    thick_tic,detector_thicksteps; // =-1;
    double fdet_vector[4]; //  = {0,0,0,1};
    double sdet_vector[4]; //  = {0,0,-1,0};
    double odet_vector[4]; //  = {0,1,0,0};
    double pix0_vector[4]; //  = {0,0,0,0};
    double detector_rotx,detector_roty,detector_rotz;
    double twotheta_axis[4]; // = {0,0,1,0};
    pivot detector_pivot; // = BEAM;
    convention beam_convention; // = MOSFLM;
    double detector_twotheta; // = 0.0;
    double airpath,omega_pixel,omega_Rsqr_pixel,omega_sum;
    bool curved_detector; // = 0;
    bool point_pixel; // = 0;
    double Xbeam,Ybeam; //=NAN;
    double Fbeam,Sbeam; //=NAN;
    double Fdet,Sdet,Odet;
    double Fdet0,Sdet0;
    double Xclose,Yclose,close_distance; //=NAN;
    double Fclose,Sclose; //=NAN;
    double ORGX,ORGY; //=NAN;
    double dials_origin[4];
    double detector_is_righthanded; //true;
    double adc_offset; // = 40.0;

    /* use these to remember "user" inputs */
    bool user_beam; //=false;
    bool user_distance; //=false;
    bool user_mosdomains; //=false;

    /* scattering vectors */
    double incident[4];
    double diffracted[4],diffracted0[4];
    double scattering[4];
    double stol,twotheta,theta;

    /* diffraction geometry stuff */
    double costwotheta,sintwotheta,psi;
    double xd,yd,zd,xd0,yd0,zd0;
    double Ewald[4],Ewald0[4],relp[4];
    double dmin; //=0;
    bool integral_form; // = 0;  experimental: use integral form to avoid need for oversampling

    /* polarization stuff */
    double polar_vector[4]; // = {0,0,0,1};
    double vert_vector[4];
    double polar; //=1.0;  // instantaneous polarization factor
    double polarization; //=0.0;  Kahn "polarization" parameter [-1:1]
    bool nopolar; // = 0;  // turn on/off polarization effect

    /* sampling */
    int steps;
    int roi_xmin,roi_xmax,roi_ymin,roi_ymax; // =-1;
    int oversample,recommended_oversample,subS,subF;
    bool user_oversample; // =False
    double subpixel_size;

    /* spindle */
    double phi,phi0,phistep,osc; // =-1.0;
    int phi_tic,phisteps; // =-1;
    double spindle_vector[4]; // = {0,0,0,1};

    /* structure factor representation */
    double phase,Fa,Fb;
    double F,Fbg,Ibg,*stol_of,*Fbg_of;  // = NULL
    double ***Fhkl;  // = NULL
    int    hkls;
    double F_latt,F_cell;
    double F000;        // to mark beam center
    double default_F;   // for spots, usually 0
    double default_Fbg; // for background, usually 0
    double Fbg_highangle,Fbg_lowangle;
    int stols,nearest; // =0;
    double stol_file_mult; // =1.0e10;  convert to meters, usually from Angstrom
    double denom;

    /* background extraction parameters */
//    double *imginfileimage;
    double *diffimage;
    double *stolimage;
    double *Fimage,pixel_F;
    int ignore_values; // =0;
    unsigned short int ignore_value[70000];
    bool *invalid_pixel;
    int valid_pixels;
    bool Fmap_pixel; // = false;

    /* radial median filter stuff */
    unsigned int bin,*pixels_in,*bin_of;
    double **bin_start;
    double median,mad,deviate,sign;
    double sum_arej,avg_arej,sumd_arej,rms_arej,rmsd_arej;

    /* pythony version of structure factors, converted by init_Fhkl */
    indices pythony_indices;
    af::shared<double> pythony_amplitudes;

    /* pythony version of amorphous structure factor table vs sin(theta)/lambda, converted by init_stolFbg */
    af::shared<vec2>  pythony_stolFbg;

    /* intensity stats */
    double I,I_bg;
    double max_I; // = 0.0;
    double max_I_x,max_I_y; // = 0.0; location of max pixel value
    double photon_scale; // = 0.0 ;arbitrary "photon scale" applied before calculating noise
    double intfile_scale; // = 0.0 ;arbitrary scale applied before rounding off to ints
    double pgm_scale; // = 0.0 ;arbitrary scale applied before rounding off to char
    double sum,sumsqr,avg,rms,rmsd;
    int sumn; // = 0;
    int overloads; // = 0;

    /* image file data */
    double *floatimage;
    /* version of image to pass back to Python */
    af::flex_double raw_pixels;
    unsigned short int *intimage;
    unsigned char *pgmimage;
//    char *byte_order; // = get_byte_order();
    /* optional input image to extract background? */
//    SMVinfo imginfile;
//    double *imginfileimage;
    /* optional mask file to speed up rendering */
//    SMVinfo maskfile;
    unsigned short int *maskimage; // = NULL;

    /* misc variables */
    int i,j,n;
    double X,Y,Z;
    double ratio,r;
    double X0,Y0,Z0,d_r;
    double test;
    double vector[4];
    double newvector[4];

    /* random number seeds */
    long seed;   // default: seed = -time((time_t *)0);
//    printf("random number seed = %u\n",seed);
    long mosaic_seed; // = 12345678;  separate seed for mosaic domains so they don't correlate with other noise
    long calib_seed;  // = different seed for calibration error, since this is the same for all images

    /* point-spread function parameters */
    shapetype psf_type;
    double psf_fwhm;
    int psf_radius;
    double photons,photons0,adu;
    double readout_noise, flicker_noise;
    double calibration_noise;  // = 0.03
    double spot_scale; // = 1
    double quantum_gain;   // = 1

    /* interpolation arrays */
    int interpolate; // = 2;
    double ***sub_Fhkl;  // = NULL
    int    h_interp[5],k_interp[5],l_interp[5];
    double h_interp_d[5],k_interp_d[5],l_interp_d[5];

    double h,k,l;
    int    h0,k0,l0,h_range,k_range,l_range,h_min,h_max,k_min,k_max,l_min,l_max;
    int    h0_flr,k0_flr,l0_flr;
    int    i1,i2,i3; // =0;

    /* unit cell stuff */
    bool user_cell; // = false;
    bool user_matrix; // = False;
    double a_A[4],b_A[4],c_A[4];  // cell vectors in Angstrom
    double a[4];                  // cell vectors in meters
    double b[4];
    double c[4];
    double a0[4],b0[4],c0[4];
    double ap[4],bp[4],cp[4];
    double alpha,beta,gamma;
    double a_star[4],b_star[4],c_star[4];
    double a_star0[4],b_star0[4],c_star0[4];
    double alpha_star,beta_star,gamma_star;
    double a_cross_b[4],b_cross_c[4],c_cross_a[4];
    double a_star_cross_b_star[4],b_star_cross_c_star[4],c_star_cross_a_star[4];
    double V_cell,V_star,skew,aavg;
    double sin_alpha,sin_beta,sin_gamma;
    double cos_alpha,cos_beta,cos_gamma;
    double sin_alpha_star,sin_beta_star,sin_gamma_star;
    double cos_alpha_star,cos_beta_star,cos_gamma_star;

    /* optional user-provided unitary rotation matrix, to be applied to provide cell or A matrix */
    bool user_umat;
    double umat[10];

    /* misseting angles, applied after any provided A and U matrices */
    double misset[4];

#ifdef NANOBRAGG_HAVE_CUDA
    int device_Id;
#endif
    /* special options */
//    bool calculate_noise; // = 1;
//    bool write_pgm; // = 1;
//    bool binary_spots; // = 0; no inter-Bragg spots, flat-top spots inside FWHM of sinc function instead

    /* the constructor that takes a DXTBX detector and beam model */
    nanoBragg(const dxtbx::model::Detector&, const dxtbx::model::Beam& beam, int verbose, int panel_id = 0);

    /* the default constructor */
//    nanoBragg();

    /* member-wise constructor, allowing all members to be initialized in various ways */
    nanoBragg(
        scitbx::vec2<int> detpixels_slowfast, // = 1024, 1024
        scitbx::vec3<int> Nabc, // 1 1 1
        cctbx::uctbx::unit_cell unitcell, // lysozyme
        vec3 misset, // 0 0 0
        vec2 beam_center, // NAN NAN
        double distance, // =100,
        double pixelsize, // =0.1,
        double wavelength, // =1,
        double divergence, // =0,
        double dispersion, // =0,
        double mosaicity, // =0
        int oversample, // =0 =auto
        int vervbose); // = 1

    inline void free_all(){
      /* Based on valgrind, these are the variables that need to be free'd within this test:
         export LIBTBX_VALGRIND=valgrind --leak-check=full --tool=memcheck --suppressions=${ROOT}/cctbx_project/libtbx/valgrind-python-cci.supp
         libtbx.valgrind libtbx.python ${ROOT}/cctbx_project/simtbx/nanoBragg/tst_nanoBragg_basic.py

         Use of the free_all() class method within the Python script is meant to be an interim
         measure to prevent memory leaks, until such time as the nanoBragg class is refactored,
         with class variables that manage memory and lifetime.  These will include std::string,
         std::map, std::vector, std::shared_ptr, and flex arrays.
       */
      if (verbose)
        printf("free all memory within nanoBragg\n");
      if(verbose>9)printf("pixels_in %p\n",pixels_in);
      free(pixels_in);
      if(verbose>9)printf("bin_start %p\n",bin_start);
      free(bin_start);
      if(verbose>9)printf("source_X %p\n",source_X);
      free(source_X);
      if(verbose>9)printf("source_Y %p\n",source_Y);
      free(source_Y);
      if(verbose>9)printf("source_Z %p\n",source_Z);
      free(source_Z);
      if(verbose>9)printf("source_I %p\n",source_I);
      free(source_I);
      if(verbose>9)printf("source_lambda %p\n",source_lambda);
      free(source_lambda);
      if(verbose>9)printf("stol_of %p\n",stol_of);
      free(stol_of);
      if(verbose>9)printf("Fbg_of %p\n",Fbg_of);
      free(Fbg_of);
      if(verbose>9)printf("mosaic_umats %p\n",mosaic_umats);
      free(mosaic_umats);
      if(verbose>9)printf("invalid_pixel %p\n",invalid_pixel);
      free(invalid_pixel);
      if(verbose>9)printf("pgmimage %p\n",pgmimage);
      free(pgmimage);
      if(verbose>9)printf("intimage %p\n",intimage);
      free(intimage);
      if(verbose>9)printf("bin_of %p\n",bin_of);
      free(bin_of);
      if(verbose>9)printf("stolimage %p\n",stolimage);
      free(stolimage);
      if(verbose>9)printf("Fimage %p\n",Fimage);
      free(Fimage);
      if(verbose>9)printf("diffimage %p\n",diffimage);
      free(diffimage);
      /* free any previous allocations */
      if(Fhkl != NULL) {
        for (h0=0; h0<=h_range;h0++) {
          for (k0=0; k0<=k_range;k0++) {
            if(verbose>6) printf("freeing %d %ld-byte double Fhkl[%d][%d] at %p\n",l_range+1,sizeof(double),h0,k0,Fhkl[h0][k0]);
            free(Fhkl[h0][k0]);
          }
          if(verbose>6) printf("freeing %d %ld-byte double* Fhkl[%d] at %p\n",k_range+1,sizeof(double*),h0,Fhkl[h0]);
          free(Fhkl[h0]);
        }
        if(verbose>6) printf("freeing %d %ld-byte double** Fhkl at %p\n",h_range+1,sizeof(double**),Fhkl);
        free(Fhkl);
      }
      hkls = 0;
      if (verbose)
        printf("finished freeing memory\n");
    }

    /* member functions to run once (might allocate memory) */
    void init_defaults();               // reset all values to defaults, nulls and NANs
    void reconcile_parameters();        // call all the below functions
    void init_detector();               // decide on pixel count and allocate raw array
    void init_beam();                   // reconcile fluence, flux and beam size
    void init_beamcenter();             // select beam center convention based on initialized values
    void init_steps();          // count up and sanitize steps across divergence, dispersion, mosaic spread and spindle rotation
    void init_interpolator();   // allocate memory for 3D spline interpolation arrays
    void init_cell();           // create A matrix from cell, misset and/or file data
    void init_Fhkl();           // copy Pythony hkl and F array in/out of internal data structure
    void init_background();     // copy Pythony Fbg vs stol array in/out of internal data structure
    void init_sources();        // generate array of sources based on divergence, dispersion and wavelength
    void init_mosaicity();      // generate mosaic domains as specified

    /* member functions for reconciling inter-related parameters */
    void update_oversample();   // automatic oversampling decision based on xtal size and pixel size
    void update_beamcenter();   // beam center, Xbeam, Fbeam, ORGX using selected convention

    /* member functions for debugging */
    void show_phisteps();       // print out everything to screen, enumerate all phi steps
    void show_detector_thicksteps();  // print out all detector layers
    void show_mosaic_blocks();  // print out individual mosaic block orientations to screen
    af::shared<mat3> get_mosaic_blocks();  // get the individual mosaic block orientations as array
    void set_mosaic_blocks(af::shared<mat3>);  // set the individual mosaic block orientations from array
    void show_params();         // print out everything to screen, just like standalone program
    void show_sources();        // print out internal source information to screen

    /* member function for randomizing crystal orientation */
    void randomize_orientation();

    /* member function for triggering spot simulation over region of interest */
    void add_nanoBragg_spots();
    void add_nanoBragg_spots_nks(boost_adaptbx::python::streambuf &);
#ifdef NANOBRAGG_HAVE_CUDA
    void add_nanoBragg_spots_cuda();
#endif

    /* member function for triggering background simulation */
    void add_background(int oversample, int source);

    /* member function for extracting background from raw image */
    void extract_background(int source);

    /* member function for applying the point-spread function */
    void apply_psf(shapetype psf_type, double fwhm_pixels, int user_psf_radius);

    /* member function for triggering noise calculation */
    void add_noise();

    /* utility function for outputting an image to examine */
    void to_smv_format(std::string const& fileout, double intfile_scale, int debug_x, int debug_y);
    void to_smv_format_streambuf(boost_adaptbx::python::streambuf &, double, int const&, int const&) const;
};


}}// namespace simtbx::nanoBragg
#endif //SIMTBX_NANOBRAGG_H
