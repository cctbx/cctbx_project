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
#include <cctbx/uctbx.h>
#include <cctbx/crystal_orientation.h>
#include <cctbx/miller.h>
#include <dxtbx/model/detector.h>
#include <boost/math/special_functions/erf.hpp>
#include <boost/math/special_functions/fpclassify.hpp>

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
#define NAN strtod("NAN",NULL)
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
typedef enum { SQUARE, ROUND, GAUSS, TOPHAT, FIBER, UNKNOWN } shapetype;
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
    int source,sources;
    double *source_X,*source_Y,*source_Z,*source_I,*source_lambda;

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

    shapetype xtal_shape;
    double hrad_sqr,fudge;
    double sample_x;            /* m */
    double sample_y;            /* m */
    double sample_z;            /* m */
    double density;             /* g/m^3 */
    double molecular_weight;    /* g/mol */
    double volume,molecules;
    /* scale factor = F^2*r_e_sqr*fluence*Avogadro*volume*density/molecular_weight
                           m^2     ph/m^2  /mol      m^3   g/m^3    g/mol   */

    /* amorphous material properties */
    double amorphous_thick;
    double amorphous_default_F;
    double amorphous_MW;
    double amorphous_density;
    double amorphous_molecules;
    /* water F = 2.57 in forward direction */

    /* detector stuff */
    double pixel_size; // = 0.1e-3;
    double pixel_pos[4];
    int fpixel,spixel,fpixels,spixels,pixels;
    double distance; // = 100.0e-3;
    double detsize_f; // = 102.4e-3;
    double detsize_s; // = 102.4e-3;
    double detector_mu; //=0.0;
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
    double adc_offset; // = 40.0;


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
    double F,F_bg,*stol_of,*F_of;
    double ***Fhkl;  // = NULL
    int    hkls;
    double F_latt,F_cell;
    double default_F;
    double F_highangle,F_lowangle;
    int stols,nearest; // =0;
    double stol_file_mult; // =1.0e10;  convert to meters, usually from Angstrom
    double denom;

    /* pythony version of structure factors, converted by init_hklF */
    indices pythony_indices;
    af::shared<double> pythony_amplitudes;

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
    af::versa<double, af::c_grid<2> > raw;
    unsigned short int *intimage;
    unsigned char *pgmimage;
//    char *byte_order; // = get_byte_order();
    /* optional input image to extract background? */
//    SMVinfo imginfile;
    double *imginfileimage;
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
    double calibration_noise;
    double quantum_gain;

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
    int user_cell; // = 0;
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

    /* misseting angles, applied after any provided matrix */
    double misset[4];


    /* special options */
//    bool calculate_noise; // = 1;
//    bool write_pgm; // = 1;
//    bool binary_spots; // = 0; no inter-Bragg spots, flat-top spots inside FWHM of sinc function instead

    /* the constructor that takes a DIALS detector model */
    nanoBragg(const dxtbx::model::Detector&);

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

    /* member functions to run once (might allocate memory) */
    void init_defaults();               // reset all values to defaults, nulls and NANs
    void reconcile_parameters();        // call all the below functions
    void init_detector();               // decide on pixel count and allocate
    void init_beam();                   // reconcile fluence, flux and beam size
    void init_beamcenter();             // select beam center convention based on initialized values
    void init_steps();          // count up and sanitize steps across divergence, dispersion, mosaic spread and spindle rotation
    void init_interpolator();
    void init_cell();
    void init_Fhkl();
    void init_background();
    void init_sources();
    void init_mosaicity();

    /* member functions for reconciling inter-related parameters */
    void update_oversample();   // automatic oversampling decision based on xtal size and pixel size
    void update_steps();        // calculate total number of steps/pixel
    void update_beamcenter();   // beam center, Xbeam, Fbeam, ORGX using selected convention

    /* member functions for debugging */
    void show_phisteps();       // print out everything to screen, enumerate all phi steps
    void show_mosaic_blocks();  // print out individual mosaic block orientations to screen
    void show_params();         // print out everything to screen, just like standalone program

    /* member function for triggering spot simulation over region of interest */
    void add_nanoBragg_spots();

    /* member function for triggering background simulation */
    void add_background();

    /* member function for applying the point-spread function */
    void apply_psf(shapetype psf_type, double fwhm_pixels, int user_psf_radius);

    /* member function for triggering noise calculation */
    void add_noise();

    void to_smv_format(std::string const& fileout, double intfile_scale, double adc_offset);

};


}}// namespace simtbx::nanoBragg
#endif //SIMTBX_NANOBRAGG_H
