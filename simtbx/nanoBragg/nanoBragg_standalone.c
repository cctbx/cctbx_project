/* perfect-lattice nanocrystal diffraction simulator            -James Holton and Ken Frankel           9-17-17

example:

gcc -O -o nanoBragg nanoBragg.c -lm

./nanoBragg -mat auto.mat -hkl P1.hkl -distance 2500

./nanoBragg -mat A.mat -hkl P1.hkl -lambda 1 -dispersion 0.1 -dispstep 3 -distance 100  -detsize 100 -pixel 0.1 \
  -hdiv 0.28 -hdivstep 0.02 -vdiv 0.28 -vdivstep 0.02 \
  -fluence 1e24 -N 0 \
  -water 0

./nanoBragg -cell 74 74 36 90 90 90 -misset 10 20 30 \
  -hkl P1.hkl -lambda 1 -dispersion 0.1 -dispstep 3 -distance 100  -detsize 100 -pixel 0.1 \
  -hdiv 0.28 -hdivstep 0.02 -vdiv 0.28 -vdivstep 0.02 \
  -fluence 1e24 -N 0 \
  -water 0

lattice positions and wavelength (lambda) should be provided in Angstrom, three numbers per line
detector distance, detsize and pixel size in mm
divergence in mrad
dispersion in percent
phi and osc are in degrees
fluence is in photons/meter^2 (integrated exposure time)
Na, Nb, Nc, are the number of unit cells along the a,b,c axes, respectively
    note that any of Na,Nb,Nc can be zero to simulate an isolated unit cell (SAXS)
water is the thickness in microns of "water" also traversed by the beam
    this generates a simplitic background: that from a material with density 1.0 and isotropic
    structure factor of 2.57 electrons (the forward-scattered structure factor of water
    more complicated backgrounds can be made in a separate run of this program using Na=Nb=Nc=0.

auto.mat can be an orientation matrix from MOSFLM, or simply a text file of the
three reciprocal lattice vector components along x,y,z:
a_star_x b_star_x c_star_x
a_star_y b_star_y c_star_y
a_star_z b_star_z c_star_z

you can also simply specify the unit cell with -cell and some miss-setting angles with -misset

P1.hkl should be a text file containing
h k l F
for EVERY spot that has an intensity (including F000).  No symmetry operators will
be imposed by this program.  Not even Friedel symmetry.

Since reading the HKL file can often be the slowest step, this program will create
a binary "dumpfile" in the current working directory that it will re-read upon
subsequent runs if -hkl is not specified.

Please note that unlike nearBragg, this program does not work in the near field,
so detector distances should always be much larger than the crystal size

 */

#define _USE_MATH_DEFINES
#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <limits.h>
#include <float.h>
#ifndef NAN
#define NAN strtod("NAN",NULL)
#endif

#define TRUE 1
#define FALSE 0
#define Avogadro 6.02214179e23

/* read in text file into double arrays at provided addresses */
size_t read_text_file(char *filename, size_t nargs, ... );

/* cubic spline interpolation functions */
void polint(double *xa, double *ya, double x, double *y);
void polin2(double *x1a, double *x2a, double **ya, double x1,double x2, double *y);
void polin3(double *x1a, double *x2a, double *x3a, double ***ya, double x1,double x2, double x3, double *y);



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


/* generate unit vector in random direction */
float uniform3Ddev(float *dx, float *dy, float *dz, long *idum);
/* generate random unitary rotation matrix within a spherical cap */
double *mosaic_rotation_umat(float mosaicity, double umat[9], long *idum);
/* convert unitary matrix into missetting angles */
double *umat2misset(double umat[9],double *missets);
/* random deviate with Poisson distribution */
float poidev(float xm, long *idum);
/* random deviate with Gaussian distribution */
float gaussdev(long *idum);
/* random deviate with Lorentzian distribution */
float lorentzdev(long *idum);
/* random deviate with triangle-shaped distribution */
float triangledev(long *idum);
/* random deviate with exponential distribution (>0) */
float expdev(long *idum);
/* random deviate with uniform distribution */
float ran1(long *idum);

/* Fourier transform of a truncated lattice */
double sincg(double x, double N);
/* Fourier transform of a sphere */
double sinc3(double x);
/* Fourier transform of a spherically-truncated lattice */
double sinc_conv_sinc3(double x);


/* file stuff */
char *matfilename = NULL;
char *hklfilename = NULL;
char *dumpfilename = "Fdump.bin\0";
char *stolfilename = NULL;
char *imginfilename = NULL;
char *maskfilename = NULL;
char *stoloutfilename = "output.stol\0";
char *sourcefilename = NULL;
char *floatfilename = "floatimage.bin\0";
//char *sinfilename = "sinimage.bin\0";
//char *cosfilename = "cosimage.bin\0";
char *intfilename = "intimage.img\0";
char *pgmfilename = "image.pgm\0";
char *noisefilename = "noiseimage.img\0";
FILE *infile = NULL;
FILE *Fdumpfile = NULL;
FILE *outfile = NULL;
FILE *stoloutfile = NULL;

typedef enum { SAMPLE, BEAM } pivot;
typedef enum { SQUARE, ROUND, GAUSS, TOPHAT } shapetype;
typedef enum { CUSTOM, ADXV, MOSFLM, XDS, DIALS, DENZO } convention;

/* frame handling routines */
typedef struct _SMVinfo
{
        char *filename;
        FILE *handle;
        int swap_bytes;
        int header_size;
        int width;
        int height;
        char *header;
        unsigned short int *mmapdata;
} SMVinfo;

/* SMV image handling routines */
SMVinfo GetFrame(char *filename);
double ValueOf( const char *keyword, SMVinfo smvfile);
char *get_byte_order();
unsigned char *read_pgm5_bytes(char *filename,unsigned int *returned_width,unsigned int *returned_height);



int main(int argc, char** argv)
{
    /* progress meter stuff */
    long progress_pixel,progress_pixels;
    int progress_meter=1;
    int babble=1;
    int printout = 0;
    int printout_spixel,printout_fpixel=-1;

    /* x-ray beam properties */
    double beam_vector[4]  = {0,1,0,0};
    int coherent = 0;
    int far_source = 1;
    int round_div = 1;
    double lambda,*lambda_of;
    double mosaic_spread=-1.0,*mosaic_umats,mosaic_missets[4];
    double umat[9];
    double dispersion=0.0,dispstep=-1,lambda0 = 1.0e-10;
    double hdiv,hdivstep=-1.0,hdivrange= -1.0;
    double vdiv,vdivstep=-1.0,vdivrange= -1.0;
    double source_path,source_distance = 10.0;
    int divsteps=-1,hdivsteps=-1,vdivsteps=-1,dispsteps=-1;
    int hdiv_tic,vdiv_tic,disp_tic,mos_tic;
    int mosaic_domains=-1;
    double weight;
    int source,sources;
    double *source_X,*source_Y,*source_Z,*source_I,*source_lambda;


    /* Thomson cross section (m^2) */
    double r_e_sqr = 7.94079248018965e-30;
    /* incident x-ray fluence in photons/m^2 */
    double fluence = 125932015286227086360700780544.0;
    double flux=0.0,exposure=1.0,beamsize=1e-4;

    /* sample size stuff */
    int    N=1;
    double Na=1.0,Nb=1.0,Nc=1.0;
    double xtalsize_max,xtalsize_a,xtalsize_b,xtalsize_c;
    double reciprocal_pixel_size;

    shapetype xtal_shape = SQUARE;
    double hrad_sqr,fudge=1;
    double sample_x   = 0;              /* m */
    double sample_y   = 0;              /* m */
    double sample_z   = 0;              /* m */
    double density    = 1.0e6;          /* g/m^3 */
    double molecular_weight = 18.0;     /* g/mol */
    double volume=0.0,molecules = 0.0;
    /* scale factor = F^2*r_e_sqr*fluence*Avogadro*volume*density/molecular_weight
                           m^2     ph/m^2  /mol      m^3   g/m^3    g/mol   */
    double water_size = 0.0;
    double water_F = 2.57;
    double water_MW = 18.0;
    /* water F = 2.57 in forward direction */

    /* detector stuff */
    double pixel_size = 0.1e-3;
    double pixel_pos[4];
    int fpixel,spixel,fpixels=0,spixels=0,pixels;
    double distance = 100.0e-3;
    double detsize_f = 102.4e-3;
    double detsize_s = 102.4e-3;
    double detector_mu=0.0,detector_thick=0.0,detector_thickstep,parallax,capture_fraction;
    int    detector_thicksteps=-1,thick_tic;
    double fdet_vector[4]  = {0,0,0,1};
    double sdet_vector[4]  = {0,0,-1,0};
    double odet_vector[4]  = {0,1,0,0};
    double pix0_vector[4]  = {0,0,0,0};
    double detector_rotx=0.0,detector_roty=0.0,detector_rotz=0.0;
    double twotheta_axis[4] = {0,0,1,0};
    pivot detector_pivot = BEAM;
    convention beam_convention = MOSFLM;
    double detector_twotheta = 0.0;
    double airpath,omega_pixel,omega_Rsqr_pixel,omega_sum;
    int curved_detector = 0;
    int point_pixel= 0;
    /* beam center value that goes into the image header */
    double Xbeam=NAN,Ybeam=NAN;
    /* direct beam coordinate on fast/slow pixel axes; used for diffraction if pivot=beam */
    double Fbeam=NAN,Sbeam=NAN;
    double Fdet,Sdet,Odet;
    double Fdet0,Sdet0;
    /* nearest point on detector for detector at rotations=0 */
    double Xclose=NAN,Yclose=NAN,close_distance=NAN;
    /* near point in fast/slow pixel units; used for diffraction if pivot=sample */
    double Fclose=NAN,Sclose=NAN;
    /* fast/slow near-point position in pixels */
    double ORGX=NAN,ORGY=NAN;
    /* similar to pix0,vector but with dials-default vectors */
    double dials_origin[4] = {0,0,0,0};
    double adc_offset = 40.0;


    /* scattering vectors */
    double incident[4];
    double diffracted[4],diffracted0[4];
    double scattering[4];
    double stol,twotheta,theta;

    /* diffraction geometry stuff */
    double costwotheta,sintwotheta,psi=0;
    double xd,yd,zd,xd0,yd0,zd0;
    double Ewald[4],Ewald0[4],relp[4];
    double dmin=0;
    int integral_form = 0;

    /* polarization stuff */
    double polar_vector[4] = {0,0,0,1};
    double vert_vector[4];
    double polar=1.0,polarization=0.0;
    int nopolar = 0;

    /* sampling */
    int steps;
    int roi_xmin=-1,roi_xmax=-1,roi_ymin=-1,roi_ymax=-1;
    int oversample = -1,recommended_oversample,subS,subF;
    double subpixel_size;

    /* spindle */
    double phi,phi0=0.0,phistep=-1.0,osc=-1.0;
    int phi_tic,phisteps=-1;
    double spindle_vector[4] = {0,0,0,1};

    /* structure factor representation */
    double phase,Fa,Fb;
    double F,F_bg,*stol_of,*F_of;
    double ***Fhkl;
    double default_F=0.0;
    int    hkls=0;
    double F_latt,F_cell;
    double F_highangle,F_lowangle;
    int stols,nearest=0;
    double stol_file_mult=1.0e10;
    double denom;


    /* intensity stats */
    double I,I_bg;
    double max_I = 0.0;
    double max_I_x = 0.0,max_I_y = 0.0;
    double intfile_scale = 0.0;
    double pgm_scale = 0.0;
    double sum,sumsqr,avg,rms,rmsd;
    int sumn = 0;
    int overloads = 0;

    /* image file data */
    float *floatimage;
    SMVinfo maskfile;
    unsigned short int *maskimage = NULL;
//    float *sinimage;
//    float *cosimage;
    unsigned short int *intimage;
    unsigned char *pgmimage;
    char *byte_order = get_byte_order();
    SMVinfo imginfile;
    float *imginfileimage;

    /* misc variables */
    int i,j,n;
    double X,Y,Z;
    double ratio,r;
    double X0,Y0,Z0,d_r;
    double RTD=180.0*M_1_PI;
    double test;
    double vector[4];
    double newvector[4];

    long seed;
    seed = -time((time_t *)0);
//    printf("random number seed = %u\n",seed);
    long mosaic_seed = -12345678;

    /* interpolation arrays */
    int interpolate = 2;
    double ***sub_Fhkl;
    int    h_interp[5],k_interp[5],l_interp[5];
    double h_interp_d[5],k_interp_d[5],l_interp_d[5];

    double h,k,l;
    int    h0,k0,l0,h_range,k_range,l_range,h_min,h_max,k_min,k_max,l_min,l_max;
    int    h0_flr,k0_flr,l0_flr;
    int    i1=0, i2=0, i3=0;


    /* unit cell stuff */
    int user_cell = 0;
    double a[4] = {0,0,0,0};
    double b[4] = {0,0,0,0};
    double c[4] = {0,0,0,0};
    double a0[4],b0[4],c0[4];
    double ap[4],bp[4],cp[4];
    double alpha=0.0,beta=0.0,gamma=0.0;
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
    double misset[4] = {0,0,0,0};


    /* special options */
    int calculate_noise = 1;
    int write_pgm = 1;



    /* check argument list */
    for(i=1; i<argc; ++i)
    {
        if(argv[i][0] == '-')
        {
            /* option specified */
            if(strstr(argv[i], "-img") && (argc > (i+1)))
            {
                imginfilename = argv[i+1];
            }
            if(strstr(argv[i], "-mask") && (argc > (i+1)))
            {
                maskfilename = argv[i+1];
            }
        }
    }



    /* read in any provided mask file */
    if(maskfilename != NULL)
    {
        /* frame handling routines */
        maskfile = GetFrame(maskfilename);
        if(maskfile.header_size > 0) {
            fpixels = maskfile.width;
            spixels = maskfile.height;
            pixels = fpixels*spixels;
            test = ValueOf("PIXEL_SIZE",maskfile);
            if(! isnan(test)) pixel_size = test/1000.0;
            detsize_f = pixel_size*fpixels;
            detsize_s = pixel_size*spixels;
            test = ValueOf("DISTANCE",maskfile);
            if(! isnan(test)) distance = test/1000.0;
            test = ValueOf("CLOSE_DISTANCE",maskfile);
            if(! isnan(test)) close_distance = test/1000.0;
            test = ValueOf("WAVELENGTH",maskfile);
            if(! isnan(test)) lambda0 = test/1e10;
            test = ValueOf("BEAM_CENTER_X",maskfile);
            if(! isnan(test)) Xbeam = test/1000.0;
            test = ValueOf("BEAM_CENTER_Y",maskfile);
            if(! isnan(test)) Ybeam = detsize_s - test/1000.0;
            test = ValueOf("ORGX",maskfile);
            if(! isnan(test)) ORGX = test;
            test = ValueOf("ORGY",maskfile);
            if(! isnan(test)) ORGY = test;
            test = ValueOf("PHI",maskfile);
            if(! isnan(test)) phi0 = test/RTD;
            test = ValueOf("OSC_RANGE",maskfile);
            if(! isnan(test)) osc = test/RTD;
            test = ValueOf("TWOTHETA",maskfile);
            if(! isnan(test)) twotheta = test/RTD;

            maskimage = (unsigned short int*) calloc(pixels+10,sizeof(unsigned short int));
            j = maskfile.header_size / sizeof(unsigned short int);
            for(i=0;i<pixels;++i){
                maskimage[i] = (float) maskfile.mmapdata[j];
                 ++j;
            }
        }
    }

    /* read in any provided img file (mostly for the header) */
    if(imginfilename != NULL)
    {
        /* frame handling routines */
        imginfile = GetFrame(imginfilename);
        if(imginfile.header_size > 0) {
            fpixels = imginfile.width;
            spixels = imginfile.height;
            pixels = fpixels*spixels;
            test = ValueOf("PIXEL_SIZE",imginfile);
            if(! isnan(test)) pixel_size = test/1000.0;
            detsize_f = pixel_size*fpixels;
            detsize_s = pixel_size*spixels;
            test = ValueOf("DISTANCE",imginfile);
            if(! isnan(test)) distance = test/1000.0;
            test = ValueOf("CLOSE_DISTANCE",imginfile);
            if(! isnan(test)) close_distance = test/1000.0;
            test = ValueOf("WAVELENGTH",imginfile);
            if(! isnan(test)) lambda0 = test/1e10;
            test = ValueOf("BEAM_CENTER_X",imginfile);
            if(! isnan(test)) Xbeam = test/1000.0;
            test = ValueOf("BEAM_CENTER_Y",imginfile);
            if(! isnan(test)) Ybeam = test/1000.0;
            test = ValueOf("ORGX",imginfile);
            if(! isnan(test)) ORGX = test;
            test = ValueOf("ORGY",imginfile);
            if(! isnan(test)) ORGY = test;
            test = ValueOf("PHI",imginfile);
            if(! isnan(test)) phi0 = test/RTD;
            test = ValueOf("OSC_RANGE",imginfile);
            if(! isnan(test)) osc = test/RTD;
            test = ValueOf("TWOTHETA",imginfile);
            if(! isnan(test)) twotheta = test/RTD;

            imginfileimage = (float *) calloc(pixels+10,sizeof(float));
            j = imginfile.header_size / sizeof(unsigned short int);
            for(i=0;i<pixels;++i){
                imginfileimage[i] = (float) imginfile.mmapdata[j];
                 ++j;
            }
        }
    }


    /* check argument list for options */
    for(i=1; i<argc; ++i)
    {
        if(argv[i][0] == '-')
        {
            /* option specified */
            if(strstr(argv[i], "-Na") && (argc > (i+1)))
            {
                Na = atoi(argv[i+1]);
            }
            if(strstr(argv[i], "-Nb") && (argc > (i+1)))
            {
                Nb = atoi(argv[i+1]);
            }
            if(strstr(argv[i], "-Nc") && (argc > (i+1)))
            {
                Nc = atoi(argv[i+1]);
            }
            if(0==strcmp(argv[i], "-N") && (argc > (i+1)))
            {
                Na = Nb = Nc = atoi(argv[i+1]);
            }
            if(strstr(argv[i], "-cell") && (argc > (i+1)))
            {
                user_cell = 1;
                if(argc <= (i+1)) continue;
                if(argv[i+1][0] == '-') continue;
                a[0] = atof(argv[i+1]);
                if(argc <= (i+2)) continue;
                if(argv[i+2][0] == '-') continue;
                b[0] = atof(argv[i+2]);
                if(argc <= (i+3)) continue;
                if(argv[i+3][0] == '-') continue;
                c[0] = atof(argv[i+3]);
                if(argc <= (i+4)) continue;
                if(argv[i+4][0] == '-') continue;
                alpha = atof(argv[i+4])/RTD;
                if(argc <= (i+5)) continue;
                if(argv[i+5][0] == '-') continue;
                beta  = atof(argv[i+5])/RTD;
                if(argc <= (i+6)) continue;
                if(argv[i+6][0] == '-') continue;
                gamma = atof(argv[i+6])/RTD;
            }
            if(strstr(argv[i], "-misset") && (argc > (i+1)))
            {
                if(strstr(argv[i+1],"rand"))
                {
                    misset[0] = -1;
                    continue;
                }
            }
            if(strstr(argv[i], "-misset") && (argc > (i+3)))
            {
                misset[0] = 1;
                misset[1] = atof(argv[i+1])/RTD;
                misset[2] = atof(argv[i+2])/RTD;
                misset[3] = atof(argv[i+3])/RTD;
            }
            if((strstr(argv[i], "-samplesize") || strstr(argv[i], "-sample_size")) && (argc > (i+1)))
            {
                sample_x = atof(argv[i+1])/1000;
                sample_y = atof(argv[i+1])/1000;
                sample_z = atof(argv[i+1])/1000;
            }
            if((strstr(argv[i], "-sample_thick") || strstr(argv[i], "-sample_x") || strstr(argv[i], "-thick")) && (argc > (i+1)))
            {
                sample_x = atof(argv[i+1])/1000;
            }
            if((strstr(argv[i], "-sample_width") || strstr(argv[i], "-sample_y")  || strstr(argv[i], "-width")) && (argc > (i+1)))
            {
                sample_y = atof(argv[i+1])/1000;
            }
            if((strstr(argv[i], "-sample_heigh") || strstr(argv[i], "-sample_z")  || strstr(argv[i], "-heigh")) && (argc > (i+1)))
            {
                sample_z = atof(argv[i+1])/1000;
            }
            if((strstr(argv[i], "-xtalsize") || strstr(argv[i], "-xtal_size")) && (argc > (i+1)))
            {
                sample_x = atof(argv[i+1])/1000;
                sample_y = atof(argv[i+1])/1000;
                sample_z = atof(argv[i+1])/1000;
            }
            if((strstr(argv[i], "-xtal_thick") || strstr(argv[i], "-xtal_x") || strstr(argv[i], "-thick")) && (argc > (i+1)))
            {
                sample_x = atof(argv[i+1])/1000;
            }
            if((strstr(argv[i], "-xtal_width") || strstr(argv[i], "-xtal_y")  || strstr(argv[i], "-width")) && (argc > (i+1)))
            {
                sample_y = atof(argv[i+1])/1000;
            }
            if((strstr(argv[i], "-xtal_heigh") || strstr(argv[i], "-xtal_z")  || strstr(argv[i], "-heigh")) && (argc > (i+1)))
            {
                sample_z = atof(argv[i+1])/1000;
            }
            if((strstr(argv[i], "-density") || strstr(argv[i], "-sample_den")) && (argc > (i+1)))
            {
                density = atof(argv[i+1])*1e6;
            }
            if((0==strcmp(argv[i], "-MW") || strstr(argv[i], "-molec")) && (argc > (i+1)))
            {
                molecular_weight = atof(argv[i+1]);
            }
            if(strstr(argv[i], "-Xbeam") && (argc > (i+1)))
            {
                Xbeam = atof(argv[i+1])/1000.0;
                detector_pivot = BEAM;
            }
            if(strstr(argv[i], "-Ybeam") && (argc > (i+1)))
            {
                Ybeam = atof(argv[i+1])/1000.0;
                detector_pivot = BEAM;
            }
            if(strstr(argv[i], "-Xclose") && (argc > (i+1)))
            {
                Xclose = atof(argv[i+1])/1000.0;
                detector_pivot = SAMPLE;
            }
            if(strstr(argv[i], "-Yclose") && (argc > (i+1)))
            {
                Yclose = atof(argv[i+1])/1000.0;
                detector_pivot = SAMPLE;
            }
            if(strstr(argv[i], "-ORGX") && (argc > (i+1)))
            {
                ORGX = atof(argv[i+1]);
                detector_pivot = SAMPLE;
            }
            if(strstr(argv[i], "-ORGY") && (argc > (i+1)))
            {
                ORGY = atof(argv[i+1]);
                detector_pivot = SAMPLE;
            }
            if(strstr(argv[i], "-pivot") && (argc > (i+1)))
            {
                if(strstr(argv[i+1], "sample")) detector_pivot = SAMPLE;
                if(strstr(argv[i+1], "beam")) detector_pivot = BEAM;
            }
            if(strstr(argv[i], "-mosflm"))
            {
                beam_convention = MOSFLM;
                detector_pivot = BEAM;
            }
            if(strstr(argv[i], "-xds"))
            {
                beam_convention = XDS;
                detector_pivot = SAMPLE;
            }
            if(strstr(argv[i], "-adxv"))
            {
                beam_convention = ADXV;
                detector_pivot = BEAM;
            }
            if(strstr(argv[i], "-denzo"))
            {
                beam_convention = DENZO;
                detector_pivot = BEAM;
            }
            if(strstr(argv[i], "-dials"))
            {
                beam_convention = DIALS;
                detector_pivot = BEAM;
            }
            if(strstr(argv[i], "-fdet_vector") && (argc > (i+3)))
            {
                beam_convention = CUSTOM;
                fdet_vector[1] = atof(argv[i+1]);
                fdet_vector[2] = atof(argv[i+2]);
                fdet_vector[3] = atof(argv[i+3]);
            }
            if(strstr(argv[i], "-sdet_vector") && (argc > (i+3)))
            {
                beam_convention = CUSTOM;
                sdet_vector[1] = atof(argv[i+1]);
                sdet_vector[2] = atof(argv[i+2]);
                sdet_vector[3] = atof(argv[i+3]);
            }
            if(strstr(argv[i], "-odet_vector") && (argc > (i+3)))
            {
                beam_convention = CUSTOM;
                odet_vector[1] = atof(argv[i+1]);
                odet_vector[2] = atof(argv[i+2]);
                odet_vector[3] = atof(argv[i+3]);
            }
            if(strstr(argv[i], "-beam_vector") && (argc > (i+3)))
            {
                beam_convention = CUSTOM;
                beam_vector[1] = atof(argv[i+1]);
                beam_vector[2] = atof(argv[i+2]);
                beam_vector[3] = atof(argv[i+3]);
            }
            if(strstr(argv[i], "-polar_vector") && (argc > (i+3)))
            {
                beam_convention = CUSTOM;
                polar_vector[1] = atof(argv[i+1]);
                polar_vector[2] = atof(argv[i+2]);
                polar_vector[3] = atof(argv[i+3]);
            }
            if(strstr(argv[i], "-spindle_axis") && (argc > (i+3)))
            {
                beam_convention = CUSTOM;
                spindle_vector[1] = atof(argv[i+1]);
                spindle_vector[2] = atof(argv[i+2]);
                spindle_vector[3] = atof(argv[i+3]);
            }
            if(strstr(argv[i], "-twotheta_axis") && (argc > (i+3)))
            {
                beam_convention = CUSTOM;
                twotheta_axis[1] = atof(argv[i+1]);
                twotheta_axis[2] = atof(argv[i+2]);
                twotheta_axis[3] = atof(argv[i+3]);
            }
            if(strstr(argv[i], "-pix0_vector") && (argc > (i+3)))
            {
                beam_convention = CUSTOM;
                pix0_vector[0] = 1.0;
                pix0_vector[1] = atof(argv[i+1]);
                pix0_vector[2] = atof(argv[i+2]);
                pix0_vector[3] = atof(argv[i+3]);
            }
            if(strstr(argv[i], "-distance") && (argc > (i+1)))
            {
                distance = atof(argv[i+1])/1000.0;
                detector_pivot = BEAM;
            }
            if(strstr(argv[i], "-close_distance") && (argc > (i+1)))
            {
                close_distance = atof(argv[i+1])/1000.0;
                detector_pivot = SAMPLE;
            }
//            if(strstr(argv[i], "-source_dist") && (argc > (i+1)))
//            {
//              source_distance = atof(argv[i+1])/1000.0;
//            }
            if(strstr(argv[i], "-detector_abs") && (argc >= (i+1)))
            {
                if(strstr(argv[i+1], "inf") || atof(argv[i+1]) == 0.0) {
                    detector_thick = 0.0;
                    detector_mu = 0.0;
                }else{
                    detector_mu = 1.0/(atof(argv[i+1])*1e-6);
                }
            }
            if(strstr(argv[i], "-detector_thick") && (strlen(argv[i]) == 15) && (argc >= (i+1)))
            {
                 detector_thick = atof(argv[i+1])*1e-6;
            }
            if(strstr(argv[i], "-detector_thicksteps") && (argc >= (i+1)))
            {
                detector_thicksteps = atoi(argv[i+1]);
            }
            if(strstr(argv[i], "-twotheta") && (argc > (i+1)))
            {
                detector_twotheta = atof(argv[i+1])/RTD;
                detector_pivot = SAMPLE;
            }
            if(strstr(argv[i], "-detector_rotx") && (argc > (i+1)))
            {
                detector_rotx = atof(argv[i+1])/RTD;
            }
            if(strstr(argv[i], "-detector_roty") && (argc > (i+1)))
            {
                detector_roty = atof(argv[i+1])/RTD;
            }
            if(strstr(argv[i], "-detector_rotz") && (argc > (i+1)))
            {
                detector_rotz = atof(argv[i+1])/RTD;
            }
            if(strstr(argv[i], "-detsize") && (strlen(argv[i]) == 8) && (argc > (i+1)))
            {
                detsize_f = atof(argv[i+1])/1000.0;
                detsize_s = atof(argv[i+1])/1000.0;
            }
            if(strstr(argv[i], "-detsize_f") && (argc > (i+1)))
            {
                detsize_f = atof(argv[i+1])/1000.0;
            }
             if(strstr(argv[i], "-detsize_s") && (argc > (i+1)))
            {
                detsize_s = atof(argv[i+1])/1000.0;
            }
            if(strstr(argv[i], "-detpixels") && (strlen(argv[i]) == 10) && (argc > (i+1)))
            {
                fpixels = spixels = atoi(argv[i+1]);
            }
            if(strstr(argv[i], "-detpixels_f") && (argc > (i+1)))
            {
                fpixels = atoi(argv[i+1]);
            }
            if(strstr(argv[i], "-detpixels_s") && (argc > (i+1)))
            {
                spixels = atoi(argv[i+1]);
            }
            if(strstr(argv[i], "-curved_det") && (argc > (i+1)))
            {
                curved_detector = 1;
            }
            if(strstr(argv[i], "-pixel") && (argc > (i+1)))
            {
                pixel_size = atof(argv[i+1])/1000.0;
            }
            if(strstr(argv[i], "-point_pixel") )
            {
                point_pixel = 1;
            }
            if(strstr(argv[i], "-polar") && (strlen(argv[i]) == 6) && (argc > (i+1)))
            {
                polarization = atof(argv[i+1]);
                nopolar = 0;
            }
            if(strstr(argv[i], "-nopolar") )
            {
                nopolar = 1;
            }
            if(strstr(argv[i], "-oversample") && (argc > (i+1)))
            {
                oversample = atoi(argv[i+1]);
            }
            if(strstr(argv[i], "-roi") && (argc > (i+4)))
            {
                roi_xmin = atoi(argv[i+1]);
                roi_xmax = atoi(argv[i+2]);
                roi_ymin = atoi(argv[i+3]);
                roi_ymax = atoi(argv[i+4]);
            }
            if((strstr(argv[i], "-lambda") || strstr(argv[i], "-wave")) && (argc > (i+1)))
            {
                lambda0 = atof(argv[i+1])/1.0e10;
            }
            if(strstr(argv[i], "-energy") && (argc > (i+1)))
            {
                lambda0 = (12398.42/atof(argv[i+1]))/1.0e10;
            }
            if(strstr(argv[i], "-fluence") && (argc > (i+1)))
            {
                fluence = atof(argv[i+1]);
            }
            if(strstr(argv[i], "-flux") && (argc > (i+1)))
            {
                flux = atof(argv[i+1]);
            }
            if(strstr(argv[i], "-exposure") && (argc > (i+1)))
            {
                exposure = atof(argv[i+1]);
            }
            if(strstr(argv[i], "-beamsize") && (argc > (i+1)))
            {
                beamsize = atof(argv[i+1])/1000;
            }
            if((strstr(argv[i], "-mosaic") && (strlen(argv[i]) == 7) || strstr(argv[i], "-mosaici") || strstr(argv[i], "-mosaic_spr")) && (argc > (i+1)))
            {
                mosaic_spread = atof(argv[i+1])/RTD;
            }
            if(strstr(argv[i], "-mosaic_dom") && (argc > (i+1)))
            {
                mosaic_domains = atoi(argv[i+1]);
            }
            if(strstr(argv[i], "-dispersion") && (argc > (i+1)))
            {
                dispersion = atof(argv[i+1])/100.0;
            }
            if(strstr(argv[i], "-dispsteps") && (argc > (i+1)))
            {
                dispsteps = atoi(argv[i+1]);
            }
            if(strstr(argv[i], "-divergence") && (argc > (i+1)))
            {
                hdivrange = vdivrange = atof(argv[i+1])/1000.0;
            }
            if(strstr(argv[i], "-hdivrange") && (argc > (i+1)))
            {
                hdivrange = atof(argv[i+1])/1000.0;
            }
            if(strstr(argv[i], "-vdivrange") && (argc > (i+1)))
            {
                vdivrange = atof(argv[i+1])/1000.0;
            }
            if(strstr(argv[i], "-hdivstep") && (strlen(argv[i]) == 9) && (argc > (i+1)))
            {
                hdivstep = atof(argv[i+1])/1000.0;
            }
            if(strstr(argv[i], "-hdivsteps") && (argc > (i+1)))
            {
                hdivsteps = atoi(argv[i+1]);
            }
            if(strstr(argv[i], "-vdivstep") && (strlen(argv[i]) == 9) && (argc > (i+1)))
            {
                vdivstep = atof(argv[i+1])/1000.0;
            }
            if(strstr(argv[i], "-vdivsteps") && (argc > (i+1)))
            {
                vdivsteps = atoi(argv[i+1]);
            }
            if(strstr(argv[i], "-divsteps") && (argc > (i+1)))
            {
                hdivsteps = vdivsteps = atoi(argv[i+1]);
            }
            if(strstr(argv[i], "-round_div") )
            {
                /* cut to circle */
                round_div = 1;
            }
            if(strstr(argv[i], "-square_div") )
            {
                /* just raster */
                round_div = 0;
            }
            if(strstr(argv[i], "-adc") && (argc > (i+1)))
            {
                adc_offset = atof(argv[i+1]);
            }
            if(strstr(argv[i], "-phi") && strlen(argv[i])==4 && (argc > (i+1)))
            {
                phi0 = atof(argv[i+1])/RTD;
            }
            if(strstr(argv[i], "-osc") && (argc > (i+1)))
            {
                osc = atof(argv[i+1])/RTD;
            }
            if(strstr(argv[i], "-phistep") && strlen(argv[i])==8 && (argc > (i+1)))
            {
                phistep = atof(argv[i+1])/RTD;
            }
            if(strstr(argv[i], "-phisteps") && (argc > (i+1)))
            {
                phisteps = atoi(argv[i+1]);
            }
            if(strstr(argv[i], "-dmin") && (argc > (i+1)))
            {
                dmin = atof(argv[i+1])*1e-10;
            }
            if(strstr(argv[i], "-mat") && (argc > (i+1)))
            {
                matfilename = argv[i+1];
            }
            if(strstr(argv[i], "-hkl") && (argc > (i+1)))
            {
                hklfilename = argv[i+1];
            }
            if(strstr(argv[i], "-default_F") && (argc > (i+1)))
            {
                default_F = atof(argv[i+1]);
            }
            if(strstr(argv[i], "-img") && (argc > (i+1)))
            {
                imginfilename = argv[i+1];
            }
            if(strstr(argv[i], "-stolout") && strlen(argv[i])>7 && (argc > (i+1)))
            {
                stoloutfilename = argv[i+1];
            }
            if(strstr(argv[i], "-stol") && strlen(argv[i])==5 && (argc > (i+1)))
            {
                stolfilename = argv[i+1];
                stol_file_mult = 1e10;
            }
            if(strstr(argv[i], "-4stol") && strlen(argv[i])==6 && (argc > (i+1)))
            {
                stolfilename = argv[i+1];
                stol_file_mult = 1e10/4;
            }
            if(strstr(argv[i], "-Q") && strlen(argv[i])==2 && (argc > (i+1)))
            {
                stolfilename = argv[i+1];
                stol_file_mult = 1e10/M_PI/4;
            }
            if(strstr(argv[i], "-sourcefile") && (argc > (i+1)))
            {
                sourcefilename = argv[i+1];
            }
            if((strstr(argv[i], "-floatfile") || strstr(argv[i], "-floatimage")) && (argc > (i+1)))
            {
                floatfilename = argv[i+1];
            }
            if((strstr(argv[i], "-intfile") || strstr(argv[i], "-intimage")) && (argc > (i+1)))
            {
                intfilename = argv[i+1];
            }
            if((strstr(argv[i], "-pgmfile") || strstr(argv[i], "-pgmimage")) && (argc > (i+1)))
            {
                pgmfilename = argv[i+1];
                write_pgm = 1;
            }
            if((strstr(argv[i], "-noisefile") || strstr(argv[i], "-noiseimage")) && (argc > (i+1)))
            {
                noisefilename = argv[i+1];
                calculate_noise = 1;
            }
            if(strstr(argv[i], "-nonoise") )
            {
                /* turn off noise */
                calculate_noise = 0;
            }
            if(strstr(argv[i], "-nopgm") )
            {
                write_pgm = 0;
            }
            if(strstr(argv[i], "-scale") && (argc > (i+1)))
            {
                /* specify the scale for the intfile */
                intfile_scale = atof(argv[i+1]);
            }
            if(strstr(argv[i], "-pgmscale") && (argc > (i+1)))
            {
                /* specify the scale for the intfile */
                pgm_scale = atof(argv[i+1]);
                write_pgm = 1;
            }
            if(strstr(argv[i], "-coherent") )
            {
                /* turn off incoherent addition */
                coherent = 1;
            }
            if(strstr(argv[i], "-printout") )
            {
                /* turn on console printing */
                printout = 1;
            }
            if(strstr(argv[i], "-noprogress") )
            {
                /* turn off progress meter */
                progress_meter = 0;
            }
            if(strstr(argv[i], "-progress") )
            {
                /* turn on progress meter */
                progress_meter = 1;
            }
            if(strstr(argv[i], "-interpolate") )
            {
                /* turn on tricubic interpolation */
                interpolate = 1;
            }
            if(strstr(argv[i], "-nointerpolate") )
            {
                /* turn off tricubic interpolation */
                interpolate = 0;
            }
            if(strstr(argv[i], "-round_xtal") )
            {
                /* use sinc3 */
                xtal_shape = ROUND;
            }
            if(strstr(argv[i], "-square_xtal") )
            {
                /* use sincg */
                xtal_shape = SQUARE;
            }
            if(strstr(argv[i], "-gauss_xtal") )
            {
                /* use Gaussian */
                xtal_shape = GAUSS;
            }
            if(strstr(argv[i], "-binary_spots") || strstr(argv[i], "-tophat_spots"))
            {
                /* top hat */
                xtal_shape = TOPHAT;
            }
            if(strstr(argv[i], "-fudge") && (argc > (i+1)))
            {
                fudge = atof(argv[i+1]);
            }
            if(strstr(argv[i], "-printout_pixel") && (argc > (i+2)))
            {
                printout_fpixel = atoi(argv[i+1]);
                printout_spixel = atoi(argv[i+2]);
            }
            if(strstr(argv[i], "-seed") && (argc > (i+1)))
            {
                seed = -atoi(argv[i+1]);
            }
            if(strstr(argv[i], "-mosaic_seed") && (argc > (i+1)))
            {
                mosaic_seed = -atoi(argv[i+1]);
            }
            if(strstr(argv[i], "-water") && (argc > (i+1)))
            {
                water_size = atof(argv[i+1])/1e6;
            }
        }
    }

    /* fill in blanks */
    if(fpixels) {

        detsize_f = pixel_size*fpixels;
    }
    if(spixels) {
        detsize_s = pixel_size*spixels;
    }
    fpixels = ceil(detsize_f/pixel_size-0.5);
    spixels = ceil(detsize_s/pixel_size-0.5);
    pixels = fpixels*spixels;

    /* get fluence from flux */
    if(flux != 0.0 && exposure > 0.0 && beamsize >= 0){
        fluence = flux*exposure/beamsize/beamsize;
    }
    if(beamsize >= 0){
        if(beamsize < sample_y){
            printf("WARNING: clipping sample (%lg m high) with beam (%lg m)\n",sample_y,beamsize);
            sample_y = beamsize;
        }
        if(beamsize < sample_z){
            printf("WARNING: clipping sample (%lg m wide) with beam (%lg m)\n",sample_z,beamsize);
            sample_z = beamsize;
        }
    }
    if(exposure > 0.0)
    {
        /* make sure flux is consistent with everything else */
        flux = fluence/exposure*beamsize*beamsize;
    }

    /* straighten up sample properties */
//    volume = sample_x*sample_y*sample_z;
//    molecules = volume*density*Avogadro/molecular_weight;


    /* defaults? */
    if(! isnan(ORGX)) Fclose = (ORGX-0.5)*pixel_size;
    if(! isnan(ORGY)) Sclose = (ORGY-0.5)*pixel_size;
    /* place beam center halfway between four middle pixels */
    /* place beam center at int(npix/2) location */
    if(isnan(Fclose)) Fclose = (detsize_f - 0*pixel_size)/2.0;
    if(isnan(Sclose)) Sclose = (detsize_s + 0*pixel_size)/2.0;
    if(isnan(Xclose)) Xclose = Fclose;
    if(isnan(Yclose)) Yclose = Sclose;
    if(isnan(Fbeam)) Fbeam = Fclose;
    if(isnan(Sbeam)) Sbeam = Sclose;
    if(roi_xmin < 0) roi_xmin = 0;
    if(roi_xmax < 0) roi_xmax = fpixels;
    if(roi_ymin < 0) roi_ymin = 0;
    if(roi_ymax < 0) roi_ymax = spixels;
    progress_pixels = (roi_xmax-roi_xmin+1)*(roi_ymax-roi_ymin+1);

    if(beam_convention == ADXV)
    {
        /* first pixel is at 0,0 pix and pixel_size,pixel_size*npixels mm */
        if(isnan(Xbeam)) Xbeam = (detsize_f + pixel_size)/2.0;
        if(isnan(Ybeam)) Ybeam = (detsize_s - pixel_size)/2.0;
           beam_vector[1]=  0;    beam_vector[2]=  0;    beam_vector[3]=  1;
           fdet_vector[1]=  1;    fdet_vector[2]=  0;    fdet_vector[3]=  0;
           sdet_vector[1]=  0;    sdet_vector[2]= -1;    sdet_vector[3]=  0;
           odet_vector[1]=  0;    odet_vector[2]=  0;    odet_vector[3]=  1;
         twotheta_axis[1]= -1;  twotheta_axis[2]=  0;  twotheta_axis[3]=  0;
          polar_vector[1]=  1;   polar_vector[2]=  0;   polar_vector[3]=  0;
        spindle_vector[1]=  1; spindle_vector[2]=  0; spindle_vector[3]=  0;
        Fbeam = Xbeam;
        Sbeam = detsize_s - Ybeam;
        detector_pivot = BEAM;
    }
    if(beam_convention == MOSFLM)
    {
        /* first pixel is at 0.5,0.5 pix and pixel_size/2,pixel_size/2 mm */
        if(isnan(Xbeam)) Xbeam = (detsize_s + pixel_size)/2.0;
        if(isnan(Ybeam)) Ybeam = (detsize_f + pixel_size)/2.0;
           beam_vector[1]=  1;    beam_vector[2]=  0;    beam_vector[3]=  0;
           odet_vector[1]=  1;    odet_vector[2]=  0;    odet_vector[3]=  0;
           fdet_vector[1]=  0;    fdet_vector[2]=  0;    fdet_vector[3]=  1;
           sdet_vector[1]=  0;    sdet_vector[2]= -1;    sdet_vector[3]=  0;
         twotheta_axis[1]=  0;  twotheta_axis[2]=  0;  twotheta_axis[3]= -1;
          polar_vector[1]=  0;   polar_vector[2]=  0;   polar_vector[3]=  1;
        spindle_vector[1]=  0; spindle_vector[2]=  0; spindle_vector[3]=  1;
        Fbeam = Ybeam + 0.5*pixel_size;
        Sbeam = Xbeam + 0.5*pixel_size;
        detector_pivot = BEAM;
    }
    if(beam_convention == DENZO)
    {
        if(isnan(Xbeam)) Xbeam = (detsize_s + pixel_size)/2.0;
        if(isnan(Ybeam)) Ybeam = (detsize_f + pixel_size)/2.0;
           beam_vector[1]=  1;    beam_vector[2]=  0;    beam_vector[3]=  0;
           odet_vector[1]=  1;    odet_vector[2]=  0;    odet_vector[3]=  0;
           fdet_vector[1]=  0;    fdet_vector[2]=  0;    fdet_vector[3]=  1;
           sdet_vector[1]=  0;    sdet_vector[2]= -1;    sdet_vector[3]=  0;
         twotheta_axis[1]=  0;  twotheta_axis[2]=  0;  twotheta_axis[3]= -1;
          polar_vector[1]=  0;   polar_vector[2]=  0;   polar_vector[3]=  1;
        spindle_vector[1]=  0; spindle_vector[2]=  0; spindle_vector[3]=  1;
        Fbeam = Ybeam + 0.0*pixel_size;
        Sbeam = Xbeam + 0.0*pixel_size;
        detector_pivot = BEAM;
    }
    if(beam_convention == XDS)
    {
        if(isnan(Xbeam)) Xbeam = Xclose;
        if(isnan(Ybeam)) Ybeam = Yclose;
           beam_vector[1]=  0;    beam_vector[2]=  0;    beam_vector[3]=  1;
           fdet_vector[1]=  1;    fdet_vector[2]=  0;    fdet_vector[3]=  0;
           sdet_vector[1]=  0;    sdet_vector[2]=  1;    sdet_vector[3]=  0;
           odet_vector[1]=  0;    odet_vector[2]=  0;    odet_vector[3]=  1;
         twotheta_axis[1]=  1;  twotheta_axis[2]=  0;  twotheta_axis[3]=  0;
          polar_vector[1]=  1;   polar_vector[2]=  0;   polar_vector[3]=  0;
        spindle_vector[1]=  1; spindle_vector[2]=  0; spindle_vector[3]=  0;
        Fbeam = Xbeam;
        Sbeam = Ybeam;
        detector_pivot = SAMPLE;
    }
    if(beam_convention == DIALS)
    {
        if(isnan(Xbeam)) Xbeam = Xclose;
        if(isnan(Ybeam)) Ybeam = Yclose;
           beam_vector[1]=  0;    beam_vector[2]=  0;    beam_vector[3]=  1;
           fdet_vector[1]=  1;    fdet_vector[2]=  0;    fdet_vector[3]=  0;
           sdet_vector[1]=  0;    sdet_vector[2]=  1;    sdet_vector[3]=  0;
           odet_vector[1]=  0;    odet_vector[2]=  0;    odet_vector[3]=  1;
         twotheta_axis[1]=  0;  twotheta_axis[2]=  1;  twotheta_axis[3]=  0;
          polar_vector[1]=  0;   polar_vector[2]=  1;   polar_vector[3]=  0;
        spindle_vector[1]=  0; spindle_vector[2]=  1; spindle_vector[3]=  0;
        Fbeam = Xbeam;
        Sbeam = Ybeam;
        detector_pivot = SAMPLE;
    }
    if(beam_convention == CUSTOM)
    {
        if(isnan(Xbeam)) Xbeam = Xclose;
        if(isnan(Ybeam)) Ybeam = Yclose;
        Fbeam = Xbeam;
        Sbeam = Ybeam;
        Fclose = Xbeam;
        Sclose = Ybeam;
    }

    /* straighten up vectors */
    unitize(beam_vector,beam_vector);
    unitize(fdet_vector,fdet_vector);
    unitize(sdet_vector,sdet_vector);
    if(unitize(odet_vector,odet_vector) != 1.0)
    {
        printf("WARNING: auto-generating odet_vector\n");
        cross_product(fdet_vector,sdet_vector,odet_vector);
        unitize(odet_vector,odet_vector);
    }
    unitize(polar_vector,polar_vector);
    unitize(spindle_vector,spindle_vector);
    cross_product(beam_vector,polar_vector,vert_vector);
    unitize(vert_vector,vert_vector);


    printf("nanoBragg nanocrystal diffraction simulator - James Holton and Ken Frankel 5-17-17\n");

    if(hklfilename == NULL)
    {
        /* see if there are Fs from a previous run */
        Fdumpfile = fopen(dumpfilename,"r");
        if(Fdumpfile == NULL && default_F == 0.0)
        {
            printf("ERROR: no hkl file and no dump file to read.");
        }
    }

    if(hklfilename == NULL && Fdumpfile == NULL && default_F == 0.0 || matfilename == NULL && a[0] == 0.0){
        printf("usage: nanoBragg -mat auto.mat -hkl Fs.hkl\n");
        printf("options:\n");\
        printf("\t-mat filename.mat\tmosflm-style matrix file containing three reciprocal unit cell vectors\n");
        printf("\t-hkl filename.hkl\ttext file containing h, k, l and F for P1 unit cell\n");
        printf("\t-distance        \tdistance from origin to detector center in mm\n");
        printf("\t-detsize         \tdetector size in mm.  may also use -detsize_f -detsize_s\n");
        printf("\t-detpixels       \tdetector size in pixels.  may also use -detpixels_x -detpixels_y\n");
        printf("\t-pixel           \tdetector pixel size in mm.\n");
        printf("\t-detector_absorb \tdetector sensor material attenuation depth (um) (default: \"inf\" to save time)\n");
        printf("\t-detector_thick  \tdetector sensor thickness (um)\n");
        printf("\t-detector_thicksteps\tnumber of layers of detector sensor material. Default: 1\n");
        printf("\t-Xbeam           \timage fast coordinate of direct-beam spot (mm). (default: center)\n");
        printf("\t-Ybeam           \timage slow coordinate of direct-beam spot (mm). (default: center)\n");
        printf("\t-mosflm          \tuse MOSFLM's direct-beam convention. (default: adxv)\n");
        printf("\t-xds             \tuse XDS detector origin convention. (default: adxv)\n");
        printf("\t-twotheta        \trotation of detector about spindle axis (deg). (default: 0)\n");
        printf("\t-N               \tnumber of unit cells in all directions. may also use -Na -Nb or -Nc\n");
        printf("\t-square_xtal     \tspecify parallelpiped crystal shape (default)\n");
        printf("\t-round_xtal      \tspecify ellipsoidal crystal shape (sort of)\n");
        printf("\t-tophat_spots    \tclip lattice transform at fwhm: no inter-Bragg maxima\n");
        printf("\t-oversample      \tnumber of sub-pixels per pixel. use this if xtalsize/lambda > distance/pixel\n");
        printf("\t-lambda          \tincident x-ray wavelength in Angstrom. may also use -energy in eV\n");
        printf("\t-mosaic          \tisotropic mosaic spread in degrees (use 90 for powder)\n");
        printf("\t-mosaic_domains  \tnumber of randomly-oriented mosaic domains to render\n");
        printf("\t-dispersion      \tspectral dispersion: delta-lambda/lambda in percent\n");
        printf("\t-dispsteps       \tnumber of wavelengths in above range\n");
        printf("\t-hdivrange       \thorizontal angular spread of source points in mrad\n");
        printf("\t-vdivrange       \tvertical angular spread of source points in mrad\n");
        printf("\t-hdivstep        \tnumber of source points in the horizontal\n");
        printf("\t-vdivstep        \tnumber of source points in the vertical\n");
        printf("\t-square_div      \tfull divergence grid (default: round off corners)\n");
        printf("\t-phi             \tstarting rotation value about spindle axis in degrees\n");
        printf("\t-osc             \trotation range about spindle axis in degrees\n");
        printf("\t-phisteps        \tnumber of rotation steps to render\n");
        printf("\t-water           \tadd contribution of x microns of water surrounding crystal\n");
        printf("\t-floatfile       \tname of binary output file (4-byte floats)\n");
        printf("\t-intfile         \tname of noiseless smv-formatted output file (not on absolute scale by default)\n");
        printf("\t-scale           \tscale factor to apply to intfile (default: autoscale)\n");
        printf("\t-noisefile       \tname of photon-scale smv-formatted output file (with Poisson noise)\n");
        printf("\t-roi             \tonly render part of the image: xmin xmax ymin ymax\n");
        printf("\t-printout        \tprint pixel values out to the screen\n");
        printf("\t-seed            \tspecify random-number seed for noisefile\n");
        printf("\t-fluence         \tincident beam intensity for photon-counting statistics (photons/m^2)\n");
        printf("\t-nonoise         \tdisable generating the noisefile\n");
        printf("\t-noprogress      \tturn off the progress meter\n");
        printf("\t-nopolar         \tturn off the polarization correction\n");
        printf("\t-nointerpolate   \tdisable inter-Bragg peak structure factor interpolation\n");
        printf("\t-interpolate     \tforce inter-Bragg peak structure factor interpolation (default: on if < 3 cells wide)\n");
        printf("\t-point_pixel     \tturn off the pixel solid angle correction\n");
        printf("\t-curved_det      \tall pixels same distance from crystal\n");
        printf("\t-fdet_vector     \tunit vector of increasing fast-axis detector pixel coordinate (default: %g %g %g)\n",fdet_vector[1],fdet_vector[2],fdet_vector[3]);
        printf("\t-sdet_vector     \tunit vector of increasing slow-axis detector pixel coordinate (default: %g %g %g)\n",sdet_vector[1],sdet_vector[2],sdet_vector[3]);
        printf("\t-odet_vector     \tunit vector of increasing detector distance (default: %g %g %g)\n",odet_vector[1],odet_vector[2],odet_vector[3]);
        printf("\t-beam_vector     \tunit vector of x-ray beam direction (default: %g %g %g)\n",beam_vector[1],beam_vector[2],beam_vector[3]);
        printf("\t-polar_vector    \tunit vector of x-ray E-vector polarization (default: %g %g %g)\n",polar_vector[1],polar_vector[2],polar_vector[3]);
        printf("\t-spindle_axis    \tunit vector of right-handed phi rotation axis (default: %g %g %g)\n",spindle_vector[1],spindle_vector[2],spindle_vector[3]);
        printf("\t-pix0_vector     \tvector from crystal to first pixel in image (default: beam centered on detector)\n");
//        printf("\t-source_distance \tdistance of x-ray source from crystal (default: 10 meters)\n");
        exit(9);
    }


    /* allocate detector memory */
    floatimage = (float*) calloc(pixels+10,sizeof(float));
    //sinimage = (float*) calloc(pixels+10,2*sizeof(float));
    //cosimage = (float*) calloc(pixels+10,2*sizeof(float));
    intimage   = (unsigned short int*) calloc(pixels+10,sizeof(unsigned short int));
    if(write_pgm) pgmimage   = (unsigned char*) calloc(pixels+10,sizeof(unsigned char));


    /* default sampling logic */
    if(phisteps < 0){
        /* auto-select number of phi steps */
        if(osc < 0.0) {
            /* auto-select osc range */
            if(phistep <= 0.0) {
                /* user doesn't care about anything */
                phisteps = 1;
                osc = 0.0;
                phistep = 0.0;
            } else {
                /* user doesn't care about osc or steps, but specified step */
                osc = phistep;
                phisteps = 2;
            }
        } else {
            /* user-speficied oscillation */
            if(phistep <= 0.0) {
                /* osc specified, but nothing else */
                phisteps = 2;
                phistep = osc/2.0;
            } else {
                /* osc and phi step specified */
                phisteps = ceil(osc/phistep);
            }
        }
    } else {
        /* user-specified number of phi steps */
        if(phisteps == 0) phisteps = 1;
        if(osc < 0.0) {
            /* auto-select osc range */
            if(phistep <= 0.0) {
                /* user cares only about number of steps */
                osc = 1.0/RTD;
                phistep = osc/phisteps;
            } else {
                /* user doesn't care about osc, but specified step */
                osc = phistep;
                phisteps = 2;
            }
        } else {
            /* user-speficied oscillation */
            if(phistep < 0.0) {
                /* osc and steps specified */
                phistep = osc/phisteps;
            } else {
                /* everything specified */
            }
        }
    }

    if(hdivsteps <= 0){
        /* auto-select number of steps */
        if(hdivrange < 0.0) {
            /* auto-select range */
            if(hdivstep <= 0.0) {
                /* user doesn't care about anything */
                hdivsteps = 1;
                hdivrange = 0.0;
                hdivstep = 0.0;
            } else {
                /* user specified stepsize and nothing else */
                hdivrange = hdivstep;
                hdivsteps = 2;
            }
        } else {
            /* user-speficied range */
            if(hdivstep <= 0.0) {
                /* range specified, but nothing else */
                hdivstep = hdivrange;
                hdivsteps = 2;
            } else {
                /* range and step specified, but not number of steps */
                hdivsteps = ceil(hdivrange/hdivstep);
            }
        }
    } else {
        /* user-specified number of steps */
        if(hdivrange < 0.0) {
            /* auto-select range */
            if(hdivstep <= 0.0) {
                /* user cares only about number of steps */
                hdivrange = 1.0;
                hdivstep = hdivrange/hdivsteps;
            } else {
                /* user doesn't care about range */
                hdivrange = hdivstep;
                hdivsteps = 2;
            }
        } else {
            /* user-speficied range */
            if(hdivstep <= 0.0) {
                /* range and steps specified */
                if(hdivsteps <=1 ) hdivsteps = 2;
                hdivstep = hdivrange/(hdivsteps-1);
            } else {
                /* everything specified */
            }
        }
    }

    if(vdivsteps <= 0){
        /* auto-select number of steps */
        if(vdivrange < 0.0) {
            /* auto-select range */
            if(vdivstep <= 0.0) {
                /* user doesn't care about anything */
                vdivsteps = 1;
                vdivrange = 0.0;
                vdivstep = 0.0;
            } else {
                /* user specified stepsize and nothing else */
                vdivrange = vdivstep;
                vdivsteps = 2;
            }
        } else {
            /* user-speficied range */
            if(vdivstep <= 0.0) {
                /* range specified, but nothing else */
                vdivstep = vdivrange;
                vdivsteps = 2;
            } else {
                /* range and step specified, but not number of steps */
                vdivsteps = ceil(vdivrange/vdivstep);
            }
        }
    } else {
        /* user-specified number of steps */
        if(vdivrange < 0.0) {
            /* auto-select range */
            if(vdivstep <= 0.0) {
                /* user cares only about number of steps */
                vdivrange = 1.0;
                vdivstep = vdivrange/vdivsteps;
            } else {
                /* user doesn't care about range */
                vdivrange = vdivstep;
                vdivsteps = 2;
            }
        } else {
            /* user-speficied range */
            if(vdivstep <= 0.0) {
                /* range and steps specified */
                if(vdivsteps <=1 ) vdivsteps = 2;
                vdivstep = vdivrange/(vdivsteps-1);
            } else {
                /* everything specified */
            }
        }
    }


    if(dispsteps <= 0){
        /* auto-select number of steps */
        if(dispersion < 0.0) {
            /* auto-select range */
            if(dispstep <= 0.0) {
                /* user doesn't care about anything */
                dispsteps = 1;
                dispersion = 0.0;
                dispstep = 0.0;
            } else {
                /* user specified stepsize and nothing else */
                dispersion = dispstep;
                dispsteps = 2;
            }
        } else {
            /* user-speficied range */
            if(dispstep <= 0.0) {
                /* range specified, but nothing else */
                dispstep = dispersion;
                dispsteps = 2;
            } else {
                /* range and step specified, but not number of steps */
                dispsteps = ceil(dispersion/dispstep);
            }
        }
    } else {
        /* user-specified number of steps */
        if(dispersion < 0.0) {
            /* auto-select range */
            if(dispstep <= 0.0) {
                /* user cares only about number of steps */
                dispersion = 1.0;
                dispstep = dispersion/dispsteps;
            } else {
                /* user doesn't care about range */
                dispersion = dispstep;
                dispsteps = 2;
            }
        } else {
            /* user-speficied range */
            if(dispstep <= 0.0) {
                /* range and steps specified */
                if(dispsteps <=1 ) dispsteps = 2;
                dispstep = dispersion/(dispsteps-1);
            } else {
                /* everything specified */
            }
        }
    }


    if(detector_thicksteps <= 0){
        /* auto-select number of steps */
        if(detector_thick < 0.0) {
            /* auto-select range */
            if(detector_thickstep <= 0.0) {
                /* user doesn't care about anything */
                detector_thicksteps = 1;
                detector_thick = 0.0;
                detector_thickstep = 0.0;
            } else {
                /* user specified stepsize and nothing else */
                detector_thick = detector_thickstep;
                detector_thicksteps = 2;
            }
        } else {
            /* user-speficied range */
            if(detector_thickstep <= 0.0) {
                /* range specified, but nothing else */
                detector_thicksteps = 2;
                detector_thickstep = detector_thick/detector_thicksteps;
            } else {
                /* range and step specified, but not number of steps */
                detector_thicksteps = ceil(detector_thick/detector_thickstep);
            }
        }
    } else {
        /* user-specified number of steps */
        if(detector_thick < 0.0) {
            /* auto-select range */
            if(detector_thickstep <= 0.0) {
                /* user cares only about number of steps */
                detector_thick = 0.5e-6;
                detector_thickstep = detector_thick/detector_thicksteps;
            } else {
                /* user doesn't care about range */
                detector_thick = detector_thickstep;
                detector_thicksteps = 2;
            }
        } else {
            /* user-speficied range */
            if(detector_thickstep <= 0.0) {
                /* range and steps specified */
                if(detector_thicksteps <=1 ) detector_thicksteps = 2;
                detector_thickstep = detector_thick/(detector_thicksteps-1);
            } else {
                /* everything specified */
            }
        }
    }

    if(mosaic_domains <= 0){
        /* auto-select number of domains */
        if(mosaic_spread < 0.0) {
            /* user doesn't care about anything */
            mosaic_domains = 1;
            mosaic_spread = 0.0;
        } else {
            /* user-speficied mosaicity, but not number of domains */
            if(mosaic_spread == 0.0)
            {
                mosaic_domains = 1;
            }
            else
            {
                printf("WARNING: finite mosaicity with only one domain! upping to 10 mosaic domains\n");
            mosaic_domains = 10;
        }
        }
    } else {
        /* user-specified number of domains */
        if(mosaic_spread < 0.0) {
            /* number of domains specified, but no spread? */
            printf("WARNING: no mosaic spread specified.  setting mosaic_domains = 1\n");
            mosaic_spread = 0.0;
            mosaic_domains = 1;
        } else {
            /* user-speficied mosaicity and number of domains */
            if(mosaic_spread == 0.0)
            {
                printf("WARNING: zero mosaic spread specified.  setting mosaic_domains = 1\n");
                mosaic_domains = 1;
            }
        }
    }


    /* sanity checks */
    if(hdivrange <= 0.0 || hdivstep <= 0.0 || hdivsteps <= 0) {
        hdivsteps = 1;
        hdivrange = 0.0;
        hdivstep = 0.0;
    }
    if(vdivrange <= 0.0 || vdivstep <= 0.0 || vdivsteps <= 0) {
        vdivsteps = 1;
        vdivrange = 0.0;
        vdivstep = 0.0;
    }
    if(dispersion <= 0.0 || dispstep <= 0.0 || dispsteps <= 0) {
        dispsteps = 1;
        dispersion = 0.0;
        dispstep = 0.0;
    }
    if(detector_thick <= 0.0 || detector_thickstep <= 0.0 || detector_thicksteps <= 0) {
        detector_thicksteps = 1;
        detector_thick = 0.0;
        detector_thickstep = 0.0;
    }


    /* initialize detector origin from a beam center and distance */
    /* there are two conventions here: mosflm and XDS */

    if(beam_convention == ADXV) printf("adxv");
    if(beam_convention == MOSFLM) printf("mosflm");
    if(beam_convention == XDS) printf("xds");
    if(beam_convention == DIALS) printf("dials");
    if(beam_convention == DENZO) printf("denzo");
    if(beam_convention == CUSTOM) printf("custom");
    printf(" convention selected.\n");

    /* first off, what is the relationship between the two "beam centers"? */
    rotate(odet_vector,vector,detector_rotx,detector_roty,detector_rotz);
    ratio = dot_product(beam_vector,vector);
    if(ratio == 0.0) { ratio = DBL_MIN; }
    if(isnan(close_distance)) close_distance = fabs(ratio*distance);
    distance = close_distance/ratio;

    if(detector_pivot == SAMPLE){
        printf("pivoting detector around sample\n");
        /* initialize detector origin before rotating detector */
        pix0_vector[1] = -Fclose*fdet_vector[1]-Sclose*sdet_vector[1]+close_distance*odet_vector[1];
        pix0_vector[2] = -Fclose*fdet_vector[2]-Sclose*sdet_vector[2]+close_distance*odet_vector[2];
        pix0_vector[3] = -Fclose*fdet_vector[3]-Sclose*sdet_vector[3]+close_distance*odet_vector[3];

        /* now swing the detector origin around */
        rotate(pix0_vector,pix0_vector,detector_rotx,detector_roty,detector_rotz);
        rotate_axis(pix0_vector,pix0_vector,twotheta_axis,detector_twotheta);
    }
    /* now orient the detector plane */
    rotate(fdet_vector,fdet_vector,detector_rotx,detector_roty,detector_rotz);
    rotate(sdet_vector,sdet_vector,detector_rotx,detector_roty,detector_rotz);
    rotate(odet_vector,odet_vector,detector_rotx,detector_roty,detector_rotz);

    /* also apply orientation part of twotheta swing */
    rotate_axis(fdet_vector,fdet_vector,twotheta_axis,detector_twotheta);
    rotate_axis(sdet_vector,sdet_vector,twotheta_axis,detector_twotheta);
    rotate_axis(odet_vector,odet_vector,twotheta_axis,detector_twotheta);

    /* make sure beam center is preserved */
    if(detector_pivot == BEAM){
        printf("pivoting detector around direct beam spot\n");
        pix0_vector[1] = -Fbeam*fdet_vector[1]-Sbeam*sdet_vector[1]+distance*beam_vector[1];
        pix0_vector[2] = -Fbeam*fdet_vector[2]-Sbeam*sdet_vector[2]+distance*beam_vector[2];
        pix0_vector[3] = -Fbeam*fdet_vector[3]-Sbeam*sdet_vector[3]+distance*beam_vector[3];
    }

    /* what is the point of closest approach between sample and detector? */
    Fclose         = -dot_product(pix0_vector,fdet_vector);
    Sclose         = -dot_product(pix0_vector,sdet_vector);
    close_distance =  dot_product(pix0_vector,odet_vector);

    /* where is the direct beam now? */
    /* difference between beam impact vector and detector origin */
    newvector[1] = close_distance/ratio*beam_vector[1]-pix0_vector[1];
    newvector[2] = close_distance/ratio*beam_vector[2]-pix0_vector[2];
    newvector[3] = close_distance/ratio*beam_vector[3]-pix0_vector[3];
    /* extract components along detector vectors */
    Fbeam = dot_product(fdet_vector,newvector);
    Sbeam = dot_product(sdet_vector,newvector);
    distance = close_distance/ratio;

    /* find origin in XDS convention */
    ORGX=Fclose/pixel_size+0.5;
    ORGY=Sclose/pixel_size+0.5;

    /* find origin in DIALS convention */
    newvector[1]=+0;newvector[2]=+0;newvector[3]=+1;
    dials_origin[1] = 1000.0*dot_product(pix0_vector,newvector);
    newvector[1]=+0;newvector[2]=+1;newvector[3]=+0;
    dials_origin[2] = 1000.0*dot_product(pix0_vector,newvector);
    newvector[1]=-1;newvector[2]=+0;newvector[3]=+0;
    dials_origin[3] = 1000.0*dot_product(pix0_vector,newvector);

    /* find the beam in the detector frame */
    newvector[1] = dot_product(beam_vector,fdet_vector);
    newvector[2] = dot_product(beam_vector,sdet_vector);
    newvector[3] = dot_product(beam_vector,odet_vector);
    printf("XDS incident beam: %g %g %g\n",newvector[1],newvector[2],newvector[3]);

    if(interpolate > 1){
        /* no user options */
        if(( Na <= 2) || (Nb <= 2) || (Nc <= 2)){
            printf("auto-selected tricubic interpolation of structure factors\n");
            interpolate = 1;
        }
        else
        {
            printf("auto-selected no interpolation\n");
            interpolate = 0;
        }
    }


    /* user-specified unit cell */
    if(user_cell)
    {
        /* a few random defaults */
        if(b[0]  <= 0.0) b[0] = a[0];
        if(c[0]  <= 0.0) c[0] = a[0];
        if(alpha <= 0.0) alpha = M_PI/2;
        if(beta  <= 0.0) beta  = M_PI/2;
        if(gamma <= 0.0) gamma = M_PI/2;

        /* get cell volume from angles */
        aavg = (alpha+beta+gamma)/2;
        skew = sin(aavg)*sin(aavg-alpha)*sin(aavg-beta)*sin(aavg-gamma);
        if(skew<0.0) skew=-skew;
        V_cell = 2.0*a[0]*b[0]*c[0]*sqrt(skew);
        if(V_cell <= 0.0)
        {
            printf("WARNING: impossible unit cell volume: %g\n",V_cell);
            V_cell = DBL_MIN;
        }
        V_star = 1.0/V_cell;

        /* now get reciprocal-cell lengths from the angles and volume */
        a_star[0] = b[0]*c[0]*sin(alpha)*V_star;
        b_star[0] = c[0]*a[0]*sin(beta)*V_star;
        c_star[0] = a[0]*b[0]*sin(gamma)*V_star;
        if(a_star[0] <= 0.0 || b_star[0] <= 0.0 || c_star[0] <= 0.0)
        {
            printf("WARNING: impossible reciprocal cell lengths: %g %g %g\n",
                a_star[0],b_star[0],c_star[0]);
            a_star[0] = fabs(a_star[0]);
            b_star[0] = fabs(b_star[0]);
            c_star[0] = fabs(c_star[0]);
            if(a_star[0] <= 0.0) a_star[0] = DBL_MIN;
            if(b_star[0] <= 0.0) b_star[0] = DBL_MIN;
            if(c_star[0] <= 0.0) c_star[0] = DBL_MIN;
        }

        /* for fun, compute the reciprocal-cell angles from direct-cell angles */
        sin_alpha_star = a[0]*V_star/b_star[0]/c_star[0];
        sin_beta_star  = b[0]*V_star/a_star[0]/c_star[0];
        sin_gamma_star = c[0]*V_star/a_star[0]/b_star[0];
        cos_alpha_star = (cos(beta)*cos(gamma)-cos(alpha))/(sin(beta)*sin(gamma));
        cos_beta_star  = (cos(gamma)*cos(alpha)-cos(beta))/(sin(gamma)*sin(alpha));
        cos_gamma_star = (cos(alpha)*cos(beta)-cos(gamma))/(sin(alpha)*sin(beta));
        if(sin_alpha_star>1.0000001 || sin_alpha_star<-1.0000001 ||
           sin_beta_star >1.0000001 || sin_beta_star <-1.0000001 ||
           sin_gamma_star>1.0000001 || sin_gamma_star<-1.0000001 ||
           cos_alpha_star>1.0000001 || cos_alpha_star<-1.0000001 ||
           cos_beta_star >1.0000001 || cos_beta_star <-1.0000001 ||
           cos_gamma_star>1.0000001 || cos_gamma_star<-1.0000001 )
        {
            printf("WARNING: oddball reciprocal cell angles:\n");
            printf("sin(alpha_star) = %.25g\n",sin_alpha_star);
            printf("cos(alpha_star) = %.25g\n",cos_alpha_star);
            printf("sin(beta_star)  = %.25g\n",sin_beta_star);
            printf("cos(beta_star)  = %.25g\n",cos_beta_star);
            printf("sin(gamma_star) = %.25g\n",sin_gamma_star);
            printf("cos9gamma_star) = %.25g\n",cos_gamma_star);
        }
        if(sin_alpha_star>1.0) sin_alpha_star=1.0;
        if(sin_beta_star >1.0) sin_beta_star =1.0;
        if(sin_gamma_star>1.0) sin_gamma_star=1.0;
        if(sin_alpha_star<-1.0) sin_alpha_star=-1.0;
        if(sin_beta_star <-1.0) sin_beta_star =-1.0;
        if(sin_gamma_star<-1.0) sin_gamma_star=-1.0;
        if(cos_alpha_star*cos_alpha_star>1.0) cos_alpha_star=1.0;
        if(cos_beta_star *cos_beta_star >1.0) cos_beta_star=1.0;
        if(cos_gamma_star*cos_gamma_star>1.0) cos_gamma_star=1.0;
        alpha_star = atan2(sin_alpha_star,cos_alpha_star);
        beta_star  = atan2(sin_beta_star ,cos_beta_star );
        gamma_star = atan2(sin_gamma_star,cos_gamma_star);


        /* construct default orientation */
        a_star[1] = a_star[0];
        b_star[1] = b_star[0]*cos_gamma_star;
        c_star[1] = c_star[0]*cos_beta_star;
        a_star[2] = 0.0;
        b_star[2] = b_star[0]*sin_gamma_star;
        c_star[2] = c_star[0]*(cos_alpha_star-cos_beta_star*cos_gamma_star)/sin_gamma_star;
        a_star[3] = 0.0;
        b_star[3] = 0.0;
        c_star[3] = c_star[0]*V_cell/(a[0]*b[0]*c[0]*sin_gamma_star);
    }

    /* load the lattice orientation (reciprocal cell vectors) from a mosflm matrix */
    if(matfilename != NULL)
    {
        infile = fopen(matfilename,"r");
        if(infile != NULL)
        {
            printf("reading %s\n",matfilename);
            if(! fscanf(infile,"%lg%lg%lg",a_star+1,b_star+1,c_star+1)) {perror("fscanf");};
            if(! fscanf(infile,"%lg%lg%lg",a_star+2,b_star+2,c_star+2)) {perror("fscanf");};
            if(! fscanf(infile,"%lg%lg%lg",a_star+3,b_star+3,c_star+3)) {perror("fscanf");};
            fclose(infile);

            /* mosflm A matrix includes the wavelength, so remove it */
            /* calculate reciprocal cell lengths, store in 0th element */
            vector_scale(a_star,a_star,1e-10/lambda0);
            vector_scale(b_star,b_star,1e-10/lambda0);
            vector_scale(c_star,c_star,1e-10/lambda0);
        }
    }

    /* check for flag to generate random missetting angle */
    if(misset[0] == -1.0)
    {
        /* use spherical cap as sphere to generate random orientation in umat */
        mosaic_rotation_umat(90.0, umat, &seed);
        /* get the missetting angles, in case we want to use them again on -misset option */
        umat2misset(umat,misset);
        printf("random orientation misset angles: %f %f %f deg\n",misset[1]*RTD,misset[2]*RTD,misset[3]*RTD);
        /* apply this orientation shift */
        //rotate_umat(a_star,a_star,umat);
        //rotate_umat(b_star,b_star,umat);
        //rotate_umat(c_star,c_star,umat);
        /* do not apply again */
        misset[0] = 1.0;
    }

    /* apply any missetting angle, if not already done */
    if(misset[0] > 0.0)
    {
        rotate(a_star,a_star,misset[1],misset[2],misset[3]);
        rotate(b_star,b_star,misset[1],misset[2],misset[3]);
        rotate(c_star,c_star,misset[1],misset[2],misset[3]);
    }

    /* various cross products */
    cross_product(a_star,b_star,a_star_cross_b_star);
    cross_product(b_star,c_star,b_star_cross_c_star);
    cross_product(c_star,a_star,c_star_cross_a_star);

    /* reciprocal lattice vector "a_star" is defined as perpendicular to both b and c, and must also preserve volume
       converse is true for direct-space lattice: a is perpendicular to both b_star and c_star
       a = ( b_star cross c_star ) / V_star    */

    /* reciprocal unit cell volume, but is it lambda-corrected? */
    V_star = dot_product(a_star,b_star_cross_c_star);

    /* make sure any user-supplied cell takes */
    if(user_cell)
    {
        /* a,b,c and V_cell were generated above */

        /* force the cross-product vectors to have proper magnitude: b_star X c_star = a*V_star */
        vector_rescale(b_star_cross_c_star,b_star_cross_c_star,a[0]/V_cell);
        vector_rescale(c_star_cross_a_star,c_star_cross_a_star,b[0]/V_cell);
        vector_rescale(a_star_cross_b_star,a_star_cross_b_star,c[0]/V_cell);
        V_star = 1.0/V_cell;
    }

    /* direct-space cell volume */
    V_cell = 1.0/V_star;

    /* generate direct-space cell vectors, also updates magnitudes */
    vector_scale(b_star_cross_c_star,a,V_cell);
    vector_scale(c_star_cross_a_star,b,V_cell);
    vector_scale(a_star_cross_b_star,c,V_cell);

    /* now that we have direct-space vectors, re-generate the reciprocal ones */
    cross_product(a,b,a_cross_b);
    cross_product(b,c,b_cross_c);
    cross_product(c,a,c_cross_a);
    vector_scale(b_cross_c,a_star,V_star);
    vector_scale(c_cross_a,b_star,V_star);
    vector_scale(a_cross_b,c_star,V_star);

    /* for fun, calculate the cell angles too */
    sin_alpha = a_star[0]*V_cell/b[0]/c[0];
    sin_beta  = b_star[0]*V_cell/a[0]/c[0];
    sin_gamma = c_star[0]*V_cell/a[0]/b[0];
    cos_alpha = dot_product(b,c)/b[0]/c[0];
    cos_beta  = dot_product(a,c)/a[0]/c[0];
    cos_gamma = dot_product(a,b)/a[0]/b[0];
    if(sin_alpha>1.0000001 || sin_alpha<-1.0000001 ||
       sin_beta >1.0000001 || sin_beta <-1.0000001 ||
       sin_gamma>1.0000001 || sin_gamma<-1.0000001 ||
       cos_alpha>1.0000001 || cos_alpha<-1.0000001 ||
       cos_beta >1.0000001 || cos_beta <-1.0000001 ||
       cos_gamma>1.0000001 || cos_gamma<-1.0000001 )
    {
        printf("WARNING: oddball cell angles:\n");
            printf("sin_alpha = %.25g\n",sin_alpha);
            printf("cos_alpha = %.25g\n",cos_alpha);
            printf("sin_beta  = %.25g\n",sin_beta);
            printf("cos_beta  = %.25g\n",cos_beta);
            printf("sin_gamma = %.25g\n",sin_gamma);
            printf("cos_gamma = %.25g\n",cos_gamma);
    }
    if(sin_alpha>1.0) sin_alpha=1.0;
    if(sin_beta >1.0) sin_beta =1.0;
    if(sin_gamma>1.0) sin_gamma=1.0;
    if(sin_alpha<-1.0) sin_alpha=-1.0;
    if(sin_beta <-1.0) sin_beta =-1.0;
    if(sin_gamma<-1.0) sin_gamma=-1.0;
    if(cos_alpha*cos_alpha>1.0) cos_alpha=1.0;
    if(cos_beta *cos_beta >1.0) cos_beta=1.0;
    if(cos_gamma*cos_gamma>1.0) cos_gamma=1.0;
    alpha = atan2(sin_alpha,cos_alpha);
    beta  = atan2(sin_beta ,cos_beta );
    gamma = atan2(sin_gamma,cos_gamma);


    /* reciprocal cell angles */
    sin_alpha_star = a[0]*V_star/b_star[0]/c_star[0];
    sin_beta_star  = b[0]*V_star/a_star[0]/c_star[0];
    sin_gamma_star = c[0]*V_star/a_star[0]/b_star[0];
    cos_alpha_star = dot_product(b_star,c_star)/b_star[0]/c_star[0];
    cos_beta_star  = dot_product(a_star,c_star)/a_star[0]/c_star[0];
    cos_gamma_star = dot_product(a_star,b_star)/a_star[0]/b_star[0];
    if(sin_alpha_star>1.0000001 || sin_alpha_star<-1.0000001 ||
       sin_beta_star >1.0000001 || sin_beta_star <-1.0000001 ||
       sin_gamma_star>1.0000001 || sin_gamma_star<-1.0000001 ||
       cos_alpha_star>1.0000001 || cos_alpha_star<-1.0000001 ||
       cos_beta_star >1.0000001 || cos_beta_star <-1.0000001 ||
       cos_gamma_star>1.0000001 || cos_gamma_star<-1.0000001 )
    {
            printf("WARNING: oddball reciprocal cell angles:\n");
            printf("sin(alpha_star) = %.25g\n",sin_alpha_star);
            printf("cos(alpha_star) = %.25g\n",cos_alpha_star);
            printf("sin(beta_star)  = %.25g\n",sin_beta_star);
            printf("cos(beta_star)  = %.25g\n",cos_beta_star);
            printf("sin(gamma_star) = %.25g\n",sin_gamma_star);
            printf("cos(gamma_star) = %.25g\n",cos_gamma_star);
    }
    if(sin_alpha_star>1.0) sin_alpha_star=1.0;
    if(sin_beta_star >1.0) sin_beta_star =1.0;
    if(sin_gamma_star>1.0) sin_gamma_star=1.0;
    if(sin_alpha_star<-1.0) sin_alpha_star=-1.0;
    if(sin_beta_star <-1.0) sin_beta_star =-1.0;
    if(sin_gamma_star<-1.0) sin_gamma_star=-1.0;
    if(cos_alpha_star*cos_alpha_star>1.0) cos_alpha_star=1.0;
    if(cos_beta_star *cos_beta_star >1.0) cos_beta_star=1.0;
    if(cos_gamma_star*cos_gamma_star>1.0) cos_gamma_star=1.0;
    alpha_star = atan2(sin_alpha_star,cos_alpha_star);
    beta_star  = atan2(sin_beta_star ,cos_beta_star );
    gamma_star = atan2(sin_gamma_star,cos_gamma_star);

    printf("Unit Cell: %g %g %g %g %g %g\n", a[0],b[0],c[0],alpha*RTD,beta*RTD,gamma*RTD);
    printf("Recp Cell: %g %g %g %g %g %g\n", a_star[0],b_star[0],c_star[0],alpha_star*RTD,beta_star*RTD,gamma_star*RTD);
    printf("volume = %g A^3\n",V_cell);

    /* print out the real-space matrix */
    printf("real-space cell vectors (Angstrom):\n");
    printf("     %-10s  %-10s  %-10s\n","a","b","c");
    printf("X: %11.8f %11.8f %11.8f\n",a[1],b[1],c[1]);
    printf("Y: %11.8f %11.8f %11.8f\n",a[2],b[2],c[2]);
    printf("Z: %11.8f %11.8f %11.8f\n",a[3],b[3],c[3]);
    printf("reciprocal-space cell vectors (Angstrom^-1):\n");
    printf("     %-10s  %-10s  %-10s\n","a_star","b_star","c_star");
    printf("X: %11.8f %11.8f %11.8f\n",a_star[1],b_star[1],c_star[1]);
    printf("Y: %11.8f %11.8f %11.8f\n",a_star[2],b_star[2],c_star[2]);
    printf("Z: %11.8f %11.8f %11.8f\n",a_star[3],b_star[3],c_star[3]);

    /* now convert these to meters */
    vector_scale(a,a,1e-10);
    vector_scale(b,b,1e-10);
    vector_scale(c,c,1e-10);

    /* define phi=0 mosaic=0 crystal orientation */
    vector_scale(a,a0,1.0);
    vector_scale(b,b0,1.0);
    vector_scale(c,c0,1.0);

    /* define phi=0 crystal orientation */
    vector_scale(a,ap,1.0);
    vector_scale(b,bp,1.0);
    vector_scale(c,cp,1.0);

    /* now we know the cell, calculate crystal size in meters */
    if(sample_x > 0) Na = ceil(sample_x/a[0]);
    if(sample_y > 0) Nb = ceil(sample_y/b[0]);
    if(sample_z > 0) Nc = ceil(sample_z/c[0]);
    if(Na <= 1.0) Na = 1.0;
    if(Nb <= 1.0) Nb = 1.0;
    if(Nc <= 1.0) Nc = 1.0;
    xtalsize_a = a[0]*Na;
    xtalsize_b = b[0]*Nb;
    xtalsize_c = c[0]*Nc;
    printf("crystal is %g x %g x %g microns\n",xtalsize_a*1e6,xtalsize_b*1e6,xtalsize_c*1e6);
    xtalsize_max = xtalsize_a;
    if(xtalsize_max < xtalsize_b) xtalsize_max = xtalsize_b;
    if(xtalsize_max < xtalsize_c) xtalsize_max = xtalsize_c;
    reciprocal_pixel_size = lambda0*distance/pixel_size;
    recommended_oversample = ceil(3.0 * xtalsize_max/reciprocal_pixel_size);
    if(recommended_oversample <= 0) recommended_oversample = 1;
    if(oversample <= 0) {
        oversample = recommended_oversample;
        printf("auto-selected %d-fold oversampling\n",oversample);
    }
    if(oversample < recommended_oversample)
    {
        printf("WARNING: maximum dimension of sample is %g A\n",xtalsize_max*1e10);
        printf("         but reciprocal pixel size is %g A\n", reciprocal_pixel_size*1e10 );
        printf("         intensity may vary significantly across a pixel!\n");
        printf("         recommend -oversample %d to work around this\n",recommended_oversample);
    }

    /* rough estimate of sample properties */
    sample_x = xtalsize_a;
    sample_y = xtalsize_b;
    sample_z = xtalsize_c;
    volume = sample_x*sample_y*sample_z;
    density = 1.2e6;
    molecules = Na*Nb*Nc;
    molecular_weight = volume*density*Avogadro/molecules;
    printf("approximate MW = %g\n",molecular_weight);

    /* load the structure factors */
    if(hklfilename == NULL)
    {
        /* try to recover Fs from a previous run */
        if(Fdumpfile != NULL)
        {
            printf("reading Fs from %s\n",dumpfilename);
//          n=0;
              if(! fscanf(Fdumpfile,"%d%d%d%d%d%d\n\f",&h_min,&h_max,&k_min,&k_max,&l_min,&l_max) ) {perror("fscanf");};
            h_range = h_max - h_min + 1;
            k_range = k_max - k_min + 1;
            l_range = l_max - l_min + 1;
            Fhkl = (double***) calloc(h_range+1,sizeof(double**));
            for (h0=0; h0<=h_range;h0++) {
                *(Fhkl +h0) = (double**) calloc(k_range+1,sizeof(double*));
                for (k0=0; k0<=k_range;k0++) {
                    *(*(Fhkl +h0)+k0) = (double*) calloc(l_range+1,sizeof(double));
                    if(! fread(*(*(Fhkl +h0)+k0),sizeof(double),l_range+1,Fdumpfile) )
                    {
                        perror("fscanf");
                    };
//                  n+=l_range;
                }
            }
            fclose(Fdumpfile);
            hkls = h_range*k_range*l_range;
        }
        else
        {
            /* no hkl file and no dumpfile */
            if(default_F == 0.0)
            {
                printf("ERROR: no hkl file and no dump file to read.");
                exit(9);
            }
        }
    }
    else
    {
        infile = fopen(hklfilename,"r");
        if(infile == NULL)
        {
            printf("ERROR: unable to open %s.",hklfilename);
            exit(9);
        }
        hkls = 0;
        h_min=k_min=l_min=1e9;
        h_max=k_max=l_max=-1e9;
        printf("counting entries in %s\n",hklfilename);
        while(4 == fscanf(infile,"%lg%lg%lg%lg",&h,&k,&l,&F_cell)){
            if(h != ceil(h-0.4)) printf("WARNING: non-integer value for h (%g) at line %d\n",h,hkls);
            if(k != ceil(k-0.4)) printf("WARNING: non-integer value for k (%g) at line %d\n",k,hkls);
            if(l != ceil(l-0.4)) printf("WARNING: non-integer value for l (%g) at line %d\n",l,hkls);
            if(h_min > h) h_min = h;
            if(k_min > k) k_min = k;
            if(l_min > l) l_min = l;
            if(h_max < h) h_max = h;
            if(k_max < k) k_max = k;
            if(l_max < l) l_max = l;
            ++hkls;
        }
        rewind(infile);
        h_range = h_max - h_min + 1;
        k_range = k_max - k_min + 1;
        l_range = l_max - l_min + 1;

        if(h_range < 0 || k_range < 0 || l_range < 0) {
            printf("h: %d - %d\n",h_min,h_max);
            printf("k: %d - %d\n",k_min,k_max);
            printf("l: %d - %d\n",l_min,l_max);
            printf("ERROR: not enough HKL indices in %s\n",hklfilename);
            exit(9);
        }

        /* allocate memory for 3d arrays */
        //printf("allocating %d %d-byte double**\n",h_range+1,sizeof(double**));
        Fhkl = (double***) calloc(h_range+1,sizeof(double**));
        if(Fhkl==NULL){perror("ERROR");exit(9);};
        for (h0=0; h0<=h_range;h0++) {
                //printf("allocating %d %d-byte double*\n",k_range+1,sizeof(double*));
                Fhkl[h0] = (double**) calloc(k_range+1,sizeof(double*));
                if(Fhkl[h0]==NULL){perror("ERROR");exit(9);};
                for (k0=0; k0<=k_range;k0++) {
                        //printf("allocating %d %d-byte double\n",k_range+1,sizeof(double));
                        Fhkl[h0][k0] = (double*) calloc(l_range+1,sizeof(double));
                        if(Fhkl[h0][k0]==NULL){perror("ERROR");exit(9);};
                }
        }
        if(default_F != 0.0) {
            printf("initializing to default_F = %g:\n",default_F);
            for (h0=0; h0<h_range;h0++) {
                for (k0=0; k0<k_range;k0++) {
                    for (l0=0; l0<l_range;l0++) {
                        Fhkl[h0][k0][l0] = default_F;
                    }
                }
            }
            printf("done initializing:\n");
        }


        printf("re-reading %s\n",hklfilename);
        while(4 == fscanf(infile,"%d%d%d%lg",&h0,&k0,&l0,&F_cell)){
            Fhkl[h0-h_min][k0-k_min][l0-l_min]=F_cell;
        }
        fclose(infile);

//      for(h0=h_min;h0<=h_max;++h0){
//          for(k0=k_min;k0<=k_max;++k0){
//              for(l0=l_min;l0<=l_max;++l0){
//                  if ( (h0<=h_max) && (h0>=h_min) && (k0<=k_max) && (k0>=k_min) && (l0<=l_max) && (l0>=l_min)  ) {
//                      /* just take nearest-neighbor */
//                      F_cell = Fhkl[h0-h_min][k0-k_min][l0-l_min];
//                  }
//                  else
//                  {
//                      F_cell = 0.0;
//                  }
//                  printf("%d %d %d = %f\n",h0,k0,l0,F_cell);
//              }
//          }
//      }

        /* make dump file */
        outfile = fopen(dumpfilename,"wb");
        if(outfile == NULL)
        {
            printf("WARNING: unable to open dump file: %s\n",dumpfilename);
        }
        else
        {
            printf("writing dump file for next time: %s\n",dumpfilename);
            fprintf(outfile,"%d %d %d %d %d %d\n\f",h_min,h_max,k_min,k_max,l_min,l_max);
            for (h0=0; h0<=h_range;h0++) {
                for (k0=0; k0<=k_range;k0++) {
                        fwrite(*(*(Fhkl +h0)+k0),sizeof(double),l_range+1,outfile);
                }
            }
            fclose(outfile);
        }
    }

    /* no point in interpolating if nothing to interpolate */
    if(hkls == 0) interpolate = 0;

    if(interpolate){
        /* allocate interpolation array */
        sub_Fhkl = (double***) calloc(6,sizeof(double**));
        for (h0=0; h0<=5;h0++) {
            *(sub_Fhkl +h0) = (double**) calloc(6,sizeof(double*));
            for (k0=0; k0<=5;k0++) {
                *(*(sub_Fhkl +h0)+k0) = (double*) calloc(6,sizeof(double));
            }
        }
    }


    /* now read in amorphous material structure factors */
    stols = 0;
    if(stolfilename != NULL)
    {
        printf("reading %s\n",stolfilename);
        stols = read_text_file(stolfilename,2,&stol_of,&F_of);
        if(stols == 0){
            perror("no data in input file");
            exit(9);
        }
    }

    if(stols == 0 && water_size != 0.0)
    {
        /* do something clever here */
    }

    if(stols > 0)
    {
        /* add two values at either end for interpolation */
        stols += 4;
        F_highangle = NAN;
        for(i=stols-3;i>1;--i){
            stol_of[i] = stol_of[i-2] * stol_file_mult;
            F_of[i]    = F_of[i-2];
            if(! isnan(F_of[i])) {
                F_lowangle = F_of[i];
                if(isnan(F_highangle)) {
                    F_highangle = F_of[i];
                }
            }
            else
            {
                /* missing values are zero */
                F_of[i] = 0.0;
            }
        }
        stol_of[0] = -1e99;
        stol_of[1] = -1e98;
        F_of[0] = F_of[1] = F_lowangle;
        stol_of[stols-2] = 1e98;
        stol_of[stols-1] = 1e99;
        F_of[stols-1] = F_of[stols-2] = F_highangle;
    }

    /* print out detector sensor thickness with sweep over all sensor layers */
    for(thick_tic=0;thick_tic<detector_thicksteps;++thick_tic){
        printf("thick%d = %g um\n",thick_tic,detector_thickstep*thick_tic*1e6);
    }

    /* show phi steps with sweep over spindle axis */
    for(phi_tic = 0; phi_tic < phisteps; ++phi_tic){
        phi = phi0 + phistep*phi_tic;
        printf("phi%d = %g\n",phi_tic,phi*RTD);
    }




    /* import sources from user file */
    sources = 0;
    if(sourcefilename != NULL) {
        sources = read_text_file(sourcefilename,5,&source_X,&source_Y,&source_Z,&source_I,&source_lambda);
        if(sources == 0) {
            perror("reading source definition file");
            exit(9);
        }
        /* apply defaults to missing values */
        for(source=0;source<sources;++source){
            if(isnan(source_X[source])) {
                source_X[source] = -source_distance*beam_vector[1];
            }
            if(isnan(source_Y[source])) {
                source_Y[source] = -source_distance*beam_vector[2];
            }
            if(isnan(source_Z[source])) {
                source_Z[source] = -source_distance*beam_vector[3];
            }
            if(isnan(source_I[source])) {
                source_I[source] = 1.0;
            }
            if(isnan(source_lambda[source])) {
                source_lambda[source] = lambda0;
            }
        }
    }


    if(sources == 0)
    {
        /* generate generic list of sources */

        /* count divsteps sweep over solid angle of beam divergence */
        divsteps = 0;
        for(hdiv_tic=0;hdiv_tic<hdivsteps;++hdiv_tic){
            for(vdiv_tic=0;vdiv_tic<vdivsteps;++vdiv_tic){
                hdiv = hdivstep * hdiv_tic - hdivrange/2.0 ;
                vdiv = vdivstep * vdiv_tic - vdivrange/2.0 ;
                /* force an elliptical divergence */
                test = (hdiv*hdiv-hdivstep*hdivstep/4.0*(1-hdivsteps%2))/hdivrange/hdivrange ;
                test += (vdiv*vdiv-vdivstep*vdivstep/4.0*(1-vdivsteps%2))/vdivrange/vdivrange ;
                if( round_div && test*4.0 > 1.1) continue;

                ++divsteps;
                printf("divergence deviation: %g %g\n",hdiv,vdiv);
            }
        }

        /* print out wavelength steps with sweep over spectral dispersion */
        for(disp_tic=0;disp_tic<dispsteps;++disp_tic){
            lambda = lambda0 * ( 1.0 + dispstep * disp_tic - dispersion/2.0 ) ;
            printf("lambda%d = %.15g\n",disp_tic,lambda);
        }

        /* allocate enough space */
        sources = divsteps*dispsteps;
        source_X = (double *) calloc(sources+10,sizeof(double));
        source_Y = (double *) calloc(sources+10,sizeof(double));
        source_Z = (double *) calloc(sources+10,sizeof(double));
        source_I = (double *) calloc(sources+10,sizeof(double));
        source_lambda = (double *) calloc(sources+10,sizeof(double));

        /* now actually create the source entries */
        weight = 1.0/sources;
        sources = 0;
        for(hdiv_tic=0;hdiv_tic<hdivsteps;++hdiv_tic){
            for(vdiv_tic=0;vdiv_tic<vdivsteps;++vdiv_tic){
                hdiv = hdivstep * hdiv_tic - hdivrange/2.0 ;
                vdiv = vdivstep * vdiv_tic - vdivrange/2.0 ;
                /* force an elliptical divergence */
                test = (hdiv*hdiv-hdivstep*hdivstep/4.0*(1-hdivsteps%2))/hdivrange/hdivrange ;
                test += (vdiv*vdiv-vdivstep*vdivstep/4.0*(1-vdivsteps%2))/vdivrange/vdivrange ;
                if( round_div && test*4.0 > 1.1) continue;

                /* construct unit vector along "beam" */
                vector[1] = -source_distance*beam_vector[1];
                vector[2] = -source_distance*beam_vector[2];
                vector[3] = -source_distance*beam_vector[3];
                /* divergence is in angle space */
                /* define "horizontal" as the E-vector of the incident beam */
                rotate_axis(vector,newvector,polar_vector,vdiv);
                rotate_axis(newvector,vector,vert_vector,hdiv);

                /* one source at each position for each wavelength */
                for(disp_tic=0;disp_tic<dispsteps;++disp_tic){
                    lambda = lambda0 * ( 1.0 + dispstep * disp_tic - dispersion/2.0 ) ;

                    source_X[sources] = vector[1];
                    source_Y[sources] = vector[2];
                    source_Z[sources] = vector[3];
                    source_I[sources] = weight;
                    source_lambda[sources] = lambda;
                    ++sources;
                }
            }
        }
    }
    printf("  created a total of %d sources:\n",sources);
    for(source=0;source<sources;++source){

        /* retrieve stuff from cache */
        X = source_X[source];
        Y = source_Y[source];
        Z = source_Z[source];
        I = source_I[source];
        lambda = source_lambda[source];

        printf("%g %g %g   %g %g\n",X,Y,Z,I,lambda);
    }


    /* allocate enough space */
    mosaic_umats = (double *) calloc(mosaic_domains+10,9*sizeof(double));

    /* now actually create the orientation of each domain */
    for(mos_tic=0;mos_tic<mosaic_domains;++mos_tic){
        mosaic_rotation_umat(mosaic_spread, mosaic_umats+9*mos_tic, &mosaic_seed);
        if(mos_tic==0)
        {
            /* force at least one domain to be "aligned"? */
            mosaic_umats[0]=1.0;mosaic_umats[1]=0.0;mosaic_umats[2]=0.0;
            mosaic_umats[3]=0.0;mosaic_umats[4]=1.0;mosaic_umats[5]=0.0;
            mosaic_umats[6]=0.0;mosaic_umats[7]=0.0;mosaic_umats[8]=1.0;
        }
//      printf("%d diagonal %f %f %f\n",mos_tic,mosaic_umats[mos_tic*9],mosaic_umats[mos_tic*9+4],mosaic_umats[mos_tic*9+8]);
        printf("%d by: %f deg\n",mos_tic,acos((mosaic_umats[mos_tic*9]+mosaic_umats[mos_tic*9+4]+mosaic_umats[mos_tic*9+8]-1)/2)*RTD);
//      umat2misset(mosaic_umats+9*mos_tic,mosaic_missets);
//      printf("%d by: %f %f %f deg\n",mos_tic,mosaic_missets[1]*RTD,mosaic_missets[2]*RTD,mosaic_missets[3]*RTD);
//      printf("%f %f %f\n",mos_tic,*(mosaic_umats+9*mos_tic+0),*(mosaic_umats+9*mos_tic+1),*(mosaic_umats+9*mos_tic+2));
//      printf("%f %f %f\n",mos_tic,*(mosaic_umats+9*mos_tic+3),*(mosaic_umats+9*mos_tic+4),*(mosaic_umats+9*mos_tic+5));
//      printf("%f %f %f\n",mos_tic,*(mosaic_umats+9*mos_tic+6),*(mosaic_umats+9*mos_tic+7),*(mosaic_umats+9*mos_tic+8));
    }

    printf("  created a total of %d mosaic domains\n",mosaic_domains);

    /* final decisions about sampling */
    if(oversample <= 0) oversample = 1;
    steps = sources*mosaic_domains*phisteps*oversample*oversample;
    subpixel_size = pixel_size/oversample;


    printf("  %d initialized hkls (all others =%g)\n",hkls,default_F);
    printf("  ");
    if(xtal_shape == ROUND)  printf("ellipsoidal");
    if(xtal_shape == SQUARE) printf("parallelpiped");
    if(xtal_shape == GAUSS ) printf("gaussian");
    if(xtal_shape == TOPHAT) printf("tophat-spot");
    printf(" xtal: %.0fx%.0fx%.0f cells\n",Na,Nb,Nc);
    printf("  wave=%g meters +/- %g%% in %d steps\n",lambda0,dispersion*100,dispsteps);
    if(nopolar) { printf("  polarization effect disabled\n"); }
           else { printf("  Kahn polarization factor: %f\n",polarization); }
    if(curved_detector) printf("  curved detector: all pixels same distance from origin\n");
    if(point_pixel) printf("  pixel obliquity effect disabled\n");
    printf("  incident fluence: %lg photons/m^2\n",fluence);
    printf("  distance=%lg detsize=%lgx%lg  pixel=%lg meters (%dx%d pixels)\n",distance,detsize_f,detsize_s,pixel_size,fpixels,spixels);
    printf("  Xbeam=%lg Ybeam=%lg\n",Xbeam,Ybeam);
    printf("  Fbeam=%lg Sbeam=%lg\n",Fbeam,Sbeam);
    printf("  Xclose=%lg Yclose=%lg\n",Xclose,Yclose);
    printf("  Fclose=%lg Sclose=%lg\n",Fclose,Sclose);
    printf("  DIRECTION_OF_DETECTOR_X-AXIS= %g %g %g\n",fdet_vector[1],fdet_vector[2],fdet_vector[3]);
    printf("  DIRECTION_OF_DETECTOR_Y-AXIS= %g %g %g\n",sdet_vector[1],sdet_vector[2],sdet_vector[3]);
    printf("  DIRECTION_OF_DETECTOR_Z-AXIS= %g %g %g\n",odet_vector[1],odet_vector[2],odet_vector[3]);
    printf("  INCIDENT_BEAM_DIRECTION= %g %g %g\n",beam_vector[1],beam_vector[2],beam_vector[3]);
    printf("  spindle ROTATION_AXIS= %g %g %g\n",spindle_vector[1],spindle_vector[2],spindle_vector[3]);
    cross_product(beam_vector,polar_vector,vector);
    printf("  POLARIZATION_PLANE_NORMAL= %g %g %g\n",vector[1],vector[2],vector[3]);
    printf("  dials origin= %g %g %g\n",dials_origin[1],dials_origin[2],dials_origin[3]);
    printf("  roi: %d < x < %d && %d < y < %d\n",roi_xmin,roi_xmax,roi_ymin,roi_ymax);
    printf("  hdivrange=%g hdivstep=%g  radians\n",hdivrange,hdivstep);
    printf("  vdivrange=%g vdivstep=%g  radians\n",vdivrange,vdivstep);
    printf("  %d divergence steps\n",divsteps);
    printf("  %d sources\n",sources);
    printf("  %d mosaic domains over mosaic spread of %g degrees\n",mosaic_domains,mosaic_spread*RTD);
    printf("  %d phi steps from %g to %g degrees\n",phisteps,phi0*RTD,(phi0+osc)*RTD);
    printf("  %dx%d pixel oversample steps\n",oversample,oversample);
    if(maskimage != NULL) printf("  skipping zero-flagged pixels in %s\n",maskfilename);
//    printf("  coherent source: %d\n",coherent);
    if(calculate_noise){
        printf("\n  noise image paramters:\n");
        printf("  seed: %ld\n",seed);
        printf("  water droplet size: %g m\n",water_size);
    }


    /* sweep over detector */
    sum = sumsqr = 0.0;
    j = sumn = 0;
    progress_pixel = 0;
    omega_sum = 0.0;
    for(spixel=0;spixel<spixels;++spixel)
    {
        for(fpixel=0;fpixel<fpixels;++fpixel)
        {

            /* allow for just one part of detector to be rendered */
            if(fpixel < roi_xmin || fpixel > roi_xmax || spixel < roi_ymin || spixel > roi_ymax)
            {
                ++j; continue;
            }
            /* allow for the use of a mask */
            if(maskimage != NULL)
            {
                /* skip any flagged pixels in the mask */
                if(maskimage[j] == 0)
                {
                    ++j; continue;
                }
            }

            /* reset photon count for this pixel */
            I = 0;

            /* add background from something amorphous */
            F_bg = water_F;
            I_bg = F_bg*F_bg*r_e_sqr*fluence*polar*water_size*water_size*water_size*1e6*Avogadro/water_MW*omega_pixel;

            /* add this now to avoid problems with skipping later */
            floatimage[j] = I_bg;

            /* loop over sub-pixels */
            for(subS=0;subS<oversample;++subS)
            {
                for(subF=0;subF<oversample;++subF)
                {
                    /* absolute mm position on detector (relative to its origin) */
                    Fdet = subpixel_size*(fpixel*oversample + subF ) + subpixel_size/2.0;
                    Sdet = subpixel_size*(spixel*oversample + subS ) + subpixel_size/2.0;
//                  Fdet = pixel_size*fpixel;
//                  Sdet = pixel_size*spixel;

                    for(thick_tic=0;thick_tic<detector_thicksteps;++thick_tic)
                    {
                        /* assume "distance" is to the front of the detector sensor layer */
                        Odet = thick_tic*detector_thickstep;

                        /* construct detector subpixel position in 3D space */
//                      pixel_X = distance;
//                      pixel_Y = Sdet-Ybeam;
//                      pixel_Z = Fdet-Xbeam;
                        pixel_pos[1] = Fdet*fdet_vector[1]+Sdet*sdet_vector[1]+Odet*odet_vector[1]+pix0_vector[1];
                        pixel_pos[2] = Fdet*fdet_vector[2]+Sdet*sdet_vector[2]+Odet*odet_vector[2]+pix0_vector[2];
                        pixel_pos[3] = Fdet*fdet_vector[3]+Sdet*sdet_vector[3]+Odet*odet_vector[3]+pix0_vector[3];
                        pixel_pos[0] = 0.0;
                        if(curved_detector) {
                            /* construct detector pixel that is always "distance" from the sample */
                            vector[1] = distance*beam_vector[1];
                            vector[2]=distance*beam_vector[2] ;
                            vector[3]=distance*beam_vector[3];
                            /* treat detector pixel coordinates as radians */
                            rotate_axis(vector,newvector,sdet_vector,pixel_pos[2]/distance);
                            rotate_axis(newvector,pixel_pos,fdet_vector,pixel_pos[3]/distance);
//                          rotate(vector,pixel_pos,0,pixel_pos[3]/distance,pixel_pos[2]/distance);
                        }
                        /* construct the diffracted-beam unit vector to this sub-pixel */
                        airpath = unitize(pixel_pos,diffracted);

                        /* solid angle subtended by a pixel: (pix/airpath)^2*cos(2theta) */
                        omega_pixel = pixel_size*pixel_size/airpath/airpath*close_distance/airpath;
                        /* option to turn off obliquity effect, inverse-square-law only */
                        if(point_pixel) omega_pixel = 1.0/airpath/airpath;
                        omega_sum += omega_pixel;

                        /* now calculate detector thickness effects */
                        if(detector_thick > 0.0)
                        {
                            /* inverse of effective thickness increase */
                            parallax = dot_product(diffracted,odet_vector);
                            capture_fraction = exp(-thick_tic*detector_thickstep*detector_mu/parallax)
                                              -exp(-(thick_tic+1)*detector_thickstep*detector_mu/parallax);
                        }
                        else
                        {
                            capture_fraction = 1.0;
                        }

                        /* loop over sources now */
                        for(source=0;source<sources;++source){

                            /* retrieve stuff from cache */
                            incident[1] = -source_X[source];
                            incident[2] = -source_Y[source];
                            incident[3] = -source_Z[source];
                            lambda = source_lambda[source];

                            /* construct the incident beam unit vector while recovering source distance */
                            source_path = unitize(incident,incident);

                            /* construct the scattering vector for this pixel */
                            scattering[1] = (diffracted[1]-incident[1])/lambda;
                            scattering[2] = (diffracted[2]-incident[2])/lambda;
                            scattering[3] = (diffracted[3]-incident[3])/lambda;

                            /* sin(theta)/lambda is half the scattering vector length */
                            stol = 0.5*magnitude(scattering);

                            /* rough cut to speed things up when we aren't using whole detector */
                            if(dmin > 0.0 && stol > 0.0)
                            {
                                if(dmin > 0.5/stol)
                                {
                                    continue;
                                }
                            }

                            /* sweep over phi angles */
                            for(phi_tic = 0; phi_tic < phisteps; ++phi_tic)
                            {
                                phi = phi0 + phistep*phi_tic;

                                if( phi != 0.0 )
                                {
                                    /* rotate about spindle if neccesary */
                                    rotate_axis(a0,ap,spindle_vector,phi);
                                    rotate_axis(b0,bp,spindle_vector,phi);
                                    rotate_axis(c0,cp,spindle_vector,phi);
                                }

                                /* enumerate mosaic domains */
                                for(mos_tic=0;mos_tic<mosaic_domains;++mos_tic)
                                {
                                    /* apply mosaic rotation after phi rotation */
                                    if( mosaic_spread > 0.0 )
                                    {
                                        rotate_umat(ap,a,&mosaic_umats[mos_tic*9]);
                                        rotate_umat(bp,b,&mosaic_umats[mos_tic*9]);
                                        rotate_umat(cp,c,&mosaic_umats[mos_tic*9]);
                                    }
                                    else
                                    {
                                        a[1]=ap[1];a[2]=ap[2];a[3]=ap[3];
                                        b[1]=bp[1];b[2]=bp[2];b[3]=bp[3];
                                        c[1]=cp[1];c[2]=cp[2];c[3]=cp[3];
                                    }
//                                  printf("%d %f %f %f\n",mos_tic,mosaic_umats[mos_tic*9+0],mosaic_umats[mos_tic*9+1],mosaic_umats[mos_tic*9+2]);
//                                  printf("%d %f %f %f\n",mos_tic,mosaic_umats[mos_tic*9+3],mosaic_umats[mos_tic*9+4],mosaic_umats[mos_tic*9+5]);
//                                  printf("%d %f %f %f\n",mos_tic,mosaic_umats[mos_tic*9+6],mosaic_umats[mos_tic*9+7],mosaic_umats[mos_tic*9+8]);

                                    /* construct fractional Miller indicies */
                                    h = dot_product(a,scattering);
                                    k = dot_product(b,scattering);
                                    l = dot_product(c,scattering);

                                    /* round off to nearest whole index */
                                    h0 = ceil(h-0.5);
                                    k0 = ceil(k-0.5);
                                    l0 = ceil(l-0.5);


                                    /* structure factor of the lattice (paralelpiped crystal)
                                        F_latt = sin(M_PI*Na*h)*sin(M_PI*Nb*k)*sin(M_PI*Nc*l)/sin(M_PI*h)/sin(M_PI*k)/sin(M_PI*l);
                                    */
                                    F_latt = 1.0;
                                    if(xtal_shape == SQUARE)
                                    {
                                        /* xtal is a paralelpiped */
                                        if(Na>1){
                                            F_latt *= sincg(M_PI*h,Na);
                                        }
                                        if(Nb>1){
                                            F_latt *= sincg(M_PI*k,Nb);
                                        }
                                        if(Nc>1){
                                            F_latt *= sincg(M_PI*l,Nc);
                                        }
                                    }
                                    else
                                    {
                                        /* handy radius in reciprocal space, squared */
                                        hrad_sqr = (h-h0)*(h-h0)*Na*Na + (k-k0)*(k-k0)*Nb*Nb + (l-l0)*(l-l0)*Nc*Nc ;
                                    }
                                    if(xtal_shape == ROUND)
                                    {
                                        /* use sinc3 for elliptical xtal shape,
                                           correcting for sqrt of volume ratio between cube and sphere */
                                        F_latt = Na*Nb*Nc*0.723601254558268*sinc3(M_PI*sqrt( hrad_sqr * fudge ) );
                                    }
                                    if(xtal_shape == GAUSS)
                                    {
                                        /* fudge the radius so that volume and FWHM are similar to square_xtal spots */
                                        F_latt = Na*Nb*Nc*exp(-( hrad_sqr / 0.63 * fudge ));
                                    }
                                    if(xtal_shape == TOPHAT)
                                    {
                                        /* make a flat-top spot of same height and volume as square_xtal spots */
                                        F_latt = Na*Nb*Nc*(hrad_sqr*fudge < 0.3969 );
                                    }
                                    /* no need to go further if result will be zero? */
                                    if(F_latt == 0.0 && water_size == 0.0) continue;


                                    /* find nearest point on Ewald sphere surface? */
                                    if( integral_form )
                                    {

                                        if( phi != 0.0 || mos_tic > 0 )
                                        {
                                            /* need to re-calculate reciprocal matrix */

                                            /* various cross products */
                                            cross_product(a,b,a_cross_b);
                                            cross_product(b,c,b_cross_c);
                                            cross_product(c,a,c_cross_a);

                                            /* new reciprocal-space cell vectors */
                                            vector_scale(b_cross_c,a_star,1e20/V_cell);
                                            vector_scale(c_cross_a,b_star,1e20/V_cell);
                                            vector_scale(a_cross_b,c_star,1e20/V_cell);
                                        }

                                        /* reciprocal-space coordinates of nearest relp */
                                        relp[1] = h0*a_star[1] + k0*b_star[1] + l0*c_star[1];
                                        relp[2] = h0*a_star[2] + k0*b_star[2] + l0*c_star[2];
                                        relp[3] = h0*a_star[3] + k0*b_star[3] + l0*c_star[3];
//                                      d_star = magnitude(relp)

                                        /* reciprocal-space coordinates of center of Ewald sphere */
                                        Ewald0[1] = -incident[1]/lambda/1e10;
                                        Ewald0[2] = -incident[2]/lambda/1e10;
                                        Ewald0[3] = -incident[3]/lambda/1e10;
//                                      1/lambda = magnitude(Ewald0)

                                        /* distance from Ewald sphere in lambda=1 units */
                                        vector[1] = relp[1]-Ewald0[1];
                                        vector[2] = relp[2]-Ewald0[2];
                                        vector[3] = relp[3]-Ewald0[3];
                                        d_r = magnitude(vector)-1.0;

                                        /* unit vector of diffracted ray through relp */
                                        unitize(vector,diffracted0);

                                        /* intersection with detector plane */
                                        xd = dot_product(fdet_vector,diffracted0);
                                        yd = dot_product(sdet_vector,diffracted0);
                                        zd = dot_product(odet_vector,diffracted0);

                                        /* where does the central direct-beam hit */
                                        xd0 = dot_product(fdet_vector,incident);
                                        yd0 = dot_product(sdet_vector,incident);
                                        zd0 = dot_product(odet_vector,incident);

                                        /* convert to mm coordinates */
                                        Fdet0 = distance*(xd/zd) + Xbeam;
                                        Sdet0 = distance*(yd/zd) + Ybeam;

                                        //printf("GOTHERE %g %g   %g %g\n",Fdet,Sdet,Fdet0,Sdet0);
                                        test = exp(-( (Fdet-Fdet0)*(Fdet-Fdet0)+(Sdet-Sdet0)*(Sdet-Sdet0) + d_r*d_r )/1e-8);
                                    } // end of integral form


                                    /* structure factor of the unit cell */
                                    if(interpolate){
                                        h0_flr = floor(h);
                                        k0_flr = floor(k);
                                        l0_flr = floor(l);


                                        if ( ((h-h_min+3)>h_range) ||
                                             (h-2<h_min)           ||
                                             ((k-k_min+3)>k_range) ||
                                             (k-2<k_min)           ||
                                             ((l-l_min+3)>l_range) ||
                                             (l-2<l_min)  ) {
                                            if(babble){
                                                babble=0;
                                                printf ("WARNING: out of range for three point interpolation: h,k,l,h0,k0,l0: %g,%g,%g,%d,%d,%d \n", h,k,l,h0,k0,l0);
                                                printf("WARNING: further warnings will not be printed! ");
                                            }
                                            F_cell = default_F;
                                            interpolate=0;
                                        }
                                    }

                                    /* only interpolate if it is safe */
                                    if(interpolate){
                                        /* integer versions of nearest HKL indicies */
                                        h_interp[0]=h0_flr-1;
                                        h_interp[1]=h0_flr;
                                        h_interp[2]=h0_flr+1;
                                        h_interp[3]=h0_flr+2;
                                        k_interp[0]=k0_flr-1;
                                        k_interp[1]=k0_flr;
                                        k_interp[2]=k0_flr+1;
                                        k_interp[3]=k0_flr+2;
                                        l_interp[0]=l0_flr-1;
                                        l_interp[1]=l0_flr;
                                        l_interp[2]=l0_flr+1;
                                        l_interp[3]=l0_flr+2;

                                        /* polin function needs doubles */
                                        h_interp_d[0] = (double) h_interp[0];
                                        h_interp_d[1] = (double) h_interp[1];
                                        h_interp_d[2] = (double) h_interp[2];
                                        h_interp_d[3] = (double) h_interp[3];
                                        k_interp_d[0] = (double) k_interp[0];
                                        k_interp_d[1] = (double) k_interp[1];
                                        k_interp_d[2] = (double) k_interp[2];
                                        k_interp_d[3] = (double) k_interp[3];
                                        l_interp_d[0] = (double) l_interp[0];
                                        l_interp_d[1] = (double) l_interp[1];
                                        l_interp_d[2] = (double) l_interp[2];
                                        l_interp_d[3] = (double) l_interp[3];

                                        /* now populate the "y" values (nearest four structure factors in each direction) */
                                        for (i1=0;i1<4;i1++) {
                                            for (i2=0;i2<4;i2++) {
                                               for (i3=0;i3<4;i3++) {
                                                      sub_Fhkl[i1][i2][i3]= Fhkl[h_interp[i1]-h_min][k_interp[i2]-k_min][l_interp[i3]-l_min];
                                               }
                                            }
                                         }


                                        /* run the tricubic polynomial interpolation */
                                        polin3(h_interp_d,k_interp_d,l_interp_d,sub_Fhkl,h,k,l,&F_cell);
                                    }

                                    if(! interpolate)
                                    {
                                        if ( hkls && (h0<=h_max) && (h0>=h_min) && (k0<=k_max) && (k0>=k_min) && (l0<=l_max) && (l0>=l_min)  ) {
                                            /* just take nearest-neighbor */
                                            F_cell = Fhkl[h0-h_min][k0-k_min][l0-l_min];
                                        }
                                        else
                                        {
                                            F_cell = default_F;  // usually zero
                                        }
                                    }

                                    /* now we have the structure factor for this pixel */

                                    /* polarization factor */
                                    if(! nopolar){
                                        /* need to compute polarization factor */
                                        polar = polarization_factor(polarization,incident,diffracted,polar_vector);
                                    }
                                    else
                                    {
                                        polar = 1.0;
                                    }

                                    /* convert amplitudes into intensity (photons per steradian) */
                                    I += F_cell*F_cell*F_latt*F_latt*capture_fraction;
                                }
                                /* end of mosaic loop */
                            }
                            /* end of phi loop */
                        }
                        /* end of source loop */
                    }
                    /* end of detector thickness loop */
                }
                /* end of sub-pixel y loop */
            }
            /* end of sub-pixel x loop */


            floatimage[j] += r_e_sqr*fluence*polar*I/steps*omega_pixel;
//          floatimage[j] = test;
            if(floatimage[j] > max_I) {
                max_I = floatimage[j];
                max_I_x = Fdet;
                max_I_y = Sdet;
            }
            sum += floatimage[j];
            sumsqr += floatimage[j]*floatimage[j];
            ++sumn;

            if( printout )
            {
                if((fpixel==printout_fpixel && spixel==printout_spixel) || printout_fpixel < 0)
                {
                    twotheta = atan2(sqrt(pixel_pos[2]*pixel_pos[2]+pixel_pos[3]*pixel_pos[3]),pixel_pos[1]);
                    test = sin(twotheta/2.0)/(lambda0*1e10);
                    printf("%4d %4d : stol = %g or %g\n", fpixel,spixel,stol,test);
                    printf("at %g %g %g\n", pixel_pos[1],pixel_pos[2],pixel_pos[3]);
                    printf("hkl= %f %f %f  hkl0= %d %d %d\n", h,k,l,h0,k0,l0);
                    printf(" F_cell=%g  F_latt=%g   I = %g\n", F_cell,F_latt,I);
                    printf("I/steps %15.10g\n", I/steps);
                    printf("polar   %15.10g\n", polar);
                    printf("omega   %15.10g\n", omega_pixel);
                    printf("capfrac %15.10g\n", capture_fraction);
                    printf("pixel   %15.10g\n", floatimage[j]);
                    printf("real-space cell vectors (Angstrom):\n");
                    printf("     %-10s  %-10s  %-10s\n","a","b","c");
                    printf("X: %11.8f %11.8f %11.8f\n",a[1]*1e10,b[1]*1e10,c[1]*1e10);
                    printf("Y: %11.8f %11.8f %11.8f\n",a[2]*1e10,b[2]*1e10,c[2]*1e10);
                    printf("Z: %11.8f %11.8f %11.8f\n",a[3]*1e10,b[3]*1e10,c[3]*1e10);
                }
            }
            else
            {
                if(progress_meter && progress_pixels/100 > 0)
                {
                    if(progress_pixel % ( progress_pixels/20 ) == 0 ||
                       ((10*progress_pixel<progress_pixels ||
                         10*progress_pixel>9*progress_pixels) &&
                        (progress_pixel % (progress_pixels/100) == 0)))
                    {
                        printf("%lu%% done\n",progress_pixel*100/progress_pixels);
                    }
                }
            }
            ++j;
            ++progress_pixel;
        }
    }
    printf("\n");

    printf("solid angle subtended by detector = %g steradian ( %g%% sphere)\n",omega_sum/steps,100*omega_sum/steps/4/M_PI);

    /* do some stats? */
    if(sumn<=0) sumn=1;
    avg = sum/sumn;
    if(sumn<=1) sumn=2;
    rms = sqrt(sumsqr/(sumn-1));
    sumsqr = 0.0;
    j = sumn = 0;
    for(spixel=0;spixel<spixels;++spixel)
    {
        for(fpixel=0;fpixel<fpixels;++fpixel)
        {
            ++j;
            if(fpixel < roi_xmin || fpixel > roi_xmax || spixel < roi_ymin || spixel > roi_ymax)
            {
                continue;
            }
            test = floatimage[j]-avg;
            sumsqr += test*test;
            ++sumn;
        }
    }
    if(sumn<=1) sumn=2;
    rmsd = sqrt(sumsqr/(sumn-1));

    printf("writing %s as %d %lu-byte floats\n",floatfilename,pixels,sizeof(float));
    outfile = fopen(floatfilename,"wb");
    if(outfile == NULL)
    {
        perror("ERROR: fopen");
        exit(9);
    }
    fwrite(floatimage,sizeof(float),pixels,outfile);
    fclose(outfile);

    /* output as ints */
    j = 0;
    printf("max_I = %g  at %g %g\n",max_I,max_I_x,max_I_y);
    printf("mean= %g rms= %g rmsd= %g\n",avg,rms,rmsd);
    if(intfile_scale <= 0.0){
        intfile_scale = 1.0;
        if(max_I > 0.0) intfile_scale = 55000.0/max_I;
    }
    printf("intfile_scale = %g\n",intfile_scale);
    for(spixel=0;spixel<spixels;++spixel)
    {
        for(fpixel=0;fpixel<fpixels;++fpixel)
        {
            if(fpixel < roi_xmin || fpixel > roi_xmax || spixel < roi_ymin || spixel > roi_ymax)
            {
               ++j; continue;
            }
            test = floatimage[j] *intfile_scale+adc_offset;
            if(test > 65535.0) test = 65535.0;
            if(test < 0.0) test = 0.0;
            intimage[j] = (unsigned short int) ( floorf(test+0.5) );
//          printf("%d %d = %d\n",fpixel,spixel,intimage[j]);
            ++j;
        }
    }

    printf("writing %s as %lu-byte integers\n",intfilename,sizeof(unsigned short int));
    outfile = fopen(intfilename,"wb");
    if(outfile == NULL)
    {
            perror("ERROR: fopen");
            exit(9);
    }
    fprintf(outfile,"{\nHEADER_BYTES=512;\nDIM=2;\nBYTE_ORDER=%s;\nTYPE=unsigned_short;\n",byte_order);
    fprintf(outfile,"SIZE1=%d;\nSIZE2=%d;\nPIXEL_SIZE=%g;\nDISTANCE=%g;\n",fpixels,spixels,pixel_size*1000.0,distance*1000.0);
    fprintf(outfile,"WAVELENGTH=%g;\n",lambda0*1e10);
    fprintf(outfile,"BEAM_CENTER_X=%g;\nBEAM_CENTER_Y=%g;\n",Xbeam*1000.0,Ybeam*1000);
    fprintf(outfile,"ADXV_CENTER_X=%g;\nADXV_CENTER_Y=%g;\n",Fbeam*1000.0,(detsize_s-Sbeam)*1000);
    fprintf(outfile,"MOSFLM_CENTER_X=%g;\nMOSFLM_CENTER_Y=%g;\n",(Sbeam-0.5*pixel_size)*1000.0,(Fbeam-0.5*pixel_size)*1000);
    fprintf(outfile,"DENZO_X_BEAM=%g;\nDENZO_Y_BEAM=%g;\n",(Sbeam-0.0*pixel_size)*1000.0,(Fbeam-0.0*pixel_size)*1000);
    fprintf(outfile,"DIALS_ORIGIN=%g,%g,%g\n",dials_origin[1],dials_origin[2],dials_origin[3]);
    fprintf(outfile,"XDS_ORGX=%g;\nXDS_ORGY=%g;\n",ORGX,ORGY);
    fprintf(outfile,"CLOSE_DISTANCE=%g;\n",close_distance*1000.0);
    fprintf(outfile,"PHI=%g;\nOSC_START=%g;\nOSC_RANGE=%g;\n",phi0*RTD,phi0*RTD,osc*RTD);
    fprintf(outfile,"TWOTHETA=%g;\n",detector_twotheta*RTD);
    fprintf(outfile,"DETECTOR_SN=000;\n");
    fprintf(outfile,"BEAMLINE=fake;\n");
    fprintf(outfile,"}\f");
    while ( ftell(outfile) < 512 ){ fprintf(outfile," "); };
    fwrite(intimage,sizeof(unsigned short int),pixels,outfile);
    fclose(outfile);


    if(write_pgm)
    {
        /* output as pgm */
        j = 0;
        if(pgm_scale <= 0.0){
            pgm_scale = intfile_scale;
            if(rmsd > 0.0) pgm_scale = 250.0/(5.0*rmsd);
        }
        printf("pgm_scale = %g\n",pgm_scale);
        j = 0;
        for(spixel=0;spixel<spixels;++spixel)
        {
            for(fpixel=0;fpixel<fpixels;++fpixel)
            {
                if(fpixel < roi_xmin || fpixel > roi_xmax || spixel < roi_ymin || spixel > roi_ymax)
                {
                    ++j; continue;
                }
                test = floatimage[j] * pgm_scale;
                if(test > 255.0) test = 255.0;
                pgmimage[j] = (unsigned char) ( test );
//              printf("%d %d = %d\n",fpixel,spixel,pgmimage[j]);
                ++j;
            }
        }

        printf("writing %s as %lu-byte integers\n",pgmfilename,sizeof(unsigned char));
        outfile = fopen(pgmfilename,"wb");
        if(outfile == NULL)
        {
                perror("ERROR: fopen");
                exit(9);
        }
        fprintf(outfile, "P5\n%d %d\n", fpixels, spixels);
        fprintf(outfile, "# pixels scaled by %lg\n", pgm_scale);
        fprintf(outfile, "255\n");
        fwrite(pgmimage,sizeof(unsigned char),pixels,outfile);
        fclose(outfile);
    }

    /* quit now if there is nothing else to do */
    if(calculate_noise == 0){
        return 0;
    }

    /* simulate Poisson noise */
    j = 0;
    sum = 0.0;
    overloads = 0;
    for(spixel=0;spixel<spixels;++spixel)
    {
        for(fpixel=0;fpixel<fpixels;++fpixel)
        {
            if(fpixel < roi_xmin || fpixel > roi_xmax || spixel < roi_ymin || spixel > roi_ymax)
            {
                ++j; continue;
            }
            test = poidev( floatimage[j], &seed );
            sum += test;
            test += adc_offset;
            if(test > 65535.0)
            {
                test = 65535.0;
                ++overloads;
            }
            intimage[j] = (unsigned short int) test;
//          printf("%d %d = %d\n",fpixel,spixel,intimage[j]);
            ++j;
        }
    }
    printf("%.0f photons on noise image (%d overloads)\n",sum,overloads);

    printf("writing %s as %lu-byte integers\n",noisefilename,sizeof(unsigned short int));
    outfile = fopen(noisefilename,"wb");
    if(outfile == NULL)
    {
            perror("ERROR: fopen");
            exit(9);
    }
    fprintf(outfile,"{\nHEADER_BYTES=512;\nDIM=2;\nBYTE_ORDER=%s;\nTYPE=unsigned_short;\n",byte_order);
    fprintf(outfile,"SIZE1=%d;\nSIZE2=%d;\nPIXEL_SIZE=%g;\nDISTANCE=%g;\n",fpixels,spixels,pixel_size*1000.0,distance*1000.0);
    fprintf(outfile,"WAVELENGTH=%g;\n",lambda0*1e10);
    fprintf(outfile,"BEAM_CENTER_X=%g;\nBEAM_CENTER_Y=%g;\n",Xbeam*1000.0,Ybeam*1000);
    fprintf(outfile,"ADXV_CENTER_X=%g;\nADXV_CENTER_Y=%g;\n",Fbeam*1000.0,(detsize_s-Sbeam)*1000);
    fprintf(outfile,"MOSFLM_CENTER_X=%g;\nMOSFLM_CENTER_Y=%g;\n",(Sbeam-0.5*pixel_size)*1000.0,(Fbeam-0.5*pixel_size)*1000);
    fprintf(outfile,"DENZO_X_BEAM=%g;\nDENZO_Y_BEAM=%g;\n",(Sbeam+0.0*pixel_size)*1000.0,(Fbeam+0.0*pixel_size)*1000);
    fprintf(outfile,"DIALS_ORIGIN=%g,%g,%g\n",dials_origin[1],dials_origin[2],dials_origin[3]);
    fprintf(outfile,"XDS_ORGX=%g;\nXDS_ORGY=%g;\n",ORGX,ORGY);
    fprintf(outfile,"CLOSE_DISTANCE=%g;\n",close_distance*1000.0);
    fprintf(outfile,"PHI=%g;\nOSC_START=%g;\nOSC_RANGE=%g;\n",phi0*RTD,phi0*RTD,osc*RTD);
    fprintf(outfile,"TWOTHETA=%g;\n",detector_twotheta*RTD);
    fprintf(outfile,"DETECTOR_SN=000;\n");
    fprintf(outfile,"BEAMLINE=fake;\n");
    fprintf(outfile,"}\f");
    while ( ftell(outfile) < 512 ){ fprintf(outfile," "); };
    fwrite(intimage,sizeof(unsigned short int),pixels,outfile);
    fclose(outfile);

    return 0;
}



/* Fourier transform of a grating */
double sincg(double x,double N) {
    if(x==0.0) return N;

    return sin(x*N)/sin(x);
}

/* Fourier transform of a sphere */
double sinc3(double x) {
    if(x==0.0) return 1.0;

    return 3.0*(sin(x)/x-cos(x))/(x*x);
}

double sinc_conv_sinc3(double x) {
    if(x==0.0) return 1.0;

    return 3.0*(sin(x)-x*cos(x))/(x*x*x);
}


double *rotate(double *v, double *newv, double phix, double phiy, double phiz) {

    double rxx,rxy,rxz,ryx,ryy,ryz,rzx,rzy,rzz;
    double new_x,new_y,new_z,rotated_x,rotated_y,rotated_z;

    new_x=v[1];
    new_y=v[2];
    new_z=v[3];

    if(phix != 0){
        /* rotate around x axis */
        //rxx= 1;         rxy= 0;         rxz= 0;
        ryx= 0;         ryy= cos(phix); ryz=-sin(phix);
        rzx= 0;         rzy= sin(phix); rzz= cos(phix);

        rotated_x = new_x;
        rotated_y = new_y*ryy + new_z*ryz;
        rotated_z = new_y*rzy + new_z*rzz;
        new_x = rotated_x; new_y = rotated_y; new_z = rotated_z;
    }

    if(phiy != 0) {
        /* rotate around y axis */
        rxx= cos(phiy); rxy= 0;         rxz= sin(phiy);
        //ryx= 0;         ryy= 1;         ryz= 0;
        rzx=-sin(phiy); rzy= 0;         rzz= cos(phiy);

        rotated_x = new_x*rxx + new_y*rxy + new_z*rxz;
        rotated_y = new_y;
        rotated_z = new_x*rzx + new_y*rzy + new_z*rzz;
        new_x = rotated_x; new_y = rotated_y; new_z = rotated_z;
    }

    if(phiz != 0){
        /* rotate around z axis */
        rxx= cos(phiz); rxy=-sin(phiz); rxz= 0;
        ryx= sin(phiz); ryy= cos(phiz); ryz= 0;
        //rzx= 0;         rzy= 0;         rzz= 1;

        rotated_x = new_x*rxx + new_y*rxy ;
        rotated_y = new_x*ryx + new_y*ryy;
        rotated_z = new_z;
        new_x = rotated_x; new_y = rotated_y; new_z = rotated_z;
    }

    newv[1]=new_x;
    newv[2]=new_y;
    newv[3]=new_z;

    return newv;
}



/* rotate a point about a unit vector axis */
double *rotate_axis(double *v, double *newv, double *axis, double phi) {

    double sinphi = sin(phi);
    double cosphi = cos(phi);
    double dot = (axis[1]*v[1]+axis[2]*v[2]+axis[3]*v[3])*(1.0-cosphi);
    double temp[4];

    temp[1] = axis[1]*dot+v[1]*cosphi+(-axis[3]*v[2]+axis[2]*v[3])*sinphi;
    temp[2] = axis[2]*dot+v[2]*cosphi+(+axis[3]*v[1]-axis[1]*v[3])*sinphi;
    temp[3] = axis[3]*dot+v[3]*cosphi+(-axis[2]*v[1]+axis[1]*v[2])*sinphi;
    newv[1]=temp[1]; newv[2]=temp[2]; newv[3]=temp[3];

    return newv;
}



/* rotate a vector using a 9-element unitary matrix */
double *rotate_umat(double *v, double *newv, double umat[9]) {

    double uxx,uxy,uxz,uyx,uyy,uyz,uzx,uzy,uzz;

    /* for convenience, assign matrix x-y coordinate */
    uxx = umat[0];
    uxy = umat[1];
    uxz = umat[2];
    uyx = umat[3];
    uyy = umat[4];
    uyz = umat[5];
    uzx = umat[6];
    uzy = umat[7];
    uzz = umat[8];

    /* rotate the vector (x=1,y=2,z=3) */
    newv[1] = uxx*v[1] + uxy*v[2] + uxz*v[3];
    newv[2] = uyx*v[1] + uyy*v[2] + uyz*v[3];
    newv[3] = uzx*v[1] + uzy*v[2] + uzz*v[3];

    return newv;
}




/* returns a unit vector in a random direction in arguments dx,dy,dz */
/* also returns a random magnitude within the unit sphere as a return value */
float uniform3Ddev(float *dx, float *dy, float *dz, long *seed)
{
    float ran1(long *idum);
    float dr;

    /* pick a random direction by cutting a sphere out of a cube */
    dr = 0;
    while(dr>1 || dr < 1e-2)
    {
        *dx = 2.1*(ran1(seed)-0.5);
        *dy = 2.1*(ran1(seed)-0.5);
        *dz = 2.1*(ran1(seed)-0.5);
        dr = sqrt(*dx**dx+*dy**dy+*dz**dz);
    }
    /* turn this into a unit vector */
    *dx/=dr;
    *dy/=dr;
    *dz/=dr;

    /* dx,dy,dz should now be a random unit vector */

    return dr;
}


/* returns a 9-element unitary matrix for a random isotropic rotation on a spherical cap of diameter "mosaicity" */
/* mosaic = 90 deg is a full sphere */
double *mosaic_rotation_umat(float mosaicity, double umat[9], long *seed)
{
    float ran1(long *idum);
    double r1,r2,r3,xyrad,rot;
    double v1,v2,v3;
    double t1,t2,t3,t6,t7,t8,t9,t11,t12,t15,t19,t20,t24;
    double uxx,uxy,uxz,uyx,uyy,uyz,uzx,uzy,uzz;

    /* make three random uniform deviates on [-1:1] */
    r1= (double) 2.0*ran1(seed)-1.0;
    r2= (double) 2.0*ran1(seed)-1.0;
    r3= (double) 2.0*ran1(seed)-1.0;

    xyrad = sqrt(1.0-r2*r2);
    rot = mosaicity*powf((1.0-r3*r3),(1.0/3.0));

    v1 = xyrad*sin(M_PI*r1);
    v2 = xyrad*cos(M_PI*r1);
    v3 = r2;

    /* commence incomprehensible quaternion calculation */
    t1 =  cos(rot);
    t2 =  1.0 - t1;
    t3 =  v1*v1;
    t6 =  t2*v1;
    t7 =  t6*v2;
    t8 =  sin(rot);
    t9 =  t8*v3;
    t11 = t6*v3;
    t12 = t8*v2;
    t15 = v2*v2;
    t19 = t2*v2*v3;
    t20 = t8*v1;
    t24 = v3*v3;

    /* populate the unitary rotation matrix */
    umat[0] = uxx = t1 + t2*t3;
    umat[1] = uxy = t7 - t9;
    umat[2] = uxz = t11 + t12;
    umat[3] = uyx = t7 + t9;
    umat[4] = uyy = t1 + t2*t15;
    umat[5] = uyz = t19 - t20;
    umat[6] = uzx = t11 - t12;
    umat[7] = uzy = t19 + t20;
    umat[8] = uzz = t1 + t2*t24;

    /* return pointer to the provided array, in case that is useful */
    return umat;
}

/* convert a unitary rotation matrix into misseting angles
   rotx roty rotz are returned as missets[1] missets[2] missets[3] */
double *umat2misset(double umat[9],double *missets)
{
    double uxx,uxy,uxz,uyx,uyy,uyz,uzx,uzy,uzz;
    double m,mx,my,mz;
    double xcy_x,xcy_y,xcy_z;
    double ycz_x,ycz_y,ycz_z;
    double zcx_x,zcx_y,zcx_z;
    double rotx,roty,rotz;

    uxx=umat[0];uxy=umat[1];uxz=umat[2];
    uyx=umat[3];uyy=umat[4];uyz=umat[5];
    uzx=umat[6];uzy=umat[7];uzz=umat[8];

    /* or transpose? */
//    uxx=umat[1];uyx=umat[2];uzx=umat[3];
//    uxy=umat[4];uyy=umat[5];uzy=umat[6];
//    uxz=umat[7];uyz=umat[8];uzz=umat[9];

    /* make sure it is unitary */
    mx = sqrt(uxx*uxx+uxy*uxy+uxz*uxz);
    my = sqrt(uyx*uyx+uyy*uyy+uyz*uyz);
    mz = sqrt(uzx*uzx+uzy*uzy+uzz*uzz);
    if(mx>0){uxx/=mx;uxy/=mx;uxz/=mx;};
    if(my>0){uyx/=my;uyy/=my;uyz/=my;};
    if(mz>0){uzx/=mz;uzy/=mz;uzz/=mz;};

    if(mx>=0 && my<=0 && mz<=0)
    {
        uyx=0;uyy=1;uyz=0;
        uzx=0;uzy=0;uzz=1;
    }
    if(mx<=0 && my>=0 && mz<=0)
    {
        uxx=1;uxy=0;uxz=0;
        uzx=0;uzy=0;uzz=1;
    }
    if(mx<=0 && my<=0 && mz>=0)
    {
        uxx=1;uxy=0;uxz=0;
        uyx=0;uyy=1;uyz=0;
    }

    /* cross products to check normality */
    xcy_x = uxy*uyz - uxz*uyy;
    xcy_y = uxz*uyx - uxx*uyz;
    xcy_z = uxx*uyy - uxy*uyx;
    m=sqrt(xcy_x*xcy_x+xcy_y*xcy_y+xcy_z*xcy_z);
    if(m>0){xcy_x/=m;xcy_y/=m;xcy_z/=m;};

    ycz_x = uyy*uzz - uyz*uzy;
    ycz_y = uyz*uzx - uyx*uzz;
    ycz_z = uyx*uzy - uyy*uzx;
    m=sqrt(ycz_x*ycz_x+ycz_y*ycz_y+ycz_z*ycz_z);
    if(m>0){ycz_x/=m;ycz_y/=m;ycz_z/=m;};

    zcx_x = uzy*uxz - uzz*uxy;
    zcx_y = uzz*uxx - uzx*uxz;
    zcx_z = uzx*uxy - uzy*uxx;
    m=sqrt(zcx_x*zcx_x+zcx_y*zcx_y+zcx_z*zcx_z);
    if(m>0){zcx_x/=m;zcx_y/=m;zcx_z/=m;};

    /* substitute any empty vectors for cross-product of other two */
    if(mx<=0){uxx=ycz_x;uxy=ycz_y;uxz=ycz_z;};
    if(my<=0){uyx=zcx_x;uyy=zcx_y;uyz=zcx_z;};
    if(mz<=0){uzx=xcy_x;uzy=xcy_y;uzz=xcy_z;};


    /* cross products to check normality */
    xcy_x = uxy*uyz - uxz*uyy;
    xcy_y = uxz*uyx - uxx*uyz;
    xcy_z = uxx*uyy - uxy*uyx;
    m=sqrt(xcy_x*xcy_x+xcy_y*xcy_y+xcy_z*xcy_z);
    if(m>0){xcy_x/=m;xcy_y/=m;xcy_z/=m;}

    ycz_x = uyy*uzz - uyz*uzy;
    ycz_y = uyz*uzx - uyx*uzz;
    ycz_z = uyx*uzy - uyy*uzx;
    m=sqrt(ycz_x*ycz_x+ycz_y*ycz_y+ycz_z*ycz_z);
    if(m>0){ycz_x/=m;ycz_y/=m;ycz_z/=m;};

    zcx_x = uzy*uxz - uzz*uxy;
    zcx_y = uzz*uxx - uzx*uxz;
    zcx_z = uzx*uxy - uzy*uxx;
    m=sqrt(zcx_x*zcx_x+zcx_y*zcx_y+zcx_z*zcx_z);
    if(m>0){zcx_x/=m;zcx_y/=m;zcx_z/=m;};


    /* substitute any empty vectors for cross-product of other two */
    if(mx<=0){uxx=ycz_x;uxy=ycz_y;uxz=ycz_z;};
    if(my<=0){uyx=zcx_x;uyy=zcx_y;uyz=zcx_z;};
    if(mz<=0){uzx=xcy_x;uzy=xcy_y;uzz=xcy_z;};



    /* make sure it is unitary */
    mx = sqrt(uxx*uxx+uxy*uxy+uxz*uxz);
    my = sqrt(uyx*uyx+uyy*uyy+uyz*uyz);
    mz = sqrt(uzx*uzx+uzy*uzy+uzz*uzz);
    if(mx>0){uxx/=mx;uxy/=mx;uxz/=mx;};
    if(my>0){uyx/=my;uyy/=my;uyz/=my;};
    if(mz>0){uzx/=mz;uzy/=mz;uzz/=mz;};

    /* see if its really orthonormal? */

    if(uzx*uzx < 1.0)
    {
        rotx = atan2(uzy,uzz);
        roty = atan2(-uzx,sqrt(uzy*uzy+uzz*uzz));
        rotz = atan2(uyx,uxx);
    }
    else
    {
        rotx = atan2(1,1)*4;
        roty = atan2(1,1)*2;
        rotz = atan2(uxy,-uyy);
    }

    missets[1] = rotx;
    missets[2] = roty;
    missets[3] = rotz;
    return missets;
}



float poidev(float xm, long *idum)
{
    float gammln(float xx);
    float ran1(long *idum);
    /* oldm is a flag for whether xm has changed since last call */
    static float sq,alxm,g,oldm=(-1.0);
    float em,t,y;

    /* routine below locks up for > 1e6 photons? */
    if (xm > 1.0e6) {
        return xm+sqrt(xm)*gaussdev(idum);
    }

    if (xm < 12.0) {
        /* use direct method: simulate exponential delays between events */
        if(xm != oldm) {
            /* xm is new, compute the exponential */
            oldm=xm;
            g=exp(-xm);
        }
        /* adding exponential deviates is equivalent to multiplying uniform deviates */
        /* final comparison is to the pre-computed exponential */
        em = -1;
        t = 1.0;
        do {
            ++em;
            t *= ran1(idum);
        } while (t > g);
    } else {
        /* Use rejection method */
        if(xm != oldm) {
            /* xm has changed, pre-compute a few things... */
            oldm=xm;
            sq=sqrt(2.0*xm);
            alxm=log(xm);
            g=xm*alxm-gammln(xm+1.0);
        }
        do {
            do {
                /* y is a deviate from a lorentzian comparison function */
                y=tan(M_PI*ran1(idum));
                /* shift and scale */
                em=sq*y+xm;
            } while (em < 0.0);         /* there are no negative Poisson deviates */
            /* round off to nearest integer */
            em=floor(em);
            /* ratio of Poisson distribution to comparison function */
            /* scale it back by 0.9 to make sure t is never > 1.0 */
            t=0.9*(1.0+y*y)*exp(em*alxm-gammln(em+1.0)-g);
        } while (ran1(idum) > t);
    }

    return em;
}


/* return gaussian deviate with rms=1 and FWHM = 2/sqrt(log(2)) */
float gaussdev(long *idum)
{
    float ran1(long *idum);
    static int iset=0;
    static float gset;
    float fac,rsq,v1,v2;

    if (iset == 0) {
        /* no extra deviats handy ... */

        /* so pick two uniform deviates on [-1:1] */
        do {
            v1=2.0*ran1(idum)-1.0;
            v2=2.0*ran1(idum)-1.0;
            rsq=v1*v1+v2*v2;
        } while (rsq >= 1.0 || rsq == 0);
        /* restrained to the unit circle */

        /* apply Box-Muller transformation to convert to a normal deviate */
        fac=sqrt(-2.0*log(rsq)/rsq);
        gset=v1*fac;
        iset=1;         /* we now have a spare deviate */
        return v2*fac;
    } else {
        /* there is an extra deviate in gset */
        iset=0;
        return gset;
    }
}


/* generate Lorentzian deviate with FWHM = 2 */
float lorentzdev(long *seed) {
    float ran1(long *idum);

    return tan(M_PI*(ran1(seed)-0.5));
}

/* return triangular deviate with FWHM = 1 */
float triangledev(long *seed) {
    float ran1(long *idum);
    float value;

    value = ran1(seed);
    if(value > 0.5){
        value = sqrt(2*(value-0.5))-1;
    }else{
        value = 1-sqrt(2*value);
    }

    return value;
}



float expdev(long *idum)
{
    float dum;

    do
    dum=ran1(idum);
    while( dum == 0.0);
    return -log(dum);
}



/* ln of the gamma function */
float gammln(float xx)
{
    double x,y,tmp,ser;
    static double cof[6]={76.18009172947146,-86.50532032941677,
    24.01409824083091,-1.231739572450155,
    0.1208650973866179e-2,-0.5395239384953e-5};
    int j;

    y=x=xx;
    tmp=x+5.5;
    tmp -= (x+0.5)*log(tmp);
    ser = 1.000000000190015;
    for(j=0;j<=5;++j) ser += cof[j]/++y;

    return -tmp+log(2.5066282746310005*ser/x);
}





/* returns a uniform random deviate between 0 and 1 */
#define IA 16807
#define IM 2147483647
#define AM (1.0/IM)
#define IQ 127773
#define IR 2836
#define NTAB 32
#define NDIV (1+(IM-1)/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)

float ran1(long *idum)
{
    int j;
    long k;
    static long iy=0;
    static long iv[NTAB];
    float temp;

    if (*idum <= 0 || !iy) {
        /* first time around.  don't want idum=0 */
        if(-(*idum) < 1) *idum=1;
        else *idum = -(*idum);

        /* load the shuffle table */
        for(j=NTAB+7;j>=0;j--) {
            k=(*idum)/IQ;
            *idum=IA*(*idum-k*IQ)-IR*k;
            if(*idum < 0) *idum += IM;
            if(j < NTAB) iv[j] = *idum;
        }
        iy=iv[0];
    }
    /* always start here after initializing */
    k=(*idum)/IQ;
    *idum=IA*(*idum-k*IQ)-IR*k;
    if (*idum < 0) *idum += IM;
    j=iy/NDIV;
    iy=iv[j];
    iv[j] = *idum;
    if((temp=AM*iy) > RNMX) return RNMX;
    else return temp;
}


void polint(double *xa, double *ya, double x, double *y)
{
        double x0,x1,x2,x3;
        x0 = (x-xa[1])*(x-xa[2])*(x-xa[3])*ya[0]/((xa[0]-xa[1])*(xa[0]-xa[2])*(xa[0]-xa[3]));
        x1 = (x-xa[0])*(x-xa[2])*(x-xa[3])*ya[1]/((xa[1]-xa[0])*(xa[1]-xa[2])*(xa[1]-xa[3]));
        x2 = (x-xa[0])*(x-xa[1])*(x-xa[3])*ya[2]/((xa[2]-xa[0])*(xa[2]-xa[1])*(xa[2]-xa[3]));
        x3 = (x-xa[0])*(x-xa[1])*(x-xa[2])*ya[3]/((xa[3]-xa[0])*(xa[3]-xa[1])*(xa[3]-xa[2]));
        *y = x0+x1+x2+x3;
}



void polin2(double *x1a, double *x2a, double **ya, double x1, double x2, double *y)
{
        void polint(double *xa, double *ya, double x, double *y);
        int j;
        double ymtmp[4];
        for (j=1;j<=4;j++) {
                polint(x2a,ya[j-1],x2,&ymtmp[j-1]);
        }
        polint(x1a,ymtmp,x1,y);
}


void polin3(double *x1a, double *x2a, double *x3a, double ***ya, double x1,
        double x2, double x3, double *y)
{
        void polint(double *xa, double ya[], double x, double *y);
        void polin2(double *x1a, double *x2a, double **ya, double x1,double x2, double *y);
        void polin1(double *x1a, double *ya, double x1, double *y);
        int j;
        double ymtmp[4];

        for (j=1;j<=4;j++) {
            polin2(x2a,x3a,&ya[j-1][0],x2,x3,&ymtmp[j-1]);
        }
        polint(x1a,ymtmp,x1,y);
}


/* FWHM = integral = 1 */
double ngauss2D(double x,double y)
{
    return log(16.)/M_PI*exp(-log(16.)*(x*x+y*y));
}
double ngauss2Dinteg(double x,double y)
{
    return 0.125*(erf(2.*x*sqrt(log(2.)))*erf(y*sqrt(log(16.)))*sqrt(log(16.)/log(2.)));
}









/* read in multi-column text file to list of double arrays */
/* provide address of undeclared arrays on command line */
size_t read_text_file(char *filename, size_t nargs, ... )
{
    /* maximum of 10240-character lines? */
    char text[10240];
    char *token;
    const char delimiters[] = " \t,;:!";
    const char numberstuf[] = "0123456789-+.EGeg";

    unsigned long line,lines;
    unsigned long i,j,m;
    double value;
    double *data;
    double **pointer;
    va_list arglist;
    FILE *infile = NULL;

    infile = fopen(filename,"r");
    if(infile == NULL) {
        perror("fopen()");
        return 0;
    }
    lines=0;
    while ( fgets ( text, sizeof text, infile ) != NULL ) {
        token = text;
        token += strspn(token,delimiters);
        if(strcmp(token,"\n")==0) {
            //printf("blank\n");
            continue;
        }
        ++lines;
    }
    rewind(infile);

    /* allocate memory for arrays */
    va_start( arglist, nargs);
    for(i=0;i<nargs;++i){
        /* allocate the array */
        data = (double*) malloc((lines+10)*sizeof(double));
        /* initialize with missing number flags */
        for(j=0;j<lines+10;++j) {
            data[j] = NAN;
        }
        /* get argument (pointer to pointer) */
        pointer = va_arg(arglist, double **);
        /* change the value of what the arg points to */
        *pointer = data;
        /* now the pointer provided as an argument points to
        something */
    }
    va_end(arglist);

    line = 0;
    while ( fgets ( text, sizeof text, infile ) != NULL ) { /* read a line */

        token = text;
        token += strspn(token,delimiters);
        if(strcmp(token,"\n")==0) {
            //printf("blank\n");
            continue;
        }
        i=0;
        va_start( arglist, nargs);
        do
        {
            value=atof(token);
            /* get argument */
            pointer = va_arg(arglist, double **);
            /* retrieve data array's address */
            data = *pointer;
            data[line] = value;

            token += strspn(token,numberstuf);
            if (strcmp(token,"\n")==0) continue;
            token += strcspn(token,delimiters);
            token += strspn(token,delimiters);
            if (strcmp(token,"\n")==0) continue;

            ++i;
            if(i>=nargs) {
                break;
            }
        }
        while (strcmp(token,"\n")!=0) ;
        va_end(arglist);

//      printf("initializing:");
//        va_start( arglist, nargs);
//        for(i=0;i<nargs;++i){
//          pointer = va_arg(arglist, double **);
//          data = *pointer;
//          printf(" %g",data[line]);
//        }
//        va_end(arglist);
//      printf("\n");

        ++line;
    }
    fclose(infile);

    return lines;
}



/* measure magnitude of provided vector */
double magnitude(double *vector) {

    /* measure the magnitude */
    vector[0] = sqrt(vector[1]*vector[1]+vector[2]*vector[2]+vector[3]*vector[3]);

    return vector[0];
}

/* make provided vector a unit vector */
double unitize(double *vector, double *new_unit_vector) {
    double mag;

    /* measure the magnitude */
    mag = magnitude(vector);

    if(mag != 0.0){
        /* normalize it */
        new_unit_vector[1]=vector[1]/mag;
        new_unit_vector[2]=vector[2]/mag;
        new_unit_vector[3]=vector[3]/mag;
    }
    else
    {
        /* can't normalize, report zero vector */
        new_unit_vector[0] = 0.0;
        new_unit_vector[1] = 0.0;
        new_unit_vector[2] = 0.0;
        new_unit_vector[3] = 0.0;
    }
    return mag;
}

/* scale magnitude of provided vector */
double vector_scale(double *vector, double *new_vector, double scale) {

    new_vector[1] = scale*vector[1];
    new_vector[2] = scale*vector[2];
    new_vector[3] = scale*vector[3];

    return magnitude(new_vector);
}

/* enforce magnitude of provided vector */
double vector_rescale(double *vector, double *new_vector, double new_magnitude) {
    double oldmag;

    oldmag = magnitude(vector);
    if(oldmag <= 0.0) oldmag = 1.0;
    new_vector[1] = new_magnitude/oldmag*vector[1];
    new_vector[2] = new_magnitude/oldmag*vector[2];
    new_vector[3] = new_magnitude/oldmag*vector[3];

    return magnitude(new_vector);
}

/* difference between two given vectors */
double vector_diff(double *vector, double *origin_vector, double *new_vector) {

    new_vector[1] = vector[1]-origin_vector[1];
    new_vector[2] = vector[2]-origin_vector[2];
    new_vector[3] = vector[3]-origin_vector[3];
    return magnitude(new_vector);
}


/* vector cross product where vector magnitude is 0th element */
double *cross_product(double *x, double *y, double *z) {
    z[1] = x[2]*y[3] - x[3]*y[2];
    z[2] = x[3]*y[1] - x[1]*y[3];
    z[3] = x[1]*y[2] - x[2]*y[1];
    z[0] = 0.0;

    return z;
}
/* vector inner product where vector magnitude is 0th element */
double dot_product(double *x, double *y) {
    return x[1]*y[1]+x[2]*y[2]+x[3]*y[3];
}


/* polarization factor */
double polarization_factor(double kahn_factor, double *incident, double *diffracted, double *axis)
{
    double cos2theta,cos2theta_sqr,sin2theta_sqr;
    double twotheta,psi;
    double E_in[4];
    double B_in[4];
    double E_out[4];
    double B_out[4];

    unitize(incident,incident);
    unitize(diffracted,diffracted);
    unitize(axis,axis);

    /* component of diffracted unit vector along incident beam unit vector */
    cos2theta = dot_product(incident,diffracted);
    cos2theta_sqr = cos2theta*cos2theta;
    sin2theta_sqr = 1-cos2theta_sqr;

    if(kahn_factor != 0.0){
        /* tricky bit here is deciding which direciton the E-vector lies in for each source
           here we assume it is closest to the "axis" defined above */

        /* cross product to get "vertical" axis that is orthogonal to the cannonical "polarization" */
        cross_product(axis,incident,B_in);
        /* make it a unit vector */
        unitize(B_in,B_in);

        /* cross product with incident beam to get E-vector direction */
        cross_product(incident,B_in,E_in);
        /* make it a unit vector */
        unitize(E_in,E_in);

        /* get components of diffracted ray projected onto the E-B plane */
        E_out[0] = dot_product(diffracted,E_in);
        B_out[0] = dot_product(diffracted,B_in);

        /* compute the angle of the diffracted ray projected onto the incident E-B plane */
        psi = -atan2(B_out[0],E_out[0]);
    }

    /* correction for polarized incident beam */
    return 0.5*(1.0 + cos2theta_sqr - kahn_factor*cos(2*psi)*sin2theta_sqr);
}


char *get_byte_order()
{
    static char *byte_order;

    typedef union
    {
        unsigned char string[2];
        unsigned short integer;
    } TWOBYTES;
    TWOBYTES twobytes;
    twobytes.integer = 24954;


    /* determine byte order on this machine */
    if(0==strncmp((const char *) twobytes.string, "az", 2))
    {
        byte_order = "big_endian";
    }
    else
    {
        byte_order = "little_endian";
    }
    return byte_order;
}


SMVinfo GetFrame(char *filename)
{
    char *string;
    SMVinfo frame;
    char *byte_order = get_byte_order();
    unsigned short int tempint;

    /* try to open the file... */
    frame.handle = fopen(filename, "rb");
    if(frame.handle != NULL)
    {
        /* just assume header will be 512 bytes?... */
        frame.header = (char *) calloc(1024,sizeof(char));
        if(! fread(frame.header, 512, 1, frame.handle))
        {
            perror("SMV file header");
            exit(9);
        }
        string = frame.header + 512;
        *string = (char) 0;

        /* remember the file name */
        frame.filename = (char *) calloc(strlen(filename)+10,sizeof(char));
        strcpy(frame.filename,filename);

        /* What kind of file is this? */
        if(0!=strncmp(frame.header, "{\nHEADER_BYTES=  512;\nDIM=2;\nBYTE_ORDER=", 12))
        {
            /* probably not an ADSC frame */

            /* inform the user */
            printf("ERROR: %s does not look like an ADSC frame!\n", filename);
            /* skip this file */
            fclose(frame.handle);

            frame.handle = NULL;
        }
        else
        {
            /* store the full header */
            frame.header_size = (int) ValueOf("HEADER_BYTES",frame);
            if(frame.header_size != 512)
            {
                free(frame.header);
                fseek(frame.handle,0,SEEK_SET);
                frame.header = (char *) calloc(2*frame.header_size,sizeof(char));
                if(! fread(frame.header, frame.header_size, 1, frame.handle))
                {
                    perror("SMV file fread");
                    exit(9);
                }
                string = frame.header + frame.header_size;
                *string = (char) 0;
            }

            /* see if we will need to swap bytes */
            string = (char *) strstr(frame.header, "BYTE_ORDER=")+11;
            /* find last instance of keyword in the header */
            while ((char *) strstr(string, "BYTE_ORDER=") != NULL)
            {
                string = (char *) strstr(string, "BYTE_ORDER=")+11;
            }
            if(0==strncmp(byte_order, string, 10))
            {
                frame.swap_bytes = FALSE;
            }
            else
            {
                frame.swap_bytes = TRUE;
            }

            /* store a couple of things */
            frame.width  = (int) ValueOf("SIZE1",frame);
            frame.height = (int) ValueOf("SIZE2",frame);

            if(frame.width == 0)
            {
                /* try other formats? */
                frame.width = frame.height = (int) ValueOf("DETECTOR_DIMENSIONS",frame);
            }

//          frame.mmapdata = mmap(NULL,2*frame.width*frame.height+frame.header_size,PROT_READ,MAP_SHARED,fileno(frame.handle),0);
            frame.mmapdata = (unsigned short int *) calloc(2,frame.width*frame.height+frame.header_size);
            if(frame.mmapdata == NULL)
            {
                perror("calloc:");
            }
            fseek(frame.handle,0,SEEK_SET);
            printf("reading %s\n",frame.filename);
            if(! fread(frame.mmapdata,1,2*frame.width*frame.height+frame.header_size,frame.handle))
            {
                perror("SMV file fread");
                exit(9);
            }

            printf("mmap(%s) = %p\n",frame.filename,frame.mmapdata);


        }
    }
    else
    {
        /* fopen() failed */
        perror(filename);
        frame.header_size=0;
    }

    return frame;
}

/* read floating-point values from keywords in an SMV header */
double ValueOf(const char *keyword, SMVinfo frame)
{
    double value;
    char *string;
    int keylen = strlen(keyword);

    /* start at the beginning */
    string = frame.header;

    /* find first instance of keyword in the header */
//    string = (char *) strstr(frame.header, keyword);
//    string = string + keylen;
    /* find last instance of keyword in the header */
    while ((char *) strstr(string, keyword) != NULL)
    {
        string = (char *) strstr(string, keyword)+keylen;
    }
    if(string == frame.header) return NAN;

    /* advance to just after the "=" sign */
    string = (char *) strstr(string, "=");
    if(string == NULL) return 0.0;
    ++string;

    value = atof(string);

    return value;
}


unsigned char *read_pgm5_bytes(char *filename,unsigned int *returned_width,unsigned int *returned_height)
{
    unsigned char test[512];
    unsigned char *array = NULL;
    FILE *handle = NULL;
    unsigned int width=0,height=0,maxvalue=0;

    handle = fopen(filename,"rb");
    if(handle)
    {
        if(! fread(test,512,1,handle))
        {
            perror("PGM fread header");
            exit(9);
        }
        if(strstr((const char *) test,"P5"))
        {
            /* PGM header: "P5<whitespace>width<whitespace>height<whitespace>maxvalue<single whitespace character>" */
            fseek(handle,3,SEEK_SET);
            if(! fscanf(handle," %u %u %u",&width,&height,&maxvalue))
            {
                perror("PGM fscanf");
                exit(9);
            }
            /* skip final single whitespsce character (first pixel could have value of "20") */
            fseek(handle,1,SEEK_CUR);
            array = (unsigned char *) calloc(sizeof(unsigned char),width*height);
            if(! fread(array,width,height,handle))
            {
                perror("PGM fread");
                exit(9);
            }
        }
        fclose(handle);
    }
    else
    {
        perror("PGM fopen");
    }

    *returned_width = width;
    *returned_height = height;
    return array;
}
