#include <simtbx/nanoBragg/nanoBragg.h>


//Contributed by James Holton, LBNL.

namespace simtbx {
namespace nanoBragg {


/* constructor that takes a DIALS detector model */
nanoBragg::nanoBragg(const dxtbx::model::Detector& detector)
{
    /* initialize to sensible default values */
    init_defaults();

    std::cout << detector << std::endl;
    spixels = detector[0].get_image_size()[1];
    fpixels = detector[0].get_image_size()[0];

    /* sensible initialization of all unititialized values */
    reconcile_parameters();
}


// constructor for the nanoBragg class that takes most any member as an argument, defaults in nanoBragg_etc.cpp
nanoBragg::nanoBragg(
        scitbx::vec2<int> detpixels_slowfast, // = 1024, 1024
        scitbx::vec3<int> Ncells_abc, // 1 1 1
        cctbx::uctbx::unit_cell unitcell, // lysozyme
        vec3 missets_deg, // 0 0 0
        vec2 beam_center_mm, // NAN NAN
        double distance_mm, // =100,
        double pixel_size_mm, // =0.1,
        double wavelength_A, // =1,
        double divergence_mrad, // =0,
        double dispersion_pct, // =0,
        double mosaic_spread_deg, // =-1
        int oversample, // =0 =auto
        int verbose)    // =0
{
    /* initialize to sensible default values */
    init_defaults();

    /* now incorporate arguments to the constructor */
    this->spixels=detpixels_slowfast[0];
    this->fpixels=detpixels_slowfast[1];
    this->Na=Ncells_abc[0]; this->Nb=Ncells_abc[1]; this->Nc=Ncells_abc[2];
    this->a_A[0]=unitcell.parameters()[0];      // units = Angstrom
    this->b_A[0]=unitcell.parameters()[1];
    this->c_A[0]=unitcell.parameters()[2];
    this->alpha=unitcell.parameters()[3]/RTD;   // input units = Degrees
    this->beta=unitcell.parameters()[4]/RTD;
    this->gamma=unitcell.parameters()[5]/RTD;
    this->user_cell=1;
    this->misset[1]=missets_deg[0];
    this->misset[2]=missets_deg[1];
    this->misset[3]=missets_deg[2];
    if(this->misset[1] != 0.0 && this->misset[2] != 0.0 && this->misset[3] != 0.0) this->misset[0]=1;
    if(! isnan(beam_center_mm[0]) && ! isnan(beam_center_mm[1])) {
        this->Xbeam = beam_center_mm[0]/1000.0;
        this->Ybeam = beam_center_mm[1]/1000.0;
    }
    this->distance=distance_mm/1000.0;
    this->pixel_size=pixel_size_mm/1000.0;
    this->lambda0=wavelength_A/1e10;
    this->hdivrange=divergence_mrad/1000.0;
    this->vdivrange=divergence_mrad/1000.0;
    this->dispersion=dispersion_pct/100.0;
    this->mosaic_spread=mosaic_spread_deg/RTD;
    this->oversample=oversample;
    if(oversample>0) this->user_oversample=true;
    this->verbose=verbose;

    /* sensible initialization of all unititialized values */
    reconcile_parameters();
}


/* sensible defaults */
void
nanoBragg::reconcile_parameters()
{
    /* now that constructor arguments are applied, fill in blanks */

    /* allocate detector pixel array */
    init_detector();

    /* initialize fluence from flux and check sample size fits inside beam size */
    init_beam();

    /* initialize interpretation of detector and beam geometry */
    init_beamcenter();

    /* sanitize options for steps across divergence, dispersion, mosaic spread and spindle rotation */
    init_steps();

    /* reconcile different conventions of beam center input and detector position */
    update_beamcenter();

    /* set up interpolation if need be */
    init_interpolator();

    /* reconcile different types of input of cell or crystal orientation */
    init_cell();

    /* decide if oversampling is needed based on crystal volume and pixel size */
    update_oversample();

    /* read in structure factors */
    init_Fhkl();

    /* read in centrosymmetric background profile (not implemented) */
    init_background();

    /* display all actual phi values */
    show_phisteps();

    /* read in or generate x-ray source properties */
    init_sources();

    /* allocate and set up mosaic domains */
    init_mosaicity();

    /* tally up the total steps/pixel */
    update_steps();

    if(verbose) printf("CONSTRUCT!!! %d   x-ray beam: %f %f %f\n",(int) raw.size(),this->beam_vector[1],this->beam_vector[2],this->beam_vector[3]);
}


/* sensible defaults */
void
nanoBragg::init_defaults()
{
    double test = NAN;
    if(! isnan(test)) {
        if(verbose) printf("ERROR! isnan(NAN) != %g\n",test);
        exit(9);
    };

    /* optional file stuff, to be removed eventually? */
    matfilename = NULL;
    hklfilename = NULL;
    stolfilename = NULL;
    imginfilename = NULL;
    maskfilename = NULL;
//    stoloutfilename = "output.stol\0";
    sourcefilename = NULL;
//    floatfilename = "floatimage.bin\0";
//    intfilename = "intimage.img\0";
//    pgmfilename = "image.pgm\0";
//    noisefilename = "noiseimage.img\0";
    infile = NULL;
//    FILE *outfile;  // = NULL;
//    FILE *stoloutfile;  // = NULL;

    /* default progress meter stuff */
    progress_pixel=progress_pixels=0;
    progress_meter = 1;
    babble = 1;
    printout = 0;
    printout_spixel=printout_fpixel=-1;
    verbose = 1;

    /* default x-ray beam properties */
//    double beam_vector[4] = {0,1,0,0};  this->beam_vector = beam_vector;
    beam_vector[0] = 0;
    beam_vector[1] = 1;
    beam_vector[2] = 0;
    beam_vector[3] = 0;
    coherent = 0;
    far_source = 1;
    round_div = 1;
    lambda=1e-10;
    lambda_of = NULL;
    mosaic_spread=-1.0;
    mosaic_umats = NULL;
    dispersion=0.0;dispstep=-1.0;
    lambda0 = 1.0e-10;
    hdiv=0;hdivstep=-1.0; hdivrange= -1.0;
    vdiv=0;vdivstep=-1.0; vdivrange= -1.0;
    source_path=0;source_distance = 10.0;
    hdiv_tic=vdiv_tic=disp_tic=mos_tic=0;
    divsteps=hdivsteps=vdivsteps=dispsteps=-1;
    mosaic_domains=-1;
    weight=0;
    source=sources=0;
    source_X=source_Y=source_Z=source_I=source_lambda=NULL;

    /* Thomson cross section (m^2) */
//    r_e_sqr = 7.94079248018965e-30;
    /* default incident x-ray fluence in photons/m^2 - exactly cancels Thomson cross section */
    fluence = 125932015286227086360700780544.0;
    flux=0.0;exposure=1.0;beamsize=1e-4;

    /* default sample size stuff */
    N=1;
    Na=Nb=Nc=1.0;
    xtalsize_max=xtalsize_a=xtalsize_b=xtalsize_c=0.0;
    reciprocal_pixel_size=0.0;

    round_xtal = 0;
    sample_x   = 0;             /* m */
    sample_y   = 0;             /* m */
    sample_z   = 0;             /* m */
    density    = 1.0e6;         /* g/m^3 */
    molecular_weight = 18.0;    /* g/mol */
    volume=0.0;molecules = 0.0;
    /* scale factor = F^2*r_e_sqr*fluence*Avogadro*volume*density/molecular_weight
                           m^2     ph/m^2  /mol      m^3   g/m^3    g/mol   */

    /* optional rudimentary background, better to use nonBragg */
    water_size = 0.0;
    water_F = 2.57;
    water_MW = 18.0;
    /* water F = 2.57 in forward direction */

    /* default detector stuff */
    pixel_size = 0.1e-3;
    fpixels=spixels=0;
    distance = 100.0e-3;
    detsize_f = 102.4e-3;
    detsize_s = 102.4e-3;
//    double fdet_vector[4]  = {0,0,0,1};  this->fdet_vector = fdet_vector;
    fdet_vector[0] = 0;
    fdet_vector[1] = 0;
    fdet_vector[2] = 0;
    fdet_vector[3] = 1;
//    double sdet_vector[4]  = {0,0,-1,0}; this->sdet_vector  = sdet_vector;
    sdet_vector[0] = 0;
    sdet_vector[1] = 0;
    sdet_vector[2] = -1;
    sdet_vector[3] = 0;
//    double odet_vector[4]  = {0,1,0,0};  this->odet_vector  = odet_vector;
    odet_vector[0] = 0;
    odet_vector[1] = 1;
    odet_vector[2] = 0;
    odet_vector[3] = 0;
//    double pix0_vector[4]  = {0,0,0,0};  this->pix0_vector  = pix0_vector;
    pix0_vector[0] = 0;
    pix0_vector[1] = 0;
    pix0_vector[2] = 0;
    pix0_vector[3] = 0;
    detector_rotx=0.0;detector_roty=0.0;detector_rotz=0.0;
//    double twotheta_axis[4] = {0,0,1,0}; this->twotheta_axis = twotheta_axis;
    twotheta_axis[0] = 0;
    twotheta_axis[1] = 0;
    twotheta_axis[2] = 1;
    twotheta_axis[3] = 0;
    detector_pivot = BEAM;
    beam_convention = MOSFLM;
    detector_twotheta = 0.0;
    airpath=omega_pixel=omega_Rsqr_pixel=omega_sum = 0.0;
    curved_detector = 0;
    point_pixel= 0;
    Xbeam=NAN;Ybeam=NAN;
    Fbeam=NAN;Sbeam=NAN;
    Fdet=Sdet=Rdet=Fdet0=Sdet0=0.0;
    Xclose=NAN;Yclose=NAN;close_distance=NAN;
    Fclose=NAN;Sclose=NAN;
    ORGX=NAN;ORGY=NAN;
    adc_offset = 40.0;

    /* scattering vectors */
//    double incident[4] = {0,0,0,0};
//    double diffracted[4],diffracted0[4] = {0,0,0,0};
//    double scattering[4];
    stol=twotheta=theta=0.0;

    /* diffraction geometry stuff */
    costwotheta=sintwotheta=psi=0;
    xd=yd=zd=xd0=yd0=zd0=0.0;
//    double Ewald[4],Ewald0[4],relp[4];
    dmin=0;
    integral_form = 0; // experimental: use integral form of spots to avoid need for oversampling

    /* polarization stuff */
//    double polar_vector[4] = {0,0,0,1}; this->polar_vector = polar_vector;
    polar_vector[0] = 0;
    polar_vector[1] = 0;
    polar_vector[2] = 0;
    polar_vector[3] = 1;
//  vert_vector[4];
    polar=1.0;
    polarization=0.0;
    nopolar = 0;

    /* sampling */
    steps=0;
    roi_xmin=-1;roi_xmax=-1;roi_ymin=-1;roi_ymax=-1;
    oversample = -1;
    recommended_oversample=subS=subF=0;
    user_oversample=false;
    subpixel_size=0.0;

    /* spindle */
    phi0=0.0;phistep=-1.0;osc=-1.0;
    phisteps=-1;
//    double spindle_vector[4] = {0,0,0,1}; this->spindle_vector = spindle_vector;
    spindle_vector[0] = 0;
    spindle_vector[1] = 0;
    spindle_vector[2] = 0;
    spindle_vector[3] = 1;

    /* structure factor representation */
    Fhkl = NULL;
    default_F = 0.0;
    nearest=0;
    stol_file_mult=1.0e10;
    pythony_indices.clear();
    pythony_amplitudes.clear();


    /* intensity stats */
    max_I = 0.0;
    max_I_x = 0.0; max_I_y = 0.0;
    photon_scale = 0.0;
    intfile_scale = 0.0;
    pgm_scale = 0.0;
    sumn = 0;
    overloads = 0;

    /* image file data */
    floatimage = NULL;
    maskimage = NULL;
    intimage = NULL;
    pgmimage = NULL;
    imginfileimage = NULL;

    /* random number seeds */
    seed = -time((time_t *)0);
//    if(verbose) printf("random number seed = %u\n",seed);
    mosaic_seed = 12345678;
    calib_seed = 123456789;

    /* interpolation defaults: auto-detect */
    interpolate = 2;
    i1=0;i2=0;i3=0;
    sub_Fhkl=NULL;

    /* unit cell stuff */
    user_cell = 0;
//    double a[4] = {0,0,0,0}; this->a = a;
//    double b[4] = {0,0,0,0}; this->b = b;
//    double c[4] = {0,0,0,0}; this->c = c;
    a_A[0]=0;a_A[1]=0;a_A[2]=0;a_A[3]=0;
    b_A[0]=0;b_A[1]=0;b_A[2]=0;b_A[3]=0;
    c_A[0]=0;c_A[1]=0;c_A[2]=0;c_A[3]=0;
    a[0]=0;a[1]=0;a[2]=0;a[3]=0;
    b[0]=0;b[1]=0;b[2]=0;b[3]=0;
    c[0]=0;c[1]=0;c[2]=0;c[3]=0;
    alpha=0.0;beta=0.0;gamma=0.0;
//    double misset[4] = {0,0,0,0}; this->misset = misset;
    misset[0] = 0;
    misset[1] = 0;
    misset[2] = 0;
    misset[3] = 0;


    /* special options */
    calculate_noise = 1;
    write_pgm = 1;
    binary_spots = false;
    }

/* initialize detector size from mm or pixel specs */
void
nanoBragg::init_detector()
    {
    /* fill in blanks */
    if(fpixels) {
            detsize_f = pixel_size*fpixels;
            }
    if(spixels) {
            detsize_s = pixel_size*spixels;
            }
    this->fpixels = ceil(detsize_f/pixel_size);
    this->spixels = ceil(detsize_s/pixel_size);
    pixels = this->fpixels*this->spixels;
    /* actually allocate memory */
    if(pixels) {
        af::versa<double, af::c_grid<2> > raw2(af::c_grid<2> (spixels,fpixels));
        this->raw=raw2;
        }
    }
// end of init_detector


/* initialize fluence from flux and check sample/beam size */
void
nanoBragg::init_beam()
            {
    /* get fluence from flux */
    if(flux != 0.0 && exposure > 0.0 && beamsize >= 0){
        fluence = flux*exposure/beamsize/beamsize;
    }
    if(beamsize >= 0){
        if(beamsize < sample_y){
            if(verbose) printf("WARNING: clipping sample (%lg m high) with beam (%lg m)\n",sample_y,beamsize);
            sample_y = beamsize;
        }
        if(beamsize < sample_z){
            if(verbose) printf("WARNING: clipping sample (%lg m wide) with beam (%lg m)\n",sample_z,beamsize);
            sample_z = beamsize;
        }
    }
    if(exposure > 0.0)
    {
        flux = fluence/exposure*beamsize*beamsize;
    }
    /* straighten up sample properties */
//    volume = sample_x*sample_y*sample_z;
//    molecules = volume*density*Avogadro/molecular_weight;
}
// end of init_beam


/* set up any uninitialized beam centers from others provided */
void
nanoBragg::init_beamcenter()
{
    /* default to center of detector */
    if(! isnan(ORGX)) Xclose = ORGX*pixel_size;
    if(! isnan(ORGY)) Yclose = ORGY*pixel_size;
    if(isnan(Xclose)) Xclose = (detsize_f - pixel_size)/2.0;
    if(isnan(Yclose)) Yclose = (detsize_s + pixel_size)/2.0;
//    if(isnan(Xbeam)) Xbeam = Xclose;
//    if(isnan(Ybeam)) Ybeam = Yclose;
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
        if(verbose) printf("beam center in adxv convention\n");
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
        if(verbose) printf("beam center in mosflm convention\n");
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
        if(verbose) printf("beam center in denzo convention\n");
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
        if(verbose) printf("beam center in XDS convention\n");
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
        if(verbose) printf("beam center in DIALS convention\n");
    }
    if(beam_convention == CUSTOM)
    {
        if(isnan(Xbeam)) Xbeam = Xclose;
        if(isnan(Ybeam)) Ybeam = Yclose;
        Fbeam = Xbeam;
        Sbeam = Ybeam;
        Fclose = Xbeam;
        Sclose = Ybeam;
        if(verbose) printf("beam center in custom convention\n");
    }

    /* straighten up vectors */
    unitize(beam_vector,beam_vector);
    unitize(fdet_vector,fdet_vector);
    unitize(sdet_vector,sdet_vector);
    if(unitize(odet_vector,odet_vector) != 1.0)
    {
        if(verbose) printf("WARNING: auto-generating odet_vector\n");
        cross_product(fdet_vector,sdet_vector,odet_vector);
        unitize(odet_vector,odet_vector);
    }
    unitize(polar_vector,polar_vector);
    unitize(spindle_vector,spindle_vector);
    cross_product(beam_vector,polar_vector,vert_vector);
    unitize(vert_vector,vert_vector);

    }
// init_beamcenter



/* count up and sanitize steps across divergence, dispersion, mosaic spread and spindle rotation */
void
nanoBragg::init_steps()
{
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
                if(verbose) printf("WARNING: finite mosaicity with only one domain! upping to 10 mosaic domains\n");
            mosaic_domains = 10;
        }
        }
    } else {
        /* user-specified number of domains */
        if(mosaic_spread < 0.0) {
            /* number of domains specified, but no spread? */
            if(verbose) printf("WARNING: no mosaic spread specified.  setting mosaic_domains = 1\n");
            mosaic_spread = 0.0;
            mosaic_domains = 1;
        } else {
            /* user-speficied mosaicity and number of domains */
            if(mosaic_spread == 0.0)
            {
                if(verbose) printf("WARNING: zero mosaic spread specified.  setting mosaic_domains = 1\n");
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
}
// end of init_steps()




/* reconcile different conventions of beam center input and detector position */
void
nanoBragg::update_beamcenter()
{
    /* initialize detector origin from a beam center and distance */
    /* there are different conventions here: mosflm, XDS etc */

    if(beam_convention == ADXV)
    {
        if(verbose) printf("adxv");
    }
    if(beam_convention == MOSFLM)
    {
        if(verbose) printf("mosflm");
    }
    if(beam_convention == XDS)
    {
        if(verbose) printf("xds");
    }
    if(beam_convention == DENZO)
    {
        if(verbose) printf("denzo");
    }
    if(beam_convention == CUSTOM)
    {
        if(verbose) printf("custom");
    }
    if(verbose) printf(" convention selected.\n");

    /* first off, what is the relationship between the two "beam centers"? */
    rotate(odet_vector,vector,detector_rotx,detector_roty,detector_rotz);
    ratio = dot_product(beam_vector,vector);
    if(ratio == 0.0) { ratio = DBL_MIN; }
    if(isnan(close_distance)) close_distance = fabs(ratio*distance);
    distance = close_distance/ratio;

    if(detector_pivot == SAMPLE){
        if(verbose) printf("pivoting detector around sample\n");
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
        if(verbose) printf("pivoting detector around direct beam spot\n");
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
    newvector[1]=+1;newvector[2]=+0;newvector[3]=+0;
    dials_origin[3] = 1000.0*dot_product(pix0_vector,newvector);

}
// end of update_beamcenter()




/* automatically decide if we are interpolating, allocate memory if yes */
void
nanoBragg::init_interpolator()
{
    /* free any previous allocations */
    if(sub_Fhkl != NULL) {
        for (h0=0; h0<=5;h0++) {
            for (k0=0; k0<=5;k0++) {
                if(verbose>6) printf("freeing %d %ld-byte double Fhkl[%d][%d] at %p\n",5,sizeof(double),h0,k0,sub_Fhkl[h0][k0]);
                free(*(*(sub_Fhkl +h0)+k0));
            }
            if(verbose>6) printf("freeing %d %ld-byte double* Fhkl[%d] at %p\n",5,sizeof(double*),h0,sub_Fhkl[h0]);
            free(*(sub_Fhkl +h0));
        }
        if(verbose>6) printf("freeing %d %ld-byte double** Fhkl at %p\n",5,sizeof(double**),sub_Fhkl);
        free(sub_Fhkl);
    }

    if(interpolate > 1){
        /* no user options */
        if(( Na <= 2) || (Nb <= 2) || (Nc <= 2)){
            if(verbose) printf("auto-selected tricubic interpolation of structure factors\n");
            interpolate = 1;
        }
        else
        {
            if(verbose) printf("auto-selected no interpolation\n");
            interpolate = 0;
        }
    }
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
}
// end of init_interpolator()




/* compute A matrix from cell, missets, mosflm matrix file, etc. */
void
nanoBragg::init_cell()
{
    /* do not run if things are not going to work */
    if(! user_cell && matfilename == NULL)
    {
        if(verbose) printf("ERROR: cannot initialize without a cell\n");
        return;
    }
    if(verbose>1) printf("raw CELL %f %f %f %f %f %f   %d\n",a_A[0],b_A[0],c_A[0],alpha,beta,gamma,user_cell);
    /* user-specified unit cell */
    if(user_cell)
    {
        /* a few random defaults */
        if(b_A[0]  <= 0.0) b_A[0] = a_A[0];
        if(c_A[0]  <= 0.0) c_A[0] = a_A[0];
        if(alpha <= 0.0) alpha = M_PI/2;
        if(beta  <= 0.0) beta  = M_PI/2;
        if(gamma <= 0.0) gamma = M_PI/2;

        /* get cell volume from angles */
        aavg = (alpha+beta+gamma)/2;
        skew = sin(aavg)*sin(aavg-alpha)*sin(aavg-beta)*sin(aavg-gamma);
        if(skew<0.0) skew=-skew;
        V_cell = 2.0*a_A[0]*b_A[0]*c_A[0]*sqrt(skew);
        if(V_cell <= 0.0)
        {
            if(verbose) printf("WARNING: impossible unit cell volume: %g\n",V_cell);
            V_cell = DBL_MIN;
        }
        V_star = 1.0/V_cell;

        /* now get reciprocal-cell lengths from the angles and volume */
        a_star[0] = b_A[0]*c_A[0]*sin(alpha)*V_star;
        b_star[0] = c_A[0]*a_A[0]*sin(beta)*V_star;
        c_star[0] = a_A[0]*b_A[0]*sin(gamma)*V_star;
        if(a_star[0] <= 0.0 || b_star[0] <= 0.0 || c_star[0] <= 0.0)
        {
            if(verbose) printf("WARNING: impossible reciprocal cell lengths: %g %g %g\n",
                a_star[0],b_star[0],c_star[0]);
            /* make up something non-zero */
            a_star[0] = fabs(a_star[0]);
            b_star[0] = fabs(b_star[0]);
            c_star[0] = fabs(c_star[0]);
            if(a_star[0] <= 0.0) a_star[0] = DBL_MIN;
            if(b_star[0] <= 0.0) b_star[0] = DBL_MIN;
            if(c_star[0] <= 0.0) c_star[0] = DBL_MIN;
        }

        /* for fun, compute the reciprocal-cell angles from direct-cell angles */
        sin_alpha_star = a_A[0]*V_star/b_star[0]/c_star[0];
        sin_beta_star  = b_A[0]*V_star/a_star[0]/c_star[0];
        sin_gamma_star = c_A[0]*V_star/a_star[0]/b_star[0];
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
            if(verbose) printf("WARNING: oddball reciprocal cell angles:\n");
            if(verbose) printf("sin(alpha_star) = %.25g\n",sin_alpha_star);
            if(verbose) printf("cos(alpha_star) = %.25g\n",cos_alpha_star);
            if(verbose) printf("sin(beta_star)  = %.25g\n",sin_beta_star);
            if(verbose) printf("cos(beta_star)  = %.25g\n",cos_beta_star);
            if(verbose) printf("sin(gamma_star) = %.25g\n",sin_gamma_star);
            if(verbose) printf("cos(gamma_star) = %.25g\n",cos_gamma_star);
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

        if(verbose>1) printf("CELL %f %f %f %f %f %f   %d\n",a_A[0],b_A[0],c_A[0],alpha,beta,gamma,user_cell);

        /* construct default orientation */
        a_star[1] = a_star[0];
        b_star[1] = b_star[0]*cos_gamma_star;
        c_star[1] = c_star[0]*cos_beta_star;
        a_star[2] = 0.0;
        b_star[2] = b_star[0]*sin_gamma_star;
        c_star[2] = c_star[0]*(cos_alpha_star-cos_beta_star*cos_gamma_star)/sin_gamma_star;
        a_star[3] = 0.0;
        b_star[3] = 0.0;
        c_star[3] = c_star[0]*V_cell/(a_A[0]*b_A[0]*c_A[0]*sin_gamma_star);
    }

    /* load the lattice orientation (reciprocal cell vectors) from a mosflm matrix */
    if(matfilename != NULL)
    {
        infile = fopen(matfilename,"r");
        if(infile != NULL)
        {
            if(verbose) printf("reading %s\n",matfilename);
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

    if(verbose>1) printf("MISSET %f %f %f \n",misset[1],misset[2],misset[3]);
    /* apply any missetting angle */
    if(misset[0] != 0.0)
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
        vector_rescale(b_star_cross_c_star,b_star_cross_c_star,a_A[0]/V_cell);
        vector_rescale(c_star_cross_a_star,c_star_cross_a_star,b_A[0]/V_cell);
        vector_rescale(a_star_cross_b_star,a_star_cross_b_star,c_A[0]/V_cell);
        V_star = 1.0/V_cell;
    }

    /* direct-space cell volume in A^3 */
    V_cell = 1.0/V_star;

    /* generate direct-space cell vectors, also updates magnitudes */
    vector_scale(b_star_cross_c_star,a_A,V_cell);
    vector_scale(c_star_cross_a_star,b_A,V_cell);
    vector_scale(a_star_cross_b_star,c_A,V_cell);

    /* now that we have direct-space vectors, re-generate the reciprocal ones */
    cross_product(a_A,b_A,a_cross_b);
    cross_product(b_A,c_A,b_cross_c);
    cross_product(c_A,a_A,c_cross_a);
    vector_scale(b_cross_c,a_star,V_star);
    vector_scale(c_cross_a,b_star,V_star);
    vector_scale(a_cross_b,c_star,V_star);

    /* for fun, calculate the cell angles too */
    sin_alpha = a_star[0]*V_cell/b_A[0]/c_A[0];
    sin_beta  = b_star[0]*V_cell/a_A[0]/c_A[0];
    sin_gamma = c_star[0]*V_cell/a_A[0]/b_A[0];
    cos_alpha = dot_product(b,c)/b_A[0]/c_A[0];
    cos_beta  = dot_product(a,c)/a_A[0]/c_A[0];
    cos_gamma = dot_product(a,b)/a_A[0]/b_A[0];
    if(sin_alpha>1.0000001 || sin_alpha<-1.0000001 ||
       sin_beta >1.0000001 || sin_beta <-1.0000001 ||
       sin_gamma>1.0000001 || sin_gamma<-1.0000001 ||
       cos_alpha>1.0000001 || cos_alpha<-1.0000001 ||
       cos_beta >1.0000001 || cos_beta <-1.0000001 ||
       cos_gamma>1.0000001 || cos_gamma<-1.0000001 )
    {
        if(verbose) printf("WARNING: oddball cell angles:\n");
            if(verbose) printf("sin_alpha = %.25g\n",sin_alpha);
            if(verbose) printf("cos_alpha = %.25g\n",cos_alpha);
            if(verbose) printf("sin_beta  = %.25g\n",sin_beta);
            if(verbose) printf("cos_beta  = %.25g\n",cos_beta);
            if(verbose) printf("sin_gamma = %.25g\n",sin_gamma);
            if(verbose) printf("cos_gamma = %.25g\n",cos_gamma);
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
    sin_alpha_star = a_A[0]*V_star/b_star[0]/c_star[0];
    sin_beta_star  = b_A[0]*V_star/a_star[0]/c_star[0];
    sin_gamma_star = c_A[0]*V_star/a_star[0]/b_star[0];
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
            if(verbose) printf("WARNING: oddball reciprocal cell angles:\n");
            if(verbose) printf("sin(alpha_star) = %.25g\n",sin_alpha_star);
            if(verbose) printf("cos(alpha_star) = %.25g\n",cos_alpha_star);
            if(verbose) printf("sin(beta_star)  = %.25g\n",sin_beta_star);
            if(verbose) printf("cos(beta_star)  = %.25g\n",cos_beta_star);
            if(verbose) printf("sin(gamma_star) = %.25g\n",sin_gamma_star);
            if(verbose) printf("cos(gamma_star) = %.25g\n",cos_gamma_star);
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

    if(verbose) printf("Unit Cell: %g %g %g %g %g %g\n", a_A[0],b_A[0],c_A[0],alpha*RTD,beta*RTD,gamma*RTD);
    if(verbose) printf("Recp Cell: %g %g %g %g %g %g\n", a_star[0],b_star[0],c_star[0],alpha_star*RTD,beta_star*RTD,gamma_star*RTD);
    if(verbose) printf("volume = %g A^3\n",V_cell);

    /* print out the real-space matrix */
    if(verbose) printf("real-space cell vectors (Angstrom):\n");
    if(verbose) printf("     %-10s  %-10s  %-10s\n","a","b","c");
    if(verbose) printf("X: %11.8f %11.8f %11.8f\n",a_A[1],b_A[1],c_A[1]);
    if(verbose) printf("Y: %11.8f %11.8f %11.8f\n",a_A[2],b_A[2],c_A[2]);
    if(verbose) printf("Z: %11.8f %11.8f %11.8f\n",a_A[3],b_A[3],c_A[3]);
    if(verbose) printf("reciprocal-space cell vectors (Angstrom^-1):\n");
    if(verbose) printf("     %-10s  %-10s  %-10s\n","a_star","b_star","c_star");
    if(verbose) printf("X: %11.8f %11.8f %11.8f\n",a_star[1],b_star[1],c_star[1]);
    if(verbose) printf("Y: %11.8f %11.8f %11.8f\n",a_star[2],b_star[2],c_star[2]);
    if(verbose) printf("Z: %11.8f %11.8f %11.8f\n",a_star[3],b_star[3],c_star[3]);

    /* now convert Angstrom to meters */
    vector_scale(a_A,a,1e-10);
    vector_scale(b_A,b,1e-10);
    vector_scale(c_A,c,1e-10);

    /* define phi=0 mosaic=0 crystal orientation */
    vector_scale(a,a0,1.0);
    vector_scale(b,b0,1.0);
    vector_scale(c,c0,1.0);

    /* define phi=0 crystal orientation */
    vector_scale(a,ap,1.0);
    vector_scale(b,bp,1.0);
    vector_scale(c,cp,1.0);
}
// end of init_cell()


/* compute sample size and required oversampling */
void
nanoBragg::update_oversample()
{
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
    if(verbose) printf("crystal is %g x %g x %g microns\n",xtalsize_a*1e6,xtalsize_b*1e6,xtalsize_c*1e6);
    xtalsize_max = xtalsize_a;
    if(xtalsize_max < xtalsize_b) xtalsize_max = xtalsize_b;
    if(xtalsize_max < xtalsize_c) xtalsize_max = xtalsize_c;
    reciprocal_pixel_size = lambda0*distance/pixel_size;
    recommended_oversample = ceil(3.0 * xtalsize_max/reciprocal_pixel_size);
    if(recommended_oversample <= 0) recommended_oversample = 1;
    if(! user_oversample) {
        oversample = recommended_oversample;
        if(verbose) printf("auto-selected %d-fold oversampling\n",oversample);
    }
    if(oversample < recommended_oversample)
    {
        if(verbose)
        {
        printf("WARNING: maximum dimension of sample is %g A\n",xtalsize_max*1e10);
        printf("         but reciprocal pixel size is %g A\n", reciprocal_pixel_size*1e10 );
        printf("         intensity may vary significantly across a pixel!\n");
            printf("         recommend oversample=%d to work around this\n",recommended_oversample);
        }
    }

    /* rough estimate of sample properties */
    sample_x = xtalsize_a;
    sample_y = xtalsize_b;
    sample_z = xtalsize_c;
    volume = sample_x*sample_y*sample_z;
    density = 1.2e6;
    molecules = Na*Nb*Nc;
    molecular_weight = volume*density*Avogadro/molecules;
    if(verbose) printf("approximate MW = %g\n",molecular_weight);

}
// end of update_oversample()




/* read in structure factors vs hkl */
void
nanoBragg::init_Fhkl()
        {
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

    /* load the structure factors */
    if(hklfilename == NULL)
    {
        /* try to recover Fs from a previous run? */
        h_min = h_max = 0;
        k_min = k_max = 0;
        l_min = l_max = 0;
    }
    else
    {
        infile = fopen(hklfilename,"r");
        if(infile == NULL)
        {
            if(verbose) printf("ERROR: unable to open %s.",hklfilename);
            exit(9);
        }
        h_min=k_min=l_min=1e9;
        h_max=k_max=l_max=-1e9;
        if(verbose) printf("counting entries in %s\n",hklfilename);
        while(4 == fscanf(infile,"%lg%lg%lg%lg",&h,&k,&l,&F_cell)){
            if(verbose && h != ceil(h-0.4)) printf("WARNING: non-integer value for h (%g) at line %d\n",h,hkls);
            if(verbose && k != ceil(k-0.4)) printf("WARNING: non-integer value for k (%g) at line %d\n",k,hkls);
            if(verbose && l != ceil(l-0.4)) printf("WARNING: non-integer value for l (%g) at line %d\n",l,hkls);
            if(h_min > h) h_min = h;
            if(k_min > k) k_min = k;
            if(l_min > l) l_min = l;
            if(h_max < h) h_max = h;
            if(k_max < k) k_max = k;
            if(l_max < l) l_max = l;
            ++hkls;
        }
        rewind(infile);
    }

    if(pythony_indices.size())
    {
        if(verbose) printf(" noticed pythony_indices.size() = %ld\n",pythony_indices.size());
        /* need to know how much memory to allocate for Fhkl array */
        h_min=k_min=l_min=1e9;
        h_max=k_max=l_max=-1e9;
        miller_t hkl;
        for (i=0; i < pythony_indices.size(); ++i)
        {
            hkl = pythony_indices[i];
            if(pythony_amplitudes.size()) F_cell = pythony_amplitudes[i];
            h0 = hkl[0];
            k0 = hkl[1];
            l0 = hkl[2];
//          if(verbose) printf("GOTHERE %d : %d %d %d = %g\n",i,h0,k0,l0,F_cell);
            if(h_min > h0) h_min = h0;
            if(k_min > k0) k_min = k0;
            if(l_min > l0) l_min = l0;
            if(h_max < h0) h_max = h0;
            if(k_max < k0) k_max = k0;
            if(l_max < l0) l_max = l0;
            ++hkls;
        }
    }

    if(1)
    {
        h_range = h_max - h_min + 1;
        k_range = k_max - k_min + 1;
        l_range = l_max - l_min + 1;

            if(verbose) printf("h: %d - %d\n",h_min,h_max);
            if(verbose) printf("k: %d - %d\n",k_min,k_max);
            if(verbose) printf("l: %d - %d\n",l_min,l_max);
        if(h_range < 0 || k_range < 0 || l_range < 0) {
            if(verbose) printf("ERROR: not enough HKL indices in %s\n",hklfilename);
            exit(9);
        }

        /* allocate memory for 3d arrays */
        if(verbose>6) printf("allocating %d %ld-byte double**\n",h_range+1,sizeof(double**));
        Fhkl = (double***) calloc(h_range+1,sizeof(double**));
        if(Fhkl==NULL){perror("ERROR");exit(9);};
        for (h0=0; h0<=h_range;h0++) {
            if(verbose>6) printf("allocating %d %ld-byte double*\n",k_range+1,sizeof(double*));
                Fhkl[h0] = (double**) calloc(k_range+1,sizeof(double*));
                if(Fhkl[h0]==NULL){perror("ERROR");exit(9);};
                for (k0=0; k0<=k_range;k0++) {
                if(verbose>6) printf("allocating %d %ld-byte double\n",l_range+1,sizeof(double));
                        Fhkl[h0][k0] = (double*) calloc(l_range+1,sizeof(double));
                        if(Fhkl[h0][k0]==NULL){perror("ERROR");exit(9);};
                }
        }
        if(verbose) printf("initializing to default_F = %g:\n",default_F);
        for (h0=0; h0<h_range;h0++) {
            for (k0=0; k0<k_range;k0++) {
                for (l0=0; l0<l_range;l0++) {
                    Fhkl[h0][k0][l0] = default_F;
                }
            }
        }
        if(verbose) printf("done initializing:\n");
    }

    if(hklfilename != NULL)
    {
        if(verbose) printf("re-reading %s\n",hklfilename);
        while(4 == fscanf(infile,"%d%d%d%lg",&h0,&k0,&l0,&F_cell)){
            Fhkl[h0-h_min][k0-k_min][l0-l_min]=F_cell;
        }
        fclose(infile);
    }

    if(pythony_indices.size() && pythony_amplitudes.size())
        {
        if(verbose) printf("initializing Fhkl with pythony indices and amplitudes\n");
        miller_t hkl;
        for (i=0; i < pythony_indices.size(); ++i)
        {
            hkl = pythony_indices[i];
            F_cell = pythony_amplitudes[i];
            h0 = hkl[0];
            k0 = hkl[1];
            l0 = hkl[2];
            Fhkl[h0-h_min][k0-k_min][l0-l_min]=F_cell;
            if(verbose>6) printf("F %d : %d %d %d = %g\n",i,h0,k0,l0,F_cell);
        }
        if(verbose) printf("done initializing:\n");
        }
    }
// end of init_Fhkl





/* read in or generate background profile */
void
nanoBragg::init_background()
{
    /* now read in amorphous material structure factors */
    stols = 0;
    if(stolfilename != NULL)
    {
        if(verbose) printf("reading %s\n",stolfilename);
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
}


/* print out actual phi values in sweep */
void
nanoBragg::show_phisteps()
{
    /* show phi steps with sweep over spindle axis */
    for(phi_tic = 0; phi_tic < phisteps; ++phi_tic){
        phi = phi0 + phistep*phi_tic;
        if(verbose) printf("phi%d = %g\n",phi_tic,phi*RTD);
    }
    }


/* read in or generate x-ray sources */
void
nanoBragg::init_sources()
{
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
                if(verbose) printf("divergence deviation: %g %g\n",hdiv,vdiv);
            }
        }

        /* print out wavelength steps with sweep over spectral dispersion */
        for(disp_tic=0;disp_tic<dispsteps;++disp_tic){
            lambda = lambda0 * ( 1.0 + dispstep * disp_tic - dispersion/2.0 ) ;
            if(verbose) printf("lambda%d = %.15g\n",disp_tic,lambda);
        }

        /* free any previous allocation */
        if(source_X != NULL) free(source_X);
        if(source_Y != NULL) free(source_Y);
        if(source_Z != NULL) free(source_Z);
        if(source_I != NULL) free(source_I);
        if(source_lambda != NULL) free(source_lambda);

        /* allocate enough space */
        sources = divsteps*dispsteps;
        if(verbose>6) printf("allocating space for %d sources\n",sources);
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
    if(verbose) printf("  created a total of %d sources:\n",sources);
    for(source=0;source<sources;++source){

        /* retrieve stuff from cache */
        X = source_X[source];
        Y = source_Y[source];
        Z = source_Z[source];
        I = source_I[source];
        lambda = source_lambda[source];

        if(verbose) printf("%g %g %g   %g %.6g\n",X,Y,Z,I,lambda);
    }

    }
// end of init_sources()






/* read in? or generate mosaic domains */
void
nanoBragg::init_mosaicity()
{
    /* free any previous allocations */
    if(mosaic_umats!=NULL) free(mosaic_umats);

    /* allocate enough space */
    if(mosaic_domains<1) mosaic_domains=1;
    if(verbose>6) printf("allocating enough space for %d mosaic domain orientation matrices\n",mosaic_domains);
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
    }

    if(verbose) printf("  created a total of %d mosaic domains\n",mosaic_domains);
}
// end of init_mosaicity()




/* print out individual mosaic domain orientations */
void
nanoBragg::show_mosaic_blocks()
{
    /* assume init_mosaicity() was already run? */
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
        umat2misset(mosaic_umats+9*mos_tic,mosaic_missets);
        printf("%d by: %f %f %f deg\n",mos_tic,mosaic_missets[1]*RTD,mosaic_missets[2]*RTD,mosaic_missets[3]*RTD);
        printf("       %f %f %f\n",*(mosaic_umats+9*mos_tic+0),*(mosaic_umats+9*mos_tic+1),*(mosaic_umats+9*mos_tic+2));
        printf("       %f %f %f\n",*(mosaic_umats+9*mos_tic+3),*(mosaic_umats+9*mos_tic+4),*(mosaic_umats+9*mos_tic+5));
        printf("       %f %f %f\n",*(mosaic_umats+9*mos_tic+6),*(mosaic_umats+9*mos_tic+7),*(mosaic_umats+9*mos_tic+8));
    }

    printf("  total of %d mosaic domains\n",mosaic_domains);
    }
// end of show_mosaic_blocks()




/* calculate total number of steps/pixel */
void
nanoBragg::update_steps()
{
    /* final decisions about sampling */
    if(oversample <= 0) oversample = 1;
    steps = sources*mosaic_domains*oversample*oversample;
    subpixel_size = pixel_size/oversample;
}
// end of update_steps()



/* reconcile different conventions of beam center input and detector position */
void
nanoBragg::show_params()
{
    printf("nanoBragg nanocrystal diffraction simulator - James Holton and Ken Frankel 3-23-16\n");

    printf("  %d initialized hkls (all others =%g)\n",hkls,default_F);
    printf("  ");
    if(round_xtal){
        printf("ellipsoidal");
    }
    else
    {
        printf("parallelpiped");
    }
    printf(" xtal: %.0fx%.0fx%.0f cells\n",Na,Nb,Nc);
    printf("Unit Cell: %g %g %g %g %g %g\n", a_A[0],b_A[0],c_A[0],alpha*RTD,beta*RTD,gamma*RTD);
    printf("Recp Cell: %g %g %g %g %g %g\n", a_star[0],b_star[0],c_star[0],alpha_star*RTD,beta_star*RTD,gamma_star*RTD);
    printf("volume = %g A^3\n",V_cell);

    printf("missets: %11.8f %11.8f %11.8f\n",misset[1]*RTD,misset[2]*RTD,misset[3]*RTD);

    /* print out the real-space matrix */
    printf("real-space cell vectors (Angstrom):\n");
    printf("     %-10s  %-10s  %-10s\n","a","b","c");
    printf("X: %11.8f %11.8f %11.8f\n",a_A[1],b_A[1],c_A[1]);
    printf("Y: %11.8f %11.8f %11.8f\n",a_A[2],b_A[2],c_A[2]);
    printf("Z: %11.8f %11.8f %11.8f\n",a_A[3],b_A[3],c_A[3]);
    printf("reciprocal-space cell vectors (Angstrom^-1):\n");
    printf("     %-10s  %-10s  %-10s\n","a_star","b_star","c_star");
    printf("X: %11.8f %11.8f %11.8f\n",a_star[1],b_star[1],c_star[1]);
    printf("Y: %11.8f %11.8f %11.8f\n",a_star[2],b_star[2],c_star[2]);
    printf("Z: %11.8f %11.8f %11.8f\n",a_star[3],b_star[3],c_star[3]);
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
}
// end of show_params()


void
nanoBragg::sweep_over_detector()
{
    max_I = 0.0;
    int j = 0;
    double* floatimage(raw.begin());
//    floatimage = (double *) calloc(spixels*fpixels+10,sizeof(double));

    if(verbose) printf("TESTING sincg(1,1)= %f\n",sincg(1,1));

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

                    /* construct detector subpixel position in 3D space */
//                  pixel_X = distance;
//                  pixel_Y = Sdet-Ybeam;
//                  pixel_Z = Fdet-Xbeam;
                    pixel_pos[1] = Fdet*fdet_vector[1]+Sdet*sdet_vector[1]+pix0_vector[1];
                    pixel_pos[2] = Fdet*fdet_vector[2]+Sdet*sdet_vector[2]+pix0_vector[2];
                    pixel_pos[3] = Fdet*fdet_vector[3]+Sdet*sdet_vector[3]+pix0_vector[3];
                    pixel_pos[0] = 0.0;
                    if(curved_detector) {
                        /* construct detector pixel that is always "distance" from the sample */
                        vector[1] = distance*beam_vector[1]; vector[2]=distance*beam_vector[2] ; vector[3]=distance*beam_vector[3];
                        /* treat detector pixel coordinates as radians */
                        rotate_axis(vector,newvector,sdet_vector,pixel_pos[2]/distance);
                        rotate_axis(newvector,pixel_pos,fdet_vector,pixel_pos[3]/distance);
//                      rotate(vector,pixel_pos,0,pixel_pos[3]/distance,pixel_pos[2]/distance);
                    }
                    /* construct the diffracted-beam unit vector to this sub-pixel */
                    airpath = unitize(pixel_pos,diffracted);

                    /* solid angle subtended by a pixel: (pix/airpath)^2*cos(2theta) */
                    omega_pixel = pixel_size*pixel_size/airpath/airpath*close_distance/airpath;
                    /* option to turn off obliquity effect, inverse-square-law only */
                    if(point_pixel) omega_pixel = 1.0/airpath/airpath;
                    omega_sum += omega_pixel;

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
//                              printf("%d %f %f %f\n",mos_tic,mosaic_umats[mos_tic*9+0],mosaic_umats[mos_tic*9+1],mosaic_umats[mos_tic*9+2]);
//                              printf("%d %f %f %f\n",mos_tic,mosaic_umats[mos_tic*9+3],mosaic_umats[mos_tic*9+4],mosaic_umats[mos_tic*9+5]);
//                              printf("%d %f %f %f\n",mos_tic,mosaic_umats[mos_tic*9+6],mosaic_umats[mos_tic*9+7],mosaic_umats[mos_tic*9+8]);

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
                                if(round_xtal){
                                    /* use sinc3 for elliptical xtal shape */
                                    F_latt = Na*Nb*Nc*sinc3(M_PI*sqrt(Na*Na*(h-h0)*(h-h0) + Nb*Nb*(k-k0)*(k-k0) + Nc*Nc*(l-l0)*(l-l0) ) );
                                }
                                else
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
                                if(binary_spots) {
                                    /* make a flat-top spot of same volume */
                                    F_latt = (fabs(F_latt) > 0.5*Na*Nb*Nc)*Na*Nb*Nc;
                                }
                                /* no need to go further if result will be zero */
                                if(F_latt == 0.0) continue;


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
//                                  d_star = magnitude(relp)

                                    /* reciprocal-space coordinates of center of Ewald sphere */
                                    Ewald0[1] = -incident[1]/lambda/1e10;
                                    Ewald0[2] = -incident[2]/lambda/1e10;
                                    Ewald0[3] = -incident[3]/lambda/1e10;
//                                  1/lambda = magnitude(Ewald0)

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
                                            if(verbose) printf ("WARNING: out of range for three point interpolation: h,k,l,h0,k0,l0: %g,%g,%g,%d,%d,%d \n", h,k,l,h0,k0,l0);
                                            if(verbose) printf("WARNING: further warnings will not be printed! ");
                                        }
                                        F_cell = default_F;
                                        continue;
                                    }

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
                                else
                                {
                                    if ( (h0<=h_max) && (h0>=h_min) && (k0<=k_max) && (k0>=k_min) && (l0<=l_max) && (l0>=l_min)  ) {
                                        /* just take nearest-neighbor */
                                        F_cell = Fhkl[h0-h_min][k0-k_min][l0-l_min];
                                    }
                                    else
                                    {
                                        F_cell = default_F; // usually zero
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
                                I += F_cell*F_cell*F_latt*F_latt;
                            }
                            /* end of mosaic loop */
                        }
                        /* end of phi loop */
                    }
                    /* end of source loop */
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
                ++j;
                ++progress_pixel;
            }
        }
    }
    if(verbose) printf("done with pixel loop\n");

    if(verbose) printf("solid angle subtended by detector = %g steradian ( %g%% sphere)\n",omega_sum/steps,100*omega_sum/steps/4/M_PI);

    }


void
nanoBragg::to_smv_format(
    std::string const& fileout, double intfile_scale, double photon_scale, bool noisify){

    int pixels = spixels * fpixels;
    floatimage = raw.begin();
    FILE* outfile;
    double max_value = (double)std::numeric_limits<unsigned short int>::max();
    double saturation = floor(max_value - 1 );
    static const char *byte_order = "little_endian";

    /* output as ints */
    af::versa<unsigned short int, af::c_grid<2> > intimage_v(
      af::c_grid<2> (spixels,fpixels));
    unsigned short int * intimage = intimage_v.begin();

    if(photon_scale <= 0) photon_scale=1.0;
    if(intfile_scale <= 0.0){
        if(verbose) printf("providing default scaling: max_I = %g at (%g %g)\n",max_I,max_I_x,max_I_y);
        intfile_scale = 55000.0/(max_I*photon_scale);
        if(verbose) printf("providing default scaling: intfile_scale = %f\n",intfile_scale);
    }

    int j = 0;
    for(int ypixel=0;ypixel<spixels;++ypixel){
      for(int xpixel=0;xpixel<fpixels;++xpixel){

        if(noisify)
            {
            /* apply photon-counting noise only */
            intimage[j] = (unsigned short int) (std::min(saturation,
                           poidev(floatimage[j]*photon_scale, &seed)*intfile_scale + adc_offset ));
            }
        else
            {
            /* no noise, just use intfile_scale */
            intimage[j] = (unsigned short int) (std::min(saturation,
                           floatimage[j]*photon_scale*intfile_scale + adc_offset ));
            }
            ++j;
        }
    }
    if (verbose){
      printf("writing %s as %d-byte integers\n",fileout.c_str(),
            (int)sizeof(unsigned short int));}
    outfile = fopen(fileout.c_str(),"w");
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
    fprintf(outfile,"TIME=%g;\n",exposure);
    fprintf(outfile,"TWOTHETA=%g;\n",detector_twotheta*RTD);
    fprintf(outfile,"DETECTOR_SN=000;\n");
    fprintf(outfile,"ADC_OFFSET=%g;\n",adc_offset);
    fprintf(outfile,"BEAMLINE=fake;\n");
    fprintf(outfile,"}\f");
    while ( ftell(outfile) < 512 ){ fprintf(outfile," "); };
    fwrite(intimage,sizeof(unsigned short int),pixels,outfile);
    fclose(outfile);

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




/* returns a 9-element unitary matrix for a random isotropic rotation on a spherical cap of diameter "mosaicity" */
/* mosaic = 90 deg is a full sphere */
double *mosaic_rotation_umat(double mosaicity, double umat[9], long *seed)
{
//    double ran1(long *idum);
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



/* random number generators */

/* Poisson deviate given expectation value of photon count (xm) */
double poidev(double xm, long *idum)
{
//    double gammln(double xx);
//    double ran1(long *idum);
    /* oldm is a flag for whether xm has changed since last call */
    static double sq,alxm,g,oldm=(-1.0);
    double em,t,y;

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
double gaussdev(long *idum)
{
//    double ran1(long *idum);
    static int iset=0;
    static double gset;
    double fac,rsq,v1,v2;

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
double lorentzdev(long *seed) {
//    double ran1(long *idum);

    return tan(M_PI*(ran1(seed)-0.5));
}

/* return triangular deviate with FWHM = 1 */
double triangledev(long *seed) {
    double ran1(long *idum);
    double value;

    value = ran1(seed);
    if(value > 0.5){
        value = sqrt(2*(value-0.5))-1;
    }else{
        value = 1-sqrt(2*value);
    }

    return value;
}



double expdev(long *idum)
{
    double dum;

    do
    dum=ran1(idum);
    while( dum == 0.0);
    return -log(dum);
}



/* ln of the gamma function */
double gammln(double xx)
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

double ran1(long *idum)
{
    int j;
    long k;
    static long iy=0;
    static long iv[NTAB];
    double temp;

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
    return log(16)/M_PI*exp(-log(16)*(x*x+y*y));
}
double ngauss2Dinteg(double x,double y)
{
    return 0.125*(erf(2*x*sqrt(log(2)))*erf(y*sqrt(log(16)))*sqrt(log(16)/log(2)));
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
    unsigned long i,j;
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
    double psi=0.0;
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





/* 2D Gaussian integral=1 */
double ngauss2D(double x, double y, double fwhm)
{
    return log(16)/M_PI*fwhm*fwhm*exp(-log(16)*((x*x+y*y)/(fwhm*fwhm) ));
}

/* integral of Gaussian fwhm=1 integral=1 */
double ngauss2D_integ(double x, double y)
    {
    return 0.125*(erf(2*x*sqrt(log(2)))*erf(y*sqrt(log(16)))*sqrt(log(16)/log(2)));
}

/* unit volume integrated over a pixel, fwhm = 1 */
double ngauss2D_pixel(double x,double y,double pix)
    {
    return ngauss2D_integ(x+pix/2.,y+pix/2.)-ngauss2D_integ(x+pix/2.,y-pix/2.)-ngauss2D_integ(x-pix/2.,y+pix/2.)+ngauss2D_integ(x-pix/2.,y-pix/2.);
    }

double integrate_gauss_over_pixel(double x, double y, double fwhm, double pix)
    {
    return ngauss2D_pixel(x/fwhm,y/fwhm,pix/fwhm);
}


double fiber2D(double x,double y,double g)
{
    /* g/(2*pi)*(g**2+x**2+y**2)**(-3/2) */
    double temp;
    temp = sqrt(g*g+x*x+y*y);
    if(temp <= 0.0) return 0.0;
    return g/2.0/M_PI/temp/temp/temp;
}
double fiber2D_integ(double x,double y,double g)
{
    return atan((x*y)/(g*sqrt(g*g + x*x + y*y)))/2.0/M_PI;
}
double fiber2D_pixel(double x,double y,double g,double pix)
    {
  return fiber2D_integ(x+pix/2.,y+pix/2.,g)-fiber2D_integ(x+pix/2.,y-pix/2.,g)-fiber2D_integ(x-pix/2.,y+pix/2.,g)+fiber2D_integ(x-pix/2.,y-pix/2.,g);
}
double integrate_fiber_over_pixel(double x, double y, double g, double pix)
        {
    return fiber2D_pixel(x,y,g,pix);
        }



/* function for applying the PSF, returns NEW image that is blurred version of input */
double *apply_psf(double *inimage, int fpixels, int spixels, psf_type psftype, double fwhm_pixels, int user_psf_radius)
        {
    double max_I;
    double *outimage=NULL;
    double *kernel;
    int x0,y0,x,y,dx,dy;
    double g,rsq;
    double photon_noise,lost_photons=0.0,total_lost_photons=0.0;
    int pixels,maxwidth,kernel_size,psf_radius;
    int i,j,k;
    double photonloss_factor = 10.0;
    int verbose=1;

    /* convert fwhm to "g" distance : fwhm = sqrt((2**(2./3)-1))/2*g */
    g = fwhm_pixels * 0.652383013252053;

    if(psftype == UNKNOWN)
    {
        if(verbose) printf("ERROR: unknown PSF type\n");
        return inimage;
        }

    pixels = fpixels*spixels;
    if(pixels == 0)
        {
        if(verbose) printf("ERROR: apply_psf image has zero size\n");
        return inimage;
    }

    if(fwhm_pixels <= 0.0)
                {
        if(verbose) printf("WARNING: apply_psf function has zero size\n");
        return inimage;
            }

    /* start with a clean slate */
    if(outimage!=NULL) free(outimage);
    outimage = (double *) calloc(pixels+10,sizeof(double));

    psf_radius = user_psf_radius;
    if(psf_radius <= 0)
            {
        /* auto-select radius */

        /* preliminary stats */
        max_I = 0.0;
        for(i=0;i<pixels;++i)
            {
            /* optionally scale the input file */
            if(max_I < inimage[i]) max_I = inimage[i];
            }
        if(verbose) printf("  maximum input photon/pixel: %g\n",max_I);

        if(max_I<=0.0)
            {
            /* nothing to blur */
            if(verbose) printf("WARNING: no photons, PSF skipped\n");
            return outimage;
            }

        /* at what level will an error in intensity be lost? */
        photon_noise = sqrt(max_I);
        lost_photons = photon_noise/photonloss_factor;

        if(psftype == GAUSS)
            {
            /* calculate the radius beyond which only 0.5 photons will fall */
            psf_radius = 1+ceil( sqrt(-log(lost_photons/max_I)/log(4)/2)*fwhm_pixels );
            if(verbose) printf("  auto-selected psf_radius = %d pixels\n",psf_radius);
            }
        if(psftype == FIBER)
            {
            /* calculate the radius r beyond which only 0.5 photons will fall */
            /* r = sqrt((g*(max_I/0.5))**2-g**2)
                 ~ 2*g*max_I */
            psf_radius = 1+ceil( g*(max_I/lost_photons)  );
            if(verbose) printf("  auto-selected psf_radius = %d pixels\n",psf_radius);
            }
        if(psf_radius == 0) psf_radius = 1;
    }
    /* limit psf kernel to be no bigger than 4x the input image */
    maxwidth = fpixels;
    if(spixels > maxwidth) maxwidth = spixels;
    if(psf_radius > maxwidth) psf_radius = maxwidth;
    kernel_size = 2*psf_radius+1;

    /* now alocate enough space to store the PSF kernel image */
    kernel = (double *) calloc(kernel_size*kernel_size,sizeof(double));
    if(kernel == NULL)
            {
        perror("apply_psf: could not allocate memory for PSF kernel");
                exit(9);
            }

    /* cache the PSF in an array */
    for(dy=-psf_radius;dy<=psf_radius;++dy)
    {
        for(dx=-psf_radius;dx<=psf_radius;++dx)
        {
            rsq = dx*dx+dy*dy;
            if(rsq > psf_radius*psf_radius) continue;

            /* this could be more efficient */
            k = kernel_size*(kernel_size/2+dy)+kernel_size/2+dx;


            if( psftype == GAUSS ) {
                kernel[k] = integrate_gauss_over_pixel(dx,dy,fwhm_pixels,1.0);
        }
            if( psftype == FIBER ) {
                kernel[k] = integrate_fiber_over_pixel(dx,dy,g,1.0);
    }
    }
}

    /* implement PSF  */
    for(i=0;i<pixels;++i)
{
        x0 = i%fpixels;
        y0 = (i-x0)/fpixels;

        /* skip if there is nothing to add */
        if(inimage[i] <= 0.0) continue;

        if(user_psf_radius != 0)
    {
            psf_radius = user_psf_radius;
    }
        else
        {
            /* at what level will an error in intensity be lost? */
            photon_noise = sqrt(inimage[i]);
            lost_photons = photon_noise/photonloss_factor;

            if(psftype == GAUSS)
            {
                /* calculate the radius beyond which only 0.5 photons will fall
                   r = sqrt(-log(lost_photons/total_photons)/log(4)/2)*fwhm */
                psf_radius = 1+ceil( sqrt(-log(lost_photons/inimage[i])/log(16))*fwhm_pixels );
//              printf("  auto-selected psf_radius = %d pixels\n",psf_radius);
}
            if(psftype == FIBER)
{
                /* calculate the radius beyond which only 0.5 photons will fall
                   r = sqrt((g*(total_photons/lost_photons))**2-g**2)
                     ~ g*total_photons/lost_photons */
                psf_radius = 1+ceil( g*(inimage[i]/lost_photons)  );
//              printf("  (%d,%d) auto-selected psf_radius = %d pixels\n",x0,y0,psf_radius);
            }
        }
        if(psf_radius == 0) psf_radius = 1;
        /* limit psf kernel to be no bigger than 4x the input image */
        maxwidth = fpixels;
        if(spixels > maxwidth) maxwidth = spixels;
        if(psf_radius > maxwidth) psf_radius = maxwidth;

        /* given the radius, how many photons will escape? */
        if(psftype == GAUSS)
    {
            /* r = sqrt(-log(lost_photons/total_photons)/log(16))*fwhm */
            /* lost_photons = total_photons*exp(-log(16)*(r^2/fwhm^2)) */
            rsq = psf_radius;
            rsq = rsq/fwhm_pixels;
            rsq = rsq*rsq;
            lost_photons = inimage[i]*exp(-log(16)*rsq);
        }
        if(psftype == FIBER)
        {
            /* r ~ g*total_photons/lost_photons
               normalized integral from r=inf to "r" :  g/sqrt(g**2+r**2) */
            lost_photons = inimage[i]*g/sqrt(g*g+psf_radius*psf_radius);
        }
        /* accumulate this so we can add it to the whole image */
        total_lost_photons += lost_photons;

        for(dx=-psf_radius;dx<=psf_radius;++dx)
        {
            for(dy=-psf_radius;dy<=psf_radius;++dy)
            {
                /* this could be more efficient */
                k = kernel_size*(kernel_size/2+dy)+kernel_size/2+dx;
                if(kernel[k] == 0.0) continue;

                rsq = dx*dx+dy*dy;
                if(rsq > psf_radius*psf_radius) continue;
                x = x0+dx;
                y = y0+dy;
                if(x<0 || x>fpixels) continue;
                if(y<0 || y>spixels) continue;

                /* index into output array */
                j = y*fpixels+x;
                /* do not wander off the output array */
                if(j<0 || j > pixels) continue;

                outimage[j] += inimage[i]*kernel[k];
            }
        }
    }
    /* now we have some lost photons, add them back "everywhere" */
    lost_photons = total_lost_photons/pixels;
    if(verbose) printf("adding back %g lost photons\n",total_lost_photons);
    for(i=0;i<pixels;++i)
    {
        outimage[i] += lost_photons;
    }

    /* don't need kernel anymore. but should we always allocate outimage? */
    free(kernel);
    return outimage;
}

}}// namespace simtbx::nanoBragg
