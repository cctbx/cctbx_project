#include <simtbx/diffBragg/src/diffBragg.h>
#include <assert.h>
namespace simtbx {
namespace nanoBragg {

// BEGIN derivative manager
derivative_manager::derivative_manager(){}

void derivative_manager::initialize(int sdim, int fdim)
{
    raw_pixels = af::flex_double(af::flex_grid<>(sdim,fdim));
    dI=0;

    // for second derivatives
    dI2=0;
    raw_pixels2 = af::flex_double(af::flex_grid<>(sdim,fdim));
}

void derivative_manager::increment_image(int idx, double value, double value2){
    double* floatimage = raw_pixels.begin();
    floatimage[idx] += value;

    // increment second derivatives
    double* floatimage2 = raw_pixels2.begin();
    floatimage2[idx] += value2;
}
// END derivative manager

//BEGIN origin manager
origin_manager::origin_manager(){}

void origin_manager::increment(
        vec3 V, mat3 N, mat3 UBO, vec3 k_diffracted, vec3 o_vec, double air_path, double wavelen,
        double Hrad, double Fcell, double Flatt, double fudge,
        double source_I, double capture_fraction, double omega_pixel, double pixel_size){

    vec3 dq = 1/air_path * (dk - k_diffracted*(1/air_path/air_path)*(dk*k_diffracted));
    dq /= wavelen;
    vec3 dV = N*(UBO*dq);
    double V_dot_dV = V*dV;
    double dHrad = 2*V_dot_dV;
    double a = 1 / 0.63 * fudge;
    double dFlatt = -1*a*Flatt*dHrad;
    double c = Fcell*Fcell*source_I*capture_fraction;
    dI += c*2*Flatt*dFlatt*omega_pixel;
    // contribution from derivative of solid angle
    //FIXME I mist be broken:
    double D = k_diffracted * o_vec;
    double pp = pixel_size*pixel_size;
    double domega  = -3*pp/pow(air_path,5)*(dk*k_diffracted)*D;
    domega += pp/pow(air_path,3)* (o_vec * dk);
    dI += c*Flatt*Flatt*domega;

    dI2 += 0;
};
//END origin manager

// BEGIN Ncells_abc manager
Ncells_manager::Ncells_manager(){}

void Ncells_manager::increment(
        vec3 V, vec3 H0_vec, vec3 H_vec,
        double Hrad, double Fcell, double Flatt, double fudge,
        double source_I, double capture_fraction, double omega_pixel){

    vec3 dV = H_vec - H0_vec;
    double V_dot_dV = V*dV;
    double dHrad = 2*V_dot_dV;
    double a = 1 / 0.63 * fudge;
    double dFlatt = Flatt * (3/value - a*dHrad); // value in this case is Ncells_abc size
    double c = Fcell*Fcell*source_I*capture_fraction*omega_pixel;
    dI += c*2*Flatt*dFlatt;

    double dFlatt2 = 6*Flatt/value/value -6*a*Flatt/value*dHrad + a*a*Flatt*dHrad*dHrad;
    dI2 += c*2*(dFlatt2*Flatt + dFlatt*dFlatt);
};

//END Ncells_abc manager

// Begin Fcell manager
Fcell_manager::Fcell_manager(){}

void Fcell_manager::increment(
    double Hrad, double Fcell, double Flatt, double fudge,
    double source_I, double capture_fraction, double omega_pixel){
    double c = Flatt*Flatt*source_I*capture_fraction*omega_pixel;
    dI += c*2*Fcell;
    dI2 += c*2;
};


//BEGIN unit cell manager
ucell_manager::ucell_manager(){}

void ucell_manager::increment(
    vec3 V,
    mat3 NABC, mat3 UR, vec3 q, mat3 Ot,
    double Hrad, double Fcell, double Flatt, double fudge,
    double source_I, double capture_fraction, double omega_pixel){

  vec3 dV = NABC*(UR*dB*Ot).transpose()*q;
  double V_dot_dV = V*dV;
  double dHrad = 2*V_dot_dV;
  double a = 1 / 0.63 * fudge;
  double dFlatt = -1*a*Flatt*dHrad;
  double c = Fcell*Fcell*source_I*capture_fraction*omega_pixel;
  dI += c*2*Flatt*dFlatt; //Fcell*Fcell*2*Flatt*source_I*capture_fraction*omega_pixel*dFlatt;

  vec3 dV2 = NABC*(UR*dB2*Ot).transpose() * q;
  double dFlatt2 = -2*a*(dFlatt * V_dot_dV + Flatt*(dV*dV) + Flatt*(V*dV2));
  dI2 += c*2*(dFlatt2*Flatt + dFlatt*dFlatt);
};


// BEGIN rotation manager begin
rot_manager::rot_manager(){}

void rot_manager::set_R(){assert (false);}

void rot_manager::increment(
        double fudge,
        mat3 X, mat3 Y, mat3 Z,
        mat3 X2, mat3 Y2, mat3 Z2,
        mat3 N, mat3 U, mat3 UBOt,
        vec3 q, vec3 V,
        double Hrad, double Fcell, double Flatt,
        double source_I, double capture_fraction, double omega_pixel)
{
  /* X,Y,Z will change depending on whether its RotX, RotY, or RotZ manager */
  XYZ = X*Y*Z;
  vec3 dV = N*(U*XYZ*UBOt).transpose() * q;
  double V_dot_dV = V*dV;
  double dHrad = 2*V_dot_dV;
  double a = 1 / 0.63 * fudge;
  double dFlatt = -1*a*Flatt*dHrad;
  double c = Fcell*Fcell*source_I*capture_fraction*omega_pixel;
  dI += c*2*Flatt*dFlatt;

  XYZ2 = X2*Y2*Z2;
  vec3 dV2 = N*(U*XYZ2*UBOt).transpose() * q;
  double dFlatt2 = -2*a*(dFlatt * V_dot_dV + Flatt*(dV*dV) + Flatt*(V*dV2));
  dI2 += c*2*(dFlatt2*Flatt + dFlatt*dFlatt);
}

rotX_manager::rotX_manager(){
    value = 0;
    set_R();
}
rotZ_manager::rotZ_manager(){
    value = 0;
    set_R();
}
rotY_manager::rotY_manager(){
    value = 0;
    set_R();
}

void rotX_manager::set_R(){
    R = mat3(1,           0,           0,
             0,  cos(value), sin(value),
             0, -sin(value), cos(value));

    dR= mat3(0,           0,           0,
             0,  -sin(value), cos(value),
             0, -cos(value), -sin(value));

    dR2= mat3(0,           0,           0,
              0,  -cos(value), -sin(value),
              0,   sin(value), -cos(value));

    //SCITBX_EXAMINE(dR2[0]);
    //SCITBX_EXAMINE(dR2[1]);
    //SCITBX_EXAMINE(dR2[2]);
    //SCITBX_EXAMINE(dR2[3]);
    //SCITBX_EXAMINE(dR2[4]);
    //SCITBX_EXAMINE(dR2[5]);
    //SCITBX_EXAMINE(dR2[6]);
    //SCITBX_EXAMINE(dR2[7]);
    //SCITBX_EXAMINE(dR2[8]);

}
void rotY_manager::set_R(){
    R= mat3(cos(value),0, -sin(value),
             0,         1,             0,
            sin(value), 0, cos(value));

    dR= mat3(-sin(value),0, -cos(value),
                0,          0,             0,
                cos(value), 0, -sin(value));

    dR2= mat3(-cos(value),0, sin(value),
              0,          0,          0,
             -sin(value), 0, -cos(value));
}
void rotZ_manager::set_R(){
    R = mat3(cos(value),  sin(value), 0,
              -sin(value), cos(value), 0,
                         0,           0, 1);

    dR = mat3(-sin(value),  cos(value), 0,
               -cos(value), -sin(value), 0,
                           0,           0, 0);

    dR2= mat3(-cos(value), -sin(value), 0,
               sin(value), -cos(value), 0,
                        0,           0, 0);
}
// END rot manager

// BEGIN diffBragg
diffBragg::diffBragg(const dxtbx::model::Detector& detector, const dxtbx::model::Beam& beam,
            int verbose, int panel_id = 0):
    nanoBragg(detector, beam, verbose, panel_id)
    {

    mat3 EYE = mat3(1,0,0,0,1,0,0,0,1);
    Omatrix = mat3(1,0,0,0,1,0,0,0,1);
    psi = 0;

    RotMats.push_back(EYE);
    RotMats.push_back(EYE);
    RotMats.push_back(EYE);

    dRotMats.push_back(EYE);
    dRotMats.push_back(EYE);
    dRotMats.push_back(EYE);
    d2RotMats.push_back(EYE);
    d2RotMats.push_back(EYE);
    d2RotMats.push_back(EYE);

    R3.push_back(EYE);
    R3.push_back(EYE);
    R3.push_back(EYE);
    R3_2.push_back(EYE);
    R3_2.push_back(EYE);
    R3_2.push_back(EYE);

    boost::shared_ptr<rot_manager> rotX = boost::shared_ptr<rot_manager>(new rotX_manager());
    boost::shared_ptr<rot_manager> rotY = boost::shared_ptr<rot_manager>(new rotY_manager());
    boost::shared_ptr<rot_manager> rotZ = boost::shared_ptr<rot_manager>(new rotZ_manager());

    boost::shared_ptr<ucell_manager> uc1 = boost::shared_ptr<ucell_manager>(new ucell_manager());
    boost::shared_ptr<ucell_manager> uc2 = boost::shared_ptr<ucell_manager>(new ucell_manager());
    boost::shared_ptr<ucell_manager> uc3 = boost::shared_ptr<ucell_manager>(new ucell_manager());
    boost::shared_ptr<ucell_manager> uc4 = boost::shared_ptr<ucell_manager>(new ucell_manager());
    boost::shared_ptr<ucell_manager> uc5 = boost::shared_ptr<ucell_manager>(new ucell_manager());
    boost::shared_ptr<ucell_manager> uc6 = boost::shared_ptr<ucell_manager>(new ucell_manager());

    boost::shared_ptr<Ncells_manager> nc = boost::shared_ptr<Ncells_manager>(new Ncells_manager());

    boost::shared_ptr<origin_manager> orig0 = boost::shared_ptr<origin_manager>(new origin_manager());

    //boost::shared_ptr<Fcell_manager> fcell_man = boost::shared_ptr<Fcell_manager>(new Fcell_manager());
    fcell_man = boost::shared_ptr<Fcell_manager>(new Fcell_manager());
    fcell_man->refine_me = false;

    rotX->refine_me = false;
    rotY->refine_me = false;
    rotZ->refine_me = false;
    uc1->refine_me = false;
    uc2->refine_me = false;
    uc3->refine_me = false;
    uc4->refine_me = false;
    uc5->refine_me = false;
    uc6->refine_me = false;

    nc->refine_me = false;

    orig0->refine_me = false;
    orig0->dk = vec3(0,0,1);

    rot_managers.push_back(rotX);
    rot_managers.push_back(rotY);
    rot_managers.push_back(rotZ);

    ucell_managers.push_back(uc1);
    ucell_managers.push_back(uc2);
    ucell_managers.push_back(uc3);
    ucell_managers.push_back(uc4);
    ucell_managers.push_back(uc5);
    ucell_managers.push_back(uc6);

    Ncells_managers.push_back(nc);

    origin_managers.push_back(orig0);


    // set ucell gradients, Bmatrix is upper triangular in diffBragg?
    // note setting these derivatives is only useful for parameter reduction code where one computes chain rule
    for (int i=0; i <6; i++){
        mat3 bb =  mat3(0,0,0,0,0,0,0,0,0);
        if (i <3)
            bb[i]= 1;
        else if (i ==3 || i==4)
            bb[i+1] = 1;
        else if (i==5)
            bb[8] = 1;

        if (verbose>5)
            printf("Param %d\nbb_real:\n%11.8f %11.8f %11.8f\n %11.8f %11.8f %11.8f\n %11.8f %11.8f %11.8f\n", i,
                bb[0], bb[1], bb[2],
                bb[3], bb[4], bb[5],
                bb[6], bb[7], bb[8]);
        ucell_managers[i]->dB = bb;
        ucell_managers[i]->dB2 = mat3(0,0,0,0,0,0,0,0,0);
        }

    init_raw_pixels_roi();
    initialize_managers();
    }


void diffBragg::update_dxtbx_geoms(
    const dxtbx::model::Detector& detector,
    const dxtbx::model::Beam& beam,
    int panel_id){

    /* BEAM properties first */

    double temp;
    vec3 xyz;
    /* direction in 3-space of beam vector */
    xyz = beam.get_unit_s0();
    beam_vector[1] = xyz[0];
    beam_vector[2] = xyz[1];
    beam_vector[3] = xyz[2];
    unitize(beam_vector,beam_vector);

    /* central wavelength, in Angstrom */
    lambda0 = beam.get_wavelength()*1e-10;

    /* divergence, what are the DXTBX units? */
    temp = beam.get_divergence();
    if(temp>0.0) hdivrange = vdivrange = temp;

    /* assume this is photons/s, unless it is zero */
    temp = beam.get_flux();
    if(temp>0.0) flux = temp;

    /* assume this is Kahn polarization parameter */
    temp = beam.get_polarization_fraction();
    if(temp>=-1.0 && temp<=1.0) polarization = temp;

    /* dxtbx polarization points down B vector, we want the E vector */
    xyz = beam.get_polarization_normal();
    vert_vector[1] = xyz[0];
    vert_vector[2] = xyz[1];
    vert_vector[3] = xyz[2];
    unitize(vert_vector,vert_vector);
    cross_product(beam_vector,vert_vector,polar_vector);
    unitize(polar_vector,polar_vector);

    /* DETECTOR properties */
    /* size of the pixels in meters, this should not vary after instantiation */
    SCITBX_ASSERT(pixel_size == detector[panel_id].get_pixel_size()[0]/1000.);

    /* pixel count in short and fast-axis directions, should not change after instantiation */
    SCITBX_ASSERT(spixels == detector[panel_id].get_image_size()[1]);
    SCITBX_ASSERT( fpixels == detector[panel_id].get_image_size()[0]);

    /* direction in 3-space of detector axes */
    SCITBX_ASSERT (beam_convention == CUSTOM);

    /* typically: 1 0 0 */
    fdet_vector[1] = detector[panel_id].get_fast_axis()[0];
    fdet_vector[2] = detector[panel_id].get_fast_axis()[1];
    fdet_vector[3] = detector[panel_id].get_fast_axis()[2];
    unitize(fdet_vector,fdet_vector);
    /* typically: 0 -1 0 */
    sdet_vector[1] = detector[panel_id].get_slow_axis()[0];
    sdet_vector[2] = detector[panel_id].get_slow_axis()[1];
    sdet_vector[3] = detector[panel_id].get_slow_axis()[2];
    unitize(sdet_vector,sdet_vector);
    /* set orthogonal vector to the detector pixel array */
    cross_product(fdet_vector,sdet_vector,odet_vector);
    unitize(odet_vector,odet_vector);

    /* dxtbx origin is location of outer corner of the first pixel */
    pix0_vector[1] = detector[panel_id].get_origin()[0]/1000.0;
    pix0_vector[2] = detector[panel_id].get_origin()[1]/1000.0;
    pix0_vector[3] = detector[panel_id].get_origin()[2]/1000.0;

    Fclose = Xclose = -dot_product(pix0_vector,fdet_vector);
    Sclose = Yclose = -dot_product(pix0_vector,sdet_vector);
    close_distance = distance =  dot_product(pix0_vector,odet_vector);

    /* set beam centre */
    scitbx::vec2<double> dials_bc = detector[panel_id].get_beam_centre(beam.get_s0());
    Xbeam = dials_bc[0]/1000.0;
    Ybeam = dials_bc[1]/1000.0;

    /* detector sensor layer properties */
    detector_thick   = detector[panel_id].get_thickness();
    temp = detector[panel_id].get_mu();        // is this really a mu? or mu/rho ?
    if(temp>0.0) detector_attnlen = 1.0/temp;

    /* quantum_gain = amp_gain * electrooptical_gain, does not include capture_fraction */
    quantum_gain = detector[panel_id].get_gain();

    //adc_offset = detector[panel_id].ADC_OFFSET;

    /* SPINDLE properties */

    /* By default align the rotation axis with the detector fast direction */
    spindle_vector[1] = fdet_vector[1];
    spindle_vector[2] = fdet_vector[2];
    spindle_vector[3] = fdet_vector[3];
    unitize(spindle_vector,spindle_vector);

    /* OMG So important otherwise center walks */
    ORGX=NAN;
    ORGY=NAN;

    init_beam();
    init_beamcenter();
    update_beamcenter();

    //SCITBX_EXAMINE(Yclose);
    //SCITBX_EXAMINE(Xclose);
    //SCITBX_EXAMINE(Ybeam);
    //SCITBX_EXAMINE(Xbeam);
    //SCITBX_EXAMINE(distance);
    //SCITBX_EXAMINE(close_distance);
    //printf("Done updating!\n");
    }

void diffBragg::init_raw_pixels_roi(){
    int fdim = roi_xmax-roi_xmin+1;
    int sdim = roi_ymax-roi_ymin+1;
    raw_pixels_roi = af::flex_double(af::flex_grid<>(sdim,fdim));
}

void diffBragg::initialize_managers()
{
    int fdim = roi_xmax-roi_xmin+1;
    int sdim = roi_ymax-roi_ymin+1;
    for (int i_rot=0; i_rot < 3; i_rot++){
        if (rot_managers[i_rot]->refine_me)
            rot_managers[i_rot]->initialize(sdim, fdim);
    }
    for (int i_uc=0; i_uc < 6; i_uc++){
        if (ucell_managers[i_uc]->refine_me)
            ucell_managers[i_uc]->initialize(sdim, fdim);
    }
    if (Ncells_managers[0]->refine_me)
        Ncells_managers[0]->initialize(sdim, fdim);
    if (origin_managers[0]->refine_me)
        origin_managers[0]->initialize(sdim, fdim);

    if (fcell_man->refine_me)
        fcell_man->initialize(sdim, fdim);
}

void diffBragg::vectorize_umats()
{
    /* vector store two copies of Umats, one unperturbed for reference */
    for(mos_tic=0;mos_tic<mosaic_domains;++mos_tic)
    {
        double uxx,uxy,uxz,uyx,uyy,uyz,uzx,uzy,uzz;
        uxx = mosaic_umats[mos_tic*9+0];
        uxy = mosaic_umats[mos_tic*9+1];
        uxz = mosaic_umats[mos_tic*9+2];
        uyx = mosaic_umats[mos_tic*9+3];
        uyy = mosaic_umats[mos_tic*9+4];
        uyz = mosaic_umats[mos_tic*9+5];
        uzx = mosaic_umats[mos_tic*9+6];
        uzy = mosaic_umats[mos_tic*9+7];
        uzz = mosaic_umats[mos_tic*9+8];
        mat3 U = mat3(uxx, uxy, uxz,
                      uyx, uyy, uyz,
                      uzx, uzy, uzz);
        UMATS.push_back(U);
        UMATS_RXYZ.push_back(U);
    }
}

void diffBragg::refine(int refine_id){
    int fdim = roi_xmax-roi_xmin+1;
    int sdim = roi_ymax-roi_ymin+1;
    if (refine_id >= 0 && refine_id < 3  ){
        // 3 possitle rotation managers (rotX, rotY, rotZ)
        rot_managers[refine_id]->refine_me=true;
        rot_managers[refine_id]->initialize(sdim, fdim);
    }
    else if (refine_id >=3 and refine_id < 9 ){
        // 6 possible unit cell managers (a,b,c,al,be,ga)
        ucell_managers[refine_id-3]->refine_me=true;
        ucell_managers[refine_id-3]->initialize(sdim, fdim);
    }
    else if (refine_id==9){
        Ncells_managers[0]->refine_me=true;
        Ncells_managers[0]->initialize(sdim, fdim);
    }
    else if (refine_id==10){
        origin_managers[0]->refine_me=true;
        origin_managers[0]->initialize(sdim, fdim);
    }
    else if(refine_id==11){
        fcell_man->refine_me=true;
        fcell_man->initialize(sdim, fdim);
    }
}

void diffBragg::set_ucell_derivative_matrix(int refine_id, af::shared<double> const& value){
    int ucell_param_idx = refine_id-3;  // its just how the API works, pass in 3 for first ucell matrix
    if (ucell_param_idx < 0 || ucell_param_idx > 5)
      printf("WARNING, passing in wrong refine_id for unit cell parameter (should be 3-8).\nNothing done.\n");
    else
        ucell_managers[ucell_param_idx]->dB = mat3(
                        value[0], value[1], value[2],
                        value[3], value[4], value[5],
                        value[6], value[7], value[8]);
}

void diffBragg::set_ucell_second_derivative_matrix(int refine_id, af::shared<double> const& value){
    int ucell_param_idx = refine_id-3;  // its just how the API works, pass in 3 for first ucell matrix
    if (ucell_param_idx < 0 || ucell_param_idx > 5)
      printf("WARNING, passing in wrong refine_id for unit cell parameter (should be 3-8).\nNothing done.\n");
    else
        ucell_managers[ucell_param_idx]->dB2 = mat3(
                        value[0], value[1], value[2],
                        value[3], value[4], value[5],
                        value[6], value[7], value[8]);
}

/* Begin parameter set/get */

// TODO : rename set_value and get_value because they dont apply to ucell derivatives...
// this function will get exeedingly complicated because it will try to ensure all the dependent parameters get
// adjusted when we update a given parameter that we are refining
// For example updating Ncells_abc should also update oversample, and should also update xtal_size
void diffBragg::set_value( int refine_id, double value ){
    if (refine_id < 3){
        rot_managers[refine_id]->value = value;
        rot_managers[refine_id]->set_R();
    }
    if (refine_id==9){
        Ncells_managers[0]->value = value;
        Na=value;
        Nb=value;
        Nc=value;
        xtal_size_x = -1;
        xtal_size_y = -1;
        xtal_size_z = -1;
        //TODO make me optional!
        update_oversample();
        NABC[0] = value;
        NABC[4] = value;
        NABC[8] = value;
    }
}

double diffBragg::get_value(int refine_id){
    double value(0);
    if (refine_id < 3)
        value = rot_managers[refine_id]->value;
    else if (refine_id ==9)
        value = Ncells_managers[0]->value;
    return value;
}
/* End parameter set/get */

af::flex_double diffBragg::get_derivative_pixels(int refine_id){
    if (refine_id>=0 and refine_id < 3){
        SCITBX_ASSERT(rot_managers[refine_id]->refine_me);
        return rot_managers[refine_id]->raw_pixels;
        }
    else if(refine_id >=3 && refine_id < 9){
        int i_uc = refine_id-3;
        SCITBX_ASSERT(i_uc >= 0);
        SCITBX_ASSERT(i_uc < 6);
        SCITBX_ASSERT(ucell_managers[i_uc]->refine_me);
        return ucell_managers[i_uc]->raw_pixels;
        }
    else if (refine_id==9)
        return Ncells_managers[0]->raw_pixels;
    else if (refine_id==10)
        return origin_managers[0]->raw_pixels;
    else
        return fcell_man->raw_pixels;
}


af::flex_double diffBragg::get_second_derivative_pixels(int refine_id){
    if (refine_id>=0 and refine_id < 3){
        SCITBX_ASSERT(rot_managers[refine_id]->refine_me);
        return rot_managers[refine_id]->raw_pixels2;
        }
    else if(refine_id >=3 && refine_id < 9){
        int i_uc = refine_id-3;
        SCITBX_ASSERT(i_uc >= 0);
        SCITBX_ASSERT(i_uc < 6);
        SCITBX_ASSERT(ucell_managers[i_uc]->refine_me);
        return ucell_managers[i_uc]->raw_pixels2;
        }
    else if (refine_id == 9)
        return Ncells_managers[0]->raw_pixels2;
    else if (refine_id==10)
        return origin_managers[0]->raw_pixels2;
    else
        return fcell_man->raw_pixels2;
}

void diffBragg::zero_raw_pixel_rois(){
    init_raw_pixels_roi();
    initialize_managers();
}

/* polarization factor */
/* override this to store variables needed for derivatives */
double diffBragg::polarization_factor(double kahn_factor, double *incident, double *diffracted, double *axis)
{
    double cos2theta,cos2theta_sqr,sin2theta_sqr;
    //double psi=0.0;
    double E_in[4];
    double B_in[4];

    unitize(incident,incident);
    unitize(diffracted,diffracted);
    unitize(axis,axis);

    /* component of diffracted unit vector along incident beam unit vector */
    cos2theta = dot_product(incident,diffracted);
    cos2theta_sqr = cos2theta*cos2theta;
    sin2theta_sqr = 1-cos2theta_sqr;
    u = cos2theta_sqr;

    if(kahn_factor != 0.0){
        //SCITBX_EXAMINE(kahn_factor);
        /* tricky bit here is deciding which direciton the E-vector lies in for each source
           here we assume it is closest to the "axis" defined above */

        /* cross product to get "vertical" axis that is orthogonal to the cannonical "polarization" */
        cross_product(axis,incident,B_in);
        /* make it a unit vector */
        unitize(B_in,B_in);
        Bi_vec[0] = B_in[1];
        Bi_vec[1] = B_in[2];
        Bi_vec[2] = B_in[3];

        /* cross product with incident beam to get E-vector direction */
        cross_product(incident,B_in,E_in);
        /* make it a unit vector */
        unitize(E_in,E_in);
        Ei_vec[0] = E_in[1];
        Ei_vec[1] = E_in[2];
        Ei_vec[2] = E_in[3];

        /* get components of diffracted ray projected onto the E-B plane */
        kEi = dot_product(diffracted,E_in);
        kBi = dot_product(diffracted,B_in);

        /* compute the angle of the diffracted ray projected onto the incident E-B plane */
        psi = -atan2(kBi,kEi);
    }

    /* correction for polarized incident beam */
    return 0.5*(1.0 + cos2theta_sqr - kahn_factor*cos(2*psi)*sin2theta_sqr);
}

// BEGIN diffBragg_add_spots
void diffBragg::add_diffBragg_spots()
{
    max_I = 0.0;
    i = 0;
    floatimage = raw_pixels.begin();
    double * floatimage_roi = raw_pixels_roi.begin();

    //floatimage_roi = raw_pixels_roi.begin();
    for (int i_rot=0; i_rot < 3; i_rot++){
        if (rot_managers[i_rot]->refine_me){
            RotMats[i_rot] = rot_managers[i_rot]->R;
            dRotMats[i_rot] = rot_managers[i_rot]->dR;
            d2RotMats[i_rot] = rot_managers[i_rot]->dR2;

            R3[i_rot] = RotMats[i_rot];
            R3_2[i_rot] = RotMats[i_rot];
        }
    }

    RXYZ = RotMats[0]*RotMats[1]*RotMats[2];

    //printf("First row: %f | %f | %f \n", RXYZ(0,0), RXYZ(0,1), RXYZ(0,2));
    //printf("Second row: %f | %f | %f \n", RXYZ(1,0), RXYZ(1,1), RXYZ(1,2));
    //printf("Third row: %f | %f | %f \n", RXYZ(2,0), RXYZ(2,1), RXYZ(2,2));

    /*  update Umats to be U*RXYZ   */
    for(mos_tic=0;mos_tic<mosaic_domains;++mos_tic)
        UMATS_RXYZ[mos_tic] = UMATS[mos_tic] * RXYZ;

    if(verbose) printf("TESTING sincg(1,1)= %f\n",sincg(1,1));

    /* make sure we are normalizing with the right number of sub-steps */
    steps = phisteps*mosaic_domains*oversample*oversample;
    subpixel_size = pixel_size/oversample;

    //int min_i= 10000000000;
    //int max_i= -1;
    int roi_fdim = roi_xmax - roi_xmin+1;
    int roi_sdim = roi_ymax - roi_ymin+1;

    sum = sumsqr = 0.0;
    i = sumn = 0;
    progress_pixel = 0;
    omega_sum = 0.0;
    //int roi_i = -1;
    for(spixel=0;spixel<spixels;++spixel)
    {
        for(fpixel=0;fpixel<fpixels;++fpixel)
        {
            /* allow for just one part of detector to be rendered */
            if(fpixel < roi_xmin || fpixel > roi_xmax || spixel < roi_ymin || spixel > roi_ymax)
            {
                ++i; continue;
            }
            //else
                //roi_i += 1;
            /* allow for the use of a mask */
            if(maskimage != NULL)
            {
                /* skip any flagged pixels in the mask */
                if(maskimage[i] == 0)
                {
                    ++i; //++roi_i;
                    continue;
                }
            }
            /* reset photon count for this pixel */
            I = 0;

            /* reset derivative photon counts for the various parameters*/
            for (int i_rot =0 ; i_rot < 3 ; i_rot++){
                if (rot_managers[i_rot]->refine_me){
                    rot_managers[i_rot]->dI =0;
                    rot_managers[i_rot]->dI2 =0;
                }
            }
            for (int i_uc =0 ; i_uc < 6 ; i_uc++){
                if (ucell_managers[i_uc]->refine_me){
                    ucell_managers[i_uc]->dI =0;
                    ucell_managers[i_uc]->dI2 =0;
                }
            }

            if (Ncells_managers[0]->refine_me){
                Ncells_managers[0]->dI =0;
                Ncells_managers[0]->dI2 =0;
            }

            if (origin_managers[0]->refine_me){
                origin_managers[0]->dI=0;
                origin_managers[0]->dI2=0;
            }

            if (fcell_man->refine_me){
                fcell_man->dI = 0;
                fcell_man->dI2 = 0;
            }

            /* loop over sub-pixels */
            for(subS=0;subS<oversample;++subS)
            {
                for(subF=0;subF<oversample;++subF)
                {
                    /* absolute mm position on detector (relative to its origin) */
                    Fdet = subpixel_size*(fpixel*oversample + subF ) + subpixel_size/2.0;
                    Sdet = subpixel_size*(spixel*oversample + subS ) + subpixel_size/2.0;

                    for(thick_tic=0;thick_tic<detector_thicksteps;++thick_tic)
                    {
                        /* assume "distance" is to the front of the detector sensor layer */
                        Odet = thick_tic*detector_thickstep;

                        /* construct detector subpixel position in 3D space */
                        pixel_pos[1] = Fdet*fdet_vector[1]+Sdet*sdet_vector[1]+Odet*odet_vector[1]+pix0_vector[1];
                        pixel_pos[2] = Fdet*fdet_vector[2]+Sdet*sdet_vector[2]+Odet*odet_vector[2]+pix0_vector[2];
                        pixel_pos[3] = Fdet*fdet_vector[3]+Sdet*sdet_vector[3]+Odet*odet_vector[3]+pix0_vector[3];
                        pixel_pos[0] = 0.0;
                        if(curved_detector) {
                            /* construct detector pixel that is always "distance" from the sample */
                                vector[1]=distance*beam_vector[1];
                                vector[2]=distance*beam_vector[2] ;
                                vector[3]=distance*beam_vector[3];
                            /* treat detector pixel coordinates as radians */
                            rotate_axis(vector,newvector,sdet_vector,pixel_pos[2]/distance);
                            rotate_axis(newvector,pixel_pos,fdet_vector,pixel_pos[3]/distance);
                        }
                        /* construct the diffracted-beam unit vector to this sub-pixel */
                        airpath = unitize(pixel_pos,diffracted);

                        /* solid angle subtended by a pixel: (pix/airpath)^2*cos(2theta) */
                        omega_pixel = pixel_size*pixel_size/airpath/airpath*close_distance/airpath;
                        /* option to turn off obliquity effect, inverse-square-law only */
                        if(point_pixel) omega_pixel = 1.0/airpath/airpath;
                        //omega_pixel = 1e-6; // FIXME
                        omega_sum += omega_pixel;

                        /* now calculate detector thickness effects */
                        if(detector_thick > 0.0 && detector_attnlen > 0.0)
                        {
                            /* inverse of effective thickness increase */
                            parallax = dot_product(diffracted,odet_vector);
                            capture_fraction = exp(-thick_tic*detector_thickstep/detector_attnlen/parallax)
                                              -exp(-(thick_tic+1)*detector_thickstep/detector_attnlen/parallax);
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
                        //vec3 incident_0 = vec3(incident[1], incident[2], incident[3]);

                        /* construct the incident beam unit vector while recovering source distance */
                        source_path = unitize(incident,incident);

                        /* construct the scattering vector for this pixel */
                        scattering[1] = (diffracted[1]-incident[1])/lambda;
                        scattering[2] = (diffracted[2]-incident[2])/lambda;
                        scattering[3] = (diffracted[3]-incident[3])/lambda;

                        k_diffracted[0] = pixel_pos[1];
                        k_diffracted[1] = pixel_pos[2];
                        k_diffracted[2] = pixel_pos[3];
                        o_vec[0] = odet_vector[1];
                        o_vec[1] = odet_vector[2];
                        o_vec[2] = odet_vector[3];

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
                                ap_vec[0] = ap[1];
                                ap_vec[1] = ap[2];
                                ap_vec[2] = ap[3];

                                bp_vec[0] = bp[1];
                                bp_vec[1] = bp[2];
                                bp_vec[2] = bp[3];

                                cp_vec[0] = cp[1];
                                cp_vec[1] = cp[2];
                                cp_vec[2] = cp[3];

                                q_vec[0] = 1e-10*scattering[1];
                                q_vec[1] = 1e-10*scattering[2];
                                q_vec[2] = 1e-10*scattering[3];

                                Bmat_realspace[0] = 1e10*ap_vec[0];
                                Bmat_realspace[3] = 1e10*ap_vec[1];
                                Bmat_realspace[6] = 1e10*ap_vec[2];

                                Bmat_realspace[1] = 1e10*bp_vec[0];
                                Bmat_realspace[4] = 1e10*bp_vec[1];
                                Bmat_realspace[7] = 1e10*bp_vec[2];

                                Bmat_realspace[2] = 1e10*cp_vec[0];
                                Bmat_realspace[5] = 1e10*cp_vec[1];
                                Bmat_realspace[8] = 1e10*cp_vec[2];

                                /* construct fractional Miller indicies */
                                UBO = (UMATS_RXYZ[mos_tic] * Umatrix*Bmat_realspace*(Omatrix.transpose())).transpose();
                                H_vec = UBO * q_vec;
                                //vec3 H_vec = (UMATS_RXYZ[mos_tic] * Bmat_realspace).transpose() * q_vec;
                                h = H_vec[0];
                                k = H_vec[1];
                                l = H_vec[2];

                                /* round off to nearest whole index */
                                h0 = static_cast<int>(ceil(h-0.5));
                                k0 = static_cast<int>(ceil(k-0.5));
                                l0 = static_cast<int>(ceil(l-0.5));

                                H0_vec[0] = h0;
                                H0_vec[1] = k0;
                                H0_vec[2] = l0;

                                NABC[0] = Na;
                                NABC[1] = 0;
                                NABC[2] = 0;
                                NABC[3] = 0;
                                NABC[4] = Nb;
                                NABC[5] = 0;
                                NABC[6] = 0;
                                NABC[7] = 0;
                                NABC[8] = Nc;

                                vec3 V = NABC*(H_vec- H0_vec);

                                //if (mos_tic==1 && fpixel==10 && spixel==10)
                                //  //printf("AAAAAAAAAAAA: %f, %f, %f \n", AA[0]*1e10, AA[1]*1e10, AA[2]*1e10);

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
                                        hrad_sqr = V*V;

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
                                /* no need to go further if result will be zero */
                                if(F_latt == 0.0) continue;

                                /* structure factor of the unit cell */
                                if(interpolate){
                                    h0_flr = static_cast<int>(floor(h));
                                    k0_flr = static_cast<int>(floor(k));
                                    l0_flr = static_cast<int>(floor(l));


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
                                        interpolate=0;
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

                                if(! interpolate)
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

                                //F_cell = Fhkl[h0-h_min][k0-k_min][l0-l_min];

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
                                //polar=1; //FIXME

                                /* convert amplitudes into intensity (photons per steradian) */
                                I += F_cell*F_cell*F_latt*F_latt*source_I[source]*capture_fraction*omega_pixel;

                                /* checkpoint for rotataion derivatives */
                                for (int i_rot =0 ; i_rot < 3 ; i_rot++){
                                    if (rot_managers[i_rot]->refine_me){
                                        R3[i_rot] = dRotMats[i_rot]; // TODO: design upgrade
                                        R3_2[i_rot] = d2RotMats[i_rot]; // TODO: design upgrade
                                        rot_managers[i_rot]->increment(
                                                            fudge,
                                                            R3[0], R3[1], R3[2],
                                                            R3_2[0], R3_2[1], R3_2[2],
                                                            NABC, UMATS[mos_tic],
                                                            Umatrix*Bmat_realspace*(Omatrix.transpose()),
                                                            q_vec, V,
                                                            hrad_sqr, F_cell, F_latt,
                                                            source_I[source], capture_fraction, omega_pixel);
                                        R3[i_rot] = RotMats[i_rot];
                                        R3_2[i_rot] = RotMats[i_rot]; // TODO: design upgrade
                                    }
                                } //* end rotation manager deriv */

                                /*Checkpoint for unit cell derivatives*/
                                for(int i_uc=0; i_uc < 6; i_uc++ ){
                                    if (ucell_managers[i_uc]->refine_me){
                                        ucell_managers[i_uc]->increment(
                                            V, NABC, UMATS_RXYZ[mos_tic]*Umatrix, q_vec,
                                            Omatrix.transpose(),
                                            hrad_sqr, F_cell, F_latt, fudge,
                                            source_I[source], capture_fraction, omega_pixel);
                                    }
                                } /*end ucell deriv */

                                /* Checkpoint for Ncells manager */
                                if (Ncells_managers[0]->refine_me){
                                    Ncells_managers[0]->increment(
                                            V, H0_vec, H_vec,
                                            hrad_sqr, F_cell, F_latt, fudge,
                                            source_I[source], capture_fraction, omega_pixel);
                                } /* end Ncells manager deriv */

                                /* Checkpoint for Origin manager */
                                if (origin_managers[0]->refine_me){
                                    origin_managers[0]->increment(
                                        V,  NABC, UBO, k_diffracted, o_vec, airpath, lambda*1e10,
                                        hrad_sqr, F_cell, F_latt, fudge,
                                        source_I[source], capture_fraction, omega_pixel, pixel_size);
                                } /* end origin manager deriv */

                                /* checkpoint for Fcell manager */
                                if (fcell_man->refine_me){
                                    fcell_man->increment(
                                            hrad_sqr, F_cell, F_latt, fudge,
                                            source_I[source], capture_fraction, omega_pixel);
                                } /* end of fcell man deriv */

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

            floatimage[i] += r_e_sqr*fluence*spot_scale*polar*I/steps;

            int roi_fs = i % fpixels - roi_xmin;
            int roi_ss = floor(i / fpixels) - roi_ymin ;
            int roi_i =  roi_ss*roi_fdim + roi_fs ;
            //SCITBX_EXAMINE(i);
            //SCITBX_EXAMINE(fpixels);
            //SCITBX_EXAMINE(spixels);
            //SCITBX_EXAMINE(roi_xmax);
            //SCITBX_EXAMINE(roi_xmin);
            //SCITBX_EXAMINE(roi_fs);
            //SCITBX_EXAMINE(roi_ss);
            //SCITBX_EXAMINE(roi_i);
            //SCITBX_ASSERT( roi_fs < roi_fdim);
            //SCITBX_ASSERT( roi_ss < roi_sdim);
            //SCITBX_ASSERT(roi_i >= 0);
            //SCITBX_ASSERT(roi_i < (roi_fdim*roi_sdim) ) ;

            floatimage_roi[roi_i] += r_e_sqr*fluence*spot_scale*polar*I/steps;

            /* udpate the rotation derivative images*/
            for (int i_rot =0 ; i_rot < 3 ; i_rot++){
                if (rot_managers[i_rot]->refine_me){
                    double value = r_e_sqr*fluence*spot_scale*polar*rot_managers[i_rot]->dI/steps;
                    double value2 = r_e_sqr*fluence*spot_scale*polar*rot_managers[i_rot]->dI2/steps;
                    rot_managers[i_rot]->increment_image(roi_i, value, value2);
                }
            } /* end rot deriv image increment */

            /*update the ucell derivative images*/
            for (int i_uc=0 ; i_uc < 6 ; i_uc++){
                if (ucell_managers[i_uc]->refine_me){
                    double value = r_e_sqr*fluence*spot_scale*polar*ucell_managers[i_uc]->dI/steps;
                    double value2 = r_e_sqr*fluence*spot_scale*polar*ucell_managers[i_uc]->dI2/steps;
                    ucell_managers[i_uc]->increment_image(roi_i, value, value2);
                }
            }/* end ucell deriv image increment */

            /*update the Ncells derivative image*/
            if (Ncells_managers[0]->refine_me){
                double value = r_e_sqr*fluence*spot_scale*polar*Ncells_managers[0]->dI/steps;
                double value2 = r_e_sqr*fluence*spot_scale*polar*Ncells_managers[0]->dI2/steps;
                Ncells_managers[0]->increment_image(roi_i, value, value2);
            }/* end Ncells deriv image increment */

            /* update Fcell derivative image */
            if(fcell_man->refine_me){
                double value = r_e_sqr*fluence*spot_scale*polar*fcell_man->dI/steps;
                double value2 = r_e_sqr*fluence*spot_scale*polar*fcell_man->dI2/steps;
                Ncells_managers[0]->increment_image(roi_i, value, value2);
            }/* end Fcell deriv image increment */

            /* update origin derivative image */
            if (origin_managers[0]->refine_me){
                double value = r_e_sqr*fluence*spot_scale*polar*  origin_managers[0]->dI  /steps;
                vec3 dk = origin_managers[0]->dk;
                double D = k_diffracted*o_vec; //origin_managers[0]->value; /* TODO check me */
                double dD = dk * o_vec;
                //SCITBX_EXAMINE(dD);
                //SCITBX_EXAMINE(dk*k_diffracted);
                double du = 2*D/airpath*(dD/airpath - D/pow(airpath,3)*(dk*k_diffracted));
                //double dpsi = 0;
                /* kahn factor is the variable called polarization */
                //if(polarization)

                double dpsi = -1/(1+kBi*kBi/kEi/kEi) * ((dk*Bi_vec)/kEi - kBi/kEi/kEi*(dk*Ei_vec));
                double dpolar = .5*(du*(1 + polarization * cos(2*psi)) + 2*dpsi*sin(2*psi)*polarization*(1-u));
                if (! nopolar)
                    value += r_e_sqr*fluence*spot_scale*  dpolar  *I/steps;

                double value2 = r_e_sqr*fluence*spot_scale*polar*origin_managers[0]->dI2/steps;
                origin_managers[0]->increment_image(roi_i, value, value2);
            } /*end origigin deriv image increment */

            if(floatimage[i] > max_I) {
                max_I = floatimage[i];
                max_I_x = Fdet;
                max_I_y = Sdet;
            }
            sum += floatimage[i];
            sumsqr += floatimage[i]*floatimage[i];
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
                    //printf("pixel   %15.10g\n", floatimage[i]);
                    printf("real-space cell vectors (Angstrom):\n");
                    printf("     %-10s  %-10s  %-10s\n","a","b","c");
                    printf("X: %11.8f %11.8f %11.8f\n",a[1]*1e10,b[1]*1e10,c[1]*1e10);
                    printf("Y: %11.8f %11.8f %11.8f\n",a[2]*1e10,b[2]*1e10,c[2]*1e10);
                    printf("Z: %11.8f %11.8f %11.8f\n",a[3]*1e10,b[3]*1e10,c[3]*1e10);
                    printf("Rot manager refine status X=%d, Y=%d, Z=%d\n",
                        rot_managers[0]->refine_me, rot_managers[1]->refine_me,
                        rot_managers[2]->refine_me);
                    printf("Ucell managers refine status a.a=%d, b.b=%d, c.c=%d, a.b=%d, a.c=%d, b.c=%d\n",
                        ucell_managers[0]->refine_me, ucell_managers[1]->refine_me, ucell_managers[2]->refine_me,
                        ucell_managers[3]->refine_me, ucell_managers[4]->refine_me, ucell_managers[5]->refine_me);
                    printf("Ncell managers refine status: %d, value=%f\n", Ncells_managers[0]->refine_me,
                            Ncells_managers[0]->value);
                    printf("Fcell manager refine status: %d\n", fcell_man->refine_me);
                    printf("Origin managers refine status: %d, value=%f\n", origin_managers[0]->refine_me,
                            origin_managers[0]->value);
                    printf("Bmatrix_real:\n%11.8f %11.8f %11.8f\n %11.8f %11.8f %11.8f\n %11.8f %11.8f %11.8f\n",
                        Bmat_realspace[0]*1e10, Bmat_realspace[1]*1e10, Bmat_realspace[2]*1e10,
                        Bmat_realspace[3]*1e10, Bmat_realspace[4]*1e10, Bmat_realspace[5]*1e10,
                        Bmat_realspace[6]*1e10, Bmat_realspace[7]*1e10, Bmat_realspace[8]*1e10);
                    printf("Umatrix_real:\n%11.8f %11.8f %11.8f\n %11.8f %11.8f %11.8f\n %11.8f %11.8f %11.8f\n",
                        Umatrix[0], Umatrix[1], Umatrix[2],
                        Umatrix[3], Umatrix[4], Umatrix[5],
                        Umatrix[6], Umatrix[7], Umatrix[8]);
                }
            }
            else
            {
                if(progress_meter && verbose && progress_pixels/100 > 0)
                {
                    if(progress_pixel % ( progress_pixels/20 ) == 0 ||
                       ((10*progress_pixel<progress_pixels ||
                         10*progress_pixel>9*progress_pixels) &&
                        (progress_pixel % (progress_pixels/100) == 0)))
                    {
                        printf("%lu%% done\n",progress_pixel*100/progress_pixels);
                    }
                }
                ++progress_pixel;
            }
            ++i;
        }
    }
    if(verbose) printf("done with pixel loop\n");

    if(verbose) printf("solid angle subtended by detector = %g steradian ( %g%% sphere)\n",omega_sum/steps,100*omega_sum/steps/4/M_PI);
    if(verbose) printf("max_I= %g sum= %g avg= %g\n",max_I,sum,sum/sumn);

} // END  of add_diffBragg_spots
// END diffBragg

} // end of namespace nanoBragg
} // end of namespace simtbx
