#include <simtbx/diffBragg/src/diffBragg.h>
#include <assert.h>
//#include <omp.h>
//#include <algorithm>
namespace simtbx {
namespace nanoBragg {

// BEGIN derivative manager
derivative_manager::derivative_manager(){}

void derivative_manager::initialize(int sdim, int fdim, bool curvatures)
{
    raw_pixels = af::flex_double(af::flex_grid<>(sdim,fdim));
    dI=0;

    // for second derivatives
    dI2=0;
    //if (curvatures)
    raw_pixels2 = af::flex_double(af::flex_grid<>(sdim,fdim));
}

void derivative_manager::increment_image(int idx, double value, double value2, bool curvatures){
    double* floatimage = raw_pixels.begin();
    floatimage[idx] += value;

    // increment second derivatives
    if (curvatures){
        double* floatimage2 = raw_pixels2.begin();
        floatimage2[idx] += value2;
    }
}
// END derivative manager

//BEGIN origin manager
origin_manager::origin_manager(){
    FF=0;
    FdF=0;
    dFdF=0;
    FdF2=0;
    }

void origin_manager::increment(
        vec3 V, mat3 N, mat3 UBO, vec3 k_diffracted, vec3 o_vec, double air_path, double wavelen,
        double Hrad, double Fcell, double Flatt, double fudge,
        double source_I, double capture_fraction, double omega_pixel, double pixel_size){

    double per_k = 1/air_path;
    double per_k3 = pow(per_k,3);
    double per_k5 = pow(per_k,5);

    double G = dk*k_diffracted;

    vec3 dQ = per_k * dk - per_k3 * G *k_diffracted;
    vec3 dQ2 = 3*per_k5*G*G*k_diffracted - per_k3 * ((dk*dk)*k_diffracted  + 2*G*dk);

    dQ /= wavelen;
    dQ2 /= wavelen;

    vec3 dV = N*(UBO*dQ);
    double V_dot_dV = V*dV;
    double dHrad = 2*V_dot_dV;
    double a = 1 / 0.63 * fudge;
    double dFlatt = -1*a*Flatt*dHrad;
    double c = Fcell*Fcell*source_I*capture_fraction;

    vec3 dV2 = N*(UBO*dQ2);
    double dFlatt2 = -2*a*(dFlatt * V_dot_dV + Flatt*(dV*dV) + Flatt*(V*dV2));

    FF += c*Flatt*Flatt;
    FdF += c*2*Flatt*dFlatt;
    dFdF += c*2*dFlatt*dFlatt;
    FdF2 += c*2*Flatt*dFlatt2;
};
//END origin manager

// BEGIN Ncells_abc manager
Ncells_manager::Ncells_manager(){}

void Ncells_manager::increment(double dI_increment, double dI2_increment){
    dI += dI_increment;
    dI2 += dI2_increment;
};

//END Ncells_abc manager

// Begin Fcell manager
Fcell_manager::Fcell_manager(){}

void Fcell_manager::increment( double value, double value2)
{
    dI += value;
    dI2 += value2;
};
// END Fcell manager

// begin lambda_manager
lambda_manager::lambda_manager(){}

void lambda_manager::increment( double value, double value2)
{
    dI += value;
    dI2 += value2;
};
//END lambda manager

// begin panel_manager
panel_manager::panel_manager(){}

void panel_manager::increment( double Iincrement, double omega_pixel, mat3 M,
    double pix2, vec3 o, vec3 k_diffracted, double per_k, double per_k3, double per_k5, vec3 V)
{
    double G = dk*k_diffracted;
    vec3 dk_hat = -per_k3*G*k_diffracted + per_k*dk;
    double coef = (M*dk_hat)*V;
    //printf("AFTER %f\n", coef);
    double coef2 = -3*pix2*per_k5*G * (o*k_diffracted);
    coef2 += pix2*per_k3*(o*dk); // + F_cross_dS*k_diffracted + dF_cross_S*k_diffracted );
    //coef2 += pix2*per_k3*(o*dk + F_cross_dS*k_diffracted + dF_cross_S*k_diffracted );
    double value = coef*Iincrement + coef2*Iincrement/omega_pixel;
    dI += value;
    dI2 += 0;
};


//END panel_manager

//BEGIN unit cell manager
ucell_manager::ucell_manager(){}

void ucell_manager::increment(double value, double value2)
{
  dI += value;
  dI2 += value2;
};

// BEGIN rotation manager begin
rot_manager::rot_manager(){}

void rot_manager::set_R(){assert (false);}

void rot_manager::increment(double value, double value2)
{
  dI += value;
  dI2 += value2;
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

    EYE = mat3(1,0,0,0,1,0,0,0,1);
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

    boost::shared_ptr<Ncells_manager> nc1 = boost::shared_ptr<Ncells_manager>(new Ncells_manager());
    boost::shared_ptr<Ncells_manager> nc2 = boost::shared_ptr<Ncells_manager>(new Ncells_manager());
    boost::shared_ptr<Ncells_manager> nc3 = boost::shared_ptr<Ncells_manager>(new Ncells_manager());

    boost::shared_ptr<lambda_manager> lam1 = boost::shared_ptr<lambda_manager>(new lambda_manager());
    boost::shared_ptr<lambda_manager> lam2 = boost::shared_ptr<lambda_manager>(new lambda_manager());

    boost::shared_ptr<panel_manager> orig0 = boost::shared_ptr<panel_manager>(new panel_manager());
    boost::shared_ptr<panel_manager> origX = boost::shared_ptr<panel_manager>(new panel_manager());
    boost::shared_ptr<panel_manager> origY = boost::shared_ptr<panel_manager>(new panel_manager());
    //boost::shared_ptr<origin_manager> orig0 = boost::shared_ptr<origin_manager>(new origin_manager());

    //boost::shared_ptr<Fcell_manager> fcell_man = boost::shared_ptr<Fcell_manager>(new Fcell_manager());
    fcell_man = boost::shared_ptr<Fcell_manager>(new Fcell_manager());
    fcell_man->refine_me = false;

    panel_rot_man = boost::shared_ptr<panel_manager>(new panel_manager());
    panel_rot_man->refine_me = false;
    panel_rot_man->F_cross_dS = vec3(0,0,0);
    panel_rot_man->dF_cross_S = vec3(0,0,0);
    panel_rot_manF = boost::shared_ptr<panel_manager>(new panel_manager());
    panel_rot_manF->refine_me = false;
    panel_rot_manF->F_cross_dS = vec3(0,0,0);
    panel_rot_manF->dF_cross_S = vec3(0,0,0);

    panel_rot_manS = boost::shared_ptr<panel_manager>(new panel_manager());
    panel_rot_manS->refine_me = false;
    panel_rot_manS->F_cross_dS = vec3(0,0,0);
    panel_rot_manS->dF_cross_S = vec3(0,0,0);

    panels.push_back(panel_rot_man);
    panels.push_back(orig0);
    panels.push_back(origX);
    panels.push_back(origY);
    panels.push_back(panel_rot_manF);
    panels.push_back(panel_rot_manS);

    rotX->refine_me = false;
    rotY->refine_me = false;
    rotZ->refine_me = false;
    uc1->refine_me = false;
    uc2->refine_me = false;
    uc3->refine_me = false;
    uc4->refine_me = false;
    uc5->refine_me = false;
    uc6->refine_me = false;

    nc1->refine_me = false;
    nc2->refine_me = false;
    nc3->refine_me = false;

    orig0->refine_me = false;
    orig0->dk = vec3(0,0,1);

    origX->refine_me = false;
    origX->dk = vec3(1,0,0);

    origY->refine_me = false;
    origY->dk = vec3(0,1,0);

    lam1->refine_me= false;
    lam2->refine_me= false;

    lambda_managers.push_back(lam1);
    lambda_managers.push_back(lam2);

    rot_managers.push_back(rotX);
    rot_managers.push_back(rotY);
    rot_managers.push_back(rotZ);

    ucell_managers.push_back(uc1);
    ucell_managers.push_back(uc2);
    ucell_managers.push_back(uc3);
    ucell_managers.push_back(uc4);
    ucell_managers.push_back(uc5);
    ucell_managers.push_back(uc6);

    Ncells_managers.push_back(nc1);
    Ncells_managers.push_back(nc2);
    Ncells_managers.push_back(nc3);
    NABC = mat3(0,0,0,0,0,0,0,0,0);
    NABC[0] = 1;
    NABC[4] = 1;
    NABC[8] = 1;

    O_reference = vec3(0,0,0);

    update_oversample_during_refinement = true;
    oversample_omega = true;
    only_save_omega_kahn = false;
    compute_curvatures = true;
    isotropic_ncells = true;

    lambda_managers[0]->value = 0;
    lambda_managers[1]->value = 1;
    use_lambda_coefficients = false;
    //source_lambda0 = 0;
    //source_lambda1 = 1;

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

    max_I_hkl = vec3(0,0,0);
    init_raw_pixels_roi();
    initialize_managers();

    Fhkl2 = NULL;
    F_cell2 = 0;
    complex_miller = false;
    pythony_amplitudes2.clear();
    }


af::flex_double diffBragg::get_panel_increment( double Iincrement, double omega_pixel, mat3 M,
    double pix2, vec3 o, vec3 k_diffracted, double per_k, double per_k3, double per_k5, vec3 V, vec3 _dk)
{
    double G = _dk*k_diffracted;
    vec3 dk_hat = -per_k3*G*k_diffracted + per_k*_dk;
    double coef = (M*dk_hat)*V;
    //printf("AFTER %f\n", coef);
    double coef2 = -3*pix2*per_k5*G * (o*k_diffracted);
    coef2 += pix2*per_k3*(o*_dk); // + F_cross_dS*k_diffracted + dF_cross_S*k_diffracted );
    //coef2 += pix2*per_k3*(o*_dk + F_cross_dS*k_diffracted + dF_cross_S*k_diffracted );
    double value = coef*Iincrement + coef2*Iincrement/omega_pixel;
    af::flex_double values = af::flex_double(2,0);
    values[0] = value;
    values[1] = 0;
    return values;
};



/* BEGIN panel Rot XYZ */
void diffBragg::rotate_fs_ss_vecs_3D(double panel_rot_angO, double panel_rot_angF, double panel_rot_angS){

    vec3 fs_vec = vec3(fdet_vector[1], fdet_vector[2], fdet_vector[3]);
    vec3 ss_vec = vec3(sdet_vector[1], sdet_vector[2], sdet_vector[3]);
    vec3 origin_vec = vec3(pix0_vector[1], pix0_vector[2], pix0_vector[3]);
    vec3 origin_diff_vec = origin_vec - O_reference; // difference between origin and a reference vector, this is the vector we will rotate

    mat3 RO = mat3(0,0,0,0,0,0,0,0,0);
    RO[1] = -odet_vector[3];
    RO[2] = odet_vector[2];
    RO[3] = odet_vector[3];
    RO[5] = -odet_vector[1];
    RO[6] = -odet_vector[2];
    RO[7] = odet_vector[1];
    mat3 RO2 = RO*RO;

    mat3 RF = mat3(0,0,0,0,0,0,0,0,0);
    RF[1] = -fdet_vector[3];
    RF[2] = fdet_vector[2];
    RF[3] = fdet_vector[3];
    RF[5] = -fdet_vector[1];
    RF[6] = -fdet_vector[2];
    RF[7] = fdet_vector[1];
    mat3 RF2 = RF*RF;

    mat3 RS = mat3(0,0,0,0,0,0,0,0,0);
    RS[1] = -sdet_vector[3];
    RS[2] = sdet_vector[2];
    RS[3] = sdet_vector[3];
    RS[5] = -sdet_vector[1];
    RS[6] = -sdet_vector[2];
    RS[7] = sdet_vector[1];
    mat3 RS2 = RS*RS;

    mat3 rotO = EYE + RO*sin(panel_rot_angO) + RO2*(1-cos(panel_rot_angO));
    mat3 rotF = EYE + RF*sin(panel_rot_angF) + RF2*(1-cos(panel_rot_angF));
    mat3 rotS = EYE + RS*sin(panel_rot_angS) + RS2*(1-cos(panel_rot_angS));

    boost::shared_ptr<panel_manager> pan_rot = boost::dynamic_pointer_cast<panel_manager>(panels[0]);
    pan_rot->dR = RO*cos(panel_rot_angO) + RO2*sin(panel_rot_angO);
    pan_rot->value = panel_rot_angO;
    mat3 ROT = (pan_rot->dR)*rotF*rotS;
    pan_rot->dF = ROT*fs_vec;
    pan_rot->dS = ROT*ss_vec;

    pan_rot = boost::dynamic_pointer_cast<panel_manager>(panels[4]);
    pan_rot->dR = RF*cos(panel_rot_angF) + RF2*sin(panel_rot_angF);
    pan_rot->value = panel_rot_angF;
    ROT = rotO*(pan_rot->dR)*rotS;
    pan_rot->dF = ROT*fs_vec;
    pan_rot->dS = ROT*ss_vec;

    pan_rot = boost::dynamic_pointer_cast<panel_manager>(panels[5]);
    pan_rot->dR = RS*cos(panel_rot_angS) + RS2*sin(panel_rot_angS);
    pan_rot->value = panel_rot_angS;
    ROT = rotO*rotF*(pan_rot->dR);
    pan_rot->dF = ROT*fs_vec;
    pan_rot->dS = ROT*ss_vec;

    ROT = rotO*rotF*rotS;
    fs_vec = ROT*fs_vec;
    ss_vec = ROT*ss_vec;
    origin_vec = O_reference + ROT*origin_diff_vec; // add the rotated difference to the reference to recover the origin of the rotated panel
    for (int i=0; i < 3; i++){
        fdet_vector[i+1] = fs_vec[i];
        sdet_vector[i+1] = ss_vec[i];
        pix0_vector[i+1] = origin_vec[i];
    }

    //int pan_mans[3] = {0,4,5};
    //for (int i_rot=0; i_rot <3; i_rot++ ){
    //    int panels_id = pan_mans[i_rot];
    //    pan_rot = boost::dynamic_pointer_cast<panel_manager>(panels[panels_id]);
    //    pan_rot->F_cross_dS = fs_vec.cross(pan_rot->dS);
    //    pan_rot->dF_cross_S = (pan_rot->dF).cross(ss_vec);
    //}


}
/* end panel Rot XYZ */


void diffBragg::update_dxtbx_geoms(
    const dxtbx::model::Detector& detector,
    const dxtbx::model::Beam& beam,
    int panel_id,
    double panel_rot_angO,
    double panel_rot_angF,
    double panel_rot_angS, double panel_offsetX, double panel_offsetY, double panel_offsetZ,
    bool force){

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
    SCITBX_ASSERT(fpixels == detector[panel_id].get_image_size()[0]);

    /* direction in 3-space of detector axes */
    SCITBX_ASSERT (beam_convention == CUSTOM);

    fdet_vector[1] = detector[panel_id].get_fast_axis()[0];
    fdet_vector[2] = detector[panel_id].get_fast_axis()[1];
    fdet_vector[3] = detector[panel_id].get_fast_axis()[2];
    unitize(fdet_vector, fdet_vector);

    sdet_vector[1] = detector[panel_id].get_slow_axis()[0];
    sdet_vector[2] = detector[panel_id].get_slow_axis()[1];
    sdet_vector[3] = detector[panel_id].get_slow_axis()[2];
    unitize(sdet_vector,sdet_vector);
    /* set orthogonal vector to the detector pixel array */
    cross_product(fdet_vector,sdet_vector,odet_vector);
    unitize(odet_vector,odet_vector);
    if (! detector_is_righthanded)
        vector_scale(odet_vector, odet_vector, -1);

    /* dxtbx origin is location of outer corner of the first pixel */
    pix0_vector[1] = detector[panel_id].get_origin()[0]/1000.0;
    pix0_vector[2] = detector[panel_id].get_origin()[1]/1000.0;
    pix0_vector[3] = detector[panel_id].get_origin()[2]/1000.0;

    //if (panel_rot_ang != 0)
    //rotate_fs_ss_vecs(panel_rot_ang);
    int n_rot_refine = 0;
    if (panels[0]->refine_me)
        n_rot_refine ++;
    if (panels[4]->refine_me)
        n_rot_refine ++;
    if (panels[5]->refine_me)
        n_rot_refine ++;
    if (n_rot_refine  > 0 || force){
        rotate_fs_ss_vecs_3D(panel_rot_angO, panel_rot_angF, panel_rot_angS);

        /* reset orthogonal vector to the detector pixel array after rotation */
        cross_product(fdet_vector,sdet_vector,odet_vector);
        unitize(odet_vector,odet_vector);
        if (! detector_is_righthanded)
            vector_scale(odet_vector, odet_vector, -1);

        boost::shared_ptr<panel_manager> pan_rot;
        vec3 fs_vec = vec3(fdet_vector[1], fdet_vector[2], fdet_vector[3]);
        vec3 ss_vec = vec3(sdet_vector[1], sdet_vector[2], sdet_vector[3]);
        int pan_mans[3] = {0,4,5};
        for (int i_rot=0; i_rot <3; i_rot++ ){
            int panels_id = pan_mans[i_rot];
            pan_rot = boost::dynamic_pointer_cast<panel_manager>(panels[panels_id]);
            if (!detector_is_righthanded){
                pan_rot->F_cross_dS = (pan_rot->dS).cross(fs_vec);
                pan_rot->dF_cross_S = ss_vec.cross(pan_rot->dF);
            }
            else{
                pan_rot->F_cross_dS = fs_vec.cross(pan_rot->dS);
                pan_rot->dF_cross_S = (pan_rot->dF).cross(ss_vec);
            }
        }
    }

    pix0_vector[1] += panel_offsetX;///1000;
    pix0_vector[2] += panel_offsetY;///1000;
    pix0_vector[3] += panel_offsetZ;///1000;

    Fclose = Xclose = -dot_product(pix0_vector, fdet_vector);
    Sclose = Yclose = -dot_product(pix0_vector, sdet_vector);
    close_distance = distance = dot_product(pix0_vector, odet_vector);

    /* set beam centre */
    mat3 dmat = mat3(fdet_vector[1], sdet_vector[1], pix0_vector[1]*1000,
                    fdet_vector[2], sdet_vector[2], pix0_vector[2]*1000,
                    fdet_vector[3], sdet_vector[3], pix0_vector[3]*1000);
    mat3 Dmat = dmat.inverse();
    vec3 s0 = vec3(beam.get_s0()[0], beam.get_s0()[1], beam.get_s0()[2]);
    vec3 dxtbx_v = Dmat*s0;
    SCITBX_ASSERT(dxtbx_v[2] > 0);

    double rotated_center_x = dxtbx_v[0] / dxtbx_v[2];
    double rotated_center_y = dxtbx_v[1] / dxtbx_v[2];

    scitbx::vec2<double> dials_bc = detector[panel_id].get_beam_centre(beam.get_s0());
    dials_bc[0] = rotated_center_x;
    dials_bc[1] = rotated_center_y;
    Xbeam = dials_bc[0]/1000.0;
    Ybeam = dials_bc[1]/1000.0;

    /* detector sensor layer properties */
    detector_thick = detector[panel_id].get_thickness();
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
    
    user_beam=true;

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
    SCITBX_ASSERT(close_distance > 0);
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
            rot_managers[i_rot]->initialize(sdim, fdim, compute_curvatures);
    }
    for (int i_uc=0; i_uc < 6; i_uc++){
        if (ucell_managers[i_uc]->refine_me)
            ucell_managers[i_uc]->initialize(sdim, fdim, compute_curvatures);
    }
    for (int i_nc=0; i_nc < 3; i_nc ++){
        if (Ncells_managers[i_nc]->refine_me)
            Ncells_managers[i_nc]->initialize(sdim, fdim, compute_curvatures);
    }

    boost::shared_ptr<panel_manager> pan_orig;
    for (int i_pan_orig=0; i_pan_orig  < 3; i_pan_orig ++){
        pan_orig = boost::dynamic_pointer_cast<panel_manager>(panels[1+i_pan_orig]);
        if (pan_orig->refine_me)
            pan_orig->initialize(sdim, fdim, compute_curvatures);
    }

    if (fcell_man->refine_me)
        fcell_man->initialize(sdim, fdim, compute_curvatures);

    for (int i_lam=0; i_lam < 2; i_lam++){
        if (lambda_managers[i_lam]->refine_me)
            lambda_managers[i_lam]->initialize(sdim, fdim, compute_curvatures);
    }

    int panrot_manager_indices[3] = {0,4,5};
    boost::shared_ptr<panel_manager> pan_rot;
    for (int i=0; i< 3; i++){
        int manager_idx = panrot_manager_indices[i];
        pan_rot = boost::dynamic_pointer_cast<panel_manager>(panels[manager_idx]);
        if (pan_rot->refine_me)
            pan_rot->initialize(sdim, fdim, compute_curvatures);
    }
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
        rot_managers[refine_id]->initialize(sdim, fdim, compute_curvatures);
    }
    else if (refine_id >=3 and refine_id < 9 ){
        // 6 possible unit cell managers (a,b,c,al,be,ga)
        ucell_managers[refine_id-3]->refine_me=true;
        ucell_managers[refine_id-3]->initialize(sdim, fdim, compute_curvatures);
    }
    else if (refine_id==9){
        for (int i_nc=0; i_nc < 3; i_nc ++){
            Ncells_managers[i_nc]->refine_me=true;
            Ncells_managers[i_nc]->initialize(sdim, fdim, compute_curvatures);
        }
    }
    else if (refine_id==10){
        boost::shared_ptr<panel_manager> pan_orig = boost::dynamic_pointer_cast<panel_manager>(panels[1]);
        pan_orig->refine_me=true;
        pan_orig->initialize(sdim, fdim, compute_curvatures);
    }
    else if(refine_id==11){
        fcell_man->refine_me=true;
        fcell_man->initialize(sdim, fdim, compute_curvatures);
    }

    else if (refine_id==12 || refine_id==13){
        use_lambda_coefficients = true;
        int i_lam = refine_id-12;
        lambda_managers[i_lam]->refine_me=true;
        lambda_managers[i_lam]->initialize(sdim, fdim, compute_curvatures);
    }

    else if (refine_id==14){
        boost::shared_ptr<panel_manager> pan_rot = boost::dynamic_pointer_cast<panel_manager>(panels[0]);
        pan_rot->refine_me=true;
        rotate_fs_ss_vecs_3D(0,0,0);
        pan_rot->initialize(sdim, fdim, compute_curvatures);
    }

    else if (refine_id==15){
        boost::shared_ptr<panel_manager> pan_orig = boost::dynamic_pointer_cast<panel_manager>(panels[2]);
        pan_orig->refine_me=true;
        pan_orig->initialize(sdim, fdim, compute_curvatures);
    }

    else if (refine_id==16){
        boost::shared_ptr<panel_manager> pan_orig = boost::dynamic_pointer_cast<panel_manager>(panels[3]);
        pan_orig->refine_me=true;
        pan_orig->initialize(sdim, fdim, compute_curvatures);
    }
    else if (refine_id==17){
        boost::shared_ptr<panel_manager> pan_rot = boost::dynamic_pointer_cast<panel_manager>(panels[4]);
        pan_rot->refine_me=true;
        rotate_fs_ss_vecs_3D(0,0,0);
        pan_rot->initialize(sdim, fdim, compute_curvatures);
    }
    else if (refine_id==18){
        boost::shared_ptr<panel_manager> pan_rot = boost::dynamic_pointer_cast<panel_manager>(panels[5]);
        pan_rot->refine_me=true;
        rotate_fs_ss_vecs_3D(0,0,0);
        pan_rot->initialize(sdim, fdim, compute_curvatures);
    }
}

void diffBragg::init_Fhkl2()
{
    /* This should only be called if init_Fhkl has already been called with python indices/amplitudes*/
    /* free any previous allocations */
    if(Fhkl2 != NULL) {
        for (h0=0; h0<=h_range;h0++) {
            for (k0=0; k0<=k_range;k0++) {
                free(Fhkl2[h0][k0]);
            }
            free(Fhkl2[h0]);
        }
        free(Fhkl2);
    }

    /* allocate memory for 3d arrays */
    Fhkl2 = (double***) calloc(h_range+1,sizeof(double**));
    if(Fhkl2==NULL){perror("ERROR");exit(9);};
    for (h0=0; h0<=h_range;h0++) {
        Fhkl2[h0] = (double**) calloc(k_range+1,sizeof(double*));
        if(Fhkl2[h0]==NULL){perror("ERROR");exit(9);};
        for (k0=0; k0<=k_range;k0++) {
            Fhkl2[h0][k0] = (double*) calloc(l_range+1,sizeof(double));
            if(Fhkl2[h0][k0]==NULL){perror("ERROR");exit(9);};
        }
    }
    for (h0=0; h0<h_range;h0++) {
        for (k0=0; k0<k_range;k0++) {
            for (l0=0; l0<l_range;l0++) {
                Fhkl2[h0][k0][l0] = 0;
            }
        }
    }
    Fhkl2[-h_min][-k_min][-l_min] = 0;

    if(verbose) printf("initializing Fhkl2 with pythony indices and amplitudes\n");
    miller_t hkl;
    for (i=0; i < pythony_indices.size(); ++i)
    {
        hkl = pythony_indices[i];
        F_cell = pythony_amplitudes2[i];
        h0 = hkl[0];
        k0 = hkl[1];
        l0 = hkl[2];
        Fhkl2[h0-h_min][k0-k_min][l0-l_min]=F_cell;
    }
    Fhkl2[-h_min][-k_min][-l_min] = 0;
    if(verbose) printf("done initializing Fhkl2:\n");
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
void diffBragg::set_ncells_values( boost::python::tuple const& values){
    Na = boost::python::extract<double>(values[0]);
    Nb = boost::python::extract<double>(values[1]);
    Nc = boost::python::extract<double>(values[2]);
    Ncells_managers[0]->value=Na;
    Ncells_managers[1]->value=Nb;
    Ncells_managers[2]->value=Nc;
    xtal_size_x = -1;
    xtal_size_y = -1;
    xtal_size_z = -1;
    if (update_oversample_during_refinement)
       update_oversample();
    NABC[0] = Na;
    NABC[4] = Nb;
    NABC[8] = Nc;
}

boost::python::tuple diffBragg::get_ncells_values(){
    boost::python::tuple values;
    values = boost::python::make_tuple( NABC[0],  NABC[4], NABC[8]);
    return values;
}


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
        if (update_oversample_during_refinement)
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

    SCITBX_ASSERT(refine_id >=0 && refine_id <= 18);

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
    else if (refine_id==10){
        boost::shared_ptr<panel_manager> pan_orig = boost::dynamic_pointer_cast<panel_manager>(panels[1]);
        return pan_orig->raw_pixels;
        }
    else if (refine_id==11)
        return fcell_man->raw_pixels;
    else if (refine_id==12)
        return lambda_managers[0]->raw_pixels;
    else if  (refine_id==13)
        return lambda_managers[1]->raw_pixels;
    else if (refine_id==14){
        boost::shared_ptr<panel_manager> pan_rot = boost::dynamic_pointer_cast<panel_manager>(panels[0]);
        return pan_rot->raw_pixels;
    }
    else if (refine_id==15){
        boost::shared_ptr<panel_manager> pan_orig = boost::dynamic_pointer_cast<panel_manager>(panels[2]);
        return pan_orig->raw_pixels;
    }
    else if(refine_id==16){
        boost::shared_ptr<panel_manager> pan_orig = boost::dynamic_pointer_cast<panel_manager>(panels[3]);
        return pan_orig->raw_pixels;
    }
    else if(refine_id==17){
        boost::shared_ptr<panel_manager> pan_rot = boost::dynamic_pointer_cast<panel_manager>(panels[4]);
        return pan_rot->raw_pixels;
    }
    else{
        boost::shared_ptr<panel_manager> pan_rot = boost::dynamic_pointer_cast<panel_manager>(panels[5]);
        return pan_rot->raw_pixels;
    }
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
    else if (refine_id==10){
        boost::shared_ptr<panel_manager> pan_orig = boost::dynamic_pointer_cast<panel_manager>(panels[1]);
        return pan_orig->raw_pixels2;}
    else if (refine_id==15){
        boost::shared_ptr<panel_manager> pan_orig = boost::dynamic_pointer_cast<panel_manager>(panels[2]);
        return pan_orig->raw_pixels2;}
    else if (refine_id==16){
        boost::shared_ptr<panel_manager> pan_orig = boost::dynamic_pointer_cast<panel_manager>(panels[3]);
        return pan_orig->raw_pixels2;}
    else
        return fcell_man->raw_pixels2;
}

boost::python::tuple diffBragg::get_ncells_derivative_pixels(){
    SCITBX_ASSERT(Ncells_managers[0]->refine_me);
    SCITBX_ASSERT(Ncells_managers[1]->refine_me);
    SCITBX_ASSERT(Ncells_managers[2]->refine_me);
    boost::python::tuple derivative_pixels;
    derivative_pixels = boost::python::make_tuple(Ncells_managers[0]->raw_pixels,
        Ncells_managers[1]->raw_pixels, Ncells_managers[2]->raw_pixels);
    return derivative_pixels;
}

boost::python::tuple diffBragg::get_lambda_derivative_pixels(){
    SCITBX_ASSERT(lambda_managers[0]->refine_me || lambda_managers[1]->refine_me);

    boost::python::tuple derivative_pixels;
    if (lambda_managers[0]->refine_me && lambda_managers[1]->refine_me){
        derivative_pixels = boost::python::make_tuple(lambda_managers[0]->raw_pixels,
            lambda_managers[1]->raw_pixels);
    }
    else{
        if (lambda_managers[0]->refine_me)
            derivative_pixels = boost::python::make_tuple(lambda_managers[0]->raw_pixels);
        else if (lambda_managers[1]->refine_me)
            derivative_pixels = boost::python::make_tuple(lambda_managers[1]->raw_pixels);
    }
    return derivative_pixels;
}

boost::python::tuple diffBragg::get_ncells_second_derivative_pixels(){
    SCITBX_ASSERT(Ncells_managers[0]->refine_me);
    SCITBX_ASSERT(Ncells_managers[1]->refine_me);
    SCITBX_ASSERT(Ncells_managers[2]->refine_me);
    boost::python::tuple second_derivative_pixels;
    second_derivative_pixels = boost::python::make_tuple(Ncells_managers[0]->raw_pixels2, Ncells_managers[1]->raw_pixels2, Ncells_managers[2]->raw_pixels2);
    return second_derivative_pixels;
}

void diffBragg::zero_raw_pixel_rois(){
    init_raw_pixels_roi();
    initialize_managers();
}

/* polarization factor */
/* override this to store variables needed for derivatives */
double diffBragg::polarization_factor(double kahn_factor, double *__incident, double *__diffracted, double *__axis)
{
    double cos2theta,cos2theta_sqr,sin2theta_sqr;
    //double psi=0.0;
    double E_in[4];
    double B_in[4];

    unitize(__incident,__incident);
    unitize(__diffracted,__diffracted);
    //unitize(__axis,__axis);

    /* component of diffracted unit vector along incident beam unit vector */
    cos2theta = dot_product(__incident,__diffracted);
    cos2theta_sqr = cos2theta*cos2theta;
    sin2theta_sqr = 1-cos2theta_sqr;
    //double _u = cos2theta_sqr;

    double _psi=0;
    if(kahn_factor != 0.0){
        //SCITBX_EXAMINE(kahn_factor);
        /* tricky bit here is deciding which direciton the E-vector lies in for each source
           here we assume it is closest to the "axis" defined above */

        /* cross product to get "vertical" axis that is orthogonal to the cannonical "polarization" */
        cross_product(__axis,__incident,B_in);
        /* make it a unit vector */
        unitize(B_in,B_in);
        //vec3 Bi_vec = vec3(B_in[1], B_in[2], B_in[3]);
        //Bi_vec[0] = B_in[1];
        //Bi_vec[1] = B_in[2];
        //Bi_vec[2] = B_in[3];

        /* cross product with incident beam to get E-vector direction */
        cross_product(__incident,B_in,E_in);
        /* make it a unit vector */
        unitize(E_in,E_in);
        //vec3 Ei_vec = vec3(E_in[1], E_in[2], E_in[3]);
        //Ei_vec[0] = E_in[1];
        //Ei_vec[1] = E_in[2];
        //Ei_vec[2] = E_in[3];

        /* get components of diffracted ray projected onto the E-B plane */
        double _kEi = dot_product(__diffracted,E_in);
        double _kBi = dot_product(__diffracted,B_in);

        /* compute the angle of the diffracted ray projected onto the incident E-B plane */
        _psi = -atan2(_kBi,_kEi);
    }
    /* correction for polarized incident beam */
    return 0.5*(1.0 + cos2theta_sqr - kahn_factor*cos(2*_psi)*sin2theta_sqr);
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
    //int roi_sdim = roi_ymax - roi_ymin+1;
    
    /* ravel the step loops into a single loop */
    int spixel_max = std::min( roi_ymax, spixels-1);
    int fpixel_max = std::min( roi_xmax, fpixels-1);
    int nspix_steps = spixel_max - roi_ymin + 1;
    int nfpix_steps = fpixel_max - roi_xmin + 1;
    //int Nsteps = nspix_steps*nfpix_steps*oversample*oversample*detector_thicksteps*sources*phisteps*mosaic_domains;
    int Nsteps = oversample*oversample*detector_thicksteps*sources*phisteps*mosaic_domains;
    int i_step = 0;
    //af::flex_int spix_pos = af::flex_int(Nsteps,0);
    //af::flex_int fpix_pos = af::flex_int(Nsteps,0);
    af::flex_int subS_pos = af::flex_int(Nsteps,0);
    af::flex_int subF_pos = af::flex_int(Nsteps,0);
    af::flex_int thick_pos = af::flex_int(Nsteps,0);
    af::flex_int source_pos = af::flex_int(Nsteps,0);
    af::flex_int phi_pos = af::flex_int(Nsteps,0);
    af::flex_int mos_pos = af::flex_int(Nsteps,0);

    //i = spixel*fpixels + fpixel;
    /* TODO theres probably a clever way to do this, but oh well */
    // TODO: time me
    //for(spixel=roi_ymin;spixel <= spixel_max;++spixel){
    //    for(fpixel=roi_xmin;fpixel <= fpixel_max;++fpixel){
            for(subS=0;subS<oversample;++subS){
                for(subF=0;subF<oversample;++subF){
                    for(thick_tic=0;thick_tic<detector_thicksteps;++thick_tic){
                        for(source=0;source<sources;++source){
                            for(phi_tic = 0; phi_tic < phisteps; ++phi_tic){
                                for(mos_tic=0;mos_tic<mosaic_domains;++mos_tic){
                                    //spix_pos[i_step] = spixel;
                                    //fpix_pos[i_step] = fpixel;
                                    subS_pos[i_step] = subS;
                                    subF_pos[i_step] = subF;
                                    thick_pos[i_step] = thick_tic;
                                    source_pos[i_step] = source;
                                    phi_pos[i_step] = phi_tic;
                                    mos_pos[i_step] = mos_tic;
                                    i_step ++;
                                }
                            }
                        }
                    }
                }
            }
        //}
    //}

    sum = sumsqr = 0.0;
    sumn = 0;
    progress_pixel = 0;
    omega_sum = 0.0;
    vec3 spindle_vec = vec3(spindle_vector[1], spindle_vector[2], spindle_vector[3]);

    //int Npix = nfpix_steps*nspix_steps;
    //af::flex_int spix_pos = af::flex_int(Npix,0);
    //af::flex_int fpix_pos = af::flex_int(Npix,0);
    //i_step = 0;
    //for(spixel=roi_ymin;spixel <= spixel_max;++spixel){
    //    for(fpixel=roi_xmin;fpixel <= fpixel_max;++fpixel){
    //        spix_pos[i_step] = spixel;
    //        fpix_pos[i_step] = fpixel;
    //        i_step += 1;
    //    }
    //}

    std::vector<vec3> dF_vecs;
    std::vector<vec3> dS_vecs;
    std::vector<vec3> panel_origin_dk_vecs;
    boost::shared_ptr<panel_manager> pan;
    int panel_rot_manager_ids[3] = {0,4,5};
    for (int i_rot=0; i_rot< 3; i_rot++){
        int i_manager = panel_rot_manager_ids[i_rot];
        pan = boost::dynamic_pointer_cast<panel_manager>(panels[i_manager]);
        dF_vecs.push_back(pan->dF);
        dS_vecs.push_back(pan->dS);
    }

    for (int i_orig=1; i_orig < 4; i_orig ++){
        pan = boost::dynamic_pointer_cast<panel_manager>(panels[i_orig]);
        panel_origin_dk_vecs.push_back(pan->dk);
    }

    /* Loop over 'em steps */
    //#pragma omp parallel for \
    //    shared(Npix, source_X, source_Y, source_Z, source_lambda, source_I, \
    //        spixel_max, fpixel_max, spixels, fpixels, \
    //        fdet_vector, sdet_vector, odet_vector, pix0_vector, \
    //        panels, lambda_managers, fcell_man, rot_managers, ucell_managers, Ncells_managers, \
    //        dmin, steps, phisteps, spindle_vec, mosaic_domains, mosaic_umats, \
    //        spot_scale, fluence, floatimage_roi, floatimage, \
    //        Omatrix, Umatrix, UMATS_RXYZ,Na, Nb, Nc,UMATS,  \
    //        fudge, h_min, k_min, l_min, h_range, k_range, l_range, verbose, babble, h_max, k_max, l_max, \
    //        default_F, polarization, polar_vector, \
    //        roi_xmax, roi_xmin, roi_ymax, roi_ymin, maskimage, Nsteps, \
    //        subS_pos, subF_pos, thick_pos, source_pos,  phi_pos, mos_pos, \
    //        fpix_pos, spix_pos, \
    //        subpixel_size, pixel_size, close_distance, oversample, detector_thickstep, detector_thick, detector_attnlen, \
    //        use_lambda_coefficients, compute_curvatures, ap, bp, cp, a0, b0, c0, \
    //        only_save_omega_kahn, complex_miller, Fhkl, Fhkl2, oversample_omega) \
    //    default(none)
    //for (int i2=0; i2 < Npix ; i2++ ){
    # pragma omp parallel for \
      collapse(2)
    for(int _spixel=roi_ymin;_spixel <= spixel_max;++_spixel)
    {
        for(int _fpixel=roi_xmin;_fpixel <= fpixel_max;++_fpixel)
        {
        //int _spixel = spix_pos[i2];
        //int _fpixel = fpix_pos[i2];
    //for (int i_step=0; i_step< Nsteps; i_step++){
    //    spixel = spix_pos[i_step];
    //    fpixel = fpix_pos[i_step];
            /* allow for just one part of detector to be rendered */
            int Ncontinue1=0;
            int Ncontinue2=0;
            int Npasses=0;
            int _i = _spixel*fpixels + _fpixel;
            //if(fpixel < roi_xmin || fpixel > roi_xmax || spixel < roi_ymin || spixel > roi_ymax)
            //{
            //    ++i; printf("EFFGS\n"); continue;
            //}
            //else
                //roi_i += 1;
            /* allow for the use of a mask */
            if(maskimage != NULL)
            {
                /* skip any flagged pixels in the mask */
                if(maskimage[_i] == 0)
                {
                    //++i; //++roi_i;
                    continue;
                }
            }

            /* reset photon count for this pixel */
            double _I = 0;

            /* reset derivative photon counts for the various parameters*/
            af::flex_double rot_manager_dI = af::flex_double(3,0);
            af::flex_double rot_manager_dI2 = af::flex_double(3,0);
            af::flex_double ucell_manager_dI = af::flex_double(6,0);
            af::flex_double ucell_manager_dI2 = af::flex_double(6,0);
            af::flex_double Ncells_manager_dI = af::flex_double(3,0);
            af::flex_double Ncells_manager_dI2 = af::flex_double(3,0);
            af::flex_double pan_orig_manager_dI = af::flex_double(3,0);
            af::flex_double pan_orig_manager_dI2 = af::flex_double(3,0);
            af::flex_double pan_rot_manager_dI = af::flex_double(3,0);
            af::flex_double pan_rot_manager_dI2 = af::flex_double(3,0);
            double fcell_manager_dI = 0;
            double fcell_manager_dI2 = 0;
            af::flex_double lambda_manager_dI = af::flex_double(2,0);
            af::flex_double lambda_manager_dI2 = af::flex_double(2,0);

        for (int _i_step=0; _i_step < Nsteps; _i_step++){

            int _subS = subS_pos[_i_step];
            int _subF = subF_pos[_i_step];
            int _thick_tic = thick_pos[_i_step];
            int _source = source_pos[_i_step];
            int _phi_tic = phi_pos[_i_step];
            int _mos_tic = mos_pos[_i_step];

            /* loop over sub-pixels */
            //for(subS=0;subS<oversample;++subS)
            //{
                //for(subF=0;subF<oversample;++subF)
                //{
                    /* absolute mm position on detector (relative to its origin) */
                    double _Fdet = subpixel_size*(_fpixel*oversample + _subF ) + subpixel_size/2.0;
                    double _Sdet = subpixel_size*(_spixel*oversample + _subS ) + subpixel_size/2.0;
                    //double _Fdet = subpixel_size*(fpixel*oversample + subF ) + subpixel_size/2.0;
                    //double _Sdet = subpixel_size*(spixel*oversample + subS ) + subpixel_size/2.0;

                    //for(thick_tic=0;thick_tic<detector_thicksteps;++thick_tic)
                    //{
                        /* assume "distance" is to the front of the detector sensor layer */
                       double _Odet = _thick_tic*detector_thickstep;
                        vec3 _o_vec;
                        _o_vec[0] = odet_vector[1];
                        _o_vec[1] = odet_vector[2];
                        _o_vec[2] = odet_vector[3];

                        /* construct detector subpixel position in 3D space */
                        //pixel_pos[1] = _Fdet*fdet_vector[1]+_Sdet*sdet_vector[1]+_Odet*odet_vector[1]+pix0_vector[1];
                        //pixel_pos[2] = _Fdet*fdet_vector[2]+_Sdet*sdet_vector[2]+_Odet*odet_vector[2]+pix0_vector[2];
                        //pixel_pos[3] = _Fdet*fdet_vector[3]+_Sdet*sdet_vector[3]+_Odet*odet_vector[3]+pix0_vector[3];
                        //pixel_pos[0] = 0.0;
                        vec3 _pixel_pos;
                        _pixel_pos[0] = _Fdet*fdet_vector[1]+_Sdet*sdet_vector[1]+_Odet*odet_vector[1]+pix0_vector[1];
                        _pixel_pos[1] = _Fdet*fdet_vector[2]+_Sdet*sdet_vector[2]+_Odet*odet_vector[2]+pix0_vector[2];
                        _pixel_pos[2] = _Fdet*fdet_vector[3]+_Sdet*sdet_vector[3]+_Odet*odet_vector[3]+pix0_vector[3];

                        int panel_rot_manager_ids[3] = {0,4,5};
                        std::vector<vec3> panel_rotation_dk_vecs;
                        for (int i_k=0;i_k < 3; i_k++){
                            int i_manager = panel_rot_manager_ids[i_k];
                            //pan_rot = boost::dynamic_pointer_cast<panel_manager>(panels[i_manager]);
                            vec3 dk = vec3(0,0,0);
                            if (panels[i_manager]->refine_me){
                                //pan_rot->dk = _Fdet*(pan_rot->dF) + _Sdet*(pan_rot->dS);
                                dk = _Fdet*(dF_vecs[i_k]) + _Sdet*(dS_vecs[i_k]);
                            }
                            panel_rotation_dk_vecs.push_back(dk);
                        }

                        double _airpath = _pixel_pos.length();
                        vec3 _diffracted = _pixel_pos/_airpath;

                        /* solid angle subtended by a pixel: (pix/airpath)^2*cos(2theta) */
                        double _omega_pixel = pixel_size*pixel_size/_airpath/_airpath*close_distance/_airpath;

                        /* option to turn off obliquity effect, inverse-square-law only */
                        if(point_pixel) _omega_pixel = 1.0/_airpath/_airpath;
                        //omega_sum += _omega_pixel;

                        /* now calculate detector thickness effects */
                        double _capture_fraction = 1;
                        if(detector_thick > 0.0 && detector_attnlen > 0.0)
                        {
                            /* inverse of effective thickness increase */
                            double _parallax = _diffracted * _o_vec ; //dot_product(diffracted,odet_vector);
                            double _capture_fraction = exp(-_thick_tic*detector_thickstep/detector_attnlen/_parallax)
                                              -exp(-(_thick_tic+1)*detector_thickstep/detector_attnlen/_parallax);
                        }

                    /* loop over sources now */
                    //for(source=0;source<sources;++source){
                        /* retrieve stuff from cache */
                        vec3 _incident;
                        _incident[0] = -source_X[_source];
                        _incident[1] = -source_Y[_source];
                        _incident[2] = -source_Z[_source];
                        double _lambda = source_lambda[_source];

                        //std::vector<vec3> = lambda_dg_vectors;
                        double lambda_ang = _lambda*1e10;
                        if (use_lambda_coefficients){
                            lambda_ang = (lambda_managers[0]->value) + (lambda_managers[1]->value)*lambda_ang;
                            _lambda = lambda_ang*1e-10;
                        }
                        //if (lambda_managers[0]->refine_me || lambda_managers[1]->refine_me){
                        //    lambda_managers[0]->dg_dlambda = 1;
                        //    lambda_managers[1]->dg_dlambda = lambda_ang;
                        //}

                        /* construct the incident beam unit vector while recovering source distance */
                        double _source_path = _incident.length();
                        _incident = _incident / _source_path;
                        //source_path = unitize(incident,incident);

                        /* construct the scattering vector for this pixel */
                        vec3 _scattering = (_diffracted - _incident) / _lambda;

                        /* sin(theta)/lambda is half the scattering vector length */
                        double _stol = 0.5*(_scattering.length()); //magnitude(scattering);

                        /* rough cut to speed things up when we aren't using whole detector */
                        if(dmin > 0.0 && _stol > 0.0)
                        {
                            if(dmin > 0.5/_stol)
                            {
                                Ncontinue1 += 1;
                                continue;
                            }
                        }

                        /* sweep over phi angles */
                        //for(phi_tic = 0; phi_tic < phisteps; ++phi_tic)
                        // {
                            vec3 ap_vec = vec3(ap[1], ap[2], ap[3]);
                            vec3 bp_vec = vec3(bp[1], bp[2], bp[3]);
                            vec3 cp_vec = vec3(cp[1], cp[2], cp[3]);

                            double _phi = phi0 + phistep*_phi_tic;
                            if( _phi != 0.0 )
                            {
                                vec3 a0_vec = vec3(a0[1], a0[2], a0[3]);
                                vec3 b0_vec = vec3(b0[1], b0[2], b0[3]);
                                vec3 c0_vec = vec3(c0[1], c0[2], c0[3]);
                                double cosphi = cos(_phi);
                                double sinphi = sin(_phi);
                                ap_vec = a0_vec*cosphi + spindle_vec.cross(a0_vec)*sinphi + spindle_vec*(spindle_vec*a0_vec)*(1-cosphi);
                                bp_vec = b0_vec*cosphi + spindle_vec.cross(b0_vec)*sinphi + spindle_vec*(spindle_vec*b0_vec)*(1-cosphi);
                                cp_vec = c0_vec*cosphi + spindle_vec.cross(c0_vec)*sinphi + spindle_vec*(spindle_vec*c0_vec)*(1-cosphi);
                                /* rotate about spindle if neccesary */
                                //rotate_axis(a0,ap,spindle_vector,_phi);
                                //rotate_axis(b0,bp,spindle_vector,_phi);
                                //rotate_axis(c0,cp,spindle_vector,_phi);
                                // TODO Legacy preserve
                                //for (int ii=0; ii <3; ii++){
                                //    ap[ii+1] = ap_vec[ii];
                                //    bp[ii+1] = bp_vec[ii];
                                //    cp[ii+1] = cp_vec[ii];
                                //}
                            }
                            /* enumerate mosaic domains */
                            //for(mos_tic=0;mos_tic<mosaic_domains;++mos_tic)
                           //{
                                vec3 q_vec = 1e-10*_scattering;

                                mat3 Bmat_realspace;
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
                                mat3 UBO = (UMATS_RXYZ[_mos_tic] * Umatrix*Bmat_realspace*(Omatrix.transpose())).transpose();
                                vec3 H_vec = UBO * q_vec;
                                double _h = H_vec[0];
                                double _k = H_vec[1];
                                double _l = H_vec[2];

                                /* round off to nearest whole index */
                                int _h0 = static_cast<int>(ceil(_h - 0.5));
                                int _k0 = static_cast<int>(ceil(_k - 0.5));
                                int _l0 = static_cast<int>(ceil(_l - 0.5));

                                vec3 H0_vec;
                                H0_vec[0] = _h0;
                                H0_vec[1] = _k0;
                                H0_vec[2] = _l0;

                                mat3 _NABC = mat3(0,0,0,0,0,0,0,0,0);
                                _NABC[0] = Na;
                                _NABC[4] = Nb;
                                _NABC[8] = Nc;

                                double C = 2 / 0.63 * fudge;
                                vec3 delta_H = (H_vec - H0_vec);
                                vec3 V = _NABC*delta_H;// (H_vec- H0_vec);

                                double _hrad_sqr = V*V;
                                double _F_latt = Na*Nb*Nc*exp(-( _hrad_sqr / 0.63 * fudge ));

                                /* no need to go further if result will be zero */
                                if(_F_latt == 0.0 && ! only_save_omega_kahn) {
                                    //Ncontinue2 +=1;
                                    continue;
                                }

                                //Npasses += 1;
                                /* structure factor of the unit cell */
                                double _F_cell = default_F;
                                double _F_cell2 = 0;

                                if ( (_h0<=h_max) && (_h0>=h_min) && (_k0<=k_max) && (_k0>=k_min) && (_l0<=l_max) && (_l0>=l_min)  ) {
                                    /* just take nearest-neighbor */
                                    _F_cell = Fhkl[_h0-h_min][_k0-k_min][_l0-l_min];
                                    if (complex_miller) _F_cell2 = Fhkl2[_h0-h_min][_k0-k_min][_l0-l_min];
                                }

                                if (complex_miller)
                                  _F_cell = sqrt(_F_cell*_F_cell + _F_cell2*_F_cell2);

                                /* now we have the structure factor for this pixel */

                                /* polarization factor */
                                //polar = 1;
                                //if(! nopolar){
                                //    /* need to compute polarization factor */
                                //    /* Note, if not oversample_omega we will compute polarization factor once for center of pixel, after looping over steps*/
                                //    if (oversample_omega)
                                //      polar = polarization_factor(polarization,incident,diffracted,polar_vector);
                                //}

                                /* convert amplitudes into intensity (photons per steradian) */
                                if (!oversample_omega)
                                    _omega_pixel = 1;

                                /* increment to intensity */
                                double Iincrement = _F_cell*_F_cell*_F_latt*_F_latt*source_I[_source]*_capture_fraction*_omega_pixel;
                                _I += Iincrement;  //  F_cell*F_cell*F_latt*F_latt*source_I[source]*capture_fraction*omega_pixel;

                                if(verbose > 3)
                                    printf("hkl= %f %f %f  hkl1= %d %d %d  Fcell=%f\n", _h,_k,_l,_h0,_k0,_l0, _F_cell);

                                double two_C = 2*C;
                                mat3 UBOt = Umatrix*Bmat_realspace*(Omatrix.transpose());
                                if (rot_managers[0]->refine_me){
                                    mat3 RyRzUBOt = RotMats[1]*RotMats[2]*UBOt;
                                    vec3 delta_H_prime = (UMATS[_mos_tic]*dRotMats[0]*RyRzUBOt).transpose()*q_vec;
                                    double V_dot_dV = V*(_NABC*delta_H_prime);
                                    double value = -two_C * V_dot_dV * Iincrement;
                                    double value2 =0;
                                    if (compute_curvatures) {
                                        vec3 delta_H_dbl_prime = (UMATS[_mos_tic]*d2RotMats[0]*RyRzUBOt).transpose()*q_vec;
                                        double dV_dot_dV = (_NABC*delta_H_prime)*(_NABC*delta_H_prime);
                                        double dV2_dot_V = (_NABC*delta_H)*(_NABC*delta_H_dbl_prime);
                                        value2 = two_C*(two_C*V_dot_dV*V_dot_dV - dV2_dot_V - dV_dot_dV)*Iincrement;
                                    }
                                    rot_manager_dI[0] += value;
                                    rot_manager_dI2[0] += value2;
                                    //rot_managers[0]->increment(value, value2);
                                }
                                if (rot_managers[1]->refine_me){
                                    mat3 UmosRx = UMATS[_mos_tic]*RotMats[0];
                                    mat3 RzUBOt = RotMats[2]*UBOt;
                                    vec3 delta_H_prime =(UmosRx*dRotMats[1]*RzUBOt).transpose()*q_vec;
                                    double V_dot_dV = V*(_NABC*delta_H_prime);
                                    double value = -two_C * V_dot_dV * Iincrement;

                                    double value2=0;
                                    if (compute_curvatures){
                                        vec3 delta_H_dbl_prime = (UmosRx*d2RotMats[1]*RzUBOt).transpose()*q_vec;
                                        double dV_dot_dV = (_NABC*delta_H_prime)*(_NABC*delta_H_prime);
                                        double dV2_dot_V = (_NABC*delta_H)*(_NABC*delta_H_dbl_prime);
                                        value2 = two_C*(two_C*V_dot_dV*V_dot_dV - dV2_dot_V - dV_dot_dV)*Iincrement;
                                    }
                                    rot_manager_dI[1] += value;
                                    rot_manager_dI2[1] += value2;
                                    //rot_managers[1]->increment(value, value2);
                                }
                                if (rot_managers[2]->refine_me){
                                    mat3 UmosRxRy = UMATS[_mos_tic]*RotMats[0]*RotMats[1];
                                    vec3 delta_H_prime = (UmosRxRy*dRotMats[2]*UBOt).transpose()*q_vec;
                                    double V_dot_dV = V*(_NABC*delta_H_prime);
                                    double value = -two_C * V_dot_dV * Iincrement;

                                    double value2=0;
                                    if (compute_curvatures){
                                        vec3 delta_H_dbl_prime = (UmosRxRy*d2RotMats[2]*UBOt).transpose()*q_vec;
                                        double dV_dot_dV = (_NABC*delta_H_prime)*(_NABC*delta_H_prime);
                                        double dV2_dot_V = (_NABC*delta_H)*(_NABC*delta_H_dbl_prime);
                                        value2 = two_C*(two_C*V_dot_dV*V_dot_dV - dV2_dot_V - dV_dot_dV)*Iincrement;
                                    }
                                    rot_manager_dI[2] += value;
                                    rot_manager_dI2[2] += value2;
                                    //rot_managers[2]->increment(value, value2);
                                }

                                /*Checkpoint for unit cell derivatives*/
                                mat3 Ot = Omatrix.transpose();
                                for(int i_uc=0; i_uc < 6; i_uc++ ){
                                    if (ucell_managers[i_uc]->refine_me){
                                        mat3 UmosRxRyRzU = UMATS_RXYZ[_mos_tic]*Umatrix;
                                        vec3 delta_H_prime = ((UmosRxRyRzU*(ucell_managers[i_uc]->dB)*Ot).transpose()*q_vec);
                                        double V_dot_dV = V*(_NABC*delta_H_prime);
                                        double value = -two_C * V_dot_dV * Iincrement;
                                        double value2 =0;
                                        if (compute_curvatures){
                                            vec3 delta_H_dbl_prime = ((UmosRxRyRzU*(ucell_managers[i_uc]->dB2)*Ot).transpose()*q_vec);
                                            double dV_dot_dV = (_NABC*delta_H_prime)*(_NABC*delta_H_prime);
                                            double dV2_dot_V = (_NABC*delta_H)*(_NABC*delta_H_dbl_prime);
                                            value2 = two_C*(two_C*V_dot_dV*V_dot_dV - dV2_dot_V - dV_dot_dV)*Iincrement;
                                        }
                                        ucell_manager_dI[i_uc] += value;
                                        ucell_manager_dI2[i_uc] += value2;

                                        //ucell_managers[i_uc]->increment(value, value2);
                                    }
                                } /*end ucell deriv */

                                /* Checkpoint for Ncells manager */
                                if (Ncells_managers[0]->refine_me){
                                    int num_ncell_deriv = 1;
                                    if (not isotropic_ncells)
                                        num_ncell_deriv = 3;
                                    for (int i_nc=0; i_nc < num_ncell_deriv; i_nc++) {
                                        int i_diag = i_nc*3 + i_nc; // diagonal term
                                        mat3 dN = mat3(0,0,0,0,0,0,0,0,0);
                                        dN[i_diag] = 1;
                                        double N_i = _NABC[i_diag];
                                        vec3 dV_dN = dN*delta_H;
                                        double deriv_coef = (1/N_i - C*dV_dN*V);
                                        double value = 2*Iincrement*deriv_coef;
                                        double value2=0;
                                        if(compute_curvatures){
                                            dN[i_diag] = 0; // TODO check maths
                                            value2 = ( -1/N_i/N_i - C*dV_dN*dV_dN)*2*Iincrement;
                                            value2 += deriv_coef*2*value;
                                        }
                                        Ncells_manager_dI[i_nc] += value;
                                        Ncells_manager_dI2[i_nc] += value2;
                                        //Ncells_managers[i_nc]->increment(value, value2);
                                    }

                                } /* end Ncells manager deriv */

                                ///* Checkpoint for Origin manager */
                                for (int i_pan_orig=0; i_pan_orig < 3; i_pan_orig++){
                                    if (panels[i_pan_orig+1]->refine_me){
                                        double per_k = 1/_airpath;
                                        double per_k3 = pow(per_k,3);
                                        double per_k5 = pow(per_k,5);
                                        double lambda_ang = _lambda*1e10;

                                        mat3 M = -two_C*(_NABC*UBO)/lambda_ang;
                                        vec3 dk = panel_origin_dk_vecs[i_pan_orig];
                                        af::flex_double dI_and_dI2 = get_panel_increment(Iincrement, _omega_pixel, M, subpixel_size*subpixel_size,
                                            _o_vec, _pixel_pos, per_k,  per_k3, per_k5, V, dk);
                                        pan_orig_manager_dI[i_pan_orig] += dI_and_dI2[0];
                                        pan_orig_manager_dI2[i_pan_orig] += dI_and_dI2[1];

                                    } /* end origin manager deriv */
                                }

                                //boost::shared_ptr<panel_manager> pan_rot;
                                int panel_manager_ids[3] = {0,4,5};
                                for (int i_pan_rot=0; i_pan_rot < 3; i_pan_rot++){
                                    int pan_manager_id = panel_manager_ids[i_pan_rot];
                                    //pan_rot = boost::dynamic_pointer_cast<panel_manager>(panels[pan_manager_id]);
                                    //if(pan_rot->refine_me){
                                    if(panels[pan_manager_id]->refine_me){
                                        double per_k = 1/_airpath;
                                        double per_k3 = pow(per_k,3);
                                        double per_k5 = pow(per_k,5);
                                        double lambda_ang = _lambda*1e10;
                                        mat3 M = -two_C*(_NABC*UBO)/lambda_ang;
                                        af::flex_double dI_and_dI2 = get_panel_increment(Iincrement, _omega_pixel, M, subpixel_size*subpixel_size,
                                            _o_vec, _pixel_pos, per_k,  per_k3, per_k5, V, panel_rotation_dk_vecs[i_pan_rot]);
                                        pan_rot_manager_dI[i_pan_rot] += dI_and_dI2[0];
                                        pan_rot_manager_dI2[i_pan_rot] += dI_and_dI2[1];
                                    }
                                }

                                /* checkpoint for Fcell manager */
                                if (fcell_man->refine_me){
                                    double value = 2*Iincrement/_F_cell ;
                                    double value2=0;
                                    if (compute_curvatures){
                                        value2 = value/_F_cell;
                                    }
                                    fcell_manager_dI += value;
                                    fcell_manager_dI2 += value2;
                                } /* end of fcell man deriv */

                                /*checkpoint for lambda manager*/
                                for(int i_lam=0; i_lam < 2; i_lam++){
                                    if (lambda_managers[i_lam]->refine_me){
                                        double lambda_ang = _lambda*1e10;
                                        double NH_dot_V = (_NABC*H_vec)*V;
                                        double dg_dlambda;
                                        if (i_lam==0)
                                            dg_dlambda = 1;
                                        else // i_lam==1
                                            dg_dlambda = lambda_ang;
                                        double coef = NH_dot_V*two_C*(dg_dlambda) / lambda_ang;
                                        double value = coef*Iincrement;
                                        double value2 = 0;
                                        //if (compute_curvatures)
                                        lambda_manager_dI[i_lam] += value;
                                        lambda_manager_dI2[i_lam] += value2;
                                        //lambda_managers[i_lam]->increment(value, value2);
                                    }
                                }
                                /*end of lambda deriv*/
                            //}
                            /* end of mosaic loop */
                        //}
                        /* end of phi loop */
                    //}
                    /* end of source loop */
                //}
                    /* end of detector thickness loop */
                //}
                /* end of sub-pixel y loop */
            //}
            /* end of sub-pixel x loop */
			} /* end of i_steps loop */
            /* absolute mm position on detector (relative to its origin) */
            if (verbose)
                printf("Pixel %d, %d: Ncont1=%d, Ncont2=%d, Npasses=%d\n" ,_fpixel, _spixel, Ncontinue1, Ncontinue2, Npasses);
            double _Fdet_ave = pixel_size*_fpixel + pixel_size/2.0;
            double _Sdet_ave = pixel_size*_spixel + pixel_size/2.0;
            double _Odet_ave = 0; //Odet; // TODO maybe make this more general for thick detectors?

            double _pixel_pos_ave[4] = {0,0,0,0};
            _pixel_pos_ave[1] = _Fdet_ave*fdet_vector[1]+_Sdet_ave*sdet_vector[1]+_Odet_ave*odet_vector[1]+pix0_vector[1];
            _pixel_pos_ave[2] = _Fdet_ave*fdet_vector[2]+_Sdet_ave*sdet_vector[2]+_Odet_ave*odet_vector[2]+pix0_vector[2];
            _pixel_pos_ave[3] = _Fdet_ave*fdet_vector[3]+_Sdet_ave*sdet_vector[3]+_Odet_ave*odet_vector[3]+pix0_vector[3];

            double _diffracted_ave[4] = {0,0,0,0};
            double _airpath_ave = unitize(_pixel_pos_ave,_diffracted_ave);
            double _omega_pixel_ave = pixel_size*pixel_size/_airpath_ave/_airpath_ave*close_distance/_airpath_ave;

            vec3 _k_diffracted_ave = vec3(pixel_pos_ave[1], pixel_pos_ave[2], pixel_pos_ave[3]);
            //k_diffracted_ave[0] = pixel_pos_ave[1];
            //k_diffracted_ave[1] = pixel_pos_ave[2];
            //k_diffracted_ave[2] = pixel_pos_ave[3];

            /* polarization factor */
            //if(!nopolar){
            //    /* TODO: what about divergence and different incident vectors ? perhaps do for averagce incident vec */
            //    if (!oversample_omega)
            //      polar = polarization_factor(polarization,incident,diffracted_ave,polar_vector);
            //}
            double _polar = 1;
            if (!nopolar){
                double incident_beam[4] = {0,0,0,0};
                // TODO what about divergence here ?
                incident_beam[1] = -source_X[0];
                incident_beam[2] = -source_Y[0];
                incident_beam[3] = -source_Z[0];
                unitize(incident_beam, incident_beam); // TODO remove ?

                _polar = polarization_factor(polarization, incident_beam, _diffracted_ave, polar_vector);

                //double cos2theta,cos2theta_sqr,sin2theta_sqr;
                //vec3 E_in[4];
                //vec3 B_in[4];

                //unitize(__incident,__incident);
                //unitize(__diffracted,__diffracted);
                //unitize(__axis,__axis);

                ///* component of diffracted unit vector along incident beam unit vector */
                //cos2theta = dot_product(__incident,__diffracted);
                //cos2theta_sqr = cos2theta*cos2theta;
                //sin2theta_sqr = 1-cos2theta_sqr;
                ////double _u = cos2theta_sqr;

                //double _psi=0;
                //if(kahn_factor != 0.0){
                //    //SCITBX_EXAMINE(kahn_factor);
                //    /* tricky bit here is deciding which direciton the E-vector lies in for each source
                //       here we assume it is closest to the "axis" defined above */

                //    /* cross product to get "vertical" axis that is orthogonal to the cannonical "polarization" */
                //    cross_product(__axis,__incident,B_in);
                //    /* make it a unit vector */
                //    unitize(B_in,B_in);
                //    //vec3 Bi_vec = vec3(B_in[1], B_in[2], B_in[3]);
                //    //Bi_vec[0] = B_in[1];
                //    //Bi_vec[1] = B_in[2];
                //    //Bi_vec[2] = B_in[3];

                //    /* cross product with incident beam to get E-vector direction */
                //    cross_product(__incident,B_in,E_in);
                //    /* make it a unit vector */
                //    unitize(E_in,E_in);
                //    //vec3 Ei_vec = vec3(E_in[1], E_in[2], E_in[3]);
                //    //Ei_vec[0] = E_in[1];
                //    //Ei_vec[1] = E_in[2];
                //    //Ei_vec[2] = E_in[3];

                //    /* get components of diffracted ray projected onto the E-B plane */
                //    double _kEi = dot_product(__diffracted,E_in);
                //    double _kBi = dot_product(__diffracted,B_in);

                //    /* compute the angle of the diffracted ray projected onto the incident E-B plane */
                //    _psi = -atan2(_kBi,_kEi);
                //}
                ///* correction for polarized incident beam */
                //return 0.5*(1.0 + cos2theta_sqr - kahn_factor*cos(2*_psi)*sin2theta_sqr);
            }

            double _om = 1;
            if (!oversample_omega)
                _om=_omega_pixel_ave;
            // final scale term to being everything to photon number units
            double _scale_term = r_e_sqr*fluence*spot_scale*_polar*_om / steps;
            //double omp_time = omp_get_wtime();
            if(only_save_omega_kahn)
                floatimage[_i] = _polar*_omega_pixel_ave;
            else
                floatimage[_i] = _scale_term*_I;

            int roi_fdim = roi_xmax - roi_xmin+1;
            int roi_fs = _i % fpixels - roi_xmin;
            int roi_ss = floor(_i / fpixels) - roi_ymin ;
            int roi_i =  roi_ss*roi_fdim + roi_fs ;
            //printf("THREAD # %d: roi_i=%d\n", omp_get_thread_num() , roi_i);

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
            floatimage_roi[roi_i] = _scale_term*_I;

            //omp_time = omp_get_wtime()- omp_time;
            //omp_time = omp_time*1000000;
            //if (verbose)
            //    printf("THREAD %d, wtime to update floatimage_roi %e microsec\n", omp_get_thread_num(), omp_time);

            /* udpate the rotation derivative images*/
            for (int i_rot =0 ; i_rot < 3 ; i_rot++){
                if (rot_managers[i_rot]->refine_me){
                    //double value = _scale_term*rot_managers[i_rot]->dI;
                    //double value2 = scale_term*rot_managers[i_rot]->dI2;
                    double value = _scale_term*rot_manager_dI[i_rot];
                    double value2 = _scale_term*rot_manager_dI2[i_rot];
                    rot_managers[i_rot]->increment_image(roi_i, value, value2, compute_curvatures);
                }
            } /* end rot deriv image increment */

            /*update the ucell derivative images*/
            for (int i_uc=0 ; i_uc < 6 ; i_uc++){
                if (ucell_managers[i_uc]->refine_me){
                    //double value = scale_term*ucell_managers[i_uc]->dI;
                    //double value2 = scale_term*ucell_managers[i_uc]->dI2;
                    double value = _scale_term*ucell_manager_dI[i_uc];
                    double value2 = _scale_term*ucell_manager_dI2[i_uc];
                    ucell_managers[i_uc]->increment_image(roi_i, value, value2, compute_curvatures);
                }
            }/* end ucell deriv image increment */

            /*update the Ncells derivative image*/
            if (Ncells_managers[0]->refine_me){
                //double value = scale_term*Ncells_managers[0]->dI;
                //double value2 = scale_term*Ncells_managers[0]->dI2;
                double value = _scale_term*Ncells_manager_dI[0];
                double value2 = _scale_term*Ncells_manager_dI2[0];
                Ncells_managers[0]->increment_image(roi_i, value, value2, compute_curvatures);
                if (! isotropic_ncells){
                    value = _scale_term*Ncells_manager_dI[1];
                    value2 = _scale_term*Ncells_manager_dI2[1];
                    Ncells_managers[1]->increment_image(roi_i, value, value2, compute_curvatures);

                    value = _scale_term*Ncells_manager_dI[2];
                    value2 = _scale_term*Ncells_manager_dI2[2];
                    Ncells_managers[2]->increment_image(roi_i, value, value2, compute_curvatures);
                }
            }/* end Ncells deriv image increment */

            /* update Fcell derivative image */
            if(fcell_man->refine_me){
                //double value = scale_term*fcell_man->dI;
                //double value2 = scale_term*fcell_man->dI2;
                double value = _scale_term*fcell_manager_dI;
                double value2 = _scale_term*fcell_manager_dI2;
                fcell_man->increment_image(roi_i, value, value2, compute_curvatures);
            }/* end Fcell deriv image increment */

            /*update the lambda derivative images*/
            for (int i_lam=0 ; i_lam < 2 ; i_lam++){
                if (lambda_managers[i_lam]->refine_me){
                    //double value = scale_term*lambda_managers[i_lam]->dI;
                    //double value2 = scale_term*lambda_managers[i_lam]->dI2;
                    double value = _scale_term*lambda_manager_dI[i_lam];
                    double value2 = _scale_term*lambda_manager_dI2[i_lam];
                    lambda_managers[i_lam]->increment_image(roi_i, value, value2, compute_curvatures);
                }
            }/* end lambda deriv image increment */

            int pan_manager_ids[3] = {0,4,5};
            for (int i_pan_rot=0; i_pan_rot < 3; i_pan_rot++){
                int pan_manager_id = pan_manager_ids[i_pan_rot];
                if(panels[pan_manager_id]->refine_me){
                    //double value = scale_term*pan_rot->dI;
                    //double value2 = scale_term*pan_rot->dI2;
                    double value = _scale_term*pan_rot_manager_dI[i_pan_rot];
                    double value2 = _scale_term*pan_rot_manager_dI2[i_pan_rot];
                    panels[pan_manager_id]->increment_image(roi_i, value, value2, compute_curvatures);
                }
            }/* end panel rot deriv image increment */

            for (int i_pan_orig=0; i_pan_orig < 3; i_pan_orig++){
                int pan_manager_id = i_pan_orig+1;
                //pan_orig = boost::dynamic_pointer_cast<panel_manager>(panels[1 + i_pan_orig]);
                if(panels[pan_manager_id]->refine_me){
                    //double value = scale_term*pan_orig->dI;
                    //double value2 = scale_term*pan_orig->dI2;
                    double value = _scale_term*pan_orig_manager_dI[i_pan_orig];
                    double value2 = _scale_term*pan_orig_manager_dI2[i_pan_orig];
                    panels[pan_manager_id]->increment_image(roi_i, value, value2, compute_curvatures);
                }/* end panel orig deriv image increment */
            }
            /* update origin derivative image */
            //if(floatimage[_i] > max_I) {
            //    max_I = floatimage[_i];
            //    max_I_x = Fdet;
            //    max_I_y = Sdet;
            //    max_I_hkl[0] = h0;
            //    max_I_hkl[1] = k0;
            //    max_I_hkl[2] = l0;
            //}
            //sum += floatimage[_i];
            //sumsqr += floatimage[_i]*floatimage[_i];
            //++sumn;

            //if( printout )
            //{
            //    if((_fpixel==printout_fpixel && _spixel==printout_spixel) || printout_fpixel < 0)
            //    {
            //        twotheta = atan2(sqrt(pixel_pos[2]*pixel_pos[2]+pixel_pos[3]*pixel_pos[3]),pixel_pos[1]);
            //        test = sin(twotheta/2.0)/(lambda0*1e10);
            //        printf("%4d %4d : stol = %g or %g\n", _fpixel,_spixel,stol,test);
            //        printf("at %g %g %g\n", pixel_pos[1],pixel_pos[2],pixel_pos[3]);
            //        printf("hkl= %f %f %f  hkl0= %d %d %d\n", h,k,l,h0,k0,l0);
            //        printf(" F_cell=%g  F_latt=%g   I = %g\n", F_cell,F_latt,I);
            //        printf("I/steps %15.10g\n", I/steps);
            //        printf("polar   %15.10g\n", polar);
            //        printf("omega   %15.10g\n", omega_pixel);
            //        /* some useful printouts for debugging purposes */
            //        //SCITBX_EXAMINE(I);
            //        //SCITBX_EXAMINE(FF);
            //        //SCITBX_EXAMINE(scale_term);
            //        //SCITBX_EXAMINE(omega_pixel_ave);
            //        //SCITBX_EXAMINE(om);
            //        //SCITBX_EXAMINE(diffracted_ave[0]);
            //        //SCITBX_EXAMINE(diffracted_ave[1]);
            //        //SCITBX_EXAMINE(diffracted_ave[2]);
            //        //SCITBX_EXAMINE(diffracted_ave[3]);
            //        //SCITBX_EXAMINE(diffracted[0]);
            //        //SCITBX_EXAMINE(diffracted[1]);
            //        //SCITBX_EXAMINE(diffracted[2]);
            //        //SCITBX_EXAMINE(diffracted[3]);
            //        //SCITBX_EXAMINE(airpath);
            //        //SCITBX_EXAMINE(airpath_ave);
            //        //SCITBX_EXAMINE(pixel_pos[0]);
            //        //SCITBX_EXAMINE(pixel_pos[1]);
            //        //SCITBX_EXAMINE(pixel_pos[2]);
            //        //SCITBX_EXAMINE(pixel_pos[3]);
            //        //SCITBX_EXAMINE(pixel_pos_ave[0]);
            //        //SCITBX_EXAMINE(pixel_pos_ave[1]);
            //        //SCITBX_EXAMINE(pixel_pos_ave[2]);
            //        //SCITBX_EXAMINE(pixel_pos_ave[3]);
            //        //SCITBX_EXAMINE(floatimage[_i]);
            //        //SCITBX_EXAMINE(u);
            //        //SCITBX_EXAMINE(du);
            //        //SCITBX_EXAMINE(v);
            //        //SCITBX_EXAMINE(dv);
            //        //SCITBX_EXAMINE(kBi);
            //        //SCITBX_EXAMINE(kEi);
            //        //SCITBX_EXAMINE(w);
            //        //SCITBX_EXAMINE(gam_sin2psi);
            //        //SCITBX_EXAMINE(gam_cos2psi);
            //        //SCITBX_EXAMINE(dpolar);
            //        //SCITBX_EXAMINE(dpolar2);
            //        //SCITBX_EXAMINE(dOmega);
            //        //SCITBX_EXAMINE(dOmega2);
            //        //SCITBX_EXAMINE(FF);
            //        //SCITBX_EXAMINE(FdF);
            //        //SCITBX_EXAMINE(dFdF);
            //        //SCITBX_EXAMINE(FdF2);
            //        //SCITBX_EXAMINE(k_diffracted_ave[0]);
            //        //SCITBX_EXAMINE(k_diffracted_ave[1]);
            //        //SCITBX_EXAMINE(k_diffracted_ave[2]);
            //        //SCITBX_EXAMINE(nopolar);
            //        //SCITBX_EXAMINE(oversample_omega);
            //        //SCITBX_EXAMINE(isotropic_ncells);
            //        //SCITBX_EXAMINE(Na);
            //        //SCITBX_EXAMINE(Nb);
            //        //SCITBX_EXAMINE(Nc);
            //        //SCITBX_EXAMINE(dN[0]);
            //        //SCITBX_EXAMINE(dN[1]);
            //        //SCITBX_EXAMINE(dN[2]);
            //        //SCITBX_EXAMINE(dN[3]);
            //        //SCITBX_EXAMINE(dN[4]);
            //        //SCITBX_EXAMINE(dN[5]);
            //        //SCITBX_EXAMINE(dN[6]);
            //        //SCITBX_EXAMINE(dN[7]);
            //        //SCITBX_EXAMINE(dN[8]);
            //        //double lambda_coef0 = lambda_managers[0]->value;
            //        //double lambda_coef1 = lambda_managers[1]->value;
            //        //SCITBX_EXAMINE(lambda_managers[0]->value);
            //        //SCITBX_EXAMINE(lambda_managers[1]->value);
            //        //SCITBX_EXAMINE(use_lambda_coefficients);
            //        printf("real-space cell vectors (Angstrom):\n");
            //        printf("     %-10s  %-10s  %-10s\n","a","b","c");
            //        printf("X: %11.8f %11.8f %11.8f\n",a[1]*1e10,b[1]*1e10,c[1]*1e10);
            //        printf("Y: %11.8f %11.8f %11.8f\n",a[2]*1e10,b[2]*1e10,c[2]*1e10);
            //        printf("Z: %11.8f %11.8f %11.8f\n",a[3]*1e10,b[3]*1e10,c[3]*1e10);
            //        printf("Rot manager refine status X=%d, Y=%d, Z=%d\n",
            //            rot_managers[0]->refine_me, rot_managers[1]->refine_me,
            //            rot_managers[2]->refine_me);
            //        printf("Ucell managers refine status a.a=%d, b.b=%d, c.c=%d, a.b=%d, a.c=%d, b.c=%d\n",
            //            ucell_managers[0]->refine_me, ucell_managers[1]->refine_me, ucell_managers[2]->refine_me,
            //            ucell_managers[3]->refine_me, ucell_managers[4]->refine_me, ucell_managers[5]->refine_me);
            //        printf("Ncell managers refine status: %d, value=%f\n", Ncells_managers[0]->refine_me,
            //                Ncells_managers[0]->value);
            //        printf("Fcell manager refine status: %d\n", fcell_man->refine_me);
            //        //boost::shared_ptr<origin_manager> pan_orig = boost::dynamic_pointer_cast<origin_manager>(panels[1]);
            //        //printf("Origin managers refine status: %d, value=%f\n", pan_orig->refine_me,
            //        //        pan_orig->value);
            //        //printf("Bmatrix_real:\n%11.8f %11.8f %11.8f\n %11.8f %11.8f %11.8f\n %11.8f %11.8f %11.8f\n",
            //        //    Bmat_realspace[0], Bmat_realspace[1], Bmat_realspace[2],
            //        //    Bmat_realspace[3], Bmat_realspace[4], Bmat_realspace[5],
            //        //    Bmat_realspace[6], Bmat_realspace[7], Bmat_realspace[8]);
            //        printf("Umatrix_real:\n%11.8f %11.8f %11.8f\n %11.8f %11.8f %11.8f\n %11.8f %11.8f %11.8f\n",
            //            Umatrix[0], Umatrix[1], Umatrix[2],
            //            Umatrix[3], Umatrix[4], Umatrix[5],
            //            Umatrix[6], Umatrix[7], Umatrix[8]);
            //    }
            //}
            //else
            //{
            //    if(progress_meter && verbose && progress_pixels/100 > 0)
            //    {
            //        if(progress_pixel % ( progress_pixels/20 ) == 0 ||
            //           ((10*progress_pixel<progress_pixels ||
            //             10*progress_pixel>9*progress_pixels) &&
            //            (progress_pixel % (progress_pixels/100) == 0)))
            //        {
            //            printf("%lu%% done\n",progress_pixel*100/progress_pixels);
            //        }
            //    }
            //    ++progress_pixel;
            //}
            //++i;
        } /* end of fpixels loop */
    } /* end of spixels loop */
    //} /* end of i2 steps */
    if(verbose) printf("done with pixel loop\n");
    //if(verbose) printf("solid angle subtended by detector = %g steradian ( %g%% sphere)\n",omega_sum/steps,100*omega_sum/steps/4/M_PI);
    //if(verbose) printf("max_I= %g sum= %g avg= %g\n",max_I,sum,sum/sumn);

} // END  of add_diffBragg_spots
// END diffBragg

} // end of namespace nanoBragg
} // end of namespace simtbx
