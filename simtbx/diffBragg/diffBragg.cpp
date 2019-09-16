#include <simtbx/diffBragg/diffBragg.h>
#include <assert.h>
namespace simtbx {
namespace nanoBragg {

// BEGIN derivative manager
derivative_manager::derivative_manager(){}

void derivative_manager::initialize(int sdim, int fdim)
{
    raw_pixels = af::flex_double(af::flex_grid<>(sdim,fdim));
    dI=0;
}

void derivative_manager::increment_image(int idx, double value){
    double* floatimage = raw_pixels.begin();
    floatimage[idx] += value;
}

// END derivative manager

//BEGIN unit cell manager
ucell_manager::ucell_manager(){}

void ucell_manager::increment(
    vec3 V,
    mat3 NABC, mat3 UR, vec3 q,
    double Hrad, double Fcell, double Flatt, double fudge,
    double source_I, double capture_fraction, double omega_pixel){

  vec3 dV = NABC*((UR*dB).transpose()) * q;
  double dHrad = V*dV + dV*V;
  double dFlatt = -1*Flatt / 0.63 * fudge * dHrad;
  dI += Fcell*Fcell*2*Flatt*source_I*capture_fraction*omega_pixel*dFlatt;
};


// BEGIN rotation manager begin
rot_manager::rot_manager(){}

void rot_manager::set_R(){assert (false);}

void rot_manager::increment(
        int Na, int Nb, int Nc,
        double hfrac, double kfrac, double lfrac,
        double fudge,
        mat3 U, mat3 X, mat3 Y, mat3 Z,
        vec3 a, vec3 b, vec3 c,  vec3 q,
        double Hrad, double Fcell, double Flatt,
        double source_I, double capture_fraction, double omega_pixel)
{
  /* X,Y,Z will change depending on whether its RotX, RotY, or RotZ manager */
  XYZ = X*Y*Z;
  double dh = U*XYZ*a*q;
  double dk = U*XYZ*b*q;
  double dl = U*XYZ*c*q;
  double dHrad = 2*hfrac*Na*Na*dh +
                 2*kfrac*Nb*Nb*dk +
                 2*lfrac*Nc*Nc*dl;

  double dFlatt = -1*Flatt / 0.63 * fudge * dHrad;

  dI += Fcell*Fcell*2*Flatt*source_I*capture_fraction*omega_pixel*dFlatt;
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
}
void rotY_manager::set_R(){
    R= mat3(cos(value),0, -sin(value),
             0,         1,             0,
            sin(value), 0, cos(value));

    dR= mat3(-sin(value),0, -cos(value),
                0,          0,             0,
                cos(value), 0, -sin(value));
}
void rotZ_manager::set_R(){
    R = mat3(cos(value),  sin(value), 0,
              -sin(value), cos(value), 0,
                         0,           0, 1);

    dR = mat3(-sin(value),  cos(value), 0,
               -cos(value), -sin(value), 0,
                           0,           0, 0);
}
// END rot manager

// BEGIN diffBragg
diffBragg::diffBragg(const dxtbx::model::Detector& detector, const dxtbx::model::Beam& beam,
            int verbose, int panel_id = 0):
    nanoBragg(detector, beam, verbose, panel_id)
    {

    mat3 EYE = mat3(1,0,0,0,1,0,0,0,1);

    RotMats.push_back(EYE);
    RotMats.push_back(EYE);
    RotMats.push_back(EYE);

    dRotMats.push_back(EYE);
    dRotMats.push_back(EYE);
    dRotMats.push_back(EYE);

    R3.push_back(EYE);
    R3.push_back(EYE);
    R3.push_back(EYE);

    boost::shared_ptr<rot_manager> rotX = boost::shared_ptr<rot_manager>(new rotX_manager());
    boost::shared_ptr<rot_manager> rotY = boost::shared_ptr<rot_manager>(new rotY_manager());
    boost::shared_ptr<rot_manager> rotZ = boost::shared_ptr<rot_manager>(new rotZ_manager());

    boost::shared_ptr<ucell_manager> uc1 = boost::shared_ptr<ucell_manager>(new ucell_manager());
    boost::shared_ptr<ucell_manager> uc2 = boost::shared_ptr<ucell_manager>(new ucell_manager());
    boost::shared_ptr<ucell_manager> uc3 = boost::shared_ptr<ucell_manager>(new ucell_manager());
    boost::shared_ptr<ucell_manager> uc4 = boost::shared_ptr<ucell_manager>(new ucell_manager());
    boost::shared_ptr<ucell_manager> uc5 = boost::shared_ptr<ucell_manager>(new ucell_manager());
    boost::shared_ptr<ucell_manager> uc6 = boost::shared_ptr<ucell_manager>(new ucell_manager());

    rotX->refine_me = false;
    rotY->refine_me = false;
    rotZ->refine_me = false;
    uc1->refine_me = false;
    uc2->refine_me = false;
    uc3->refine_me = false;
    uc4->refine_me = false;
    uc5->refine_me = false;
    uc6->refine_me = false;

    rot_managers.push_back(rotX);
    rot_managers.push_back(rotY);
    rot_managers.push_back(rotZ);

    ucell_managers.push_back(uc1);
    ucell_managers.push_back(uc2);
    ucell_managers.push_back(uc3);
    ucell_managers.push_back(uc4);
    ucell_managers.push_back(uc5);
    ucell_managers.push_back(uc6);

    init_raw_pixels_roi();
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
    rot_managers[refine_id]->refine_me=true;
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
}

void diffBragg::set_value( int refine_id, double value ){
    rot_managers[refine_id]->value = value;
    if (refine_id < 3)
        rot_managers[refine_id]->set_R();
}

double diffBragg::get_value( int refine_id){
    return rot_managers[refine_id]->value;
}

af::flex_double diffBragg::get_derivative_pixels(int refine_id){
    return rot_managers[refine_id]->raw_pixels;
}

void diffBragg::zero_raw_pixel_rois(){
    init_raw_pixels_roi();
    initialize_managers();
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
            R3[i_rot] = RotMats[i_rot];
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
    int roi_fdim = roi_xmax - roi_xmin;

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

            /* reset derivative photon count for this pixel */
            for (int i_rot =0 ; i_rot < 3 ; i_rot++){
                if (rot_managers[i_rot]->refine_me){
                    rot_managers[i_rot]->dI =0;
                }
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
                                vector[1] = distance*beam_vector[1];
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
                                ap_vec[0] = ap[1];
                                ap_vec[1] = ap[2];
                                ap_vec[2] = ap[3];

                                bp_vec[0] = bp[1];
                                bp_vec[1] = bp[2];
                                bp_vec[2] = bp[3];

                                cp_vec[0] = cp[1];
                                cp_vec[1] = cp[2];
                                cp_vec[2] = cp[3];

                                q_vec[0] = scattering[1];
                                q_vec[1] = scattering[2];
                                q_vec[2] = scattering[3];


                                mat3 Bmat_realspace;
                                Bmat_realspace[0] = ap_vec[0];
                                Bmat_realspace[3] = ap_vec[1];
                                Bmat_realspace[6] = ap_vec[2];

                                Bmat_realspace[1] = bp_vec[0];
                                Bmat_realspace[4] = bp_vec[1];
                                Bmat_realspace[7] = bp_vec[2];

                                Bmat_realspace[2] = cp_vec[0];
                                Bmat_realspace[5] = cp_vec[1];
                                Bmat_realspace[8] = cp_vec[2];

                                /* construct fractional Miller indicies */
                                vec3 H_vec = (UMATS_RXYZ[mos_tic] * Bmat_realspace).transpose() * q_vec;
                                h = H_vec[0];
                                k = H_vec[1];
                                l = H_vec[2];

                                /* round off to nearest whole index */
                                h0 = static_cast<int>(ceil(h-0.5));
                                k0 = static_cast<int>(ceil(k-0.5));
                                l0 = static_cast<int>(ceil(l-0.5));

                                vec3 H0_vec;
                                H0_vec[0] = h0;
                                H0_vec[1] = k0;
                                H0_vec[2] = l0;

                                mat3 NABC;
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

                                /* convert amplitudes into intensity (photons per steradian) */
                                I += F_cell*F_cell*F_latt*F_latt*source_I[source]*capture_fraction*omega_pixel;

                                /* checkpoint for rotataion derivatives */
                                for (int i_rot =0 ; i_rot < 3 ; i_rot++){
                                    if (rot_managers[i_rot]->refine_me){
                                        R3[i_rot] = dRotMats[i_rot]; // TODO: design upgrade
                                        rot_managers[i_rot]->increment(
                                                            Na, Nb, Nc,
                                                            h-h0, k-k0, l-l0,
                                                            fudge, UMATS[mos_tic],
                                                            R3[0], R3[1], R3[2],
                                                            ap_vec, bp_vec, cp_vec, q_vec,
                                                            hrad_sqr, F_cell, F_latt,
                                                            source_I[source], capture_fraction, omega_pixel);
                                        R3[i_rot] = RotMats[i_rot];
                                    }
                                }
                                /*Checkpoint for unit cell derivatives*/
                                for(int i_uc=0; i_uc < 6; i_uc++ ){
                                    if (ucell_managers[i_uc]->refine_me){
                                        ucell_managers[i_uc]->increment(
                                            V, NABC, UMATS_RXYZ[mos_tic], q_vec,
                                            hrad_sqr, F_cell, F_latt, fudge,
                                            source_I[source], capture_fraction, omega_pixel);
                                    }
                                }
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
            int roi_ss = floor(i / spixels) - roi_ymin ;
            int roi_i =  roi_ss*(roi_fdim+1) + roi_fs ;
            //int iii = i % fpixels;
            //int jjj = floor(i / spixels);
            //if ( roi_fs >= (roi_xmax - roi_xmin))
            //  printf("OOPS:  %d, %d, %d, %d, %d\n", i, fpixels, roi_xmin, roi_i, iii);
            //SCITBX_ASSERT( roi_fs < (roi_xmax - roi_xmin));
            //if ( roi_ss >= (roi_ymax - roi_ymin))
            //  printf("OOPS:  %d, %d, %d, %d, %d\n", i, spixels, roi_ymin, roi_i, jjj);
            //SCITBX_ASSERT( roi_ss < (roi_ymax - roi_ymin));
            //min_i = std::min(roi_i, min_i);
            //max_i = std::max(roi_i, max_i);
            //if (roi_i < 0)
            // printf("OOPS: %d %d %d\n", i, iii, jjj);

            //SCITBX_ASSERT(roi_i >= 0);
            // SCITBX_ASSERT(roi_i < (roi_xmax - roi_xmin + 1)*(roi_ymax - roi_ymin+1) ) ;
            floatimage_roi[roi_i] += r_e_sqr*fluence*spot_scale*polar*I/steps;

            /* udpate the derivative images*/
            for (int i_rot =0 ; i_rot < 3 ; i_rot++){
                if (rot_managers[i_rot]->refine_me){
                    //if (mos_tic==0 && fpixel==0 && spixel==0)
                    //  printf("Updating image for parameter %d <><><<><>>\n", i_rot);
                    double value = r_e_sqr*fluence*spot_scale*polar*rot_managers[i_rot]->dI/steps;
                    rot_managers[i_rot]->increment_image(roi_i, value);
                }
            }

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
