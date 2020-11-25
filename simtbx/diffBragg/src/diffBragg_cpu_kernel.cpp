#include <simtbx/diffBragg/src/diffBragg.h>
#include <assert.h>
#include <stdbool.h>

namespace simtbx { namespace nanoBragg { // BEGIN namespace simtbx::nanoBragg

void diffBragg::diffBragg_sum_over_steps(
        int Npix_to_model, std::vector<unsigned int>& panels_fasts_slows,
        image_type& floatimage,
        image_type& d_Umat_images, image_type& d2_Umat_images,
        image_type& d_Bmat_images, image_type& d2_Bmat_images,
        image_type& d_Ncells_images, image_type& d2_Ncells_images,
        image_type& d_fcell_images, image_type& d2_fcell_images,
        image_type& d_eta_images,
        image_type& d_lambda_images, image_type& d2_lambda_images,
        image_type& d_panel_rot_images, image_type& d2_panel_rot_images,
        image_type& d_panel_orig_images, image_type& d2_panel_orig_images,
        image_type& d_sausage_XYZ_scale_images,
        int* subS_pos, int* subF_pos, int* thick_pos,
        int* source_pos, int* phi_pos, int* mos_pos, int* sausage_pos,
        const int Nsteps, int _printout_fpixel, int _printout_spixel, bool _printout, double _default_F,
        int oversample, bool _oversample_omega, double subpixel_size, double pixel_size,
        double detector_thickstep, double _detector_thick, double close_distance, double detector_attnlen,
        bool use_lambda_coefficients, double lambda0, double lambda1,
        Eigen::Matrix3d& eig_U, Eigen::Matrix3d& eig_O, Eigen::Matrix3d& eig_B, Eigen::Matrix3d& RXYZ,
        std::vector<Eigen::Vector3d,Eigen::aligned_allocator<Eigen::Vector3d> >& dF_vecs,
        std::vector<Eigen::Vector3d,Eigen::aligned_allocator<Eigen::Vector3d> >& dS_vecs,
        std::vector<Eigen::Matrix3d,Eigen::aligned_allocator<Eigen::Matrix3d> >& UMATS_RXYZ,
        std::vector<Eigen::Matrix3d,Eigen::aligned_allocator<Eigen::Matrix3d> >& UMATS_RXYZ_prime,
        std::vector<Eigen::Matrix3d,Eigen::aligned_allocator<Eigen::Matrix3d> >& RotMats,
        std::vector<Eigen::Matrix3d,Eigen::aligned_allocator<Eigen::Matrix3d> >& dRotMats,
        std::vector<Eigen::Matrix3d,Eigen::aligned_allocator<Eigen::Matrix3d> >& d2RotMats,
        std::vector<Eigen::Matrix3d,Eigen::aligned_allocator<Eigen::Matrix3d> >& UMATS,
        std::vector<Eigen::Matrix3d,Eigen::aligned_allocator<Eigen::Matrix3d> >& dB_Mats,
        std::vector<Eigen::Matrix3d,Eigen::aligned_allocator<Eigen::Matrix3d> >& dB2_Mats,
        std::vector<Eigen::Matrix3d,Eigen::aligned_allocator<Eigen::Matrix3d> >& sausages_RXYZ,
        std::vector<Eigen::Matrix3d,Eigen::aligned_allocator<Eigen::Matrix3d> >& d_sausages_RXYZ,
        std::vector<Eigen::Matrix3d,Eigen::aligned_allocator<Eigen::Matrix3d> >& sausages_U,
        std::vector<double>& sausages_scale,
        double* source_X, double* source_Y, double* source_Z, double* source_lambda, double* source_I,
        double kahn_factor,
        double Na, double Nb, double Nc,
        double phi0, double phistep,
        Eigen::Vector3d& spindle_vec, Eigen::Vector3d& _polarization_axis,
        int h_range, int k_range, int l_range,
        int h_max, int h_min, int k_max, int k_min, int l_max, int l_min, double dmin,
        double fudge, bool complex_miller, int verbose, bool only_save_omega_kahn,
        bool isotropic_ncells, bool compute_curvatures,
        std::vector<double>& _FhklLinear, std::vector<double>& _Fhkl2Linear,
        std::vector<bool>& refine_Bmat, std::vector<bool>& refine_Ncells, std::vector<bool>& refine_panel_origin, std::vector<bool>& refine_panel_rot,
        bool refine_fcell, std::vector<bool>& refine_lambda, bool refine_eta, std::vector<bool>& refine_Umat,
        bool refine_sausages,
        int num_sausages,
        std::vector<double>& fdet_vectors, std::vector<double>& sdet_vectors,
        std::vector<double>& odet_vectors, std::vector<double>& pix0_vectors,
        bool _nopolar, bool _point_pixel, double _fluence, double _r_e_sqr, double _spot_scale ) {

    #pragma omp parallel for
    for (int i_pix=0; i_pix < Npix_to_model; i_pix++){
        int _pid = panels_fasts_slows[i_pix*3];
        int _fpixel = panels_fasts_slows[i_pix*3+1];
        int _spixel = panels_fasts_slows[i_pix*3+2];

        // reset photon count for this pixel
        double _I=0;

        // reset derivative photon counts for the various parameters
        double rot_manager_dI[3] = {0,0,0};
        double rot_manager_dI2[3] = {0,0,0};
        double ucell_manager_dI[6]= {0,0,0,0,0,0};
        double ucell_manager_dI2[6]= {0,0,0,0,0,0};
        double Ncells_manager_dI[3]= {0,0,0};
        double Ncells_manager_dI2[3]= {0,0,0};
        double pan_orig_manager_dI[3]= {0,0,0};
        double pan_orig_manager_dI2[3]= {0,0,0};
        double pan_rot_manager_dI[3]= {0,0,0};
        double pan_rot_manager_dI2[3]= {0,0,0};
        double fcell_manager_dI=0;
        double fcell_manager_dI2=0;
        double eta_manager_dI = 0;
        double lambda_manager_dI[2] = {0,0};
        double lambda_manager_dI2[2] = {0,0};
        std::vector<double> sausage_manager_dI(num_sausages*4,0);

        for (int _i_step=0; _i_step < Nsteps; _i_step++){

            int _subS = subS_pos[_i_step];
            int _subF = subF_pos[_i_step];
            int _thick_tic = thick_pos[_i_step];
            int _source = source_pos[_i_step];
            int _phi_tic = phi_pos[_i_step];
            int _mos_tic = mos_pos[_i_step];
            int _sausage_tic = sausage_pos[_i_step];

            /* absolute mm position on detector (relative to its origin) */
            double _Fdet = subpixel_size*(_fpixel*oversample + _subF ) + subpixel_size/2.0;
            double _Sdet = subpixel_size*(_spixel*oversample + _subS ) + subpixel_size/2.0;

            /* assume "distance" is to the front of the detector sensor layer */
            double _Odet = _thick_tic*detector_thickstep;
            int pid_x = _pid*3;
            int pid_y = _pid*3+1;
            int pid_z = _pid*3+2;
            Eigen::Vector3d _o_vec(odet_vectors[pid_x], odet_vectors[pid_y], odet_vectors[pid_z]);

            double pixposX = _Fdet*fdet_vectors[pid_x]+_Sdet*sdet_vectors[pid_x]+_Odet*odet_vectors[pid_x]+pix0_vectors[pid_x];
            double pixposY = _Fdet*fdet_vectors[pid_y]+_Sdet*sdet_vectors[pid_y]+_Odet*odet_vectors[pid_y]+pix0_vectors[pid_y];
            double pixposZ = _Fdet*fdet_vectors[pid_z]+_Sdet*sdet_vectors[pid_z]+_Odet*odet_vectors[pid_z]+pix0_vectors[pid_z];
            Eigen::Vector3d _pixel_pos(pixposX, pixposY, pixposZ);

            double _airpath = _pixel_pos.norm();
            Eigen::Vector3d _diffracted = _pixel_pos/_airpath;

            double _omega_pixel = pixel_size*pixel_size/_airpath/_airpath*close_distance/_airpath;

            if(_point_pixel) _omega_pixel = 1.0/_airpath/_airpath;

            double _capture_fraction = 1;
            if(_detector_thick > 0.0 && detector_attnlen > 0.0)
            {
                double _parallax = _diffracted.dot(_o_vec) ;
                _capture_fraction = exp(-_thick_tic*detector_thickstep/detector_attnlen/_parallax)
                                  -exp(-(_thick_tic+1)*detector_thickstep/detector_attnlen/_parallax);
            }
            Eigen::Vector3d _incident(-source_X[_source], -source_Y[_source], -source_Z[_source]);
            double _lambda = source_lambda[_source];
            double lambda_ang = _lambda*1e10;
            if (use_lambda_coefficients){
                lambda_ang = lambda0 + lambda1*lambda_ang;
                _lambda = lambda_ang*1e-10;
            }

            /* construct the incident beam unit vector while recovering source distance */
            double _source_path = _incident.norm();
            _incident /= _source_path;

            /* construct the scattering vector for this pixel */
            //vec3 _scattering = (_diffracted - _incident) / _lambda;
            Eigen::Vector3d _scattering = (_diffracted - _incident) / _lambda;

            /* sin(theta)/lambda is half the scattering vector length */
            double _stol = 0.5*(_scattering.norm()); //magnitude(scattering);

            /* rough cut to speed things up when we aren't using whole detector */
            if(dmin > 0.0 && _stol > 0.0)
            {
                if(dmin > 0.5/_stol)
                {
                    continue;
                }
            }

            double _phi = phi0 + phistep*_phi_tic;
            Eigen::Matrix3d Bmat_realspace = eig_B;
            if( _phi != 0.0 )
            {
                double cosphi = cos(_phi);
                double sinphi = sin(_phi);
                Eigen::Vector3d ap_vec(eig_B(0,0), eig_B(1,0), eig_B(2,0));
                Eigen::Vector3d bp_vec(eig_B(0,1), eig_B(1,1), eig_B(2,1));
                Eigen::Vector3d cp_vec(eig_B(0,2), eig_B(1,2), eig_B(2,2));

                ap_vec = ap_vec*cosphi + spindle_vec.cross(ap_vec)*sinphi + spindle_vec*(spindle_vec.dot(ap_vec))*(1-cosphi);
                bp_vec = bp_vec*cosphi + spindle_vec.cross(bp_vec)*sinphi + spindle_vec*(spindle_vec.dot(bp_vec))*(1-cosphi);
                cp_vec = cp_vec*cosphi + spindle_vec.cross(cp_vec)*sinphi + spindle_vec*(spindle_vec.dot(cp_vec))*(1-cosphi);

                Bmat_realspace << ap_vec[0], bp_vec[0], cp_vec[0],
                                    ap_vec[1], bp_vec[1], cp_vec[1],
                                    ap_vec[2], bp_vec[2], cp_vec[2];
            }
            Bmat_realspace *= 1e10;

            Eigen::Matrix3d U = sausages_U[_sausage_tic] * eig_U;
            Eigen::Matrix3d UBO = (UMATS_RXYZ[_mos_tic] * U*Bmat_realspace*(eig_O.transpose())).transpose();

            Eigen::Vector3d q_vec(_scattering[0], _scattering[1], _scattering[2]);
            q_vec *= 1e-10;
            Eigen::Vector3d H_vec = UBO*q_vec;

            double _h = H_vec[0];
            double _k = H_vec[1];
            double _l = H_vec[2];

            /* round off to nearest whole index */
            int _h0 = static_cast<int>(ceil(_h - 0.5));
            int _k0 = static_cast<int>(ceil(_k - 0.5));
            int _l0 = static_cast<int>(ceil(_l - 0.5));

            Eigen::Vector3d H0(_h0, _k0, _l0);
            Eigen::Matrix3d _NABC;
            _NABC << Na,0,0,
                    0,Nb,0,
                    0,0,Nc;

            double C = 2 / 0.63 * fudge;
            Eigen::Vector3d delta_H = H_vec - H0;
            Eigen::Vector3d V = _NABC*delta_H;
            double _hrad_sqr = V.dot(V);
            double _F_latt = Na*Nb*Nc*exp(-( _hrad_sqr / 0.63 * fudge ));

            /* no need to go further if result will be zero */
            if(_F_latt == 0.0 && ! only_save_omega_kahn) {
                continue;
            }

            /* structure factor of the unit cell */
            double _F_cell = _default_F;
            double _F_cell2 = 0;

            if ( (_h0<=h_max) && (_h0>=h_min) && (_k0<=k_max) && (_k0>=k_min) && (_l0<=l_max) && (_l0>=l_min)  ) {
                /* just take nearest-neighbor */
                int Fhkl_linear_index = (_h0-h_min) * k_range * l_range + (_k0-k_min) * l_range + (_l0-l_min);
                _F_cell = _FhklLinear[Fhkl_linear_index];
                if (complex_miller) _F_cell2 = _Fhkl2Linear[Fhkl_linear_index];
            }
            //else{
            // _F_cell = _default_F;
            //}

            if (complex_miller)
              _F_cell = sqrt(_F_cell*_F_cell + _F_cell2*_F_cell2);

            /* convert amplitudes into intensity (photons per steradian) */
            if (!_oversample_omega)
                _omega_pixel = 1;

            /* increment to intensity */
            double Iincrement = _F_cell*_F_cell*_F_latt*_F_latt*source_I[_source]*_capture_fraction*_omega_pixel;
            Iincrement *= sausages_scale[_sausage_tic]*sausages_scale[_sausage_tic];
            _I += Iincrement;

            if(verbose > 3)
                printf("hkl= %f %f %f  hkl1= %d %d %d  Fcell=%f\n", _h,_k,_l,_h0,_k0,_l0, _F_cell);

            double two_C = 2*C;
            Eigen::Matrix3d UBOt = U*Bmat_realspace*(eig_O.transpose());
            if (refine_Umat[0]){
                Eigen::Matrix3d RyRzUBOt = RotMats[1]*RotMats[2]*UBOt;
                Eigen::Vector3d delta_H_prime = (UMATS[_mos_tic]*dRotMats[0]*RyRzUBOt).transpose()*q_vec;
                double V_dot_dV = V.dot(_NABC*delta_H_prime);
                double value = -two_C * V_dot_dV * Iincrement;
                double value2 =0;
                if (compute_curvatures) {
                    Eigen::Vector3d delta_H_dbl_prime = (UMATS[_mos_tic]*d2RotMats[0]*RyRzUBOt).transpose()*q_vec;
                    double dV_dot_dV = (_NABC*delta_H_prime).dot(_NABC*delta_H_prime);
                    double dV2_dot_V = (_NABC*delta_H).dot(_NABC*delta_H_dbl_prime);
                    value2 = two_C*(two_C*V_dot_dV*V_dot_dV - dV2_dot_V - dV_dot_dV)*Iincrement;
                }
                rot_manager_dI[0] += value;
                rot_manager_dI2[0] += value2;
            }
            if (refine_Umat[1]){
                Eigen::Matrix3d UmosRx = UMATS[_mos_tic]*RotMats[0];
                Eigen::Matrix3d RzUBOt = RotMats[2]*UBOt;
                Eigen::Vector3d delta_H_prime =(UmosRx*dRotMats[1]*RzUBOt).transpose()*q_vec;
                double V_dot_dV = V.dot(_NABC*delta_H_prime);
                double value = -two_C * V_dot_dV * Iincrement;

                double value2=0;
                if (compute_curvatures){
                    Eigen::Vector3d delta_H_dbl_prime = (UmosRx*d2RotMats[1]*RzUBOt).transpose()*q_vec;
                    double dV_dot_dV = (_NABC*delta_H_prime).dot(_NABC*delta_H_prime);
                    double dV2_dot_V = (_NABC*delta_H).dot(_NABC*delta_H_dbl_prime);
                    value2 = two_C*(two_C*V_dot_dV*V_dot_dV - dV2_dot_V - dV_dot_dV)*Iincrement;
                }
                rot_manager_dI[1] += value;
                rot_manager_dI2[1] += value2;
            }
            if (refine_Umat[2]){
                Eigen::Matrix3d UmosRxRy = UMATS[_mos_tic]*RotMats[0]*RotMats[1];
                Eigen::Vector3d delta_H_prime = (UmosRxRy*dRotMats[2]*UBOt).transpose()*q_vec;
                double V_dot_dV = V.dot(_NABC*delta_H_prime);
                double value = -two_C * V_dot_dV * Iincrement;

                double value2=0;
                if (compute_curvatures){
                    Eigen::Vector3d delta_H_dbl_prime = (UmosRxRy*d2RotMats[2]*UBOt).transpose()*q_vec;
                    double dV_dot_dV = (_NABC*delta_H_prime).dot(_NABC*delta_H_prime);
                    double dV2_dot_V = (_NABC*delta_H).dot(_NABC*delta_H_dbl_prime);
                    value2 = two_C*(two_C*V_dot_dV*V_dot_dV - dV2_dot_V - dV_dot_dV)*Iincrement;
                }
                rot_manager_dI[2] += value;
                rot_manager_dI2[2] += value2;
            }
            /*Checkpoint for unit cell derivatives*/
            Eigen::Matrix3d Ot = eig_O.transpose();
            Eigen::Matrix3d UmosRxRyRzU ;
            Eigen::Vector3d delta_H_prime;
            for(int i_uc=0; i_uc < 6; i_uc++ ){
                if (refine_Bmat[i_uc]){
                    UmosRxRyRzU = UMATS_RXYZ[_mos_tic]*U;
                    delta_H_prime = ((UmosRxRyRzU*dB_Mats[i_uc]*Ot).transpose()*q_vec);
                    double V_dot_dV = V.dot(_NABC*delta_H_prime);
                    double value = -two_C * V_dot_dV * Iincrement;
                    double value2 =0;
                    if (compute_curvatures){
                        Eigen::Vector3d delta_H_dbl_prime = ((UmosRxRyRzU*dB2_Mats[i_uc]*Ot).transpose()*q_vec);
                        double dV_dot_dV = (_NABC*delta_H_prime).dot(_NABC*delta_H_prime);
                        double dV2_dot_V = (_NABC*delta_H).dot(_NABC*delta_H_dbl_prime);
                        value2 = two_C*(two_C*V_dot_dV*V_dot_dV - dV2_dot_V - dV_dot_dV)*Iincrement;
                    }
                    ucell_manager_dI[i_uc] += value;
                    ucell_manager_dI2[i_uc] += value2;
                }
            } /*end ucell deriv */

            /* Checkpoint for Ncells manager */
            if (refine_Ncells[0]){
                int num_ncell_deriv = 1;
                if (not isotropic_ncells)
                    num_ncell_deriv = 3;
                for (int i_nc=0; i_nc < num_ncell_deriv; i_nc++) {
                    Eigen::Matrix3d dN;
                    dN << 0,0,0,0,0,0,0,0,0;
                    dN(i_nc, i_nc) = 1;
                    double N_i = _NABC(i_nc, i_nc);
                    Eigen::Vector3d dV_dN = dN*delta_H;
                    double deriv_coef = 1/N_i - C* ( dV_dN.dot(V));
                    double value = 2*Iincrement*deriv_coef;
                    double value2=0;
                    if(compute_curvatures){
                        dN(i_nc, i_nc) = 0; // TODO check maths
                        value2 = ( -1/N_i/N_i - C*(dV_dN.dot(dV_dN))) *2*Iincrement;
                        value2 += deriv_coef*2*value;
                    }
                    Ncells_manager_dI[i_nc] += value;
                    Ncells_manager_dI2[i_nc] += value2;
                }

            } /* end Ncells manager deriv */

            ///* Checkpoint for Origin manager */
            for (int i_pan_orig=0; i_pan_orig < 3; i_pan_orig++){
                if (refine_panel_origin[i_pan_orig]){
                    double per_k = 1/_airpath;
                    double per_k3 = pow(per_k,3);
                    double per_k5 = pow(per_k,5);
                    double lambda_ang = _lambda*1e10;

                    Eigen::Matrix3d M = -two_C*(_NABC*UBO)/lambda_ang;
                    Eigen::Vector3d dk;
                    if (i_pan_orig == 0)
                        dk << 0,0,1;
                    else if (i_pan_orig == 1)
                        dk << 1,0,0;
                    else
                        dk << 0,1,0;
                    af::flex_double dI_and_dI2 = get_panel_increment(Iincrement, _omega_pixel, M, subpixel_size*subpixel_size,
                        _o_vec, _pixel_pos, per_k,  per_k3, per_k5, V, dk);
                    pan_orig_manager_dI[i_pan_orig] += dI_and_dI2[0];
                    pan_orig_manager_dI2[i_pan_orig] += dI_and_dI2[1];

                } /* end origin manager deriv */
            }

            for (int i_pan_rot=0; i_pan_rot < 3; i_pan_rot++){
                if(refine_panel_rot[i_pan_rot]){
                    double per_k = 1/_airpath;
                    double per_k3 = pow(per_k,3);
                    double per_k5 = pow(per_k,5);
                    double lambda_ang = _lambda*1e10;
                    Eigen::Matrix3d M = -two_C*(_NABC*UBO)/lambda_ang;

                    Eigen::Vector3d dk = _Fdet*(dF_vecs[_pid*3 + i_pan_rot]) + _Sdet*(dS_vecs[_pid*3 + i_pan_rot]);
                    af::flex_double dI_and_dI2 = get_panel_increment(Iincrement, _omega_pixel, M, subpixel_size*subpixel_size,
                        _o_vec, _pixel_pos, per_k,  per_k3, per_k5, V, dk);
                    pan_rot_manager_dI[i_pan_rot] += dI_and_dI2[0];
                    pan_rot_manager_dI2[i_pan_rot] += dI_and_dI2[1];
                }
            }

            /* checkpoint for Fcell manager */
            if (refine_fcell){
                double value = 2*Iincrement/_F_cell ;
                double value2=0;
                if (compute_curvatures){
                    value2 = value/_F_cell;
                }
                fcell_manager_dI += value;
                fcell_manager_dI2 += value2;
            } /* end of fcell man deriv */

            /* checkpoint for eta manager */
            if (refine_eta){
                Eigen::Vector3d DeltaH_deriv = (UMATS_RXYZ_prime[_mos_tic]*UBOt).transpose()*q_vec;
                // vector V is _Nabc*Delta_H
                double value = -two_C*(V.dot(_NABC*DeltaH_deriv))*Iincrement;
                eta_manager_dI += value;
            } /* end of eta man deriv */

            // sausage deriv
            if (refine_sausages){
                Eigen::Matrix3d UBOt = eig_U*Bmat_realspace*(eig_O.transpose());
                int x = _sausage_tic*3;
                int y = _sausage_tic*3+1;
                int z = _sausage_tic*3+2;
                double value=0;
                for (int i=0;i<3; i++){
                    Eigen::Matrix3d UprimeBOt;
                    if (i==0)
                        UprimeBOt = d_sausages_RXYZ[x] * sausages_RXYZ[y] * sausages_RXYZ[z] * UBOt;
                    else if (i==1)
                        UprimeBOt = sausages_RXYZ[x] * d_sausages_RXYZ[y] * sausages_RXYZ[z] * UBOt;
                    else
                        UprimeBOt = sausages_RXYZ[x] * sausages_RXYZ[y] * d_sausages_RXYZ[z] * UBOt;

                    Eigen::Vector3d DeltaH_deriv = (UMATS_RXYZ[_mos_tic]*UprimeBOt).transpose()*q_vec;
                    value = -two_C*(V.dot(_NABC*DeltaH_deriv))*Iincrement;
                    sausage_manager_dI[_sausage_tic*4 + i] += value;
                }
                // sausage scale derivative
                value = 2* Iincrement / sausages_scale[_sausage_tic];
                sausage_manager_dI[_sausage_tic*4 + 3] += value;
            }
            // end of sausage deriv

            /*checkpoint for lambda manager*/
            for(int i_lam=0; i_lam < 2; i_lam++){
                if (refine_lambda[i_lam]){
                    double lambda_ang = _lambda*1e10;
                    double NH_dot_V = (_NABC*H_vec).dot(V);
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

            if( _printout && _i_step==0 ){
                if((_fpixel==_printout_fpixel && _spixel==_printout_spixel) || _printout_fpixel < 0)
                {
                    printf("%4d %4d : stol = %g, lambda = %g\n", _fpixel,_spixel,_stol, _lambda);
                    printf("at %g %g %g\n", _pixel_pos[0],_pixel_pos[1],_pixel_pos[2]);
                    printf("hkl= %f %f %f  hkl0= %d %d %d\n", _h,_k,_l,_h0,_k0,_l0);
                    printf(" F_cell=%g  F_latt=%g   I = %g\n", _F_cell,_F_latt,_I);
                    printf("I/steps %15.10g\n", _I/steps);
                    printf("Fdet= %g; Sdet= %g ; Odet= %g\n", _Fdet, _Sdet, _Odet);
                    printf("omega   %15.10g\n", _omega_pixel);
                    printf("default_F= %f\n", _default_F);
                    printf("PIX0: %f %f %f\n" , pix0_vectors[pid_x], pix0_vectors[pid_y], pix0_vectors[pid_z]);
                    printf("F: %f %f %f\n" , fdet_vectors[pid_x], fdet_vectors[pid_y], fdet_vectors[pid_z]);
                    printf("S: %f %f %f\n" , sdet_vectors[pid_x], sdet_vectors[pid_y], sdet_vectors[pid_z]);
                    printf("O: %f %f %f\n" , odet_vectors[pid_x], odet_vectors[pid_y], odet_vectors[pid_z]);
                    printf("QVECTOR: %f %f %f\n" , q_vec[0], q_vec[1], q_vec[2]);
                    printf("MOSAIC UMAT RXYZ\n");
                    std::cout << UMATS_RXYZ[_mos_tic] << std::endl;
                    printf("Bmat_realspace\n");
                    std::cout << Bmat_realspace << std::endl;
                    printf("UBO\n");
                    std::cout << UBO << std::endl;
                    printf("UBOt\n");
                    std::cout << UBOt << std::endl;
                    printf("UmosRxRyRzU\n");
                    std::cout << UmosRxRyRzU << std::endl;
                    printf("deltaHprime\n");
                    std::cout << delta_H_prime << std::endl;
                    printf("Iincrement: %f\n", Iincrement);
                    printf("pid_x=%d, pid_y=%d; pid_z=%d\n", pid_x, pid_y, pid_z);
                    SCITBX_EXAMINE(fdet_vectors[0]);
                    SCITBX_EXAMINE(fdet_vectors[1]);
                    SCITBX_EXAMINE(fdet_vectors[2]);
                    SCITBX_EXAMINE(pix0_vectors[0]);
                    SCITBX_EXAMINE(pix0_vectors[1]);
                    SCITBX_EXAMINE(pix0_vectors[2]);
                    for (int i=0; i < num_sausages; i++){
                        printf("Matrix for sausage %d (scale=%f):\n", i, sausages_scale[i]);
                        std::cout << sausages_U[i] << std::endl;
                    }
                }
            }
        } /* end of i_steps loop */

        double _Fdet_ave = pixel_size*_fpixel + pixel_size/2.0;
        double _Sdet_ave = pixel_size*_spixel + pixel_size/2.0;
        double _Odet_ave = 0; //Odet; // TODO maybe make this more general for thick detectors?

        Eigen::Vector3d _pixel_pos_ave(0,0,0);
        int pid_x = _pid*3;
        int pid_y = _pid*3+1;
        int pid_z = _pid*3+2;
        _pixel_pos_ave[0] = _Fdet_ave * fdet_vectors[pid_x]+_Sdet_ave*sdet_vectors[pid_x]+_Odet_ave*odet_vectors[pid_x]+pix0_vectors[pid_x];
        _pixel_pos_ave[1] = _Fdet_ave * fdet_vectors[pid_y]+_Sdet_ave*sdet_vectors[pid_y]+_Odet_ave*odet_vectors[pid_y]+pix0_vectors[pid_y];
        _pixel_pos_ave[2] = _Fdet_ave * fdet_vectors[pid_z]+_Sdet_ave*sdet_vectors[pid_z]+_Odet_ave*odet_vectors[pid_z]+pix0_vectors[pid_z];

        double _airpath_ave = _pixel_pos_ave.norm();
        Eigen::Vector3d _diffracted_ave = _pixel_pos_ave/_airpath_ave;
        double _omega_pixel_ave = pixel_size*pixel_size/_airpath_ave/_airpath_ave*close_distance/_airpath_ave;

        double _polar = 1;
        if (!_nopolar){
            Eigen::Vector3d _incident(-source_X[0], -source_Y[0], -source_Z[0]);
            _incident = _incident / _incident.norm();
            // component of diffracted unit vector along incident beam unit vector
            double cos2theta = _incident.dot(_diffracted_ave);
            double cos2theta_sqr = cos2theta*cos2theta;
            double sin2theta_sqr = 1-cos2theta_sqr;

            double _psi=0;
            if(kahn_factor != 0.0){
                // cross product to get "vertical" axis that is orthogonal to the cannonical "polarization"
                Eigen::Vector3d B_in = _polarization_axis.cross(_incident);
                // cross product with incident beam to get E-vector direction
                Eigen::Vector3d E_in = _incident.cross(B_in);
                // get components of diffracted ray projected onto the E-B plane
                double _kEi = _diffracted_ave.dot(E_in);
                double _kBi = _diffracted_ave.dot(B_in);
                // compute the angle of the diffracted ray projected onto the incident E-B plane
                _psi = -atan2(_kBi,_kEi);
            }
            // correction for polarized incident beam
            _polar = 0.5*(1.0 + cos2theta_sqr - kahn_factor*cos(2*_psi)*sin2theta_sqr);
        }

        double _om = 1;
        if (!_oversample_omega)
            _om=_omega_pixel_ave;

        // final scale term to being everything to photon number units
        double _scale_term = _r_e_sqr*_fluence*_spot_scale*_polar*_om/Nsteps*num_sausages;

        int roi_i = i_pix; // TODO replace roi_i with i_pix

        floatimage[i_pix] = _scale_term*_I;

        /* udpate the rotation derivative images*/
        for (int i_rot =0 ; i_rot < 3 ; i_rot++){
            if (refine_Umat[i_rot]){
                double value = _scale_term*rot_manager_dI[i_rot];
                double value2 = _scale_term*rot_manager_dI2[i_rot];
                int idx = i_rot*Npix_to_model + i_pix;
                d_Umat_images[idx] = value;
                d2_Umat_images[idx] = value2;
            }
        } /* end rot deriv image increment */

        /*update the ucell derivative images*/
        for (int i_uc=0 ; i_uc < 6 ; i_uc++){
            if (refine_Bmat[i_uc]){
                double value = _scale_term*ucell_manager_dI[i_uc];
                double value2 = _scale_term*ucell_manager_dI2[i_uc];
                int idx= i_uc*Npix_to_model + i_pix;
                d_Bmat_images[idx] = value;
                d2_Bmat_images[idx] = value2;
            }
        }/* end ucell deriv image increment */

        /*update the Ncells derivative image*/
        if (refine_Ncells[0]){
            double value = _scale_term*Ncells_manager_dI[0];
            double value2 = _scale_term*Ncells_manager_dI2[0];
            double idx = i_pix;
            d_Ncells_images[idx] = value;
            d2_Ncells_images[idx] = value2;

            if (! isotropic_ncells){
                value = _scale_term*Ncells_manager_dI[1];
                value2 = _scale_term*Ncells_manager_dI2[1];
                idx = Npix_to_model + i_pix;
                d_Ncells_images[idx] = value;
                d2_Ncells_images[idx] = value2;

                value = _scale_term*Ncells_manager_dI[2];
                value2 = _scale_term*Ncells_manager_dI2[2];
                idx = Npix_to_model*2 + i_pix;
                d_Ncells_images[idx] = value;
                d2_Ncells_images[idx] = value2;
            }
        }/* end Ncells deriv image increment */

        /* update Fcell derivative image */
        if(refine_fcell){
            double value = _scale_term*fcell_manager_dI;
            double value2 = _scale_term*fcell_manager_dI2;
            d_fcell_images[i_pix] = value;
            d2_fcell_images[i_pix] = value2;
        }/* end Fcell deriv image increment */

        /* update eta derivative image */
        if(refine_eta){
            double value = _scale_term*eta_manager_dI;
            double value2 = 0;
            d_eta_images[i_pix] = value;
        }/* end eta deriv image increment */

        /*update the lambda derivative images*/
        for (int i_lam=0 ; i_lam < 2 ; i_lam++){
            if (refine_lambda[i_lam]){
                //double value = scale_term*lambda_managers[i_lam]->dI;
                //double value2 = scale_term*lambda_managers[i_lam]->dI2;
                double value = _scale_term*lambda_manager_dI[i_lam];
                double value2 = _scale_term*lambda_manager_dI2[i_lam];
                int idx = i_lam*Npix_to_model + i_pix;
                d_lambda_images[idx] = value;
                d2_lambda_images[idx] = value2;
            }
        }/* end lambda deriv image increment */

        // sausage increment
        if (refine_sausages){
            for (int i_sausage=0; i_sausage<num_sausages; i_sausage++){
                for (int i=0; i < 4; i++){
                    int sausage_parameter_i = i_sausage*4+i;
                    double value = _scale_term*sausage_manager_dI[sausage_parameter_i];
                    int idx = sausage_parameter_i*Npix_to_model + i_pix;
                    d_sausage_XYZ_scale_images[idx] = value;
                }
            }
        }
        // end sausage

        // panel rotation
        for (int i_pan_rot=0; i_pan_rot < 3; i_pan_rot++){
            if(refine_panel_rot[i_pan_rot]){
                double value = _scale_term*pan_rot_manager_dI[i_pan_rot];
                double value2 = _scale_term*pan_rot_manager_dI2[i_pan_rot];
                int idx = i_pan_rot*Npix_to_model + i_pix;
                d_panel_rot_images[idx] = value;
                d2_panel_rot_images[idx] = value2;
            }
        }/* end panel rot deriv image increment */

        // panel origin
        for (int i_pan_orig=0; i_pan_orig < 3; i_pan_orig++){
            if(refine_panel_origin[i_pan_orig]){
                double value = _scale_term*pan_orig_manager_dI[i_pan_orig];
                double value2 = _scale_term*pan_orig_manager_dI2[i_pan_orig];
                int idx = i_pan_orig*Npix_to_model + i_pix;
                d_panel_orig_images[idx] = value;
                d2_panel_orig_images[idx] = value2;
            }/* end panel orig deriv image increment */
        }
    } // end i_pix loop
} // END of CPU kernel

}} // END namespace simtbx::nanoBragg
