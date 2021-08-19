#include <simtbx/diffBragg/src/diffBragg.h>
#include <assert.h>
#include <stdbool.h>
#include<unordered_map>
#include<string>

namespace simtbx { namespace nanoBragg { // BEGIN namespace simtbx::nanoBragg

void diffBragg::diffBragg_sum_over_steps(
        int Npix_to_model, std::vector<unsigned int>& panels_fasts_slows,
        image_type& floatimage,
        images& d_image,
        images& d2_image,
        int* subS_pos, int* subF_pos, int* thick_pos,
        int* source_pos, int* phi_pos, int* mos_pos, int* sausage_pos,
        const int Nsteps, int _printout_fpixel, int _printout_spixel, bool _printout, double _default_F,
        int oversample, bool _oversample_omega, double subpixel_size, double pixel_size,
        double detector_thickstep, double _detector_thick, std::vector<double>& close_distances, double detector_attnlen,
        bool use_lambda_coefficients, double lambda0, double lambda1,
        Eigen::Matrix3d& eig_U, Eigen::Matrix3d& eig_O, Eigen::Matrix3d& eig_B, Eigen::Matrix3d& RXYZ,
        std::vector<Eigen::Vector3d,Eigen::aligned_allocator<Eigen::Vector3d> >& dF_vecs,
        std::vector<Eigen::Vector3d,Eigen::aligned_allocator<Eigen::Vector3d> >& dS_vecs,
        std::vector<Eigen::Matrix3d,Eigen::aligned_allocator<Eigen::Matrix3d> >& UMATS_RXYZ,
        std::vector<Eigen::Matrix3d,Eigen::aligned_allocator<Eigen::Matrix3d> >& UMATS_RXYZ_prime,
        std::vector<Eigen::Matrix3d,Eigen::aligned_allocator<Eigen::Matrix3d> >& UMATS_RXYZ_dbl_prime,
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
        double Nd, double Ne, double Nf,
        double phi0, double phistep,
        Eigen::Vector3d& spindle_vec, Eigen::Vector3d& _polarization_axis,
        int h_range, int k_range, int l_range,
        int h_max, int h_min, int k_max, int k_min, int l_max, int l_min, double dmin,
        double fudge, bool complex_miller, int verbose, bool only_save_omega_kahn,
        bool isotropic_ncells, bool compute_curvatures,
        std::vector<double>& _FhklLinear, std::vector<double>& _Fhkl2Linear,
        std::vector<bool>& refine_Bmat, std::vector<bool>& refine_Ncells, bool refine_Ncells_def,  std::vector<bool>& refine_panel_origin, std::vector<bool>& refine_panel_rot,
        bool refine_fcell, std::vector<bool>& refine_lambda, bool refine_eta, std::vector<bool>& refine_Umat,
        bool refine_sausages,
        int num_sausages,
        bool refine_fp_fdp,
        std::vector<double>& fdet_vectors, std::vector<double>& sdet_vectors,
        std::vector<double>& odet_vectors, std::vector<double>& pix0_vectors,
        bool _nopolar, bool _point_pixel, double _fluence, double _r_e_sqr, double _spot_scale ,
        bool no_Nabc_scale,
        std::vector<double>& fpfdp, std::vector<double>& fpfdp_derivs,
        std::vector<double>& atom_data, bool track_Fhkl, std::vector<int>& nominal_hkl) {

    #pragma omp parallel for
    for (int i_pix=0; i_pix < Npix_to_model; i_pix++){
        int _pid = panels_fasts_slows[i_pix*3];
        int _fpixel = panels_fasts_slows[i_pix*3+1];
        int _spixel = panels_fasts_slows[i_pix*3+2];
        double _close_distance = close_distances[_pid];
        //std::unordered_map<int, int> Fhkl_tracker;
        std::unordered_map<std::string, int> Fhkl_tracker;
        bool use_nominal_hkl = false;
        if (!nominal_hkl.empty())
            use_nominal_hkl = true;
        // reset photon count for this pixel
        double _I=0;
        double II_max = -1;
        double max_stats[11] = {0,0,0,0,0,
                               0,0,0,0,0,0};
        double max_stats2[11] = {0,0,0,0,0,
                               0,0,0,0,0,0};
        // reset derivative photon counts for the various parameters
        double rot_manager_dI[3] = {0,0,0};
        double rot_manager_dI2[3] = {0,0,0};
        double ucell_manager_dI[6]= {0,0,0,0,0,0};
        double ucell_manager_dI2[6]= {0,0,0,0,0,0};
        double Ncells_manager_dI[6]= {0,0,0,0,0,0};
        double Ncells_manager_dI2[6]= {0,0,0,0,0,0};
        double pan_orig_manager_dI[3]= {0,0,0};
        double pan_orig_manager_dI2[3]= {0,0,0};
        double pan_rot_manager_dI[3]= {0,0,0};
        double pan_rot_manager_dI2[3]= {0,0,0};
        double fcell_manager_dI[3] = {0,0,0};
        double fcell_manager_dI2[3] = {0,0,0};
        double eta_manager_dI[3] = {0,0,0};
        double eta_manager_dI2[3] = {0,0,0};
        double lambda_manager_dI[2] = {0,0};
        double lambda_manager_dI2[2] = {0,0};
        double fp_fdp_manager_dI[2] = {0,0};
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

            double _omega_pixel = pixel_size*pixel_size/_airpath/_airpath*_close_distance/_airpath;

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
            _NABC << Na,Nd,Nf,
                    Nd,Nb,Ne,
                    Nf,Ne,Nc;

            double C = 2 / 0.63 * fudge;
            Eigen::Vector3d delta_H = H_vec - H0;
            Eigen::Vector3d V = _NABC*delta_H;
            double _hrad_sqr = V.dot(V);
            double _F_latt;
            if (no_Nabc_scale)
                _F_latt = exp(-( _hrad_sqr / 0.63 * fudge ));
            else
                _F_latt = Na*Nb*Nc*exp(-( _hrad_sqr / 0.63 * fudge ));

            /* no need to go further if result will be zero */
            //if(_F_latt == 0.0 && ! only_save_omega_kahn) {
            //    continue;
            //}
            /* convert amplitudes into intensity (photons per steradian) */
            if (!_oversample_omega)
                _omega_pixel = 1;

            /* increment to intensity */
            double sauce = pow(sausages_scale[_sausage_tic],2);
            double I_noFcell = _F_latt*_F_latt*source_I[_source]*_capture_fraction*_omega_pixel*sauce;

            /* structure factor of the unit cell */
            double _F_cell = _default_F;
            double _F_cell2 = 0;

            if ( (_h0<=h_max) && (_h0>=h_min) && (_k0<=k_max) && (_k0>=k_min) && (_l0<=l_max) && (_l0>=l_min)  ) {
                /* just take nearest-neighbor */
                int Fhkl_linear_index = (_h0-h_min) * k_range * l_range + (_k0-k_min) * l_range + (_l0-l_min);
                _F_cell = _FhklLinear[Fhkl_linear_index];
                if (track_Fhkl){
                    std::string hkl_s ;
                    hkl_s = std::to_string(_h0) + ","+ std::to_string(_k0) + "," + std::to_string(_l0);
                    if (Fhkl_tracker.count(hkl_s))
                        Fhkl_tracker[hkl_s] += 1;
                    else
                        Fhkl_tracker[hkl_s] = 0;
                    continue;
                }
                if (complex_miller) _F_cell2 = _Fhkl2Linear[Fhkl_linear_index];
            }
            //else{
            // _F_cell = _default_F;
            //}
            bool do_max = false;
            if(verbose > 3){//} && _i_step==0){
                double F2 = sqrt(_F_cell*_F_cell + _F_cell2*_F_cell2);
                double II = I_noFcell*F2*F2;
                if (I_noFcell > II_max){
                    do_max = true;
                    max_stats[0] = _F_cell;
                    max_stats[1] = _F_cell2;
                    max_stats[2] = I_noFcell;
                    max_stats[3] = _F_latt;
                    max_stats[4] = II;
                    max_stats[5] = _omega_pixel;
                    max_stats[6] = source_I[_source];
                    max_stats[7]= source_lambda[_source];
                    max_stats[8] = _h;
                    max_stats[9] = _k;
                    max_stats[10] = _l;
                    II_max = I_noFcell;
                }
                else{
                do_max=false;}
                //printf("hkl= %f %f %f  hkl1= %d %d %d  |Fcell| before=%f, Ino=%f, I=%f, Flatt=%f, cap_frac=%f, omega_pix=%f, sourceI=%f\n", _h,_k,_l,_h0,_k0,_l0, F2, I_noFcell, II, _F_latt,
                //_capture_fraction, _omega_pixel, source_I[_source] );
            }
            double c_deriv_Fcell_real = 0;
            double c_deriv_Fcell_imag = 0;
            double d_deriv_Fcell_real = 0;
            double d_deriv_Fcell_imag = 0;
            double c_deriv_Fcell = 0;
            double d_deriv_Fcell = 0;
            if (complex_miller){
            // TODO shouldnt this be constant for each HKL?
              if (fpfdp.size() > 0){
                   double S_2 = 1.e-20*(_scattering[0]*_scattering[0]+_scattering[1]*_scattering[1]+_scattering[2]*_scattering[2]);

                    // fp is always followed by the fdp value
                   double val_fp = fpfdp[2*_source];
                   double val_fdp = fpfdp[2*_source+1];

                   double c_deriv_prime=0;
                   double c_deriv_dblprime=0;
                   double d_deriv_prime = 0;
                   double d_deriv_dblprime = 0;
                   if (refine_fp_fdp){
                   //   currently only supports two parameter model
                       int nsources_times_two = fpfdp.size();
                       int d_idx = 2*_source;
                       c_deriv_prime = fpfdp_derivs[d_idx];
                       c_deriv_dblprime = fpfdp_derivs[d_idx+1];
                       d_deriv_prime = fpfdp_derivs[d_idx+nsources_times_two];
                       d_deriv_dblprime = fpfdp_derivs[d_idx+1+nsources_times_two];
                   }
                   // 5 valeus per atom: x,y,z,B,occupancy
                   int num_atoms = atom_data.size()/5;
                   for (int  i_atom=0; i_atom < num_atoms; i_atom++){
                       if (verbose>5)
                         printf("Processing atom %d, _F_cell=%10.3f\n", i_atom, _F_cell);
                       // fractional atomic coordinates
                       double atom_x = atom_data[i_atom*5];
                       double atom_y = atom_data[i_atom*5+1];
                       double atom_z = atom_data[i_atom*5+2];
                       double B = atom_data[i_atom*5+3]; // B factor
                       B = exp(-B*S_2/4.0); // TODO: speed me up?
                       double occ = atom_data[i_atom*5+4]; // occupancy
                       double r_dot_h = _h0*atom_x + _k0*atom_y + _l0*atom_z;
                       double phase = 2*M_PI*r_dot_h;
                       double s_rdoth = sin(phase);
                       double c_rdoth = cos(phase);
                       double Bocc = B*occ;
                       double BC = Bocc*c_rdoth;
                       double BS = Bocc*s_rdoth;
                       double real_part = BC*val_fp - BS*val_fdp;
                       double imag_part = BS*val_fp + BC*val_fdp;
                       _F_cell += real_part;
                       _F_cell2 += imag_part;
                       if (refine_fp_fdp){
                            c_deriv_Fcell_real += BC*c_deriv_prime - BS*c_deriv_dblprime;
                            c_deriv_Fcell_imag += BS*c_deriv_prime + BC*c_deriv_dblprime;

                            d_deriv_Fcell_real += BC*d_deriv_prime - BS*d_deriv_dblprime;
                            d_deriv_Fcell_imag += BS*d_deriv_prime + BC*d_deriv_dblprime;
                       }
                   }
               }
               double Freal = _F_cell;
               double Fimag = _F_cell2;
               _F_cell = sqrt(pow(Freal,2) + pow(Fimag,2));// TODO work around the sqrt ?
               //if (verbose> 3 && _i_step==0){
               //     double II = I_noFcell*_F_cell*_F_cell;
               //     printf("hkl= %f %f %f  hkl1= %d %d %d  |Fcell| after=%f, Ino=%f, I=%f\n", _h,_k,_l,_h0,_k0,_l0, _F_cell, I_noFcell, II);
               // }
               if (refine_fp_fdp){
                   c_deriv_Fcell = Freal*c_deriv_Fcell_real + Fimag*c_deriv_Fcell_imag;
                   d_deriv_Fcell = Freal*d_deriv_Fcell_real + Fimag*d_deriv_Fcell_imag;
               }
            }

            if(verbose > 3){
                double F2 = sqrt(_F_cell*_F_cell + _F_cell2*_F_cell2);
                double II = I_noFcell*F2*F2;
                if (do_max){
                    max_stats2[0] = F2 ;
                    max_stats2[1] = 0 ; //_F_cell2;
                    max_stats2[2] = I_noFcell;
                    max_stats2[3] = _F_latt;
                    max_stats2[4] = II;
                    max_stats2[5] = _omega_pixel;
                    max_stats2[6] = source_I[_source];
                    max_stats2[7]= source_lambda[_source];
                    max_stats2[8] = _h;
                    max_stats2[9] = _k;
                    max_stats2[10] = _l;
                }
            }
            double Iincrement = _F_cell*_F_cell*I_noFcell;
            //Iincrement *= sausages_scale[_sausage_tic]*sausages_scale[_sausage_tic];
            _I += Iincrement;

            // TODO if test?
            if (refine_fp_fdp){
                fp_fdp_manager_dI[0] += 2*I_noFcell * (c_deriv_Fcell);
                fp_fdp_manager_dI[1] += 2*I_noFcell * (d_deriv_Fcell);
            }

            //if(verbose > 3)
            //    printf("hkl= %f %f %f  hkl1= %d %d %d  Fcell=%f\n", _h,_k,_l,_h0,_k0,_l0, _F_cell);

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
                if (! isotropic_ncells)
                    num_ncell_deriv = 3;
                for (int i_nc=0; i_nc < num_ncell_deriv; i_nc++) {
                    Eigen::Matrix3d dN;
                    if (!isotropic_ncells){
                        dN << 0,0,0,0,0,0,0,0,0;
                        dN(i_nc, i_nc) = 1;
                        }
                    else
                        dN << 1,0,0,0,1,0,0,0,1;

                    double N_i = _NABC(i_nc, i_nc);
                    Eigen::Vector3d dV_dN = dN*delta_H;
                    //if (_i_step==0 && _fpixel ==0 && _spixel == 0){
                    //    printf("dN matrix: %f %f %f\n %f %f %f\n %f %f %f\n",
                    //        dN(0,0), dN(0,1), dN(0,2),
                    //        dN(1,0), dN(1,1), dN(1,2),
                    //        dN(2,0), dN(2,1), dN(2,2)
                    //        );
                    //    printf("NABC matrix: %f %f %f\n %f %f %f\n %f %f %f\n",
                    //        _NABC(0,0), _NABC(0,1), _NABC(0,2),
                    //        _NABC(1,0), _NABC(1,1), _NABC(1,2),
                    //        _NABC(2,0), _NABC(2,1), _NABC(2,2)
                    //        );
                    //}
                    double deriv_coef;
                    if (isotropic_ncells)
                        deriv_coef= 3/N_i - C* ( dV_dN.dot(V));
                    else
                        deriv_coef= 1/N_i - C* ( dV_dN.dot(V));
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

            if (refine_Ncells_def){
                for (int i_nc =3; i_nc < 6; i_nc++ ){
                    Eigen::Matrix3d dN;
                    if (i_nc ==3)
                        dN << 0,1,0,1,0,0,0,0,0;
                    else if (i_nc == 4)
                        dN << 0,0,0,0,0,1,0,1,0;
                    else
                        dN << 0,0,1,0,0,0,1,0,0;
                    Eigen::Vector3d dV_dN = dN*delta_H;
                    double deriv_coef = -C* (2* dV_dN.dot(V));
                    double value = Iincrement*deriv_coef;
                    Ncells_manager_dI[i_nc] += value;
                    double value2 =0;
                    if (compute_curvatures){
                        value2 = deriv_coef*value;
                        value2 +=  -2*C*Iincrement*(dV_dN.dot(dV_dN));
                    Ncells_manager_dI2[i_nc] += value2;
                    }
                }
            }

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
            //    TODO rewrite so no divide by Fcell
                int nom_h=_h0;
                int nom_k=_k0;
                int nom_l=_l0;
                //int f_cell_idx = 1;
                if (use_nominal_hkl){
                    nom_h = nominal_hkl[i_pix*3];
                    nom_k = nominal_hkl[i_pix*3+1];
                    nom_l = nominal_hkl[i_pix*3+2];
                    //f_cell_idx = _l0 - nom_l + 1;
                }
                double value = 2*I_noFcell*_F_cell; //2*Iincrement/_F_cell ;
                double value2=0;
                if (compute_curvatures){
                    if (_F_cell > 0)
                        value2 = value/_F_cell;
                }
                //if (f_cell_idx >= 0 && f_cell_idx <= 2){
                // NOTE use nominal hkl to regulate when gradients are compputed, as h,k,l can drift within a shoebox
                if (use_nominal_hkl){
                    if (nom_h==_h0 && nom_k==_k0 && nom_l==_l0 ){
                        fcell_manager_dI[1] += value;
                        fcell_manager_dI2[1] += value2;
                    }
                }
                else{
                    fcell_manager_dI[1] += value;
                    fcell_manager_dI2[1] += value2;
                }
            } /* end of fcell man deriv */

            /* checkpoint for eta manager */
            if (refine_eta){
                for (int i_eta=0; i_eta < 3; i_eta++){
                    if (i_eta != 0 && UMATS_RXYZ_prime.size()==UMATS_RXYZ.size()){
                        continue;
                    }
                    int mtic2 = _mos_tic  + i_eta*mosaic_domains;
                    Eigen::Vector3d DeltaH_deriv = (UMATS_RXYZ_prime[mtic2]*UBOt).transpose()*q_vec;
                    // vector V is _Nabc*Delta_H
                    Eigen::Vector3d dV = _NABC*DeltaH_deriv;
                    double V_dot_dV = V.dot(dV);
                    double Iprime = -two_C*(V_dot_dV)*Iincrement;
                    eta_manager_dI[i_eta] += Iprime;

                    double Idbl_prime = 0;
                    if (compute_curvatures){
                        Eigen::Vector3d DeltaH_second_deriv = (UMATS_RXYZ_dbl_prime[mtic2]*UBOt).transpose()*q_vec;
                        Eigen::Vector3d dV2 = _NABC*DeltaH_second_deriv;
                        Idbl_prime = -two_C*(dV.dot(dV) + V.dot(dV2))*Iincrement;
                        Idbl_prime += -two_C*(V_dot_dV)*Iprime;
                    }
                    eta_manager_dI2[i_eta] += Idbl_prime;
                }
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
                    printf(" F_cell=%g  F_cell_2=%g F_latt=%g   I = %g\n", _F_cell,_F_cell2,_F_latt,_I);
                    printf("I/steps %15.10g\n", _I/steps);
                    printf("cap frac   %f\n", _capture_fraction);
                    printf("sauce   %f\n", sauce);
                    printf("Fdet= %g; Sdet= %g ; Odet= %g\n", _Fdet, _Sdet, _Odet);
                    printf("omega   %15.10g\n", _omega_pixel);
                    printf("close_distance    %15.10g\n", _close_distance);
                    printf("fluence    %15.10g\n", _fluence);
                    printf("spot scale   %15.10g\n", _spot_scale);
                    printf("sourceI[0]:   %15.10g\n", source_I[0]);
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

        if (verbose >3){
            printf("hkl= %f %f %f  |Fcell| before=(%f,%f), Ino=%f, I=%f, Flatt=%f, omega_pix=%2.3e, sourceI=%f, sourceLambda=%f\n",
             max_stats[8],max_stats[9],max_stats[10],
             max_stats[0], max_stats[1], max_stats[2], max_stats[4], max_stats[3],
            max_stats[5], max_stats[6], max_stats[7]*1e10 );

            printf("hkl= %f %f %f  |Fcell| after=(%f,%f), Ino=%f, I=%f, Flatt=%f, omega_pix=%2.3e, sourceI=%f, sourceLambda=%f\n",
             max_stats2[8],max_stats2[9],max_stats2[10],
             max_stats2[0], max_stats2[1], max_stats2[2], max_stats2[4], max_stats2[3],
            max_stats2[5], max_stats2[6], max_stats2[7]*1e10 );

        }

        double _scale_term=1;
        if (!track_Fhkl){

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
            double _omega_pixel_ave = pixel_size*pixel_size/_airpath_ave/_airpath_ave*_close_distance/_airpath_ave;

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
            _scale_term = _r_e_sqr*_fluence*_spot_scale*_polar*_om/Nsteps*num_sausages;

            //int roi_i = i_pix; // TODO replace roi_i with i_pix

            floatimage[i_pix] = _scale_term*_I;

        }
        /* udpate the rotation derivative images*/
        for (int i_rot =0 ; i_rot < 3 ; i_rot++){
            if (refine_Umat[i_rot]){
                double value = _scale_term*rot_manager_dI[i_rot];
                double value2 = _scale_term*rot_manager_dI2[i_rot];
                int idx = i_rot*Npix_to_model + i_pix;
                d_image.Umat[idx] = value;
                d2_image.Umat[idx] = value2;
            }
        } /* end rot deriv image increment */

        /*update the ucell derivative images*/
        for (int i_uc=0 ; i_uc < 6 ; i_uc++){
            if (refine_Bmat[i_uc]){
                double value = _scale_term*ucell_manager_dI[i_uc];
                double value2 = _scale_term*ucell_manager_dI2[i_uc];
                int idx= i_uc*Npix_to_model + i_pix;
                d_image.Bmat[idx] = value;
                d2_image.Bmat[idx] = value2;
            }
        }/* end ucell deriv image increment */

        /*update the Ncells derivative image*/
        if (refine_Ncells[0]){
            double value = _scale_term*Ncells_manager_dI[0];
            double value2 = _scale_term*Ncells_manager_dI2[0];
            double idx = i_pix;
            d_image.Ncells[idx] = value;
            d2_image.Ncells[idx] = value2;
            //d_Ncells_images[idx] = value;
            //d2_Ncells_images[idx] = value2;

            if (! isotropic_ncells){
                value = _scale_term*Ncells_manager_dI[1];
                value2 = _scale_term*Ncells_manager_dI2[1];
                idx = Npix_to_model + i_pix;
                d_image.Ncells[idx] = value;
                d2_image.Ncells[idx] = value2;
                //d_Ncells_images[idx] = value;
                //d2_Ncells_images[idx] = value2;

                value = _scale_term*Ncells_manager_dI[2];
                value2 = _scale_term*Ncells_manager_dI2[2];
                idx = Npix_to_model*2 + i_pix;
                d_image.Ncells[idx] = value;
                d2_image.Ncells[idx] = value2;
                //d_Ncells_images[idx] = value;
                //d2_Ncells_images[idx] = value2;
            }
        }/* end Ncells deriv image increment */
        if (refine_Ncells_def){
            for (int i_nc=3; i_nc<6; i_nc++){
                double value = _scale_term*Ncells_manager_dI[i_nc];
                double value2 = _scale_term*Ncells_manager_dI2[i_nc];
                int idx = i_nc* Npix_to_model + i_pix;
                //d_Ncells_images[idx] = value;
                //d2_Ncells_images[idx] = value2;
                d_image.Ncells[idx] = value;
                d2_image.Ncells[idx] = value2;
            }
        }

        /* update Fcell derivative image */
        if(refine_fcell){
            double value = _scale_term*fcell_manager_dI[1];
            double value2 = _scale_term*fcell_manager_dI2[1];
            d_image.fcell[Npix_to_model+ i_pix] = value;
            d2_image.fcell[Npix_to_model + i_pix] = value2;
            //for (int i_fcell=0; i_fcell < 3; i_fcell++){
            //    double value = _scale_term*fcell_manager_dI[i_fcell];
            //    double value2 = _scale_term*fcell_manager_dI2[i_fcell];
            //    d_fcell_images[i_fcell*Npix_to_model+ i_pix] = value;
            //    d2_fcell_images[i_fcell*Npix_to_model + i_pix] = value2;
            //}
        }/* end Fcell deriv image increment */

        if (refine_fp_fdp){
            // c derivative
            double value = _scale_term*fp_fdp_manager_dI[0];
            d_image.fp_fdp[i_pix] = value;
            // d derivative
            value = _scale_term*fp_fdp_manager_dI[1];
            d_image.fp_fdp[Npix_to_model + i_pix] = value;
        }

        /* update eta derivative image */
        if(refine_eta){
            //double value = _scale_term*eta_manager_dI[0];
            //double value2 = _scale_term*eta_manager_dI2[0];
            //d_eta_images[i_pix] = value;
            //d2_eta_images[i_pix] = value2;
            //if (UMATS_RXYZ_prime.size() >= 3*mosaic_domains){
            for(int i_eta=0; i_eta<3; i_eta++){
                if (i_eta != 0 && UMATS_RXYZ.size() == UMATS_RXYZ_prime.size())
                    continue;
                int idx = i_pix + Npix_to_model*i_eta;
                double value = _scale_term*eta_manager_dI[i_eta];
                double value2 = _scale_term*eta_manager_dI2[i_eta];
                d_image.fp_fdp[idx] = value;
                d2_image.fp_fdp[idx] = value2;
            }
            //}
        }/* end eta deriv image increment */

        /*update the lambda derivative images*/
        for (int i_lam=0 ; i_lam < 2 ; i_lam++){
            if (refine_lambda[i_lam]){
                //double value = scale_term*lambda_managers[i_lam]->dI;
                //double value2 = scale_term*lambda_managers[i_lam]->dI2;
                double value = _scale_term*lambda_manager_dI[i_lam];
                double value2 = _scale_term*lambda_manager_dI2[i_lam];
                int idx = i_lam*Npix_to_model + i_pix;
                d_image.lambda[idx] = value;
                d2_image.lambda[idx] = value2;
            }
        }/* end lambda deriv image increment */

        // sausage increment
        if (refine_sausages){
            for (int i_sausage=0; i_sausage<num_sausages; i_sausage++){
                for (int i=0; i < 4; i++){
                    int sausage_parameter_i = i_sausage*4+i;
                    double value = _scale_term*sausage_manager_dI[sausage_parameter_i];
                    int idx = sausage_parameter_i*Npix_to_model + i_pix;
                    d_image.sausage[idx] = value;
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
                d_image.panel_rot[idx] = value;
                d2_image.panel_rot[idx] = value2;
            }
        }/* end panel rot deriv image increment */

        // panel origin
        for (int i_pan_orig=0; i_pan_orig < 3; i_pan_orig++){
            if(refine_panel_origin[i_pan_orig]){
                double value = _scale_term*pan_orig_manager_dI[i_pan_orig];
                double value2 = _scale_term*pan_orig_manager_dI2[i_pan_orig];
                int idx = i_pan_orig*Npix_to_model + i_pix;
                d_image.panel_orig[idx] = value;
                d2_image.panel_orig[idx] = value2;
            }/* end panel orig deriv image increment */
        }
        if (track_Fhkl){
            for(auto &x: Fhkl_tracker)
                printf("Pixel %d: Fhkl linear index %s came up %d times\n", i_pix, x.first.c_str(), x.second);
        }
    } // end i_pix loop
} // END of CPU kernel

}} // END namespace simtbx::nanoBragg
