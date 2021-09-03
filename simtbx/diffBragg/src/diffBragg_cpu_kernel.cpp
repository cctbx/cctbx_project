#include <simtbx/diffBragg/src/diffBragg.h>
#include <assert.h>
#include <stdbool.h>
#include<unordered_map>
#include<string>

namespace simtbx { namespace nanoBragg { // BEGIN namespace simtbx::nanoBragg

void diffBragg_sum_over_steps(
        int Npix_to_model, std::vector<unsigned int>& panels_fasts_slows,
        image_type& floatimage,
        images& d_image,
        images& d2_image,
        step_arrays& db_steps,
        detector& db_det,
        beam& db_beam,
        crystal& db_cryst,
        flags& db_flags){

    #pragma omp parallel for
    for (int i_pix=0; i_pix < Npix_to_model; i_pix++){
        int pid = panels_fasts_slows[i_pix*3];
        int fpixel = panels_fasts_slows[i_pix*3+1];
        int spixel = panels_fasts_slows[i_pix*3+2];
        double close_distance = db_det.close_distances[pid];
        //std::unordered_map<int, int> Fhkl_tracker;
        std::unordered_map<std::string, int> Fhkl_tracker;
        bool use_nominal_hkl = false;
        if (!db_cryst.nominal_hkl.empty())
            use_nominal_hkl = true;
        // reset photon count for this pixel
        double I=0;
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
        double fcell_manager_dI=0;
        double fcell_manager_dI2=0;
        double eta_manager_dI[3] = {0,0,0};
        double eta_manager_dI2[3] = {0,0,0};
        double lambda_manager_dI[2] = {0,0};
        double lambda_manager_dI2[2] = {0,0};
        double fp_fdp_manager_dI[2] = {0,0};

        for (int i_step=0; i_step < db_steps.Nsteps; i_step++){

            int subS = db_steps.subS_pos[i_step];
            int subF = db_steps.subF_pos[i_step];
            int thick_tic = db_steps.thick_pos[i_step];
            int source = db_steps.source_pos[i_step];
            int phi_tic = db_steps.phi_pos[i_step];
            int mos_tic = db_steps.mos_pos[i_step];

            /* absolute mm position on detector (relative to its origin) */
            double Fdet = db_det.subpixel_size*(fpixel*db_det.oversample + subF ) + db_det.subpixel_size/2.0;
            double Sdet = db_det.subpixel_size*(spixel*db_det.oversample + subS ) + db_det.subpixel_size/2.0;

            /* assume "distance" is to the front of the detector sensor layer */
            double Odet = thick_tic*db_det.detector_thickstep;
            int pid_x = pid*3;
            int pid_y = pid*3+1;
            int pid_z = pid*3+2;
            Eigen::Vector3d o_vec(db_det.odet_vectors[pid_x], db_det.odet_vectors[pid_y], db_det.odet_vectors[pid_z]);

            double pixposX = Fdet*db_det.fdet_vectors[pid_x]+Sdet*db_det.sdet_vectors[pid_x]+Odet*db_det.odet_vectors[pid_x]+db_det.pix0_vectors[pid_x];
            double pixposY = Fdet*db_det.fdet_vectors[pid_y]+Sdet*db_det.sdet_vectors[pid_y]+Odet*db_det.odet_vectors[pid_y]+db_det.pix0_vectors[pid_y];
            double pixposZ = Fdet*db_det.fdet_vectors[pid_z]+Sdet*db_det.sdet_vectors[pid_z]+Odet*db_det.odet_vectors[pid_z]+db_det.pix0_vectors[pid_z];
            Eigen::Vector3d pixel_pos(pixposX, pixposY, pixposZ);

            double airpath = pixel_pos.norm();
            Eigen::Vector3d diffracted = pixel_pos/airpath;

            double omega_pixel = db_det.pixel_size*db_det.pixel_size/airpath/airpath*close_distance/airpath;

            if(db_flags.point_pixel) omega_pixel = 1.0/airpath/airpath;

            double capture_fraction = 1;
            if(db_det.detector_thick > 0.0 && db_det.detector_attnlen > 0.0)
            {
                double parallax = diffracted.dot(o_vec) ;
                capture_fraction = exp(-thick_tic*db_det.detector_thickstep/db_det.detector_attnlen/parallax)
                                  -exp(-(thick_tic+1)*db_det.detector_thickstep/db_det.detector_attnlen/parallax);
            }
            Eigen::Vector3d incident(-db_beam.source_X[source], -db_beam.source_Y[source], -db_beam.source_Z[source]);
            double lambda = db_beam.source_lambda[source];
            double lambda_ang = lambda*1e10;
            if (db_flags.use_lambda_coefficients){
                lambda_ang = db_beam.lambda0 + db_beam.lambda1*lambda_ang;
                lambda = lambda_ang*1e-10;
            }

            /* construct the incident beam unit vector while recovering source distance */
            double source_path = incident.norm();
            incident /= source_path;

            /* construct the scattering vector for this pixel */
            //vec3 scattering = (diffracted - incident) / lambda;
            Eigen::Vector3d scattering = (diffracted - incident) / lambda;

            /* sin(theta)/lambda is half the scattering vector length */
            double stol = 0.5*(scattering.norm()); //magnitude(scattering);

            /* rough cut to speed things up when we aren't using whole detector */
            if(db_cryst.dmin > 0.0 && stol > 0.0)
            {
                if(db_cryst.dmin > 0.5/stol)
                {
                    continue;
                }
            }

            double phi = db_cryst.phi0 + db_cryst.phistep*phi_tic;
            Eigen::Matrix3d Bmat_realspace = db_cryst.eig_B;
            if( phi != 0.0 )
            {
                double cosphi = cos(phi);
                double sinphi = sin(phi);
                Eigen::Vector3d ap_vec(db_cryst.eig_B(0,0), db_cryst.eig_B(1,0), db_cryst.eig_B(2,0));
                Eigen::Vector3d bp_vec(db_cryst.eig_B(0,1), db_cryst.eig_B(1,1), db_cryst.eig_B(2,1));
                Eigen::Vector3d cp_vec(db_cryst.eig_B(0,2), db_cryst.eig_B(1,2), db_cryst.eig_B(2,2));

                ap_vec = ap_vec*cosphi + db_cryst.spindle_vec.cross(ap_vec)*sinphi + db_cryst.spindle_vec*(db_cryst.spindle_vec.dot(ap_vec))*(1-cosphi);
                bp_vec = bp_vec*cosphi + db_cryst.spindle_vec.cross(bp_vec)*sinphi + db_cryst.spindle_vec*(db_cryst.spindle_vec.dot(bp_vec))*(1-cosphi);
                cp_vec = cp_vec*cosphi + db_cryst.spindle_vec.cross(cp_vec)*sinphi + db_cryst.spindle_vec*(db_cryst.spindle_vec.dot(cp_vec))*(1-cosphi);

                Bmat_realspace << ap_vec[0], bp_vec[0], cp_vec[0],
                                    ap_vec[1], bp_vec[1], cp_vec[1],
                                    ap_vec[2], bp_vec[2], cp_vec[2];
            }
            Bmat_realspace *= 1e10;

            Eigen::Matrix3d U = db_cryst.eig_U;
            Eigen::Matrix3d UBO = (db_cryst.UMATS_RXYZ[mos_tic] * U*Bmat_realspace*(db_cryst.eig_O.transpose())).transpose();

            Eigen::Vector3d q_vec(scattering[0], scattering[1], scattering[2]);
            q_vec *= 1e-10;
            Eigen::Vector3d H_vec = UBO*q_vec;

            double h = H_vec[0];
            double k = H_vec[1];
            double l = H_vec[2];

            /* round off to nearest whole index */
            int h0 = static_cast<int>(ceil(h - 0.5));
            int k0 = static_cast<int>(ceil(k - 0.5));
            int l0 = static_cast<int>(ceil(l - 0.5));

            Eigen::Vector3d H0(h0, k0, l0);
            Eigen::Matrix3d NABC;
            NABC << db_cryst.Na, db_cryst.Nd, db_cryst.Nf,
                     db_cryst.Nd, db_cryst.Nb, db_cryst.Ne,
                     db_cryst.Nf, db_cryst.Ne, db_cryst.Nc;

            double C = 2 / 0.63 * db_cryst.fudge;
            Eigen::Vector3d delta_H = H_vec - H0;
            Eigen::Vector3d V = NABC*delta_H;
            double hrad_sqr = V.dot(V);
            double F_latt;
            if (db_flags.no_Nabc_scale)
                F_latt = exp(-( hrad_sqr / 0.63 * db_cryst.fudge ));
            else
                F_latt = db_cryst.Na*db_cryst.Nb*db_cryst.Nc*exp(-( hrad_sqr / 0.63 * db_cryst.fudge ));

            if (db_flags.use_diffuse){
                Eigen::Vector3d delta_q_vec = UBO.inverse()*delta_H;
                Eigen::Vector3d bragg_q_vec = UBO.inverse()*H0;

                Eigen::Matrix3d anisoU; // TODO make this matrix outside the loop
                anisoU << db_cryst.this_sigma,0,0,
                          0,db_cryst.this_sigma,0,
                          0,0,db_cryst.this_sigma;
                double exparg = 4*M_PI*M_PI*bragg_q_vec.dot(anisoU*bragg_q_vec);
                double dwf = exp(-exparg);

                //F_latt_diffuse = 8.*M_PI*db_cryst.this_gamma*db_cryst.this_gamma*db_cryst.this_gamma/ pow((1.+db_cryst.this_gamma*db_cryst.this_gamma*2.*M_PI*2.*M_PI* rsqr),2.);
                Eigen::Matrix3d anisoG; anisoG << db_cryst.this_gamma,0,0, // TODO make this matrix outside loop
                                          0,db_cryst.this_gamma,0,
                                          0,0,db_cryst.this_gamma;
                Eigen::Vector3d anisoG_q = anisoG*delta_q_vec;
                double F_latt_diffuse = 4.*M_PI*anisoG.determinant() /
                        (1.+ anisoG_q.dot(anisoG_q)* 4*M_PI*M_PI);
                F_latt_diffuse *= (dwf*exparg);
                F_latt += F_latt_diffuse;
            }

            /* no need to go further if result will be zero */
            //if(F_latt == 0.0 && ! only_save_omega_kahn) {
            //    continue;
            //}
            /* convert amplitudes into intensity (photons per steradian) */
            if (!db_flags.oversample_omega)
                omega_pixel = 1;

            /* increment to intensity */
            double I_noFcell = F_latt*F_latt*db_beam.source_I[source]*capture_fraction*omega_pixel;

            /* structure factor of the unit cell */
            double F_cell = db_cryst.default_F;
            double F_cell2 = 0;

            if ( (h0<=db_cryst.h_max) && (h0>=db_cryst.h_min) && (k0<=db_cryst.k_max) && (k0>=db_cryst.k_min) && (l0<=db_cryst.l_max) && (l0>=db_cryst.l_min)  ) {
                /* just take nearest-neighbor */
                int Fhkl_linear_index = (h0-db_cryst.h_min) * db_cryst.k_range * db_cryst.l_range +
                                (k0- db_cryst.k_min) * db_cryst.l_range +
                                (l0-db_cryst.l_min);
                F_cell = db_cryst.FhklLinear[Fhkl_linear_index];
                if (db_flags.track_Fhkl){
                    std::string hkl_s ;
                    hkl_s = std::to_string(h0) + ","+ std::to_string(k0) + "," + std::to_string(l0);
                    if (Fhkl_tracker.count(hkl_s))
                        Fhkl_tracker[hkl_s] += 1;
                    else
                        Fhkl_tracker[hkl_s] = 0;
                    continue;
                }
                if (db_flags.complex_miller) F_cell2 = db_cryst.Fhkl2Linear[Fhkl_linear_index];
            }
            //else{
            // F_cell = default_F;
            //}
            bool do_max = false;
            if(db_flags.verbose > 3){//} && i_step==0){
                double F2 = sqrt(F_cell*F_cell + F_cell2*F_cell2);
                double II = I_noFcell*F2*F2;
                if (I_noFcell > II_max){
                    do_max = true;
                    max_stats[0] = F_cell;
                    max_stats[1] = F_cell2;
                    max_stats[2] = I_noFcell;
                    max_stats[3] = F_latt;
                    max_stats[4] = II;
                    max_stats[5] = omega_pixel;
                    max_stats[6] = db_beam.source_I[source];
                    max_stats[7]= db_beam.source_lambda[source];
                    max_stats[8] = h;
                    max_stats[9] = k;
                    max_stats[10] = l;
                    II_max = I_noFcell;
                }
                else{
                do_max=false;}
                //printf("hkl= %f %f %f  hkl1= %d %d %d  |Fcell| before=%f, Ino=%f, I=%f, Flatt=%f, cap_frac=%f, omega_pix=%f, sourceI=%f\n", h,k,l,h0,k0,l0, F2, I_noFcell, II, F_latt,
                //capture_fraction, omega_pixel, source_I[source] );
            }
            double c_deriv_Fcell_real = 0;
            double c_deriv_Fcell_imag = 0;
            double d_deriv_Fcell_real = 0;
            double d_deriv_Fcell_imag = 0;
            double c_deriv_Fcell = 0;
            double d_deriv_Fcell = 0;
            if (db_flags.complex_miller){
            // TODO shouldnt this be constant for each HKL?
              if (db_cryst.fpfdp.size() > 0){
                   double S_2 = 1.e-20*(scattering[0]*scattering[0]+scattering[1]*scattering[1]+scattering[2]*scattering[2]);

                    // fp is always followed by the fdp value
                   double val_fp = db_cryst.fpfdp[2*source];
                   double val_fdp = db_cryst.fpfdp[2*source+1];

                   double c_deriv_prime=0;
                   double c_deriv_dblprime=0;
                   double d_deriv_prime = 0;
                   double d_deriv_dblprime = 0;
                   if (db_flags.refine_fp_fdp){
                   //   currently only supports two parameter model
                       int nsources_times_two = db_cryst.fpfdp.size();
                       int d_idx = 2*source;
                       c_deriv_prime = db_cryst.fpfdp_derivs[d_idx];
                       c_deriv_dblprime = db_cryst.fpfdp_derivs[d_idx+1];
                       d_deriv_prime = db_cryst.fpfdp_derivs[d_idx+nsources_times_two];
                       d_deriv_dblprime = db_cryst.fpfdp_derivs[d_idx+1+nsources_times_two];
                   }
                   // 5 valeus per atom: x,y,z,B,occupancy
                   int num_atoms = db_cryst.atom_data.size()/5;
                   for (int  i_atom=0; i_atom < num_atoms; i_atom++){
                       if (db_flags.verbose>5)
                         printf("Processing atom %d, F_cell=%10.3f\n", i_atom, F_cell);
                       // fractional atomic coordinates
                       double atom_x = db_cryst.atom_data[i_atom*5];
                       double atom_y = db_cryst.atom_data[i_atom*5+1];
                       double atom_z = db_cryst.atom_data[i_atom*5+2];
                       double B = db_cryst.atom_data[i_atom*5+3]; // B factor
                       B = exp(-B*S_2/4.0); // TODO: speed me up?
                       double occ = db_cryst.atom_data[i_atom*5+4]; // occupancy
                       double r_dot_h = h0*atom_x + k0*atom_y + l0*atom_z;
                       double phase = 2*M_PI*r_dot_h;
                       double s_rdoth = sin(phase);
                       double c_rdoth = cos(phase);
                       double Bocc = B*occ;
                       double BC = Bocc*c_rdoth;
                       double BS = Bocc*s_rdoth;
                       double real_part = BC*val_fp - BS*val_fdp;
                       double imag_part = BS*val_fp + BC*val_fdp;
                       F_cell += real_part;
                       F_cell2 += imag_part;
                       if (db_flags.refine_fp_fdp){
                            c_deriv_Fcell_real += BC*c_deriv_prime - BS*c_deriv_dblprime;
                            c_deriv_Fcell_imag += BS*c_deriv_prime + BC*c_deriv_dblprime;

                            d_deriv_Fcell_real += BC*d_deriv_prime - BS*d_deriv_dblprime;
                            d_deriv_Fcell_imag += BS*d_deriv_prime + BC*d_deriv_dblprime;
                       }
                   }
               }
               double Freal = F_cell;
               double Fimag = F_cell2;
               F_cell = sqrt(pow(Freal,2) + pow(Fimag,2));// TODO work around the sqrt ?
               //if (verbose> 3 && i_step==0){
               //     double II = I_noFcell*F_cell*F_cell;
               //     printf("hkl= %f %f %f  hkl1= %d %d %d  |Fcell| after=%f, Ino=%f, I=%f\n", h,k,l,h0,k0,l0, F_cell, I_noFcell, II);
               // }
               if (db_flags.refine_fp_fdp){
                   c_deriv_Fcell = Freal*c_deriv_Fcell_real + Fimag*c_deriv_Fcell_imag;
                   d_deriv_Fcell = Freal*d_deriv_Fcell_real + Fimag*d_deriv_Fcell_imag;
               }
            }

            if(db_flags.verbose > 3){
                double F2 = sqrt(F_cell*F_cell + F_cell2*F_cell2);
                double II = I_noFcell*F2*F2;
                if (do_max){
                    max_stats2[0] = F2 ;
                    max_stats2[1] = 0 ; //F_cell2;
                    max_stats2[2] = I_noFcell;
                    max_stats2[3] = F_latt;
                    max_stats2[4] = II;
                    max_stats2[5] = omega_pixel;
                    max_stats2[6] = db_beam.source_I[source];
                    max_stats2[7]= db_beam.source_lambda[source];
                    max_stats2[8] = h;
                    max_stats2[9] = k;
                    max_stats2[10] = l;
                }
            }
            double Iincrement = F_cell*F_cell*I_noFcell;
            I += Iincrement;

            // TODO if test?
            if (db_flags.refine_fp_fdp){
                fp_fdp_manager_dI[0] += 2*I_noFcell * (c_deriv_Fcell);
                fp_fdp_manager_dI[1] += 2*I_noFcell * (d_deriv_Fcell);
            }

            //if(verbose > 3)
            //    printf("hkl= %f %f %f  hkl1= %d %d %d  Fcell=%f\n", h,k,l,h0,k0,l0, F_cell);

            double two_C = 2*C;
            Eigen::Matrix3d UBOt = U*Bmat_realspace*(db_cryst.eig_O.transpose());
            if (db_flags.refine_Umat[0]){
                Eigen::Matrix3d RyRzUBOt = db_cryst.RotMats[1]*db_cryst.RotMats[2]*UBOt;
                Eigen::Vector3d delta_H_prime = (db_cryst.UMATS[mos_tic]*db_cryst.dRotMats[0]*RyRzUBOt).transpose()*q_vec;
                double V_dot_dV = V.dot(NABC*delta_H_prime);
                double value = -two_C * V_dot_dV * Iincrement;
                double value2 =0;
                if (db_flags.compute_curvatures) {
                    Eigen::Vector3d delta_H_dbl_prime = (db_cryst.UMATS[mos_tic]*db_cryst.d2RotMats[0]*RyRzUBOt).transpose()*q_vec;
                    double dV_dot_dV = (NABC*delta_H_prime).dot(NABC*delta_H_prime);
                    double dV2_dot_V = (NABC*delta_H).dot(NABC*delta_H_dbl_prime);
                    value2 = two_C*(two_C*V_dot_dV*V_dot_dV - dV2_dot_V - dV_dot_dV)*Iincrement;
                }
                rot_manager_dI[0] += value;
                rot_manager_dI2[0] += value2;
            }
            if (db_flags.refine_Umat[1]){
                Eigen::Matrix3d UmosRx = db_cryst.UMATS[mos_tic]*db_cryst.RotMats[0];
                Eigen::Matrix3d RzUBOt = db_cryst.RotMats[2]*UBOt;
                Eigen::Vector3d delta_H_prime =(UmosRx*db_cryst.dRotMats[1]*RzUBOt).transpose()*q_vec;
                double V_dot_dV = V.dot(NABC*delta_H_prime);
                double value = -two_C * V_dot_dV * Iincrement;

                double value2=0;
                if (db_flags.compute_curvatures){
                    Eigen::Vector3d delta_H_dbl_prime = (UmosRx*db_cryst.d2RotMats[1]*RzUBOt).transpose()*q_vec;
                    double dV_dot_dV = (NABC*delta_H_prime).dot(NABC*delta_H_prime);
                    double dV2_dot_V = (NABC*delta_H).dot(NABC*delta_H_dbl_prime);
                    value2 = two_C*(two_C*V_dot_dV*V_dot_dV - dV2_dot_V - dV_dot_dV)*Iincrement;
                }
                rot_manager_dI[1] += value;
                rot_manager_dI2[1] += value2;
            }
            if (db_flags.refine_Umat[2]){
                Eigen::Matrix3d UmosRxRy = db_cryst.UMATS[mos_tic]*db_cryst.RotMats[0]*db_cryst.RotMats[1];
                Eigen::Vector3d delta_H_prime = (UmosRxRy*db_cryst.dRotMats[2]*UBOt).transpose()*q_vec;
                double V_dot_dV = V.dot(NABC*delta_H_prime);
                double value = -two_C * V_dot_dV * Iincrement;

                double value2=0;
                if (db_flags.compute_curvatures){
                    Eigen::Vector3d delta_H_dbl_prime = (UmosRxRy*db_cryst.d2RotMats[2]*UBOt).transpose()*q_vec;
                    double dV_dot_dV = (NABC*delta_H_prime).dot(NABC*delta_H_prime);
                    double dV2_dot_V = (NABC*delta_H).dot(NABC*delta_H_dbl_prime);
                    value2 = two_C*(two_C*V_dot_dV*V_dot_dV - dV2_dot_V - dV_dot_dV)*Iincrement;
                }
                rot_manager_dI[2] += value;
                rot_manager_dI2[2] += value2;
            }
            /*Checkpoint for unit cell derivatives*/
            Eigen::Matrix3d Ot = db_cryst.eig_O.transpose();
            Eigen::Matrix3d UmosRxRyRzU ;
            Eigen::Vector3d delta_H_prime;
            for(int i_uc=0; i_uc < 6; i_uc++ ){
                if (db_flags.refine_Bmat[i_uc]){
                    UmosRxRyRzU = db_cryst.UMATS_RXYZ[mos_tic]*U;
                    delta_H_prime = ((UmosRxRyRzU*db_cryst.dB_Mats[i_uc]*Ot).transpose()*q_vec);
                    double V_dot_dV = V.dot(NABC*delta_H_prime);
                    double value = -two_C * V_dot_dV * Iincrement;
                    double value2 =0;
                    if (db_flags.compute_curvatures){
                        Eigen::Vector3d delta_H_dbl_prime = ((UmosRxRyRzU*db_cryst.dB2_Mats[i_uc]*Ot).transpose()*q_vec);
                        double dV_dot_dV = (NABC*delta_H_prime).dot(NABC*delta_H_prime);
                        double dV2_dot_V = (NABC*delta_H).dot(NABC*delta_H_dbl_prime);
                        value2 = two_C*(two_C*V_dot_dV*V_dot_dV - dV2_dot_V - dV_dot_dV)*Iincrement;
                    }
                    ucell_manager_dI[i_uc] += value;
                    ucell_manager_dI2[i_uc] += value2;
                }
            } /*end ucell deriv */

            /* Checkpoint for Ncells manager */
            if (db_flags.refine_Ncells[0]){
                int num_ncell_deriv = 1;
                if (! db_flags.isotropic_ncells)
                    num_ncell_deriv = 3;
                for (int i_nc=0; i_nc < num_ncell_deriv; i_nc++) {
                    Eigen::Matrix3d dN;
                    if (!db_flags.isotropic_ncells){
                        dN << 0,0,0,0,0,0,0,0,0;
                        dN(i_nc, i_nc) = 1;
                        }
                    else
                        dN << 1,0,0,0,1,0,0,0,1;

                    double N_i = NABC(i_nc, i_nc);
                    Eigen::Vector3d dV_dN = dN*delta_H;
                    //if (i_step==0 && fpixel ==0 && spixel == 0){
                    //    printf("dN matrix: %f %f %f\n %f %f %f\n %f %f %f\n",
                    //        dN(0,0), dN(0,1), dN(0,2),
                    //        dN(1,0), dN(1,1), dN(1,2),
                    //        dN(2,0), dN(2,1), dN(2,2)
                    //        );
                    //    printf("NABC matrix: %f %f %f\n %f %f %f\n %f %f %f\n",
                    //        NABC(0,0), NABC(0,1), NABC(0,2),
                    //        NABC(1,0), NABC(1,1), NABC(1,2),
                    //        NABC(2,0), NABC(2,1), NABC(2,2)
                    //        );
                    //}
                    double deriv_coef;
                    if (db_flags.isotropic_ncells)
                        deriv_coef= 3/N_i - C* ( dV_dN.dot(V));
                    else
                        deriv_coef= 1/N_i - C* ( dV_dN.dot(V));
                    double value = 2*Iincrement*deriv_coef;
                    double value2=0;
                    if(db_flags.compute_curvatures){
                        dN(i_nc, i_nc) = 0; // TODO check maths
                        value2 = ( -1/N_i/N_i - C*(dV_dN.dot(dV_dN))) *2*Iincrement;
                        value2 += deriv_coef*2*value;
                    }
                    Ncells_manager_dI[i_nc] += value;
                    Ncells_manager_dI2[i_nc] += value2;
                }

            } /* end Ncells manager deriv */

            if (db_flags.refine_Ncells_def){
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
                    if (db_flags.compute_curvatures){
                        value2 = deriv_coef*value;
                        value2 +=  -2*C*Iincrement*(dV_dN.dot(dV_dN));
                    Ncells_manager_dI2[i_nc] += value2;
                    }
                }
            }

            ///* Checkpoint for Origin manager */
            for (int i_pan_orig=0; i_pan_orig < 3; i_pan_orig++){
                if (db_flags.refine_panel_origin[i_pan_orig]){
                    double per_k = 1/airpath;
                    double per_k3 = pow(per_k,3);
                    double per_k5 = pow(per_k,5);
                    double lambda_ang = lambda*1e10;

                    Eigen::Matrix3d M = -two_C*(NABC*UBO)/lambda_ang;
                    Eigen::Vector3d dk;
                    if (i_pan_orig == 0)
                        dk << 0,0,1;
                    else if (i_pan_orig == 1)
                        dk << 1,0,0;
                    else
                        dk << 0,1,0;
                    af::flex_double dI_and_dI2 = get_panel_increment(Iincrement, omega_pixel, M, db_det.subpixel_size*db_det.subpixel_size,
                        o_vec, pixel_pos, per_k,  per_k3, per_k5, V, dk);
                    pan_orig_manager_dI[i_pan_orig] += dI_and_dI2[0];
                    pan_orig_manager_dI2[i_pan_orig] += dI_and_dI2[1];

                } /* end origin manager deriv */
            }

            for (int i_pan_rot=0; i_pan_rot < 3; i_pan_rot++){
                if(db_flags.refine_panel_rot[i_pan_rot]){
                    double per_k = 1/airpath;
                    double per_k3 = pow(per_k,3);
                    double per_k5 = pow(per_k,5);
                    double lambda_ang = lambda*1e10;
                    Eigen::Matrix3d M = -two_C*(NABC*UBO)/lambda_ang;

                    Eigen::Vector3d dk = Fdet*(db_det.dF_vecs[pid*3 + i_pan_rot]) + Sdet*(db_det.dS_vecs[pid*3 + i_pan_rot]);
                    af::flex_double dI_and_dI2 = get_panel_increment(Iincrement, omega_pixel, M, db_det.subpixel_size*db_det.subpixel_size,
                        o_vec, pixel_pos, per_k,  per_k3, per_k5, V, dk);
                    pan_rot_manager_dI[i_pan_rot] += dI_and_dI2[0];
                    pan_rot_manager_dI2[i_pan_rot] += dI_and_dI2[1];
                }
            }

            /* checkpoint for Fcell manager */
            if (db_flags.refine_fcell){
            //    TODO rewrite so no divide by Fcell
                int nom_h=h0;
                int nom_k=k0;
                int nom_l=l0;
                //int f_cell_idx = 1;
                if (use_nominal_hkl){
                    nom_h = db_cryst.nominal_hkl[i_pix*3];
                    nom_k = db_cryst.nominal_hkl[i_pix*3+1];
                    nom_l = db_cryst.nominal_hkl[i_pix*3+2];
                    //f_cell_idx = l0 - nom_l + 1;
                }
                double value = 2*I_noFcell*F_cell; //2*Iincrement/F_cell ;
                double value2=0;
                if (db_flags.compute_curvatures){
                    if (F_cell > 0)
                        value2 =2*I_noFcell;
                }
                //if (f_cell_idx >= 0 && f_cell_idx <= 2){
                // NOTE use nominal hkl to regulate when gradients are compputed, as h,k,l can drift within a shoebox
                if (use_nominal_hkl){
                    if (nom_h==h0 && nom_k==k0 && nom_l==l0 ){
                        fcell_manager_dI += value;
                        fcell_manager_dI2 += value2;
                    }
                }
                else{
                    fcell_manager_dI += value;
                    fcell_manager_dI2 += value2;
                }
            } /* end of fcell man deriv */

            /* checkpoint for eta manager */
            if (db_flags.refine_eta){
                for (int i_eta=0; i_eta < 3; i_eta++){
                    if (i_eta != 0 && db_cryst.UMATS_RXYZ_prime.size()==db_cryst.UMATS_RXYZ.size()){
                        continue;
                    }
                    int mtic2 = mos_tic  + i_eta*db_cryst.mosaic_domains;
                    Eigen::Vector3d DeltaH_deriv = (db_cryst.UMATS_RXYZ_prime[mtic2]*UBOt).transpose()*q_vec;
                    // vector V is Nabc*Delta_H
                    Eigen::Vector3d dV = NABC*DeltaH_deriv;
                    double V_dot_dV = V.dot(dV);
                    double Iprime = -two_C*(V_dot_dV)*Iincrement;
                    eta_manager_dI[i_eta] += Iprime;

                    double Idbl_prime = 0;
                    if (db_flags.compute_curvatures){
                        Eigen::Vector3d DeltaH_second_deriv = (db_cryst.UMATS_RXYZ_dbl_prime[mtic2]*UBOt).transpose()*q_vec;
                        Eigen::Vector3d dV2 = NABC*DeltaH_second_deriv;
                        Idbl_prime = -two_C*(dV.dot(dV) + V.dot(dV2))*Iincrement;
                        Idbl_prime += -two_C*(V_dot_dV)*Iprime;
                    }
                    eta_manager_dI2[i_eta] += Idbl_prime;
                }
            } /* end of eta man deriv */

            /*checkpoint for lambda manager*/
            for(int i_lam=0; i_lam < 2; i_lam++){
                if (db_flags.refine_lambda[i_lam]){
                    double lambda_ang = lambda*1e10;
                    double NH_dot_V = (NABC*H_vec).dot(V);
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

            if( db_flags.printout && i_step==0 ){
                if((fpixel==db_flags.printout_fpixel && spixel==db_flags.printout_spixel) || db_flags.printout_fpixel < 0)
                {
                    printf("%4d %4d : stol = %g, lambda = %g\n", fpixel,spixel,stol, lambda);
                    printf("at %g %g %g\n", pixel_pos[0],pixel_pos[1],pixel_pos[2]);
                    printf("hkl= %f %f %f  hkl0= %d %d %d\n", h,k,l,h0,k0,l0);
                    printf(" F_cell=%g  F_cell_2=%g F_latt=%g   I = %g\n", F_cell,F_cell2,F_latt,I);
                    printf("I/steps %15.10g\n", I/db_steps.Nsteps);
                    printf("cap frac   %f\n", capture_fraction);
                    printf("Fdet= %g; Sdet= %g ; Odet= %g\n", Fdet, Sdet, Odet);
                    printf("omega   %15.10g\n", omega_pixel);
                    printf("close_distance    %15.10g\n", close_distance);
                    printf("fluence    %15.10g\n", db_beam.fluence);
                    printf("spot scale   %15.10g\n", db_cryst.spot_scale);
                    printf("sourceI[0]:   %15.10g\n", db_beam.source_I[0]);
                    printf("default_F= %f\n", db_cryst.default_F);
                    printf("PIX0: %f %f %f\n" , db_det.pix0_vectors[pid_x], db_det.pix0_vectors[pid_y], db_det.pix0_vectors[pid_z]);
                    printf("F: %f %f %f\n" , db_det.fdet_vectors[pid_x], db_det.fdet_vectors[pid_y], db_det.fdet_vectors[pid_z]);
                    printf("S: %f %f %f\n" , db_det.sdet_vectors[pid_x], db_det.sdet_vectors[pid_y], db_det.sdet_vectors[pid_z]);
                    printf("O: %f %f %f\n" , db_det.odet_vectors[pid_x], db_det.odet_vectors[pid_y], db_det.odet_vectors[pid_z]);
                    printf("QVECTOR: %f %f %f\n" , q_vec[0], q_vec[1], q_vec[2]);
                    printf("MOSAIC UMAT RXYZ\n");
                    std::cout << db_cryst.UMATS_RXYZ[mos_tic] << std::endl;
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
                    SCITBX_EXAMINE(db_det.fdet_vectors[0]);
                    SCITBX_EXAMINE(db_det.fdet_vectors[1]);
                    SCITBX_EXAMINE(db_det.fdet_vectors[2]);
                    SCITBX_EXAMINE(db_det.pix0_vectors[0]);
                    SCITBX_EXAMINE(db_det.pix0_vectors[1]);
                    SCITBX_EXAMINE(db_det.pix0_vectors[2]);
                }
            }
        } /* end of i_steps loop */

        if (db_flags.verbose >3){
            printf("hkl= %f %f %f  |Fcell| before=(%f,%f), Ino=%f, I=%f, Flatt=%f, omega_pix=%2.3e, sourceI=%f, sourceLambda=%f\n",
             max_stats[8],max_stats[9],max_stats[10],
             max_stats[0], max_stats[1], max_stats[2], max_stats[4], max_stats[3],
            max_stats[5], max_stats[6], max_stats[7]*1e10 );

            printf("hkl= %f %f %f  |Fcell| after=(%f,%f), Ino=%f, I=%f, Flatt=%f, omega_pix=%2.3e, sourceI=%f, sourceLambda=%f\n",
             max_stats2[8],max_stats2[9],max_stats2[10],
             max_stats2[0], max_stats2[1], max_stats2[2], max_stats2[4], max_stats2[3],
            max_stats2[5], max_stats2[6], max_stats2[7]*1e10 );

        }

        double scale_term=1;
        if (!db_flags.track_Fhkl){

            double Fdet_ave = db_det.pixel_size*fpixel + db_det.pixel_size/2.0;
            double Sdet_ave = db_det.pixel_size*spixel + db_det.pixel_size/2.0;
            double Odet_ave = 0; //Odet; // TODO maybe make this more general for thick detectors?

            Eigen::Vector3d pixel_pos_ave(0,0,0);
            int pid_x = pid*3;
            int pid_y = pid*3+1;
            int pid_z = pid*3+2;
            pixel_pos_ave[0] = Fdet_ave * db_det.fdet_vectors[pid_x]+Sdet_ave*db_det.sdet_vectors[pid_x]+Odet_ave*db_det.odet_vectors[pid_x]+db_det.pix0_vectors[pid_x];
            pixel_pos_ave[1] = Fdet_ave * db_det.fdet_vectors[pid_y]+Sdet_ave*db_det.sdet_vectors[pid_y]+Odet_ave*db_det.odet_vectors[pid_y]+db_det.pix0_vectors[pid_y];
            pixel_pos_ave[2] = Fdet_ave * db_det.fdet_vectors[pid_z]+Sdet_ave*db_det.sdet_vectors[pid_z]+Odet_ave*db_det.odet_vectors[pid_z]+db_det.pix0_vectors[pid_z];

            double airpath_ave = pixel_pos_ave.norm();
            Eigen::Vector3d diffracted_ave = pixel_pos_ave/airpath_ave;
            double omega_pixel_ave = db_det.pixel_size*db_det.pixel_size/airpath_ave/airpath_ave*close_distance/airpath_ave;

            double polar = 1;
            if (!db_flags.nopolar){
                Eigen::Vector3d incident(-db_beam.source_X[0], -db_beam.source_Y[0], -db_beam.source_Z[0]);
                incident = incident / incident.norm();
                // component of diffracted unit vector along incident beam unit vector
                double cos2theta = incident.dot(diffracted_ave);
                double cos2theta_sqr = cos2theta*cos2theta;
                double sin2theta_sqr = 1-cos2theta_sqr;

                double psi=0;
                if(db_beam.kahn_factor != 0.0){
                    // cross product to get "vertical" axis that is orthogonal to the cannonical "polarization"
                    Eigen::Vector3d B_in = db_beam.polarization_axis.cross(incident);
                    // cross product with incident beam to get E-vector direction
                    Eigen::Vector3d E_in = incident.cross(B_in);
                    // get components of diffracted ray projected onto the E-B plane
                    double kEi = diffracted_ave.dot(E_in);
                    double kBi = diffracted_ave.dot(B_in);
                    // compute the angle of the diffracted ray projected onto the incident E-B plane
                    psi = -atan2(kBi,kEi);
                }
                // correction for polarized incident beam
                polar = 0.5*(1.0 + cos2theta_sqr - db_beam.kahn_factor*cos(2*psi)*sin2theta_sqr);
            }

            double om = 1;
            if (!db_flags.oversample_omega)
                om=omega_pixel_ave;

            // final scale term to being everything to photon number units
            scale_term = db_cryst.r_e_sqr*db_beam.fluence*db_cryst.spot_scale*polar*om/db_steps.Nsteps;

            //int roi_i = i_pix; // TODO replace roi_i with i_pix

            floatimage[i_pix] = scale_term*I;

        }
        /* udpate the rotation derivative images*/
        for (int i_rot =0 ; i_rot < 3 ; i_rot++){
            if (db_flags.refine_Umat[i_rot]){
                double value = scale_term*rot_manager_dI[i_rot];
                double value2 = scale_term*rot_manager_dI2[i_rot];
                int idx = i_rot*Npix_to_model + i_pix;
                d_image.Umat[idx] = value;
                d2_image.Umat[idx] = value2;
            }
        } /* end rot deriv image increment */

        /*update the ucell derivative images*/
        for (int i_uc=0 ; i_uc < 6 ; i_uc++){
            if (db_flags.refine_Bmat[i_uc]){
                double value = scale_term*ucell_manager_dI[i_uc];
                double value2 = scale_term*ucell_manager_dI2[i_uc];
                int idx= i_uc*Npix_to_model + i_pix;
                d_image.Bmat[idx] = value;
                d2_image.Bmat[idx] = value2;
            }
        }/* end ucell deriv image increment */

        /*update the Ncells derivative image*/
        if (db_flags.refine_Ncells[0]){
            double value = scale_term*Ncells_manager_dI[0];
            double value2 = scale_term*Ncells_manager_dI2[0];
            double idx = i_pix;
            d_image.Ncells[idx] = value;
            d2_image.Ncells[idx] = value2;
            //d_Ncells_images[idx] = value;
            //d2_Ncells_images[idx] = value2;

            if (! db_flags.isotropic_ncells){
                value = scale_term*Ncells_manager_dI[1];
                value2 = scale_term*Ncells_manager_dI2[1];
                idx = Npix_to_model + i_pix;
                d_image.Ncells[idx] = value;
                d2_image.Ncells[idx] = value2;
                //d_Ncells_images[idx] = value;
                //d2_Ncells_images[idx] = value2;

                value = scale_term*Ncells_manager_dI[2];
                value2 = scale_term*Ncells_manager_dI2[2];
                idx = Npix_to_model*2 + i_pix;
                d_image.Ncells[idx] = value;
                d2_image.Ncells[idx] = value2;
                //d_Ncells_images[idx] = value;
                //d2_Ncells_images[idx] = value2;
            }
        }/* end Ncells deriv image increment */
        if (db_flags.refine_Ncells_def){
            for (int i_nc=3; i_nc<6; i_nc++){
                double value = scale_term*Ncells_manager_dI[i_nc];
                double value2 = scale_term*Ncells_manager_dI2[i_nc];
                int idx = i_nc* Npix_to_model + i_pix;
                //d_Ncells_images[idx] = value;
                //d2_Ncells_images[idx] = value2;
                d_image.Ncells[idx] = value;
                d2_image.Ncells[idx] = value2;
            }
        }

        /* update Fcell derivative image */
        if(db_flags.refine_fcell){
            double value = scale_term*fcell_manager_dI;
            double value2 = scale_term*fcell_manager_dI2;
            d_image.fcell[i_pix] = value;
            d2_image.fcell[i_pix] = value2;
        }/* end Fcell deriv image increment */

        if (db_flags.refine_fp_fdp){
            // c derivative
            double value = scale_term*fp_fdp_manager_dI[0];
            d_image.fp_fdp[i_pix] = value;
            // d derivative
            value = scale_term*fp_fdp_manager_dI[1];
            d_image.fp_fdp[Npix_to_model + i_pix] = value;
        }

        /* update eta derivative image */
        if(db_flags.refine_eta){
            //double value = scale_term*eta_manager_dI[0];
            //double value2 = scale_term*eta_manager_dI2[0];
            //d_eta_images[i_pix] = value;
            //d2_eta_images[i_pix] = value2;
            //if (db_cryst.UMATS_RXYZ_prime.size() >= 3*mosaic_domains){
            for(int i_eta=0; i_eta<3; i_eta++){
                if (i_eta != 0 && db_cryst.UMATS_RXYZ.size() == db_cryst.UMATS_RXYZ_prime.size())
                    continue;
                int idx = i_pix + Npix_to_model*i_eta;
                double value = scale_term*eta_manager_dI[i_eta];
                double value2 = scale_term*eta_manager_dI2[i_eta];
                d_image.eta[idx] = value;
                d2_image.eta[idx] = value2;
            }
            //}
        }/* end eta deriv image increment */

        /*update the lambda derivative images*/
        for (int i_lam=0 ; i_lam < 2 ; i_lam++){
            if (db_flags.refine_lambda[i_lam]){
                //double value = scale_term*lambda_managers[i_lam]->dI;
                //double value2 = scale_term*lambda_managers[i_lam]->dI2;
                double value = scale_term*lambda_manager_dI[i_lam];
                double value2 = scale_term*lambda_manager_dI2[i_lam];
                int idx = i_lam*Npix_to_model + i_pix;
                d_image.lambda[idx] = value;
                d2_image.lambda[idx] = value2;
            }
        }/* end lambda deriv image increment */

        // panel rotation
        for (int i_pan_rot=0; i_pan_rot < 3; i_pan_rot++){
            if(db_flags.refine_panel_rot[i_pan_rot]){
                double value = scale_term*pan_rot_manager_dI[i_pan_rot];
                double value2 = scale_term*pan_rot_manager_dI2[i_pan_rot];
                int idx = i_pan_rot*Npix_to_model + i_pix;
                d_image.panel_rot[idx] = value;
                d2_image.panel_rot[idx] = value2;
            }
        }/* end panel rot deriv image increment */

        // panel origin
        for (int i_pan_orig=0; i_pan_orig < 3; i_pan_orig++){
            if(db_flags.refine_panel_origin[i_pan_orig]){
                double value = scale_term*pan_orig_manager_dI[i_pan_orig];
                double value2 = scale_term*pan_orig_manager_dI2[i_pan_orig];
                int idx = i_pan_orig*Npix_to_model + i_pix;
                d_image.panel_orig[idx] = value;
                d2_image.panel_orig[idx] = value2;
            }/* end panel orig deriv image increment */
        }
        if (db_flags.track_Fhkl){
            for(auto &x: Fhkl_tracker)
                printf("Pixel %d: Fhkl linear index %s came up %d times\n", i_pix, x.first.c_str(), x.second);
        }
    } // end i_pix loop
} // END of CPU kernel

}} // END namespace simtbx::nanoBragg
