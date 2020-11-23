#include "diffBraggCUDA.h"

__global__
void gpu_sum_over_steps(
        int Npix_to_model, unsigned int* panels_fasts_slows,
        CUDAREAL* floatimage,
        CUDAREAL* d_Umat_images, CUDAREAL* d2_Umat_images,
        CUDAREAL* d_Bmat_images, CUDAREAL* d2_Bmat_images,
        CUDAREAL* d_Ncells_images, CUDAREAL* d2_Ncells_images,
        CUDAREAL* d_fcell_images, CUDAREAL* d2_fcell_images,
        CUDAREAL* d_eta_images,
        CUDAREAL* d_lambda_images, CUDAREAL* d2_lambda_images,
        CUDAREAL* d_panel_rot_images, CUDAREAL* d2_panel_rot_images,
        CUDAREAL* d_panel_orig_images, CUDAREAL* d2_panel_orig_images,
        CUDAREAL* d_sausage_XYZ_scale_images,
        CUDAREAL* d_sausage_XYZ_scale_images2,
        CUDAREAL* d_sausage_XYZ_scale_images3,
        CUDAREAL* d_sausage_XYZ_scale_images4,
        CUDAREAL* d_sausage_XYZ_scale_images5,
        CUDAREAL* d_sausage_XYZ_scale_images6,
        int* subS_pos, int* subF_pos, int* thick_pos,
        int* source_pos, int* phi_pos, int* mos_pos, int* sausage_pos,
        const int Nsteps, int _printout_fpixel, int _printout_spixel, bool _printout, CUDAREAL _default_F,
        int oversample, bool _oversample_omega, CUDAREAL subpixel_size, CUDAREAL pixel_size,
        CUDAREAL detector_thickstep, CUDAREAL _detector_thick, CUDAREAL close_distance, CUDAREAL detector_attnlen,
        bool use_lambda_coefficients, CUDAREAL lambda0, CUDAREAL lambda1,
        MAT3 eig_U, MAT3 eig_O, MAT3 eig_B, MAT3 RXYZ,
        VEC3* dF_vecs,
        VEC3* dS_vecs,
        MAT3* UMATS_RXYZ,
        MAT3* UMATS_RXYZ_prime,
        MAT3* RotMats,
        MAT3* dRotMats,
        MAT3* d2RotMats,
        MAT3* UMATS,
        MAT3* dB_mats,
        MAT3* dB2_mats,
        MAT3* sausages_RXYZ, MAT3* d_sausages_RXYZ, MAT3* sausages_U, CUDAREAL* sausages_scale,
        CUDAREAL* source_X, CUDAREAL* source_Y, CUDAREAL* source_Z, CUDAREAL* source_lambda, CUDAREAL* source_I,
        CUDAREAL kahn_factor,
        CUDAREAL Na, CUDAREAL Nb, CUDAREAL Nc,
        CUDAREAL phi0, CUDAREAL phistep,
        VEC3 spindle_vec, VEC3 _polarization_axis,
        int h_range, int k_range, int l_range,
        int h_max, int h_min, int k_max, int k_min, int l_max, int l_min, CUDAREAL dmin,
        CUDAREAL fudge, bool complex_miller, int verbose, bool only_save_omega_kahn,
        bool isotropic_ncells, bool compute_curvatures,
        CUDAREAL* _FhklLinear, CUDAREAL* _Fhkl2Linear,
        bool* refine_Bmat, bool* refine_Ncells, bool* refine_panel_origin, bool* refine_panel_rot,
        bool refine_fcell, bool* refine_lambda, bool refine_eta, bool* refine_Umat,
        bool refine_sausages, int num_sausages,
        CUDAREAL* fdet_vectors, CUDAREAL* sdet_vectors,
        CUDAREAL* odet_vectors, CUDAREAL* pix0_vectors,
        bool _nopolar, bool _point_pixel, CUDAREAL _fluence, CUDAREAL _r_e_sqr, CUDAREAL _spot_scale)
{ // BEGIN GPU kernel

    int tid = blockIdx.x * blockDim.x + threadIdx.x;
    int thread_stride = blockDim.x * gridDim.x;
    for (int i_pix=tid; i_pix < Npix_to_model; i_pix+= thread_stride){
       int _pid = panels_fasts_slows[i_pix*3];
       int _fpixel = panels_fasts_slows[i_pix*3+1];
       int _spixel = panels_fasts_slows[i_pix*3+2];

       // reset photon count for this pixel
       CUDAREAL _I=0;

       // reset derivative photon counts for the various parameters
       CUDAREAL rot_manager_dI[3] = {0,0,0};
       CUDAREAL rot_manager_dI2[3] = {0,0,0};
       CUDAREAL ucell_manager_dI[6]= {0,0,0,0,0,0};
       CUDAREAL ucell_manager_dI2[6]= {0,0,0,0,0,0};
       CUDAREAL Ncells_manager_dI[3]= {0,0,0};
       CUDAREAL Ncells_manager_dI2[3]= {0,0,0};
       CUDAREAL pan_orig_manager_dI[3]= {0,0,0};
       CUDAREAL pan_orig_manager_dI2[3]= {0,0,0};
       CUDAREAL pan_rot_manager_dI[3]= {0,0,0};
       CUDAREAL pan_rot_manager_dI2[3]= {0,0,0};
       CUDAREAL fcell_manager_dI=0;
       CUDAREAL fcell_manager_dI2=0;
       CUDAREAL eta_manager_dI = 0;
       CUDAREAL lambda_manager_dI[2] = {0,0};
       CUDAREAL lambda_manager_dI2[2] = {0,0};

       CUDAREAL sausage_manager_dI[24] = {0,0,0,0,0, // TODO use shared memory determined at runtime to increase max sausages
                                          0,0,0,0,0,
                                          0,0,0,0,0,
                                          0,0,0,0,0,
                                          0,0,0,0}; // maximum of 6 sausages!

       for (int _i_step=0; _i_step < Nsteps; _i_step++){

           int _subS = subS_pos[_i_step];
           int _subF = subF_pos[_i_step];
           int _thick_tic = thick_pos[_i_step];
           int _source = source_pos[_i_step];
           int _phi_tic = phi_pos[_i_step];
           int _mos_tic = mos_pos[_i_step];
           int _sausage_tic = sausage_pos[_i_step];

           // absolute mm position on detector (relative to its origin)
           CUDAREAL _Fdet = subpixel_size*(_fpixel*oversample + _subF ) + subpixel_size/2.0;
           CUDAREAL _Sdet = subpixel_size*(_spixel*oversample + _subS ) + subpixel_size/2.0;

           // assume "distance" is to the front of the detector sensor layer
           CUDAREAL _Odet = _thick_tic*detector_thickstep;
           int pid_x = _pid*3;
           int pid_y = _pid*3+1;
           int pid_z = _pid*3+2;
           VEC3 _o_vec(odet_vectors[pid_x], odet_vectors[pid_y], odet_vectors[pid_z]);

           CUDAREAL pixposX = _Fdet*fdet_vectors[pid_x]+_Sdet*sdet_vectors[pid_x]+_Odet*odet_vectors[pid_x]+pix0_vectors[pid_x];
           CUDAREAL pixposY = _Fdet*fdet_vectors[pid_y]+_Sdet*sdet_vectors[pid_y]+_Odet*odet_vectors[pid_y]+pix0_vectors[pid_y];
           CUDAREAL pixposZ = _Fdet*fdet_vectors[pid_z]+_Sdet*sdet_vectors[pid_z]+_Odet*odet_vectors[pid_z]+pix0_vectors[pid_z];
           VEC3 _pixel_pos(pixposX, pixposY, pixposZ);

           CUDAREAL _airpath = _pixel_pos.norm();
           VEC3 _diffracted = _pixel_pos/_airpath;

           // solid angle subtended by a pixel: (pix/airpath)^2*cos(2theta)
           CUDAREAL _omega_pixel = pixel_size*pixel_size/_airpath/_airpath*close_distance/_airpath;

           // option to turn off obliquity effect, inverse-square-law only
           if(_point_pixel) _omega_pixel = 1.0/_airpath/_airpath;

           // now calculate detector thickness effects
           CUDAREAL _capture_fraction = 1;
           if(_detector_thick > 0.0 && detector_attnlen > 0.0)
           {
               // inverse of effective thickness increase
               CUDAREAL _parallax = _diffracted.dot(_o_vec) ; //dot_product(diffracted,odet_vector);
               _capture_fraction = exp(-_thick_tic*detector_thickstep/detector_attnlen/_parallax)
                                 -exp(-(_thick_tic+1)*detector_thickstep/detector_attnlen/_parallax);
           }
           VEC3 _incident(-source_X[_source], -source_Y[_source], -source_Z[_source]);
           CUDAREAL _lambda = source_lambda[_source];
           CUDAREAL lambda_ang = _lambda*1e10;
           if (use_lambda_coefficients){
               lambda_ang = lambda0 + lambda1*lambda_ang;
               _lambda = lambda_ang*1e-10;
           }

           CUDAREAL _source_path = _incident.norm();
           _incident /= _source_path;

           VEC3 _scattering = (_diffracted - _incident) / _lambda;

           CUDAREAL _stol = 0.5*(_scattering.norm()); //magnitude(scattering);

           //if(dmin > 0.0 && _stol > 0.0)
           //{
           //    if(dmin > 0.5/_stol)
           //    {
           //        continue;
           //    }
           //}

          CUDAREAL _phi = phi0 + phistep*_phi_tic;
          MAT3 Bmat_realspace = eig_B;
          if( _phi != 0.0 )
          {
              CUDAREAL cosphi = cos(_phi);
              CUDAREAL sinphi = sin(_phi);
              VEC3 ap_vec(eig_B(0,0), eig_B(1,0), eig_B(2,0));
              VEC3 bp_vec(eig_B(0,1), eig_B(1,1), eig_B(2,1));
              VEC3 cp_vec(eig_B(0,2), eig_B(1,2), eig_B(2,2));

              ap_vec = ap_vec*cosphi + spindle_vec.cross(ap_vec)*sinphi + spindle_vec*(spindle_vec.dot(ap_vec))*(1-cosphi);
              bp_vec = bp_vec*cosphi + spindle_vec.cross(bp_vec)*sinphi + spindle_vec*(spindle_vec.dot(bp_vec))*(1-cosphi);
              cp_vec = cp_vec*cosphi + spindle_vec.cross(cp_vec)*sinphi + spindle_vec*(spindle_vec.dot(cp_vec))*(1-cosphi);

              Bmat_realspace << ap_vec[0], bp_vec[0], cp_vec[0],
                                  ap_vec[1], bp_vec[1], cp_vec[1],
                                  ap_vec[2], bp_vec[2], cp_vec[2];
          }
          Bmat_realspace *= 1e10;

          MAT3 U = sausages_U[_sausage_tic] * eig_U;
          MAT3 UBO = (UMATS_RXYZ[_mos_tic] * U*Bmat_realspace*(eig_O.transpose())).transpose();

          VEC3 q_vec(_scattering[0], _scattering[1], _scattering[2]);
          q_vec *= 1e-10;
          VEC3 H_vec = UBO*q_vec;

          CUDAREAL _h = H_vec[0];
          CUDAREAL _k = H_vec[1];
          CUDAREAL _l = H_vec[2];

          int _h0 = ceil(_h - 0.5);
          int _k0 = ceil(_k - 0.5);
          int _l0 = ceil(_l - 0.5);

          VEC3 H0(_h0, _k0, _l0);
          MAT3 _NABC;
          _NABC << Na,0,0,
                  0,Nb,0,
                  0,0,Nc;

          CUDAREAL C = 2 / 0.63 * fudge;
          VEC3 delta_H = H_vec - H0;
          VEC3 V = _NABC*delta_H;
          CUDAREAL _hrad_sqr = V.dot(V);
          CUDAREAL _F_latt = Na*Nb*Nc*exp(-( _hrad_sqr / 0.63 * fudge ));

          //if(_F_latt == 0.0 && ! only_save_omega_kahn) {
          //    continue;
          //}

          CUDAREAL _F_cell = _default_F;
          CUDAREAL _F_cell2 = 0;

          if ( (_h0<=h_max) && (_h0>=h_min) && (_k0<=k_max) && (_k0>=k_min) && (_l0<=l_max) && (_l0>=l_min)  ) {
              int Fhkl_linear_index = (_h0-h_min) * k_range * l_range + (_k0-k_min) * l_range + (_l0-l_min);
              _F_cell = _FhklLinear[Fhkl_linear_index];
              if (complex_miller) _F_cell2 = _Fhkl2Linear[Fhkl_linear_index];
          }

          if (complex_miller)
            _F_cell = sqrt(_F_cell*_F_cell + _F_cell2*_F_cell2);

          if (!_oversample_omega)
              _omega_pixel = 1;

          CUDAREAL Iincrement = _F_cell*_F_cell*_F_latt*_F_latt*source_I[_source]*_capture_fraction*_omega_pixel;
          Iincrement *= sausages_scale[_sausage_tic]*sausages_scale[_sausage_tic];
          _I += Iincrement;

          if(verbose > 3)
              printf("hkl= %f %f %f  hkl1= %d %d %d  Fcell=%f\n", _h,_k,_l,_h0,_k0,_l0, _F_cell);

          CUDAREAL two_C = 2*C;
          MAT3 UBOt = U*Bmat_realspace*(eig_O.transpose());
          if (refine_Umat[0]){
              MAT3 RyRzUBOt = RotMats[1]*RotMats[2]*UBOt;
              VEC3 delta_H_prime = (UMATS[_mos_tic]*dRotMats[0]*RyRzUBOt).transpose()*q_vec;
              CUDAREAL V_dot_dV = V.dot(_NABC*delta_H_prime);
              CUDAREAL value = -two_C * V_dot_dV * Iincrement;
              CUDAREAL value2 =0;
              if (compute_curvatures) {
                  VEC3 delta_H_dbl_prime = (UMATS[_mos_tic]*d2RotMats[0]*RyRzUBOt).transpose()*q_vec;
                  CUDAREAL dV_dot_dV = (_NABC*delta_H_prime).dot(_NABC*delta_H_prime);
                  CUDAREAL dV2_dot_V = (_NABC*delta_H).dot(_NABC*delta_H_dbl_prime);
                  value2 = two_C*(two_C*V_dot_dV*V_dot_dV - dV2_dot_V - dV_dot_dV)*Iincrement;
              }
              rot_manager_dI[0] += value;
              rot_manager_dI2[0] += value2;
          }
          if (refine_Umat[1]){
              MAT3 UmosRx = UMATS[_mos_tic]*RotMats[0];
              MAT3 RzUBOt = RotMats[2]*UBOt;
              VEC3 delta_H_prime =(UmosRx*dRotMats[1]*RzUBOt).transpose()*q_vec;
              CUDAREAL V_dot_dV = V.dot(_NABC*delta_H_prime);
              CUDAREAL value = -two_C * V_dot_dV * Iincrement;

              CUDAREAL value2=0;
              if (compute_curvatures){
                  VEC3 delta_H_dbl_prime = (UmosRx*d2RotMats[1]*RzUBOt).transpose()*q_vec;
                  CUDAREAL dV_dot_dV = (_NABC*delta_H_prime).dot(_NABC*delta_H_prime);
                  CUDAREAL dV2_dot_V = (_NABC*delta_H).dot(_NABC*delta_H_dbl_prime);
                  value2 = two_C*(two_C*V_dot_dV*V_dot_dV - dV2_dot_V - dV_dot_dV)*Iincrement;
              }
              rot_manager_dI[1] += value;
              rot_manager_dI2[1] += value2;
          }
          if (refine_Umat[2]){
              MAT3 UmosRxRy = UMATS[_mos_tic]*RotMats[0]*RotMats[1];
              VEC3 delta_H_prime = (UmosRxRy*dRotMats[2]*UBOt).transpose()*q_vec;
              CUDAREAL V_dot_dV = V.dot(_NABC*delta_H_prime);
              CUDAREAL value = -two_C * V_dot_dV * Iincrement;

              CUDAREAL value2=0;
              if (compute_curvatures){
                  VEC3 delta_H_dbl_prime = (UmosRxRy*d2RotMats[2]*UBOt).transpose()*q_vec;
                  CUDAREAL dV_dot_dV = (_NABC*delta_H_prime).dot(_NABC*delta_H_prime);
                  CUDAREAL dV2_dot_V = (_NABC*delta_H).dot(_NABC*delta_H_dbl_prime);
                  value2 = two_C*(two_C*V_dot_dV*V_dot_dV - dV2_dot_V - dV_dot_dV)*Iincrement;
              }
              rot_manager_dI[2] += value;
              rot_manager_dI2[2] += value2;
          }
          //Checkpoint for unit cell derivatives
          MAT3 Ot = eig_O.transpose();
          for(int i_uc=0; i_uc < 6; i_uc++ ){
              if (refine_Bmat[i_uc]){
                  MAT3 UmosRxRyRzU = UMATS_RXYZ[_mos_tic]*U;
                  VEC3 delta_H_prime = ((UmosRxRyRzU*(dB_mats[i_uc])*Ot).transpose()*q_vec);
                  CUDAREAL V_dot_dV = V.dot(_NABC*delta_H_prime);
                  CUDAREAL value = -two_C * V_dot_dV * Iincrement;
                  CUDAREAL value2 =0;
                  if (compute_curvatures){
                      VEC3 delta_H_dbl_prime = ((UmosRxRyRzU*(dB2_mats[i_uc])*Ot).transpose()*q_vec);
                      CUDAREAL dV_dot_dV = (_NABC*delta_H_prime).dot(_NABC*delta_H_prime);
                      CUDAREAL dV2_dot_V = (_NABC*delta_H).dot(_NABC*delta_H_dbl_prime);
                      value2 = two_C*(two_C*V_dot_dV*V_dot_dV - dV2_dot_V - dV_dot_dV)*Iincrement;
                  }
                  ucell_manager_dI[i_uc] += value;
                  ucell_manager_dI2[i_uc] += value2;
              }
          } //end ucell deriv

          // Checkpoint for Ncells manager
          if (refine_Ncells[0]){
              int num_ncell_deriv = 1;
              if (not isotropic_ncells)
                  num_ncell_deriv = 3;
              for (int i_nc=0; i_nc < num_ncell_deriv; i_nc++) {
                  MAT3 dN;
                  dN << 0,0,0,0,0,0,0,0,0;
                  dN(i_nc, i_nc) = 1;
                  CUDAREAL N_i = _NABC(i_nc, i_nc);
                  VEC3 dV_dN = dN*delta_H;
                  CUDAREAL deriv_coef = 1/N_i - C* ( dV_dN.dot(V));
                  CUDAREAL value = 2*Iincrement*deriv_coef;
                  CUDAREAL value2=0;
                  if(compute_curvatures){
                      dN(i_nc, i_nc) = 0; // TODO check maths
                      value2 = ( -1/N_i/N_i - C*(dV_dN.dot(dV_dN))) *2*Iincrement;
                      value2 += deriv_coef*2*value;
                  }
                  Ncells_manager_dI[i_nc] += value;
                  Ncells_manager_dI2[i_nc] += value2;
              }

          } // end Ncells manager deriv

          // Checkpoint for Origin manager
          for (int i_pan_orig=0; i_pan_orig < 3; i_pan_orig++){
              if (refine_panel_origin[i_pan_orig]){
                  CUDAREAL per_k = 1/_airpath;
                  CUDAREAL per_k3 = pow(per_k,3);
                  CUDAREAL per_k5 = pow(per_k,5);
                  CUDAREAL lambda_ang = _lambda*1e10;

                  MAT3 M = -two_C*(_NABC*UBO)/lambda_ang;
                  VEC3 dk;
                  if (i_pan_orig == 0)
                      dk << 0,0,1;
                  else if (i_pan_orig == 1)
                      dk << 1,0,0;
                  else
                      dk << 0,1,0;

                  CUDAREAL G = dk.dot(_pixel_pos);
                  CUDAREAL pix2 = subpixel_size*subpixel_size;
                  VEC3 dk_hat = -per_k3*G*_pixel_pos + per_k*dk;
                  CUDAREAL coef = (M*dk_hat).dot(V);
                  CUDAREAL coef2 = -3*pix2*per_k5*G * (_o_vec.dot(_pixel_pos));
                  coef2 += pix2*per_k3*(_o_vec.dot(dk));
                  CUDAREAL value = coef*Iincrement + coef2*Iincrement/_omega_pixel;

                  pan_orig_manager_dI[i_pan_orig] += value;
                  pan_orig_manager_dI2[i_pan_orig] += 0;

              } // end origin manager deriv
          }

          for (int i_pan_rot=0; i_pan_rot < 3; i_pan_rot++){
              if(refine_panel_rot[i_pan_rot]){
                  CUDAREAL per_k = 1/_airpath;
                  CUDAREAL per_k3 = pow(per_k,3);
                  CUDAREAL per_k5 = pow(per_k,5);
                  CUDAREAL lambda_ang = _lambda*1e10;
                  MAT3 M = -two_C*(_NABC*UBO)/lambda_ang;
                  VEC3 dk = _Fdet*(dF_vecs[_pid*3 + i_pan_rot]) + _Sdet*(dS_vecs[_pid*3 + i_pan_rot]);
                  CUDAREAL G = dk.dot(_pixel_pos);
                  CUDAREAL pix2 = subpixel_size*subpixel_size;
                  VEC3 dk_hat = -per_k3*G*_pixel_pos + per_k*dk;
                  CUDAREAL coef = (M*dk_hat).dot(V);
                  CUDAREAL coef2 = -3*pix2*per_k5*G * (_o_vec.dot(_pixel_pos));
                  coef2 += pix2*per_k3*(_o_vec.dot(dk));
                  CUDAREAL value = coef*Iincrement + coef2*Iincrement/_omega_pixel;

                  pan_rot_manager_dI[i_pan_rot] += value;
                  pan_rot_manager_dI2[i_pan_rot] += 0;
              }
          }

          // checkpoint for Fcell manager
          if (refine_fcell){
              CUDAREAL value = 2*Iincrement/_F_cell ;
              CUDAREAL value2=0;
              if (compute_curvatures){
                  value2 = value/_F_cell;
              }
              fcell_manager_dI += value;
              fcell_manager_dI2 += value2;
          } // end of fcell man deriv

          // checkpoint for eta manager
          if (refine_eta){
              VEC3 DeltaH_deriv = (UMATS_RXYZ_prime[_mos_tic]*UBOt).transpose()*q_vec;
              // vector V is _Nabc*Delta_H
              CUDAREAL value = -two_C*(V.dot(_NABC*DeltaH_deriv))*Iincrement;
              eta_manager_dI += value;
          } // end of eta man deriv

            // sausage deriv
          if (refine_sausages){
              MAT3 UBOt = eig_U*Bmat_realspace*(eig_O.transpose());
              int x = _sausage_tic*3;
              int y = _sausage_tic*3+1;
              int z = _sausage_tic*3+2;
              double value=0;
              for (int i=0;i<3; i++){
                  MAT3 UprimeBOt;
                  if (i==0)
                      UprimeBOt = d_sausages_RXYZ[x] * sausages_RXYZ[y] * sausages_RXYZ[z] * UBOt;
                  else if (i==1)
                      UprimeBOt = sausages_RXYZ[x] * d_sausages_RXYZ[y] * sausages_RXYZ[z] * UBOt;
                  else
                      UprimeBOt = sausages_RXYZ[x] * sausages_RXYZ[y] * d_sausages_RXYZ[z] * UBOt;

                  VEC3 DeltaH_deriv = (UMATS_RXYZ[_mos_tic]*UprimeBOt).transpose()*q_vec;
                  value = -two_C*(V.dot(_NABC*DeltaH_deriv))*Iincrement;
                  sausage_manager_dI[_sausage_tic*4 + i] += value;
              }
              // sausage scale derivative
              value = 2* Iincrement / sausages_scale[_sausage_tic];
              sausage_manager_dI[_sausage_tic*4 + 3] += value;
          }
          // end of sausage deriv

          // checkpoint for lambda manager
          for(int i_lam=0; i_lam < 2; i_lam++){
              if (refine_lambda[i_lam]){
                  CUDAREAL lambda_ang = _lambda*1e10;
                  CUDAREAL NH_dot_V = (_NABC*H_vec).dot(V);
                  CUDAREAL dg_dlambda;
                  if (i_lam==0)
                      dg_dlambda = 1;
                  else // i_lam==1
                      dg_dlambda = lambda_ang;
                  CUDAREAL coef = NH_dot_V*two_C*(dg_dlambda) / lambda_ang;
                  CUDAREAL value = coef*Iincrement;
                  CUDAREAL value2 = 0;
                  lambda_manager_dI[i_lam] += value;
                  lambda_manager_dI2[i_lam] += value2;
              }
          }
          //end of lambda deriv
          if( _printout && _i_step==0 ){
              if((_fpixel==_printout_fpixel && _spixel==_printout_spixel) || _printout_fpixel < 0)
              {
                 printf("%4d %4d : stol = %g, lambda = %g\n", _fpixel,_spixel,_stol, _lambda);
                 printf("at %g %g %g\n", _pixel_pos[0],_pixel_pos[1],_pixel_pos[2]);
                 printf("source XYZ %g %g %g\n", source_X[0],source_Y[0],source_Z[0]);
                 printf("hkl= %f %f %f  hkl0= %d %d %d\n", _h,_k,_l,_h0,_k0,_l0);
                 printf(" F_cell=%g  F_latt=%g   I = %g\n", _F_cell,_F_latt,_I);
                 printf("I/steps %15.10g\n", _I/Nsteps);
                 printf("omega   %15.10g\n", _omega_pixel);
                 printf("default_F= %f\n", _default_F);
                 printf("Incident[0]=%g, Incident[1]=%g, Incident[2]=%g\n", _incident[0], _incident[1], _incident[2]);
                 printf("source_path %g\n", _source_path);
                 for (int i_saus=0; i_saus<num_sausages; i_saus++){
                   printf("Sausages U (i_sausage=%d, scale=%f) :\n%f  %f  %f\n%f  %f  %f\n%f  %f  %f\n",
                    i_saus,sausages_scale[i_saus],
                    sausages_U[i_saus](0,0),  sausages_U[i_saus](0,1), sausages_U[i_saus](0,2),
                    sausages_U[i_saus](1,0),  sausages_U[i_saus](1,1), sausages_U[i_saus](1,2),
                    sausages_U[i_saus](2,0),  sausages_U[i_saus](2,1), sausages_U[i_saus](2,2));
                 }
              }
          }

       } // end of i_steps loop

       CUDAREAL _Fdet_ave = pixel_size*_fpixel + pixel_size/2.0;
       CUDAREAL _Sdet_ave = pixel_size*_spixel + pixel_size/2.0;
       CUDAREAL _Odet_ave = 0; //Odet; // TODO maybe make this more general for thick detectors?

       VEC3 _pixel_pos_ave(0,0,0);
       int pid_x = _pid*3;
       int pid_y = _pid*3+1;
       int pid_z = _pid*3+2;
       _pixel_pos_ave[0] = _Fdet_ave * fdet_vectors[pid_x]+_Sdet_ave*sdet_vectors[pid_x]+_Odet_ave*odet_vectors[pid_x]+pix0_vectors[pid_x];
       _pixel_pos_ave[1] = _Fdet_ave * fdet_vectors[pid_y]+_Sdet_ave*sdet_vectors[pid_y]+_Odet_ave*odet_vectors[pid_y]+pix0_vectors[pid_y];
       _pixel_pos_ave[2] = _Fdet_ave * fdet_vectors[pid_z]+_Sdet_ave*sdet_vectors[pid_z]+_Odet_ave*odet_vectors[pid_z]+pix0_vectors[pid_z];

       CUDAREAL _airpath_ave = _pixel_pos_ave.norm();
       VEC3 _diffracted_ave = _pixel_pos_ave/_airpath_ave;
       CUDAREAL _omega_pixel_ave = pixel_size*pixel_size/_airpath_ave/_airpath_ave*close_distance/_airpath_ave;

       CUDAREAL _polar = 1;
       if (!_nopolar){
           VEC3 _incident(-source_X[0], -source_Y[0], -source_Z[0]);
           _incident = _incident / _incident.norm();
           // component of diffracted unit vector along incident beam unit vector
           CUDAREAL cos2theta = _incident.dot(_diffracted_ave);
           CUDAREAL cos2theta_sqr = cos2theta*cos2theta;
           CUDAREAL sin2theta_sqr = 1-cos2theta_sqr;

           CUDAREAL _psi=0;
           if(kahn_factor != 0.0){
               // cross product to get "vertical" axis that is orthogonal to the cannonical "polarization"
               VEC3 B_in = _polarization_axis.cross(_incident);
               // cross product with incident beam to get E-vector direction
               VEC3 E_in = _incident.cross(B_in);
               // get components of diffracted ray projected onto the E-B plane
               CUDAREAL _kEi = _diffracted_ave.dot(E_in);
               CUDAREAL _kBi = _diffracted_ave.dot(B_in);
               // compute the angle of the diffracted ray projected onto the incident E-B plane
               _psi = -atan2(_kBi,_kEi);
           }
           // correction for polarized incident beam
           _polar = 0.5*(1.0 + cos2theta_sqr - kahn_factor*cos(2*_psi)*sin2theta_sqr);
       }

       CUDAREAL _om = 1;
       if (!_oversample_omega)
           _om=_omega_pixel_ave;
       // final scale term to being everything to photon number units
       CUDAREAL _scale_term = _r_e_sqr*_fluence*_spot_scale*_polar*_om / Nsteps*num_sausages;

       floatimage[i_pix] = _scale_term*_I;

       // udpate the rotation derivative images*
       for (int i_rot =0 ; i_rot < 3 ; i_rot++){
           if (refine_Umat[i_rot]){
               CUDAREAL value = _scale_term*rot_manager_dI[i_rot];
               CUDAREAL value2 = _scale_term*rot_manager_dI2[i_rot];
               int idx = i_rot*Npix_to_model + i_pix;
               d_Umat_images[idx] = value;
               //d2_Umat_images[idx] = value2;
           }
       } // end rot deriv image increment

       //update the ucell derivative images
       for (int i_uc=0 ; i_uc < 6 ; i_uc++){
           if (refine_Bmat[i_uc]){
               CUDAREAL value = _scale_term*ucell_manager_dI[i_uc];
               CUDAREAL value2 = _scale_term*ucell_manager_dI2[i_uc];
               int idx= i_uc*Npix_to_model + i_pix;
               d_Bmat_images[idx] = value;
               //d2_Bmat_images[idx] = value2;
           }
       }// end ucell deriv image increment

       //update the Ncells derivative image
       if (refine_Ncells[0]){
           CUDAREAL value = _scale_term*Ncells_manager_dI[0];
           CUDAREAL value2 = _scale_term*Ncells_manager_dI2[0];
           int idx = i_pix;
           d_Ncells_images[idx] = value;
           //d2_Ncells_images[idx] = value2;

           if (! isotropic_ncells){
               value = _scale_term*Ncells_manager_dI[1];
               value2 = _scale_term*Ncells_manager_dI2[1];
               idx = Npix_to_model + i_pix;
               d_Ncells_images[idx] = value;
               //d2_Ncells_images[idx] = value2;

               value = _scale_term*Ncells_manager_dI[2];
               value2 = _scale_term*Ncells_manager_dI2[2];
               idx = Npix_to_model*2 + i_pix;
               d_Ncells_images[idx] = value;
               //d2_Ncells_images[idx] = value2;
           }
       }// end Ncells deriv image increment

       // update Fcell derivative image
       if(refine_fcell){
           CUDAREAL value = _scale_term*fcell_manager_dI;
           CUDAREAL value2 = _scale_term*fcell_manager_dI2;
           d_fcell_images[i_pix] = value;
           //d2_fcell_images[i_pix] = value2;
       }// end Fcell deriv image increment

       // update eta derivative image
       if(refine_eta){
           CUDAREAL value = _scale_term*eta_manager_dI;
           CUDAREAL value2 = 0;
           d_eta_images[i_pix] = value;
       }// end eta deriv image increment

       //update the lambda derivative images
       for (int i_lam=0 ; i_lam < 2 ; i_lam++){
           if (refine_lambda[i_lam]){
               CUDAREAL value = _scale_term*lambda_manager_dI[i_lam];
               CUDAREAL value2 = _scale_term*lambda_manager_dI2[i_lam];
               int idx = i_lam*Npix_to_model + i_pix;
               d_lambda_images[idx] = value;
               //d2_lambda_images[idx] = value2;
           }
       }// end lambda deriv image increment

       // sausage increment
       if (refine_sausages){
           for (int i_sausage=0; i_sausage<num_sausages; i_sausage++){
               for (int i=0; i < 4; i++){
                   int sausage_parameter_i = i_sausage*4+i;
                   double value = _scale_term*sausage_manager_dI[sausage_parameter_i];
                   //int idx = sausage_parameter_i*Npix_to_model + i_pix;
                   int idx = i*Npix_to_model + i_pix;
                   if (i_sausage==0)
                       d_sausage_XYZ_scale_images[idx] = value;
                   else if (i_sausage==1)
                       d_sausage_XYZ_scale_images2[idx] = value;
                   else if (i_sausage==2)
                       d_sausage_XYZ_scale_images3[idx] = value;
                   else if (i_sausage==3)
                       d_sausage_XYZ_scale_images4[idx] = value;
                   else if (i_sausage==4)
                       d_sausage_XYZ_scale_images5[idx] = value;
                   else if (i_sausage==5)
                       d_sausage_XYZ_scale_images6[idx] = value;
                   }
           }
       }
       // end sausage

       for (int i_pan_rot=0; i_pan_rot < 3; i_pan_rot++){
           if(refine_panel_rot[i_pan_rot]){
               CUDAREAL value = _scale_term*pan_rot_manager_dI[i_pan_rot];
               CUDAREAL value2 = _scale_term*pan_rot_manager_dI2[i_pan_rot];
               int idx = i_pan_rot*Npix_to_model + i_pix;
               d_panel_rot_images[idx] = value;
               //d2_panel_rot_images[idx] = value2;
           }
       }// end panel rot deriv image increment

       for (int i_pan_orig=0; i_pan_orig < 3; i_pan_orig++){
           if(refine_panel_origin[i_pan_orig]){
               CUDAREAL value = _scale_term*pan_orig_manager_dI[i_pan_orig];
               CUDAREAL value2 = _scale_term*pan_orig_manager_dI2[i_pan_orig];
               int idx = i_pan_orig*Npix_to_model + i_pix;
               d_panel_orig_images[idx] = value;
               //d2_panel_orig_images[idx] = value2;
           }//end panel orig deriv image increment
       }
    } // end i_pix loop
}  // END of GPU kernel

