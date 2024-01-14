#ifndef SIMTBX_DIFFBRAGG_DIFFUSE_UTIL_KOKKOS
#define SIMTBX_DIFFBRAGG_DIFFUSE_UTIL_KOKKOS

#include <simtbx/diffBragg/src/util_kokkos.h>

KOKKOS_INLINE_FUNCTION
int gen_laue_mats(int laue_group_num, vector_mat3_t lmats, KOKKOS_MAT3 rpa) {

  assert(laue_group_num>0);
  assert(laue_group_num<15);
  int num_mats = 0;

  const double one_over_root2 = 1./sqrt(2.);

  if ( laue_group_num == 1 ) {
  // P -1

              // x,y,z
    lmats(0) <<  1, 0, 0,
                 0, 1, 0,
                 0, 0, 1;

    num_mats= 1;
  }
  if ( laue_group_num == 2 ) {
  // P 1 1 2/m

              // x,y,z
    lmats(0) <<  1, 0, 0,
                 0, 1, 0,
                 0, 0, 1;

              // -x,-y,z
    lmats(1) << -1, 0, 0,
                 0,-1, 0,
                 0, 0, 1;

    num_mats= 2;
  }
  if ( laue_group_num == 3 ) {
  // P 1 2/m 1

              // x,y,z
    lmats(0) <<  1, 0, 0,
                 0, 1, 0,
                 0, 0, 1;

              // -x,y,-z
    lmats(1) << -1, 0, 0,
                 0, 1, 0,
                 0, 0,-1;

    num_mats= 2;
  }
  if ( laue_group_num == 4 ) {
  // P 2/m 1 1

              // x,y,z
    lmats(0) <<  1, 0, 0,
                 0, 1, 0,
                 0, 0, 1;

              // x,-y,-z
    lmats(1) <<  1, 0, 0,
                 0,-1, 0,
                 0, 0,-1;

    num_mats= 2;
  }
  if ( laue_group_num == 5 ) {
  // P m m m

              // x,y,z
    lmats(0) <<  1, 0, 0,
                 0, 1, 0,
                 0, 0, 1;

              // x,-y,-z
    lmats(1) <<  1, 0, 0,
                 0,-1, 0,
                 0, 0,-1;

              // -x,y,-z
    lmats(2) << -1, 0, 0,
                 0, 1, 0,
                 0, 0,-1;

              // -x,-y,z
    lmats(3) << -1, 0, 0,
                 0,-1, 0,
                 0, 0, 1;

    num_mats= 4;
  }
  if ( laue_group_num == 6 ) {
  // P 4/m

              // x,y,z
    lmats(0) <<  1, 0, 0,
                 0, 1, 0,
                 0, 0, 1;

              // -y,x,z
    lmats(1) <<  0,-1, 0,
                 1, 0, 0,
                 0, 0, 1;

              // y,-x,z
    lmats(2) <<  0, 1, 0,
                -1, 0, 0,
                 0, 0, 1;

              // -x,-y,z
    lmats(3) << -1, 0, 0,
                 0,-1, 0,
                 0, 0, 1;

    num_mats= 4;
  }
  if ( laue_group_num == 7 ) {
  // P 4/m m m

              // x,y,z
    lmats(0) <<  1, 0, 0,
                 0, 1, 0,
                 0, 0, 1;

              // -y,x,z
    lmats(1) <<  0,-1, 0,
                 1, 0, 0,
                 0, 0, 1;

              // y,-x,z
    lmats(2) <<  0, 1, 0,
                -1, 0, 0,
                 0, 0, 1;

              // x,-y,-z
    lmats(3) <<  1, 0, 0,
                 0,-1, 0,
                 0, 0,-1;

              // -x,y,-z
    lmats(4) << -1, 0, 0,
                 0, 1, 0,
                 0, 0,-1;

              // -x,-y,z
    lmats(5) << -1, 0, 0,
                 0,-1, 0,
                 0, 0, 1;

              // y,x,-z
    lmats(6) <<  0, 1, 0,
                 1, 0, 0,
                 0, 0,-1;

              // -y,-x,-z
    lmats(7) <<  0,-1, 0,
                -1, 0, 0,
                 0, 0,-1;

    num_mats= 8;
  }
  if ( laue_group_num == 8 ) {
  // P -3

              // x,y,z
    lmats(0) <<  1, 0, 0,
                 0, 1, 0,
                 0, 0, 1;

              // -y,x-y,z
    lmats(1) <<  0,-1, 0,
                 one_over_root2,-one_over_root2, 0,
                 0, 0, 1;

              // -x+y,-x,z
    lmats(2) << -one_over_root2, one_over_root2, 0,
                -1, 0, 0,
                 0, 0, 1;

    num_mats= 3;
  }
  if ( laue_group_num == 9 ) {
  // P -3 m 1

              // x,y,z
    lmats(0) <<  1, 0, 0,
                 0, 1, 0,
                 0, 0, 1;

              // -y,x-y,z
    lmats(1) <<  0,-1, 0,
                 one_over_root2,-one_over_root2, 0,
                 0, 0, 1;

              // -x+y,-x,z
    lmats(2) << -one_over_root2, one_over_root2, 0,
                -1, 0, 0,
                 0, 0, 1;

              // x-y,-y,-z
    lmats(3) <<  one_over_root2,-one_over_root2, 0,
                 0,-1, 0,
                 0, 0,-1;

              // -x,-x+y,-z
    lmats(4) << -1, 0, 0,
                -one_over_root2, one_over_root2, 0,
                 0, 0,-1;

              // y,x,-z
    lmats(5) <<  0, 1, 0,
                 1, 0, 0,
                 0, 0,-1;

    num_mats= 6;
  }
  if ( laue_group_num == 10 ) {
  // P -3 1 m

              // x,y,z
    lmats(0) <<  1, 0, 0,
                 0, 1, 0,
                 0, 0, 1;

              // -y,x-y,z
    lmats(1) <<  0,-1, 0,
                 one_over_root2,-one_over_root2, 0,
                 0, 0, 1;

              // -x+y,-x,z
    lmats(2) << -one_over_root2, one_over_root2, 0,
                -1, 0, 0,
                 0, 0, 1;

              // -y,-x,-z
    lmats(3) <<  0,-1, 0,
                -1, 0, 0,
                 0, 0,-1;

              // -x+y,y,-z
    lmats(4) << -one_over_root2, one_over_root2, 0,
                 0, 1, 0,
                 0, 0,-1;

              // x,x-y,-z
    lmats(5) <<  1, 0, 0,
                 one_over_root2,-one_over_root2, 0,
                 0, 0,-1;

    num_mats= 6;
  }
  if ( laue_group_num == 11 ) {
  // P 6/m

              // x,y,z
    lmats(0) <<  1, 0, 0,
                 0, 1, 0,
                 0, 0, 1;

              // x-y,x,z
    lmats(1) <<  one_over_root2,-one_over_root2, 0,
                 1, 0, 0,
                 0, 0, 1;

              // y,-x+y,z
    lmats(2) <<  0, 1, 0,
                -one_over_root2, one_over_root2, 0,
                 0, 0, 1;

              // -y,x-y,z
    lmats(3) <<  0,-1, 0,
                 one_over_root2,-one_over_root2, 0,
                 0, 0, 1;

              // -x+y,-x,z
    lmats(4) << -one_over_root2, one_over_root2, 0,
                -1, 0, 0,
                 0, 0, 1;

              // -x,-y,z
    lmats(5) << -1, 0, 0,
                 0,-1, 0,
                 0, 0, 1;

    num_mats= 6;
  }
  if ( laue_group_num == 12 ) {
  // P 6/m m m

              // x,y,z
    lmats(0) <<  1, 0, 0,
                 0, 1, 0,
                 0, 0, 1;

              // x-y,x,z
    lmats(1) <<  one_over_root2,-one_over_root2, 0,
                 1, 0, 0,
                 0, 0, 1;

              // y,-x+y,z
    lmats(2) <<  0, 1, 0,
                -one_over_root2, one_over_root2, 0,
                 0, 0, 1;

              // -y,x-y,z
    lmats(3) <<  0,-1, 0,
                 one_over_root2,-one_over_root2, 0,
                 0, 0, 1;

              // -x+y,-x,z
    lmats(4) << -one_over_root2, one_over_root2, 0,
                -1, 0, 0,
                 0, 0, 1;

              // x-y,-y,-z
    lmats(5) <<  one_over_root2,-one_over_root2, 0,
                 0,-1, 0,
                 0, 0,-1;

              // -x,-x+y,-z
    lmats(6) << -1, 0, 0,
                -one_over_root2, one_over_root2, 0,
                 0, 0,-1;

              // -x,-y,z
    lmats(7) << -1, 0, 0,
                 0,-1, 0,
                 0, 0, 1;

              // y,x,-z
    lmats(8) <<  0, 1, 0,
                 1, 0, 0,
                 0, 0,-1;

              // -y,-x,-z
    lmats(9) <<  0,-1, 0,
                -1, 0, 0,
                 0, 0,-1;

              // -x+y,y,-z
    lmats(10) << -one_over_root2, one_over_root2, 0,
                  0, 1, 0,
                  0, 0,-1;

              // x,x-y,-z
    lmats(11) <<  1, 0, 0,
                  one_over_root2,-one_over_root2, 0,
                  0, 0,-1;

    num_mats= 12;
  }
  if ( laue_group_num == 13 ) {
  // P m -3

              // x,y,z
    lmats(0) <<  1, 0, 0,
                 0, 1, 0,
                 0, 0, 1;

              // z,x,y
    lmats(1) <<  0, 0, 1,
                 1, 0, 0,
                 0, 1, 0;

              // y,z,x
    lmats(2) <<  0, 1, 0,
                 0, 0, 1,
                 1, 0, 0;

              // -y,-z,x
    lmats(3) <<  0,-1, 0,
                 0, 0,-1,
                 1, 0, 0;

              // z,-x,-y
    lmats(4) <<  0, 0, 1,
                -1, 0, 0,
                 0,-1, 0;

              // -y,z,-x
    lmats(5) <<  0,-1, 0,
                 0, 0, 1,
                -1, 0, 0;

              // -z,-x,y
    lmats(6) <<  0, 0,-1,
                -1, 0, 0,
                 0, 1, 0;

              // -z,x,-y
    lmats(7) <<  0, 0,-1,
                 1, 0, 0,
                 0,-1, 0;

              // y,-z,-x
    lmats(8) <<  0, 1, 0,
                 0, 0,-1,
                -1, 0, 0;

              // x,-y,-z
    lmats(9) <<  1, 0, 0,
                 0,-1, 0,
                 0, 0,-1;

              // -x,y,-z
    lmats(10) << -1, 0, 0,
                  0, 1, 0,
                  0, 0,-1;

              // -x,-y,z
    lmats(11) << -1, 0, 0,
                  0,-1, 0,
                  0, 0, 1;

    num_mats= 12;
  }
  if ( laue_group_num == 14 ) {
  // P m -3 m

              // x,y,z
    lmats(0) <<  1, 0, 0,
                 0, 1, 0,
                 0, 0, 1;

              // x,-z,y
    lmats(1) <<  1, 0, 0,
                 0, 0,-1,
                 0, 1, 0;

              // x,z,-y
    lmats(2) <<  1, 0, 0,
                 0, 0, 1,
                 0,-1, 0;

              // z,y,-x
    lmats(3) <<  0, 0, 1,
                 0, 1, 0,
                -1, 0, 0;

              // -z,y,x
    lmats(4) <<  0, 0,-1,
                 0, 1, 0,
                 1, 0, 0;

              // -y,x,z
    lmats(5) <<  0,-1, 0,
                 1, 0, 0,
                 0, 0, 1;

              // y,-x,z
    lmats(6) <<  0, 1, 0,
                -1, 0, 0,
                 0, 0, 1;

              // z,x,y
    lmats(7) <<  0, 0, 1,
                 1, 0, 0,
                 0, 1, 0;

              // y,z,x
    lmats(8) <<  0, 1, 0,
                 0, 0, 1,
                 1, 0, 0;

              // -y,-z,x
    lmats(9) <<  0,-1, 0,
                 0, 0,-1,
                 1, 0, 0;

              // z,-x,-y
    lmats(10) <<  0, 0, 1,
                 -1, 0, 0,
                  0,-1, 0;

              // -y,z,-x
    lmats(11) <<  0,-1, 0,
                  0, 0, 1,
                 -1, 0, 0;

              // -z,-x,y
    lmats(12) <<  0, 0,-1,
                 -1, 0, 0,
                  0, 1, 0;

              // -z,x,-y
    lmats(13) <<  0, 0,-1,
                  1, 0, 0,
                  0,-1, 0;

              // y,-z,-x
    lmats(14) <<  0, 1, 0,
                  0, 0,-1,
                 -1, 0, 0;

              // x,-y,-z
    lmats(15) <<  1, 0, 0,
                  0,-1, 0,
                  0, 0,-1;

              // -x,y,-z
    lmats(16) << -1, 0, 0,
                  0, 1, 0,
                  0, 0,-1;

              // -x,-y,z
    lmats(17) << -1, 0, 0,
                  0,-1, 0,
                  0, 0, 1;

              // y,x,-z
    lmats(18) <<  0, 1, 0,
                  1, 0, 0,
                  0, 0,-1;

              // -y,-x,-z
    lmats(19) <<  0,-1, 0,
                 -1, 0, 0,
                  0, 0,-1;

              // z,-y,x
    lmats(20) <<  0, 0, 1,
                  0,-1, 0,
                  1, 0, 0;

              // -z,-y,-x
    lmats(21) <<  0, 0,-1,
                  0,-1, 0,
                 -1, 0, 0;

              // -x,z,y
    lmats(22) << -1, 0, 0,
                  0, 0, 1,
                  0, 1, 0;

              // -x,-z,-y
    lmats(23) << -1, 0, 0,
                  0, 0,-1,
                  0,-1, 0;

    num_mats = 24;
  }

  for (int i_mat=0; i_mat < num_mats; i_mat ++){
    lmats(i_mat) = lmats(i_mat) * rpa;
  }
  return num_mats;
};

// ***NEEDS UPDATE: use legacy API for passing diffuse scale as KOKKOS_MAT3
KOKKOS_FUNCTION
void calc_diffuse_at_hkl(KOKKOS_VEC3 H_vec, KOKKOS_VEC3 H0, KOKKOS_VEC3 dHH, int h_min, int k_min, int l_min, int h_max, int k_max, int l_max, int h_range, int k_range, int l_range, const KOKKOS_MAT3 diffuse_scale_mat3, const vector_cudareal_t FhklLinear, int num_laue_mats, const vector_mat3_t laue_mats, KOKKOS_MAT3 anisoG_local, vector_cudareal_t dG_trace, CUDAREAL anisoG_determ, KOKKOS_MAT3 anisoU_local, const vector_vec3_t dG_dgam, bool refine_diffuse, CUDAREAL *I0, CUDAREAL *step_diffuse_param){
   constexpr CUDAREAL four_mpi_sq = 4.*M_PI*M_PI;
   // loop over laue matrices
   int num_stencil_points = (2*dHH[0] + 1) * (2*dHH[1] + 1) * (2*dHH[2] + 1);
   bool h_bounded= (H0[0]+dHH[0]<=h_max) && (H0[0]-dHH[0]>=h_min) ;
   bool k_bounded= (H0[1]+dHH[1]<=k_max) && (H0[1]-dHH[1]>=k_min) ;
   bool l_bounded= (H0[2]+dHH[2]<=l_max) && (H0[2]-dHH[2]>=l_min) ;
   if (h_bounded && k_bounded && l_bounded) {
         int Fhkl_linear_index_0 = (H0[0]-h_min) * k_range * l_range + (H0[1]-k_min) * l_range + (H0[2]-l_min);
         const CUDAREAL _F_cell_0 = FhklLinear(Fhkl_linear_index_0);
         for (int hh=-dHH[0]; hh <= dHH[0]; hh++){
            for (int kk=-dHH[1]; kk <= dHH[1]; kk++){
               for (int ll=-dHH[2]; ll <= dHH[2]; ll++){
                     CUDAREAL ID_this = 0;
                     CUDAREAL step_diffuse_param_this[6]  = {0,0,0,0,0,0};
                     const int Fhkl_linear_index_this = hh * k_range * l_range + kk * l_range + ll + Fhkl_linear_index_0;
                     const CUDAREAL _F_cell_this = FhklLinear(Fhkl_linear_index_this);
                     CUDAREAL _this_diffuse_scale;
                     if (_F_cell_0 != 0.0)
                        _this_diffuse_scale = _F_cell_this/_F_cell_0;
                     else
                        _this_diffuse_scale = 1.0;

                     _this_diffuse_scale *= _this_diffuse_scale * diffuse_scale_mat3(0,0) / (CUDAREAL)num_laue_mats; // legacy API for passing diffuse_scale

                     const KOKKOS_VEC3 H0_offset(H0[0]+hh, H0[1]+kk, H0[2]+ll);
                     const KOKKOS_VEC3 delta_H_offset = H_vec - H0_offset;

                     for ( int iL = 0; iL < num_laue_mats; iL++ ){
                        const KOKKOS_VEC3 Q0 = laue_mats(iL)*H0;
                        const CUDAREAL exparg = four_mpi_sq*Q0.dot(anisoU_local*Q0);
                        const CUDAREAL dwf = exp(-exparg);
                        const KOKKOS_VEC3 delta_Q = laue_mats(iL)*delta_H_offset;
                        const KOKKOS_VEC3 anisoG_q = anisoG_local*delta_Q;

                        const CUDAREAL V_dot_V = anisoG_q.length_sqr();
                        const CUDAREAL gamma_portion_denom = 1 / (1.+ V_dot_V* four_mpi_sq);
                        const CUDAREAL gamma_portion = 8.*M_PI*anisoG_determ * gamma_portion_denom * gamma_portion_denom;
                        const CUDAREAL this_I_latt_diffuse = dwf*exparg*gamma_portion;

                        ID_this += this_I_latt_diffuse;
                        if (refine_diffuse){ // add the contributions to diffuse scattering gradients here
                           for (int i_gam=0; i_gam<3; i_gam++){
                                 const CUDAREAL dV = dG_dgam(i_gam).dot(delta_Q);
                                 const CUDAREAL V_dot_dV = anisoG_q[i_gam] * dV;
                                 const CUDAREAL deriv = dG_trace(i_gam) - 4.*four_mpi_sq*V_dot_dV*gamma_portion_denom;
                                 step_diffuse_param_this[i_gam] += gamma_portion*deriv*dwf*exparg;
                           }
                           for (int i_sig = 0;i_sig<3; i_sig++){
                                 const CUDAREAL dexparg = 2 * four_mpi_sq * sqrt(anisoU_local(i_sig,i_sig)) * Q0[i_sig] * Q0[i_sig];
                                 step_diffuse_param_this[i_sig+3] += gamma_portion*dwf*dexparg*(1. - exparg);
                           }
                        }
                     } // end loop over iL (laue group mats)
                     // Update the lattice interference term here to include diffuse scattering (F_latt squared)
                     *I0 += ID_this * _this_diffuse_scale;

                     if (refine_diffuse) {
                        for (int idp=0; idp < 6; idp++) {
                           step_diffuse_param[idp] += step_diffuse_param_this[idp]*_this_diffuse_scale;
                        }
                     }
               } // end ll loop
            } // end kk loop
         } // end hh loop
   } // end if bounded
}

#endif
