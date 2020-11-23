#include <Eigen/Dense>
#include <iostream>
#include <math.h>
#include <vector>
#ifndef CUDAREAL
#define CUDAREAL double
#endif

typedef std::vector<CUDAREAL> image_type;
typedef Eigen::Matrix<CUDAREAL,3,1> VEC3;
typedef Eigen::Matrix<CUDAREAL,3,3> MAT3;
typedef std::vector<MAT3,Eigen::aligned_allocator<MAT3> > eigMat3_vec;
typedef std::vector<VEC3,Eigen::aligned_allocator<VEC3> > eigVec3_vec;

__global__ void gpu_sum_over_steps(
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
        bool _nopolar, bool _point_pixel, CUDAREAL _fluence, CUDAREAL _r_e_sqr, CUDAREAL _spot_scale);
