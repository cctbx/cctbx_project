#include <Eigen/Dense>
#include <iostream>
#include <math.h>
#include <vector>
#include <simtbx/diffBragg/src/util.h>
#ifndef CUDAREAL
#define CUDAREAL double
#endif

typedef std::vector<CUDAREAL> image_type;
typedef Eigen::Matrix<CUDAREAL,3,1> VEC3;
typedef Eigen::Matrix<CUDAREAL,3,3> MAT3;
typedef std::vector<MAT3,Eigen::aligned_allocator<MAT3> > eigMat3_vec;
typedef std::vector<VEC3,Eigen::aligned_allocator<VEC3> > eigVec3_vec;

//#define CUDA_CHECK_RETURN(value) CheckCudaErrorAux(__FILE__,__LINE__, #value, value)

struct diffBragg_cudaPointers {

  bool device_is_allocated = false;
  int npix_allocated=0;

  unsigned int* cu_panels_fasts_slows;

  CUDAREAL* cu_floatimage;
  CUDAREAL* cu_d_Umat_images=NULL;
  CUDAREAL* cu_d_Bmat_images=NULL;
  CUDAREAL* cu_d_Ncells_images=NULL;
  CUDAREAL* cu_d_fcell_images=NULL;
  CUDAREAL* cu_d_eta_images=NULL;
  CUDAREAL* cu_d2_eta_images=NULL;
  CUDAREAL* cu_d_lambda_images=NULL;
  CUDAREAL* cu_d_panel_rot_images=NULL;
  CUDAREAL* cu_d_panel_orig_images=NULL;

  CUDAREAL* cu_d2_Umat_images=NULL;
  CUDAREAL* cu_d2_Bmat_images=NULL;
  CUDAREAL* cu_d2_Ncells_images=NULL;
  CUDAREAL* cu_d2_fcell_images=NULL;
  CUDAREAL* cu_d2_lambda_images=NULL;
  CUDAREAL* cu_d2_panel_rot_images=NULL;
  CUDAREAL* cu_d2_panel_orig_images=NULL;

  CUDAREAL* cu_d_sausage_XYZ_scale_images=NULL;
  CUDAREAL* cu_d_fp_fdp_images=NULL;

  int* cu_subS_pos;
  int* cu_subF_pos;
  int* cu_thick_pos;
  int* cu_source_pos;
  int* cu_mos_pos;
  int* cu_phi_pos;
  int* cu_sausage_pos;

  CUDAREAL * cu_Fhkl;
  CUDAREAL * cu_Fhkl2=NULL;

  CUDAREAL * cu_fdet_vectors;
  CUDAREAL * cu_sdet_vectors;
  CUDAREAL * cu_odet_vectors;
  CUDAREAL * cu_pix0_vectors;
  CUDAREAL * cu_close_distances;

  int * cu_nominal_hkl=NULL;
  CUDAREAL * cu_fpfdp=NULL;
  CUDAREAL * cu_fpfdp_derivs=NULL;
  CUDAREAL * cu_atom_data=NULL;

  CUDAREAL * cu_source_X, * cu_source_Y, * cu_source_Z, * cu_source_I, * cu_source_lambda;
  int cu_sources;
  bool sources_are_allocated = false;
  bool sources_recopy = false;

  Eigen::Matrix3d* cu_UMATS;
  Eigen::Matrix3d* cu_dB_Mats;
  Eigen::Matrix3d* cu_dB2_Mats;
  Eigen::Matrix3d* cu_UMATS_RXYZ;
  Eigen::Matrix3d* cu_UMATS_RXYZ_prime=NULL;
  Eigen::Matrix3d* cu_UMATS_RXYZ_dbl_prime=NULL;
  Eigen::Matrix3d* cu_RotMats;
  Eigen::Matrix3d* cu_dRotMats;
  Eigen::Matrix3d* cu_d2RotMats;

  Eigen::Matrix3d* cu_AMATS;

  Eigen::Vector3d* cu_dF_vecs;
  Eigen::Vector3d* cu_dS_vecs;

  Eigen::Matrix3d* cu_sausages_RXYZ;
  Eigen::Matrix3d* cu_d_sausages_RXYZ;
  Eigen::Matrix3d* cu_sausages_U;
  CUDAREAL* cu_sausages_scale;

  bool* cu_refine_Bmat;
  bool* cu_refine_Umat;
  bool* cu_refine_Ncells;
  bool* cu_refine_lambda;
  bool* cu_refine_panel_origin;
  bool* cu_refine_panel_rot;

};

//TODO can loopy take a diffBragg object as argument ?
void diffBragg_loopy(
        int Npix_to_model, std::vector<unsigned int>& panels_fasts_slows,
        image_type& floatimage,
        images& d_image, images& d2_image,
        const int Nsteps, int _printout_fpixel, int _printout_spixel, bool _printout, CUDAREAL _default_F,
        int oversample, bool _oversample_omega, CUDAREAL subpixel_size, CUDAREAL pixel_size,
        CUDAREAL detector_thickstep, CUDAREAL _detector_thick, std::vector<CUDAREAL>& close_distances, CUDAREAL detector_attnlen,
        bool use_lambda_coefficients, CUDAREAL lambda0, CUDAREAL lambda1,
        MAT3& eig_U, MAT3& eig_O, MAT3& eig_B, MAT3& RXYZ,
        eigVec3_vec& dF_vecs,
        eigVec3_vec& dS_vecs,
        eigMat3_vec& UMATS_RXYZ,
        eigMat3_vec& UMATS_RXYZ_prime,
        eigMat3_vec& UMATS_RXYZ_dbl_prime,
        eigMat3_vec& RotMats,
        eigMat3_vec& dRotMats,
        eigMat3_vec& d2RotMats,
        eigMat3_vec& UMATS,
        eigMat3_vec& dB_Mats,
        eigMat3_vec& dB2_Mats,
        eigMat3_vec& sausages_RXYZ, eigMat3_vec& d_sausages_RXYZ,
        eigMat3_vec& sausages_U, image_type& sausages_scale,
        CUDAREAL* source_X, CUDAREAL* source_Y, CUDAREAL* source_Z, CUDAREAL* source_lambda, CUDAREAL* source_I,
        CUDAREAL kahn_factor,
        CUDAREAL Na, CUDAREAL Nb, CUDAREAL Nc,
        CUDAREAL Nd, CUDAREAL Ne, CUDAREAL Nf,
        CUDAREAL phi0, CUDAREAL phistep,
        VEC3& spindle_vec, VEC3 _polarization_axis,
        int h_range, int k_range, int l_range,
        int h_max, int h_min, int k_max, int k_min, int l_max, int l_min, CUDAREAL dmin,
        CUDAREAL fudge, bool complex_miller, int verbose, bool only_save_omega_kahn,
        bool isotropic_ncells, bool compute_curvatures,
        std::vector<CUDAREAL>& _FhklLinear, std::vector<CUDAREAL>& _Fhkl2Linear,
        std::vector<bool>& refine_Bmat, std::vector<bool>& refine_Ncells, bool refine_Ncells_def, std::vector<bool>& refine_panel_origin, std::vector<bool>& refine_panel_rot,
        bool refine_fcell, std::vector<bool>& refine_lambda, bool refine_eta, std::vector<bool>& refine_Umat,
        bool refine_sausages, int num_sausages,
        bool refine_fp_fdp,
        std::vector<CUDAREAL>& fdet_vectors, std::vector<CUDAREAL>& sdet_vectors,
        std::vector<CUDAREAL>& odet_vectors, std::vector<CUDAREAL>& pix0_vectors,
        bool _nopolar, bool _point_pixel, CUDAREAL _fluence, CUDAREAL _r_e_sqr, CUDAREAL _spot_scale,
        int number_of_sources, int device_Id,
        diffBragg_cudaPointers& cp,
        bool update_step_positions, bool update_panels_fasts_slows, bool update_sources, bool update_umats,
        bool update_dB_mats, bool update_rotmats, bool update_Fhkl, bool update_detector, bool update_refine_flags,
        bool update_panel_deriv_vecs, bool update_sausages_on_device,
        int detector_thicksteps, int phisteps, int Npix_to_allocate, bool no_Nabc_scale,
        std::vector<double>& fpfdp,
        std::vector<double>& fpfdp_derivs,
        std::vector<double>& atom_data, std::vector<int>& nominal_hkl , timer_variables& TIMERS);


void freedom(diffBragg_cudaPointers& cp);

