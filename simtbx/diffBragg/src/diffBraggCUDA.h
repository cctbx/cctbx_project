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

void diffBragg_sum_over_steps_cuda(
        int Npix_to_model,
        std::vector<unsigned int>& panels_fasts_slows,
        image_type& floatimage,
        images& d_image,
        images& d2_image,
        step_arrays& db_steps,
        detector& db_det,
        beam& db_beam,
        crystal& db_cryst,
        flags& db_flags,
        cuda_flags& db_cu_flags,
        diffBragg_cudaPointers& cp,
        timer_variables& TIMERS);


void freedom(diffBragg_cudaPointers& cp);

