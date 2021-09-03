
#ifndef SIMTBX_DIFFBRAGG_UTIL
#define SIMTBX_DIFFBRAGG_UTIL

#include <Eigen/Dense>
#include<Eigen/StdVector>
#include<vector>

typedef std::vector<double> image_type;

struct timer_variables{
    double add_spots_pre=0; // times the initializations for add spots kernel
    double add_spots_post=0; // times the copies that occur after add spots kernel
    double add_spots_kernel_wrapper=0; // times the add spots kernel overall, either CPU or GPU
    double cuda_alloc=0; // times the allocation of the device
    double cuda_copy_to_dev=0; // times the copying from host to device
    double cuda_copy_from_dev=0; // times the copying back from device to host
    double cuda_kernel=0; // times the GPU kernel
    int timings=0; // how many times these variables were incremented
    bool recording=true;
  };


// CONTAINERS
struct images{
    image_type Umat;
    image_type Bmat;
    image_type Ncells;
    image_type fcell;
    image_type eta;
    image_type lambda;
    image_type panel_rot;
    image_type panel_orig;
    image_type sausage;
    image_type fp_fdp;
};


struct step_arrays{
   int* subS_pos;
   int* subF_pos;
   int* thick_pos;
   int* source_pos;
   int* phi_pos;
   int* mos_pos;
   int Nsteps;
};


struct cuda_flags{
    int device_Id=0;
    bool update_step_positions;
    bool update_panels_fasts_slows;
    bool update_sources;
    bool update_umats;
    bool update_dB_mats;
    bool update_rotmats;
    bool update_Fhkl;
    bool update_detector;
    bool update_refine_flags;
    bool update_panel_deriv_vecs;
    int Npix_to_allocate;
};

struct flags{
    bool track_Fhkl;
    bool printout;
    bool nopolar;
    bool point_pixel;
    bool only_save_omega_kahn;
    bool compute_curvatures;
    bool isotropic_ncells;
    bool complex_miller;
    bool no_Nabc_scale;
    std::vector<bool> refine_Bmat;
    std::vector<bool> refine_Ncells;
    bool refine_Ncells_def;
    std::vector<bool> refine_panel_origin;
    std::vector<bool> refine_panel_rot;
    bool refine_fcell;
    std::vector<bool> refine_lambda;
    bool refine_eta;
    std::vector<bool> refine_Umat;
    bool refine_fp_fdp;
    bool use_lambda_coefficients;
    bool oversample_omega;
    int printout_fpixel, printout_spixel;
    int verbose;
    bool use_diffuse = false;
};

struct crystal{
    double this_gamma=50;
    double this_sigma=10;


    int mosaic_domains;
    double Na, Nb, Nc, Nd, Ne, Nf;
    double phi0;
    double phistep;
    double phisteps;
    double fudge;
    double spot_scale;
    int h_range, k_range, l_range;
    int h_max, h_min, k_max, k_min, l_max, l_min;
    double dmin;
    std::vector<double> FhklLinear, Fhkl2Linear;
    std::vector<double> fpfdp;
    std::vector<double> fpfdp_derivs;
    std::vector<double> atom_data;
    std::vector<int> nominal_hkl;
    double default_F;
    double r_e_sqr;

    Eigen::Matrix3d eig_U;
    Eigen::Matrix3d eig_O;
    Eigen::Matrix3d eig_B;
    Eigen::Matrix3d RXYZ;
    Eigen::Vector3d spindle_vec;

    std::vector<Eigen::Matrix3d,Eigen::aligned_allocator<Eigen::Matrix3d> > UMATS_RXYZ;
    std::vector<Eigen::Matrix3d,Eigen::aligned_allocator<Eigen::Matrix3d> > UMATS_RXYZ_prime;
    std::vector<Eigen::Matrix3d,Eigen::aligned_allocator<Eigen::Matrix3d> > UMATS_RXYZ_dbl_prime;
    std::vector<Eigen::Matrix3d,Eigen::aligned_allocator<Eigen::Matrix3d> > RotMats;
    std::vector<Eigen::Matrix3d,Eigen::aligned_allocator<Eigen::Matrix3d> > dRotMats;
    std::vector<Eigen::Matrix3d,Eigen::aligned_allocator<Eigen::Matrix3d> > d2RotMats;
    std::vector<Eigen::Matrix3d,Eigen::aligned_allocator<Eigen::Matrix3d> > UMATS;
    std::vector<Eigen::Matrix3d,Eigen::aligned_allocator<Eigen::Matrix3d> > UMATS_prime;
    std::vector<Eigen::Matrix3d,Eigen::aligned_allocator<Eigen::Matrix3d> > UMATS_dbl_prime;
    std::vector<Eigen::Matrix3d,Eigen::aligned_allocator<Eigen::Matrix3d> > dB_Mats;
    std::vector<Eigen::Matrix3d,Eigen::aligned_allocator<Eigen::Matrix3d> > dB2_Mats;

};

struct beam{
    Eigen::Vector3d polarization_axis;
    double fluence;
    double kahn_factor;
    double *source_X, *source_Y, *source_Z, *source_lambda, *source_I;
    double lambda0,lambda1;
    int number_of_sources;
};

struct detector{
    std::vector<Eigen::Vector3d,Eigen::aligned_allocator<Eigen::Vector3d> > dF_vecs;
    std::vector<Eigen::Vector3d,Eigen::aligned_allocator<Eigen::Vector3d> > dS_vecs;
    double detector_thickstep, detector_thicksteps, detector_thick, detector_attnlen;
    std::vector<double> close_distances;
    int oversample;
    double subpixel_size, pixel_size;
    std::vector<double> fdet_vectors, sdet_vectors, odet_vectors, pix0_vectors;
};

#endif
