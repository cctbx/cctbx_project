
#ifndef SIMTBX_DIFFBRAGG_UTIL
#define SIMTBX_DIFFBRAGG_UTIL

#include <Eigen/Dense>
#include<Eigen/StdVector>


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
};

struct eigen_objects{
   Eigen::Matrix3d eig_U;
   Eigen::Matrix3d eig_O;
   Eigen::Matrix3d eig_B;
   Eigen::Matrix3d RXYZ;

   std::vector<Eigen::Vector3d,Eigen::aligned_allocator<Eigen::Vector3d> > dF_vecs;
   std::vector<Eigen::Vector3d,Eigen::aligned_allocator<Eigen::Vector3d> > dS_vecs;

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

   Eigen::Vector3d spindle_vec;
   Eigen::Vector3d _polarization_axis;
};

#endif
