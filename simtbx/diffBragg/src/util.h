
#ifndef SIMTBX_DIFFBRAGG_UTIL
#define SIMTBX_DIFFBRAGG_UTIL


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

#endif
