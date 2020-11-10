#include <cudatbx/cuda_base.cuh>
#include <simtbx/gpu/detector.h>

namespace simtbx {
namespace gpu {

  gpu_detector::gpu_detector(int const& arg_device_id, 
                             dxtbx::model::Detector const & arg_detector):
    h_deviceID(arg_device_id),
    detector(arg_detector),
    cu_accumulate_floatimage(NULL) {
    cudaSetDevice(arg_device_id);

    //1) determine the size
    cu_n_panels = detector.size();
    SCITBX_ASSERT( cu_n_panels >= 1);

    //2) confirm that array dimensions are similar for each size
    cu_slow_pixels = detector[0].get_image_size()[0];
    cu_fast_pixels = detector[0].get_image_size()[1];
    for (int ipanel=1; ipanel < detector.size(); ++ipanel){
      SCITBX_ASSERT(detector[ipanel].get_image_size()[0] == cu_slow_pixels);
      SCITBX_ASSERT(detector[ipanel].get_image_size()[1] == cu_fast_pixels);
    }
    _image_size = cu_n_panels * cu_slow_pixels * cu_fast_pixels;

    //3) allocate a cuda array with these dimensions
    /* separate accumulator image outside the usual nanoBragg data structure.
           1. accumulate contributions from a sequence of source energy channels computed separately
           2. represent multiple panels, all same rectangular shape; slowest dimension = n_panels */
    cudaSafeCall(cudaMalloc((void ** )&cu_accumulate_floatimage, 
                            sizeof(*cu_accumulate_floatimage) * _image_size));
    cudaSafeCall(cudaMemset((void *)cu_accumulate_floatimage, 0, 
                            sizeof(*cu_accumulate_floatimage) * _image_size));
    std::cout<<"allocated cuda array of "<<_image_size<<std::endl;
  };

  void gpu_detector::free_detail(){
    cudaSetDevice(h_deviceID);
    //4) make sure we can deallocate cuda array later on
    if (cu_accumulate_floatimage != NULL) {
      cudaSafeCall(cudaFree(cu_accumulate_floatimage));
    }
  };

  void
  gpu_detector::each_image_allocate_cuda(){
    cudaSetDevice(h_deviceID);
    /*allocate and zero reductions */
    bool * rangemap = (bool*) calloc(_image_size, sizeof(bool));
    float * omega_reduction = (float*) calloc(_image_size, sizeof(float));
    float * max_I_x_reduction = (float*) calloc(_image_size, sizeof(float));
    float * max_I_y_reduction = (float*) calloc(_image_size, sizeof(float));
    //It is not quite clear why we must zero them on CPU, why not just on GPU?
    
    cu_omega_reduction = NULL;
    cudaSafeCall(cudaMalloc((void ** )&cu_omega_reduction, sizeof(*cu_omega_reduction) * _image_size));
    printf("Allocated %d bytes to \n",sizeof(*cu_omega_reduction) * _image_size);
    std::cout<< cu_omega_reduction<<std::endl;
    cudaSafeCall(cudaMemcpy(cu_omega_reduction, 
                 omega_reduction, sizeof(*cu_omega_reduction) * _image_size, 
                 cudaMemcpyHostToDevice));

    cu_max_I_x_reduction = NULL;
    cudaSafeCall(cudaMalloc((void ** )&cu_max_I_x_reduction, sizeof(*cu_max_I_x_reduction) * _image_size));
    cudaSafeCall(cudaMemcpy(cu_max_I_x_reduction, 
                 max_I_x_reduction, sizeof(*cu_max_I_x_reduction) * _image_size,
                 cudaMemcpyHostToDevice));

    cu_max_I_y_reduction = NULL;
    cudaSafeCall(cudaMalloc((void ** )&cu_max_I_y_reduction, sizeof(*cu_max_I_y_reduction) * _image_size));
    cudaSafeCall(cudaMemcpy(cu_max_I_y_reduction, max_I_y_reduction, sizeof(*cu_max_I_y_reduction) * _image_size,
                 cudaMemcpyHostToDevice));

    cu_rangemap = NULL;
    cudaSafeCall(cudaMalloc((void ** )&cu_rangemap, sizeof(*cu_rangemap) * _image_size));
    cudaSafeCall(cudaMemcpy(cu_rangemap, 
                 rangemap, sizeof(*cu_rangemap) * _image_size,
                 cudaMemcpyHostToDevice));

    // deallocate host arrays
    // potential memory leaks
    free(rangemap);
    free(omega_reduction);
    free(max_I_x_reduction);
    free(max_I_y_reduction);

    cu_maskimage = NULL;
    int unsigned short * maskimage = NULL; //default case, must implement non-trivial initializer elsewhere
    if (maskimage != NULL) {
      cudaSafeCall(cudaMalloc((void ** )&cu_maskimage, sizeof(*cu_maskimage) * _image_size));
      cudaSafeCall(cudaMemcpy(cu_maskimage, maskimage, sizeof(*cu_maskimage) * _image_size,
                   cudaMemcpyHostToDevice));
    }

    // In contrast to old API, new API initializes its own accumulator, does not take values from CPU
    cu_floatimage = NULL;
    cudaSafeCall(cudaMalloc((void ** )&cu_floatimage, sizeof(*cu_floatimage) * _image_size));

  }

  void
  gpu_detector::each_image_free_cuda(){
    cudaSetDevice(h_deviceID);
    SCITBX_EXAMINE(h_deviceID);
    SCITBX_EXAMINE(cu_omega_reduction);
    SCITBX_EXAMINE(cu_max_I_x_reduction);
    SCITBX_EXAMINE(cu_max_I_y_reduction);
    SCITBX_EXAMINE(cu_rangemap);
    SCITBX_EXAMINE(cu_maskimage);
    SCITBX_EXAMINE(cu_floatimage);
    cudaSafeCall(cudaDeviceSynchronize());
Line:  112
Error: misaligned address
    cudaSafeCall(cudaFree(cu_omega_reduction));
    cudaSafeCall(cudaFree(cu_max_I_x_reduction));
    cudaSafeCall(cudaFree(cu_max_I_y_reduction));
    cudaSafeCall(cudaFree(cu_rangemap));
    cudaSafeCall(cudaFree(cu_maskimage));
    cudaSafeCall(cudaFree(cu_floatimage));
  }

} // gpu
} // simtbx
