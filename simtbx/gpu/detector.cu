#include <scitbx/array_family/boost_python/flex_fwd.h>
#include <cudatbx/cuda_base.cuh>
#include <simtbx/gpu/detector.h>
#include <simtbx/gpu/detector.cuh>
#include <scitbx/vec3.h>
#include <scitbx/vec2.h>
#define THREADS_PER_BLOCK_X 128
#define THREADS_PER_BLOCK_Y 1
#define THREADS_PER_BLOCK_TOTAL (THREADS_PER_BLOCK_X * THREADS_PER_BLOCK_Y)

namespace simtbx {
namespace gpu {

//refactor later into helper file
  static cudaError_t detMemcpyVectorDoubleToDevice(CUDAREAL *dst, const double *src, size_t vector_items) {
        CUDAREAL * temp = new CUDAREAL[vector_items];
        for (size_t i = 0; i < vector_items; i++) {
                temp[i] = src[i];
        }
        cudaError_t ret = cudaMemcpy(dst, temp, sizeof(*dst) * vector_items, cudaMemcpyHostToDevice);
        delete temp;
        return ret;
  }

  packed_metrology::packed_metrology(dxtbx::model::Detector const & arg_detector,
                                   dxtbx::model::Beam const & arg_beam) {

    for (std::size_t panel_id = 0; panel_id < arg_detector.size(); panel_id++){
          // helper code arising from the nanoBragg constructor, with user_beam=True
      typedef scitbx::vec3<double> vec3;

      /* DETECTOR properties */
      /* typically: 1 0 0 */
      vec3 fdet_vector = arg_detector[panel_id].get_fast_axis();
      fdet_vector = fdet_vector.normalize();

      /* typically: 0 -1 0 */
      vec3 sdet_vector = arg_detector[panel_id].get_slow_axis();
      sdet_vector = sdet_vector.normalize();

      /* set orthogonal vector to the detector pixel array */
      vec3 odet_vector = fdet_vector.cross(sdet_vector);
      odet_vector = odet_vector.normalize();

      /* dxtbx origin is location of outer corner of the first pixel */
      vec3 pix0_vector = arg_detector[panel_id].get_origin()/1000.0;

      /* what is the point of closest approach between sample and detector? */
      double close_distance = pix0_vector * odet_vector;
      if (close_distance < 0){
        bool verbose = false;
        if(verbose)printf("WARNING: dxtbx model is lefthanded. Inverting odet_vector.\n");
        odet_vector = -1. * odet_vector;
        close_distance = -1*close_distance;
      }

      sdet.push_back(sdet_vector.length());
      fdet.push_back(fdet_vector.length());
      odet.push_back(odet_vector.length());
      pix0.push_back(0.);
      for (std::size_t idx_vec = 0; idx_vec < 3; idx_vec++){
            sdet.push_back(sdet_vector[idx_vec]);
            fdet.push_back(fdet_vector[idx_vec]);
            odet.push_back(odet_vector[idx_vec]);
            pix0.push_back(pix0_vector[idx_vec]);
      }
      /* set beam centre */
      scitbx::vec2<double> dials_bc=arg_detector[panel_id].get_beam_centre(arg_beam.get_s0());
      dists.push_back(close_distance);
      Xbeam.push_back(dials_bc[0]/1000.0);
      Ybeam.push_back(dials_bc[1]/1000.0);
    }
  };

  packed_metrology::packed_metrology(const simtbx::nanoBragg::nanoBragg& nB){
      for (std::size_t idx_vec = 0; idx_vec < 4; idx_vec++){
            sdet.push_back(nB.sdet_vector[idx_vec]);
            fdet.push_back(nB.fdet_vector[idx_vec]);
            odet.push_back(nB.odet_vector[idx_vec]);
            pix0.push_back(nB.pix0_vector[idx_vec]);
      }
      dists.push_back(nB.close_distance);
      Xbeam.push_back(nB.Xbeam);
      Ybeam.push_back(nB.Ybeam);
  }

  void
  packed_metrology::show() const {
    for (std::size_t idx_p = 0; idx_p < Xbeam.size(); idx_p++){
      printf(" Panel %3d\n",idx_p);
      printf(" Panel %3d sdet %9.6f %9.6f %9.6f %9.6f fdet %9.6f %9.6f %9.6f %9.6f\n",
             idx_p,sdet[4*idx_p+0],sdet[4*idx_p+1],sdet[4*idx_p+2],sdet[4*idx_p+3],
                          fdet[4*idx_p+0],fdet[4*idx_p+1],fdet[4*idx_p+2],fdet[4*idx_p+3]
      );
      printf(" Panel %3d odet %9.6f %9.6f %9.6f %9.6f pix0 %9.6f %9.6f %9.6f %9.6f\n",
             idx_p,odet[4*idx_p+0],odet[4*idx_p+1],odet[4*idx_p+2],odet[4*idx_p+3],
                          pix0[4*idx_p+0],pix0[4*idx_p+1],pix0[4*idx_p+2],pix0[4*idx_p+3]
      );
      printf(" Panel %3d beam %11.8f %11.8f\n",idx_p,Xbeam[idx_p],Ybeam[idx_p]);
    }
  }

  void
  gpu_detector::construct_detail(int const& arg_device_id,
                                 dxtbx::model::Detector const & arg_detector){
    cudaSetDevice(arg_device_id);

    //1) determine the size
    cu_n_panels = detector.size();
    SCITBX_ASSERT( cu_n_panels >= 1);

    //2) confirm that array dimensions are similar for each size
    cu_slow_pixels = detector[0].get_image_size()[1];
    cu_fast_pixels = detector[0].get_image_size()[0];
    for (int ipanel=1; ipanel < detector.size(); ++ipanel){
      SCITBX_ASSERT(detector[ipanel].get_image_size()[1] == cu_slow_pixels);
      SCITBX_ASSERT(detector[ipanel].get_image_size()[0] == cu_fast_pixels);
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
  };

  gpu_detector::gpu_detector(int const& arg_device_id,
                             dxtbx::model::Detector const & arg_detector,
                             dxtbx::model::Beam const& arg_beam):
    h_deviceID(arg_device_id),
    detector(arg_detector),
    cu_active_pixel_list(NULL),
    cu_accumulate_floatimage(NULL),
    metrology(arg_detector, arg_beam){
    construct_detail(arg_device_id, arg_detector);
  }

  gpu_detector::gpu_detector(int const& arg_device_id,
                             const simtbx::nanoBragg::nanoBragg& nB):
    h_deviceID(arg_device_id),
    metrology(nB),
    cu_active_pixel_list(NULL),
    cu_accumulate_floatimage(NULL){
    cudaSetDevice(arg_device_id);

    //1) determine the size
    cu_n_panels = 1;

    //2) confirm that array dimensions are similar for each size
    cu_slow_pixels = nB.spixels;
    cu_fast_pixels = nB.fpixels;
    _image_size = cu_n_panels * cu_slow_pixels * cu_fast_pixels;

    //3) allocate a cuda array with these dimensions
    /* separate accumulator image outside the usual nanoBragg data structure.
           1. accumulate contributions from a sequence of source energy channels computed separately
           2. represent multiple panels, all same rectangular shape; slowest dimension = n_panels */
    cudaSafeCall(cudaMalloc((void ** )&cu_accumulate_floatimage,
                            sizeof(*cu_accumulate_floatimage) * _image_size));
    cudaSafeCall(cudaMemset((void *)cu_accumulate_floatimage, 0,
                            sizeof(*cu_accumulate_floatimage) * _image_size));
  }

  void gpu_detector::free_detail(){
    cudaSetDevice(h_deviceID);
    //4) make sure we can deallocate cuda array later on
    if (cu_accumulate_floatimage != NULL) {
      cudaSafeCall(cudaFree(cu_accumulate_floatimage));
    }
  };

  void
  gpu_detector::scale_in_place(const double& factor){
    cudaSafeCall(cudaSetDevice(h_deviceID));
    cudaDeviceProp deviceProps = { 0 };
    cudaSafeCall(cudaGetDeviceProperties(&deviceProps, h_deviceID));
  int smCount = deviceProps.multiProcessorCount;
  dim3 threadsPerBlock(THREADS_PER_BLOCK_X, THREADS_PER_BLOCK_Y);
  dim3 numBlocks(smCount * 8, 1);
  int total_pixels = _image_size;
  scale_array_CUDAKernel<<<numBlocks, threadsPerBlock>>>(
    factor, cu_accumulate_floatimage, total_pixels);
  }

  void
  gpu_detector::write_raw_pixels(simtbx::nanoBragg::nanoBragg& nB){
    //only implement the monolithic detector case, one panel
    SCITBX_ASSERT(nB.spixels == cu_slow_pixels);
    SCITBX_ASSERT(nB.fpixels == cu_fast_pixels);
    SCITBX_ASSERT(cu_n_panels == 1);
    /* nB.raw_pixels = af::flex_double(af::flex_grid<>(nB.spixels,nB.fpixels));
       do not reallocate CPU memory for the data write, as it is not needed
     */
    double * double_floatimage = nB.raw_pixels.begin();
    cudaSafeCall(cudaSetDevice(nB.device_Id));
    cudaSafeCall(cudaMemcpy(
     double_floatimage,
     cu_accumulate_floatimage,
     sizeof(*cu_accumulate_floatimage) * _image_size,
     cudaMemcpyDeviceToHost));
  }

  af::flex_double
  gpu_detector::get_raw_pixels(){
    //return the data array for the multipanel detector case
    af::flex_double z(af::flex_grid<>(cu_n_panels,cu_slow_pixels,cu_fast_pixels), af::init_functor_null<double>());
    double* begin = z.begin();
    cudaSafeCall(cudaSetDevice(h_deviceID));
    cudaSafeCall(cudaMemcpy(
     begin,
     cu_accumulate_floatimage,
     sizeof(*cu_accumulate_floatimage) * _image_size,
     cudaMemcpyDeviceToHost));
    return z;
  }

  void
  gpu_detector::set_active_pixels_on_GPU(af::shared<std::size_t> active_pixel_list_value){
    active_pixel_list = active_pixel_list_value;
    cudaSafeCall(cudaSetDevice(h_deviceID));
    std::size_t * ptr_active_pixel_list = active_pixel_list.begin();
    cudaSafeCall(cudaMalloc((void ** )&cu_active_pixel_list, sizeof(*cu_active_pixel_list) * active_pixel_list.size() ));
    cudaSafeCall(cudaMemcpy(cu_active_pixel_list,
                            ptr_active_pixel_list,
                            sizeof(*cu_active_pixel_list) * active_pixel_list.size(),
                            cudaMemcpyHostToDevice));
  }

  af::shared<double>
  gpu_detector::get_whitelist_raw_pixels(af::shared<std::size_t> selection
  ){
    //return the data array for the multipanel detector case, but only for whitelist pixels
    af::shared<double> z(selection.size(), af::init_functor_null<double>());
    double* begin = z.begin();
    cudaSafeCall(cudaSetDevice(h_deviceID));
    CUDAREAL * cu_active_pixel_results;
    std::size_t * cu_active_pixel_selection;

    cudaSafeCall(cudaMalloc((void ** )&cu_active_pixel_results, sizeof(*cu_active_pixel_results) * z.size() ));
    cudaSafeCall(cudaMalloc((void ** )&cu_active_pixel_selection, sizeof(*cu_active_pixel_selection) * z.size() ));
    cudaSafeCall(cudaMemcpy(cu_active_pixel_selection,
                 selection.begin(), sizeof(*cu_active_pixel_selection) * z.size(),
                 cudaMemcpyHostToDevice));

    cudaDeviceProp deviceProps = { 0 };
    cudaSafeCall(cudaGetDeviceProperties(&deviceProps, h_deviceID));
    int smCount = deviceProps.multiProcessorCount;
    dim3 threadsPerBlock(THREADS_PER_BLOCK_X, THREADS_PER_BLOCK_Y);
    dim3 numBlocks(smCount * 8, 1);
    int total_pixels = z.size();
    get_active_pixel_selection_CUDAKernel<<<numBlocks, threadsPerBlock>>>(
      cu_active_pixel_results, cu_active_pixel_selection, cu_accumulate_floatimage, z.size());

    cudaSafeCall(cudaMemcpy(
      begin,
      cu_active_pixel_results,
      sizeof(*cu_active_pixel_results) * z.size(),
      cudaMemcpyDeviceToHost));
    cudaSafeCall(cudaFree(cu_active_pixel_selection));
    cudaSafeCall(cudaFree(cu_active_pixel_results));
    return z;
  }

  void
  gpu_detector::each_image_allocate(){
    cudaSetDevice(h_deviceID);
    /*allocate but do not initialize (set to 0) the reductions (the code was too inefficient and was removed as the reductions
      are not utilized in practice.  Should they be needed in the future a faster zeroing API must be found*/
    cu_omega_reduction = NULL;
    cudaSafeCall(cudaMalloc((void ** )&cu_omega_reduction, sizeof(*cu_omega_reduction) * _image_size));

    cu_max_I_x_reduction = NULL;
    cudaSafeCall(cudaMalloc((void ** )&cu_max_I_x_reduction, sizeof(*cu_max_I_x_reduction) * _image_size));

    cu_max_I_y_reduction = NULL;
    cudaSafeCall(cudaMalloc((void ** )&cu_max_I_y_reduction, sizeof(*cu_max_I_y_reduction) * _image_size));

    cu_rangemap = NULL;
    cudaSafeCall(cudaMalloc((void ** )&cu_rangemap, sizeof(*cu_rangemap) * _image_size));

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

        const int met_length = metrology.sdet.size();
        cudaSafeCall(cudaMalloc((void ** )&cu_sdet_vector, sizeof(*cu_sdet_vector) * met_length));
        cudaSafeCall(detMemcpyVectorDoubleToDevice(cu_sdet_vector, metrology.sdet.begin(), met_length));

        cudaSafeCall(cudaMalloc((void ** )&cu_fdet_vector, sizeof(*cu_fdet_vector) * met_length));
        cudaSafeCall(detMemcpyVectorDoubleToDevice(cu_fdet_vector, metrology.fdet.begin(), met_length));

        cudaSafeCall(cudaMalloc((void ** )&cu_odet_vector, sizeof(*cu_odet_vector) * met_length));
        cudaSafeCall(detMemcpyVectorDoubleToDevice(cu_odet_vector, metrology.odet.begin(), met_length));

        cudaSafeCall(cudaMalloc((void ** )&cu_pix0_vector, sizeof(*cu_pix0_vector) * met_length));
        cudaSafeCall(detMemcpyVectorDoubleToDevice(cu_pix0_vector, metrology.pix0.begin(), met_length));

        cudaSafeCall(cudaMalloc((void ** )&cu_distance, sizeof(*cu_distance) * metrology.dists.size()));
        cudaSafeCall(detMemcpyVectorDoubleToDevice(cu_distance, metrology.dists.begin(), metrology.dists.size()));

        cudaSafeCall(cudaMalloc((void ** )&cu_Xbeam,    sizeof(*cu_Xbeam) * metrology.Xbeam.size()));
        cudaSafeCall(detMemcpyVectorDoubleToDevice(cu_Xbeam,    metrology.Xbeam.begin(), metrology.Xbeam.size()));

        cudaSafeCall(cudaMalloc((void ** )&cu_Ybeam,    sizeof(*cu_Ybeam) * metrology.Ybeam.size()));
        cudaSafeCall(detMemcpyVectorDoubleToDevice(cu_Ybeam,    metrology.Ybeam.begin(), metrology.Ybeam.size()));
  }

  void
  gpu_detector::each_image_free(){
    cudaSetDevice(h_deviceID);
    cudaSafeCall(cudaDeviceSynchronize());
    cudaSafeCall(cudaFree(cu_omega_reduction));
    cudaSafeCall(cudaFree(cu_max_I_x_reduction));
    cudaSafeCall(cudaFree(cu_max_I_y_reduction));
    cudaSafeCall(cudaFree(cu_rangemap));
    cudaSafeCall(cudaFree(cu_maskimage));
    cudaSafeCall(cudaFree(cu_floatimage));
    cudaSafeCall(cudaFree(cu_sdet_vector));
    cudaSafeCall(cudaFree(cu_fdet_vector));
    cudaSafeCall(cudaFree(cu_odet_vector));
    cudaSafeCall(cudaFree(cu_pix0_vector));
    cudaSafeCall(cudaFree(cu_distance));
    cudaSafeCall(cudaFree(cu_Xbeam));
    cudaSafeCall(cudaFree(cu_Ybeam));
    cudaSafeCall(cudaFree(cu_active_pixel_list));
  }

} // gpu
} // simtbx
