#include <scitbx/array_family/boost_python/flex_fwd.h>
#include <cudatbx/cuda_base.cuh>
#include <simtbx/gpu/simulation.h>
#include <simtbx/gpu/simulation.cuh>
#include <scitbx/array_family/flex_types.h>
#define THREADS_PER_BLOCK_X 128
#define THREADS_PER_BLOCK_Y 1
#define THREADS_PER_BLOCK_TOTAL (THREADS_PER_BLOCK_X * THREADS_PER_BLOCK_Y)

namespace simtbx {
namespace gpu {

namespace af = scitbx::af;
//refactor later into helper file
  static cudaError_t cudaMemcpyVectorDoubleToDevice(CUDAREAL *dst, const double *src, size_t vector_items) {
        CUDAREAL * temp = new CUDAREAL[vector_items];
        for (size_t i = 0; i < vector_items; i++) {
                temp[i] = src[i];
        }
        cudaError_t ret = cudaMemcpy(dst, temp, sizeof(*dst) * vector_items, cudaMemcpyHostToDevice);
        delete temp;
        return ret;
  }

/* make a unit vector pointing in same direction and report magnitude (both args can be same vector) */
  double cpu_unitize(const double * vector, double * new_unit_vector) {

        double v1 = vector[1];
        double v2 = vector[2];
        double v3 = vector[3];

        double mag = sqrt(v1 * v1 + v2 * v2 + v3 * v3);

        if (mag != 0.0) {
                /* normalize it */
                new_unit_vector[0] = mag;
                new_unit_vector[1] = v1 / mag;
                new_unit_vector[2] = v2 / mag;
                new_unit_vector[3] = v3 / mag;
        } else {
                /* can't normalize, report zero vector */
                new_unit_vector[0] = 0.0;
                new_unit_vector[1] = 0.0;
                new_unit_vector[2] = 0.0;
                new_unit_vector[3] = 0.0;
        }
        return mag;
  }

  void
  exascale_api::show(){
    SCITBX_EXAMINE(SIM.roi_xmin);
    SCITBX_EXAMINE(SIM.roi_xmax);
    SCITBX_EXAMINE(SIM.roi_ymin);
    SCITBX_EXAMINE(SIM.roi_ymax);
    SCITBX_EXAMINE(SIM.oversample);
    SCITBX_EXAMINE(SIM.point_pixel);
    SCITBX_EXAMINE(SIM.pixel_size);
    SCITBX_EXAMINE(cu_subpixel_size);
    SCITBX_EXAMINE(cu_steps);
    SCITBX_EXAMINE(SIM.detector_thickstep);
    SCITBX_EXAMINE(SIM.detector_thicksteps);
    SCITBX_EXAMINE(SIM.detector_thick);
    SCITBX_EXAMINE(SIM.detector_attnlen);
    SCITBX_EXAMINE(SIM.curved_detector);
    SCITBX_EXAMINE(SIM.distance);
    SCITBX_EXAMINE(SIM.close_distance);
    SCITBX_EXAMINE(SIM.dmin);
    SCITBX_EXAMINE(SIM.phi0);
    SCITBX_EXAMINE(SIM.phistep);
    SCITBX_EXAMINE(SIM.phisteps);
    SCITBX_EXAMINE(SIM.sources);
    SCITBX_EXAMINE(SIM.mosaic_spread);
    SCITBX_EXAMINE(SIM.mosaic_domains);
    SCITBX_EXAMINE(SIM.Na);
    SCITBX_EXAMINE(SIM.Nb);
    SCITBX_EXAMINE(SIM.Nc);
    SCITBX_EXAMINE(SIM.fluence);
    SCITBX_EXAMINE(SIM.spot_scale);
    SCITBX_EXAMINE(SIM.integral_form);
    SCITBX_EXAMINE(SIM.default_F);
    SCITBX_EXAMINE(SIM.interpolate);
    SCITBX_EXAMINE(SIM.nopolar);
    SCITBX_EXAMINE(SIM.polarization);
    SCITBX_EXAMINE(SIM.fudge);
  }

  void
  exascale_api::add_energy_channel_from_gpu_amplitudes(
    int const& ichannel,
    simtbx::gpu::gpu_energy_channels & gec,
    simtbx::gpu::gpu_detector & gdt,
    double const& weight
  ){
        SCITBX_ASSERT(SIM.device_Id == stash_device_Id);
        cudaSafeCall(cudaSetDevice(stash_device_Id));

        // transfer source_I, source_lambda
        // the int arguments are for sizes of the arrays
        int source_count = SIM.sources;
        af::shared<double> weighted_sources_I = af::shared<double>(source_count);
        double* wptr = weighted_sources_I.begin();
        for (std::size_t iwt = 0; iwt < source_count; iwt++){wptr[iwt] = weight*(SIM.source_I[iwt]);}
        cudaSafeCall(cudaMemcpyVectorDoubleToDevice(cu_source_I, wptr, SIM.sources));
        cudaSafeCall(cudaMemcpyVectorDoubleToDevice(cu_source_lambda, SIM.source_lambda, SIM.sources));

        // magic happens here: take pointer from singleton, temporarily use it for add Bragg iteration:
        cu_current_channel_Fhkl = gec.d_channel_Fhkl[ichannel];

        cudaDeviceProp deviceProps = { 0 };
        cudaSafeCall(cudaGetDeviceProperties(&deviceProps, stash_device_Id));
        int smCount = deviceProps.multiProcessorCount;
        dim3 threadsPerBlock(THREADS_PER_BLOCK_X, THREADS_PER_BLOCK_Y);
        dim3 numBlocks(smCount * 8, 1);

        std::size_t panel_size = gdt.cu_slow_pixels * gdt.cu_fast_pixels;
        const int vec_len = 4;

        // the for loop around panels.  Offsets given.
        for (std::size_t idx_p = 0; idx_p < gdt.cu_n_panels; idx_p++){
          // loop thru panels and increment the array ptrs
          nanoBraggSpotsCUDAKernel<<<numBlocks, threadsPerBlock>>>(
          gdt.cu_slow_pixels, gdt.cu_fast_pixels, SIM.roi_xmin,
          SIM.roi_xmax, SIM.roi_ymin, SIM.roi_ymax, SIM.oversample, SIM.point_pixel,
          SIM.pixel_size, cu_subpixel_size, cu_steps, SIM.detector_thickstep, SIM.detector_thicksteps,
          SIM.detector_thick, SIM.detector_attnlen,
          &(gdt.cu_sdet_vector[vec_len * idx_p]),
          &(gdt.cu_fdet_vector[vec_len * idx_p]),
          &(gdt.cu_odet_vector[vec_len * idx_p]),
          &(gdt.cu_pix0_vector[vec_len * idx_p]),
          SIM.curved_detector, gdt.metrology.dists[idx_p], gdt.metrology.dists[idx_p], cu_beam_vector,
          gdt.metrology.Xbeam[idx_p], gdt.metrology.Ybeam[idx_p],
          SIM.dmin, SIM.phi0, SIM.phistep, SIM.phisteps, cu_spindle_vector,
          SIM.sources, cu_source_X, cu_source_Y, cu_source_Z,
          cu_source_I, cu_source_lambda, cu_a0, cu_b0,
          cu_c0, SIM.xtal_shape, SIM.mosaic_spread, SIM.mosaic_domains, cu_mosaic_umats,
          SIM.Na, SIM.Nb, SIM.Nc, SIM.V_cell,
          cu_water_size, cu_water_F, cu_water_MW, simtbx::nanoBragg::r_e_sqr, SIM.fluence,
          simtbx::nanoBragg::Avogadro, SIM.spot_scale, SIM.integral_form, SIM.default_F,
          SIM.interpolate, cu_current_channel_Fhkl, gec.cu_FhklParams, SIM.nopolar,
          cu_polar_vector, SIM.polarization, SIM.fudge,
          /* &(gdt.cu_maskimage[panel_size * idx_p]), */
          NULL,
          &(gdt.cu_floatimage[panel_size * idx_p]) /*out*/,
          &(gdt.cu_omega_reduction[panel_size * idx_p]) /*out*/,
          &(gdt.cu_max_I_x_reduction[panel_size * idx_p]) /*out*/,
          &(gdt.cu_max_I_y_reduction[panel_size * idx_p]) /*out*/,
          &(gdt.cu_rangemap[panel_size * idx_p]) /*out*/);

          cudaSafeCall(cudaPeekAtLastError());
        }
        cudaSafeCall(cudaDeviceSynchronize());

        //don't want to free the gec data when the nanoBragg goes out of scope, so switch the pointer
        cu_current_channel_Fhkl = NULL;

        add_array_CUDAKernel<<<numBlocks, threadsPerBlock>>>(gdt.cu_accumulate_floatimage,
          gdt.cu_floatimage,
          gdt.cu_n_panels * gdt.cu_slow_pixels * gdt.cu_fast_pixels);
  }

  void
  exascale_api::add_energy_channel_mask_allpanel(
    int const& ichannel,
    simtbx::gpu::gpu_energy_channels & gec,
    simtbx::gpu::gpu_detector & gdt,
    af::shared<bool> all_panel_mask
  ){
        // here or there, need to convert the all_panel_mask (3D map) into a 1D list of accepted pixels
        // coordinates for the active pixel list are absolute offsets into the detector array
        af::shared<std::size_t> active_pixel_list;
        const bool* jptr = all_panel_mask.begin();
        for (int j=0; j < all_panel_mask.size(); ++j){
          if (jptr[j]) {
            active_pixel_list.push_back(j);
          }
        }
        add_energy_channel_mask_allpanel(
          ichannel, gec, gdt, active_pixel_list);
  }

  void
  exascale_api::add_energy_channel_mask_allpanel(
    int const& ichannel,
    simtbx::gpu::gpu_energy_channels & gec,
    simtbx::gpu::gpu_detector & gdt,
    af::shared<std::size_t> const active_pixel_list
  ){
        SCITBX_ASSERT(SIM.device_Id == stash_device_Id);
        cudaSafeCall(cudaSetDevice(stash_device_Id));

        gdt.set_active_pixels_on_GPU(active_pixel_list);

        // transfer source_I, source_lambda
        // the int arguments are for sizes of the arrays
        cudaSafeCall(cudaMemcpyVectorDoubleToDevice(cu_source_I, SIM.source_I, SIM.sources));
        cudaSafeCall(cudaMemcpyVectorDoubleToDevice(cu_source_lambda, SIM.source_lambda, SIM.sources));

        // magic happens here: take pointer from singleton, temporarily use it for add Bragg iteration:
        cu_current_channel_Fhkl = gec.d_channel_Fhkl[ichannel];

        cudaDeviceProp deviceProps = { 0 };
        cudaSafeCall(cudaGetDeviceProperties(&deviceProps, stash_device_Id));
        int smCount = deviceProps.multiProcessorCount;
        dim3 threadsPerBlock(THREADS_PER_BLOCK_X, THREADS_PER_BLOCK_Y);
        dim3 numBlocks(smCount * 8, 1);

        const int vec_len = 4;
        // for call for all panels at the same time

          debranch_maskall_CUDAKernel<<<numBlocks, threadsPerBlock>>>(
          gdt.cu_n_panels, gdt.cu_slow_pixels, gdt.cu_fast_pixels, active_pixel_list.size(),
          SIM.oversample, SIM.point_pixel,
          SIM.pixel_size, cu_subpixel_size, cu_steps,
          SIM.detector_thickstep, SIM.detector_thicksteps,
          SIM.detector_thick, SIM.detector_attnlen,
          vec_len,
          gdt.cu_sdet_vector,
          gdt.cu_fdet_vector,
          gdt.cu_odet_vector,
          gdt.cu_pix0_vector,
          gdt.cu_distance, gdt.cu_distance, cu_beam_vector,
          gdt.cu_Xbeam, gdt.cu_Ybeam,
          SIM.dmin, SIM.phi0, SIM.phistep, SIM.phisteps, cu_spindle_vector,
          SIM.sources, cu_source_X, cu_source_Y, cu_source_Z,
          cu_source_I, cu_source_lambda, cu_a0, cu_b0,
          cu_c0, SIM.xtal_shape, SIM.mosaic_domains, cu_mosaic_umats,
          SIM.Na, SIM.Nb, SIM.Nc, SIM.V_cell,
          cu_water_size, cu_water_F, cu_water_MW, simtbx::nanoBragg::r_e_sqr, SIM.fluence,
          simtbx::nanoBragg::Avogadro, SIM.spot_scale, SIM.integral_form, SIM.default_F,
          cu_current_channel_Fhkl, gec.cu_FhklParams, SIM.nopolar,
          cu_polar_vector, SIM.polarization, SIM.fudge,
          gdt.cu_active_pixel_list,
          gdt.cu_floatimage /*out*/,
          gdt.cu_omega_reduction /*out*/,
          gdt.cu_max_I_x_reduction /*out*/,
          gdt.cu_max_I_y_reduction /*out*/,
          gdt.cu_rangemap /*out*/);

          cudaSafeCall(cudaPeekAtLastError());
        cudaSafeCall(cudaDeviceSynchronize());

        //don't want to free the gec data when the nanoBragg goes out of scope, so switch the pointer
        cu_current_channel_Fhkl = NULL;

        add_array_CUDAKernel<<<numBlocks, threadsPerBlock>>>(gdt.cu_accumulate_floatimage,
          gdt.cu_floatimage,
          gdt.cu_n_panels * gdt.cu_slow_pixels * gdt.cu_fast_pixels);
  }


  void
  exascale_api::add_background(simtbx::gpu::gpu_detector & gdt, int const& override_source){
        SCITBX_ASSERT(SIM.device_Id == stash_device_Id);
        cudaSafeCall(cudaSetDevice(stash_device_Id));

        // transfer source_I, source_lambda
        // the int arguments are for sizes of the arrays
        cudaSafeCall(cudaMemcpyVectorDoubleToDevice(cu_source_I, SIM.source_I, SIM.sources));
        cudaSafeCall(cudaMemcpyVectorDoubleToDevice(cu_source_lambda, SIM.source_lambda, SIM.sources));

        CUDAREAL * cu_stol_of;
        cudaSafeCall(cudaMalloc((void ** )&cu_stol_of, sizeof(*cu_stol_of) * SIM.stols));
        cudaSafeCall(cudaMemcpyVectorDoubleToDevice(cu_stol_of, SIM.stol_of, SIM.stols));

        CUDAREAL * cu_Fbg_of;
        cudaSafeCall(cudaMalloc((void ** )&cu_Fbg_of, sizeof(*cu_Fbg_of) * SIM.stols));
        cudaSafeCall(cudaMemcpyVectorDoubleToDevice(cu_Fbg_of, SIM.Fbg_of, SIM.stols));

        cudaDeviceProp deviceProps = { 0 };
        cudaSafeCall(cudaGetDeviceProperties(&deviceProps, stash_device_Id));
        int smCount = deviceProps.multiProcessorCount;
        dim3 threadsPerBlock(THREADS_PER_BLOCK_X, THREADS_PER_BLOCK_Y);
        dim3 numBlocks(smCount * 8, 1);

        //  initialize the device memory within a kernel.
        //  modify the arguments to initialize multipanel detector.
        nanoBraggSpotsInitCUDAKernel<<<numBlocks, threadsPerBlock>>>(
          gdt.cu_n_panels * gdt.cu_slow_pixels, gdt.cu_fast_pixels,
          gdt.cu_floatimage, gdt.cu_omega_reduction,
          gdt.cu_max_I_x_reduction, gdt.cu_max_I_y_reduction,
          gdt.cu_rangemap);
        cudaSafeCall(cudaPeekAtLastError());
        cudaSafeCall(cudaDeviceSynchronize());

        std::size_t panel_size = gdt.cu_slow_pixels * gdt.cu_fast_pixels;
        const int vec_len = 4;

        // the for loop around panels.  Offsets given.
        for (std::size_t idx_p = 0; idx_p < gdt.cu_n_panels; idx_p++){
          add_background_CUDAKernel<<<numBlocks, threadsPerBlock>>>(SIM.sources,
          SIM.oversample, override_source,
          SIM.pixel_size, gdt.cu_slow_pixels, gdt.cu_fast_pixels, SIM.detector_thicksteps,
          SIM.detector_thickstep, SIM.detector_attnlen,
          &(gdt.cu_sdet_vector[vec_len * idx_p]),
          &(gdt.cu_fdet_vector[vec_len * idx_p]),
          &(gdt.cu_odet_vector[vec_len * idx_p]),
          &(gdt.cu_pix0_vector[vec_len * idx_p]),
          gdt.metrology.dists[idx_p], SIM.point_pixel, SIM.detector_thick,
          cu_source_X, cu_source_Y, cu_source_Z,
          cu_source_lambda, cu_source_I,
          SIM.stols, cu_stol_of, cu_Fbg_of,
          SIM.nopolar, SIM.polarization, cu_polar_vector,
          simtbx::nanoBragg::r_e_sqr, SIM.fluence, SIM.amorphous_molecules,
          &(gdt.cu_floatimage[panel_size * idx_p]) /*out*/);

          cudaSafeCall(cudaPeekAtLastError());
        }
        cudaSafeCall(cudaDeviceSynchronize());
        add_array_CUDAKernel<<<numBlocks, threadsPerBlock>>>(gdt.cu_accumulate_floatimage,
          gdt.cu_floatimage,
          gdt.cu_n_panels * gdt.cu_slow_pixels * gdt.cu_fast_pixels);

        cudaSafeCall(cudaFree(cu_stol_of));
        cudaSafeCall(cudaFree(cu_Fbg_of));
}

  void
  exascale_api::allocate(){
    SCITBX_ASSERT(SIM.device_Id == stash_device_Id);
    cudaSafeCall(cudaSetDevice(stash_device_Id));

    /* water_size not defined in class, CLI argument, defaults to 0 */
    double water_size = 0.0;
    /* missing constants */
    double water_F = 2.57;
    double water_MW = 18.0;

    /* make sure we are normalizing with the right number of sub-steps */
    int nb_steps = SIM.phisteps*SIM.mosaic_domains*SIM.oversample*SIM.oversample;
    double nb_subpixel_size = SIM.pixel_size/SIM.oversample;

        /*create transfer arguments to device space*/
        cu_subpixel_size = nb_subpixel_size; //check for conflict?
        cu_steps = nb_steps; //check for conflict?

        /* presumably thickness and attenuation can be migrated to the gpu detector class XXX FIXME*/
        //cu_detector_thick = SIM.detector_thick;
        //cu_detector_mu = SIM.detector_attnlen; // synonyms
        //cu_distance = SIM.distance; /* distance and close distance, detector properties? XXX FIXME */
        //cu_close_distance = SIM.close_distance;

        cu_water_size = water_size;
        cu_water_F = water_F;
        cu_water_MW = water_MW;

        const int vector_length = 4;
        int cu_sources = SIM.sources;
        int cu_mosaic_domains = SIM.mosaic_domains;

        cudaSafeCall(cudaMalloc((void ** )&cu_beam_vector, sizeof(*cu_beam_vector) * vector_length));
        cudaSafeCall(cudaMemcpyVectorDoubleToDevice(cu_beam_vector, SIM.beam_vector, vector_length));

        cudaSafeCall(cudaMalloc((void ** )&cu_spindle_vector, sizeof(*cu_spindle_vector) * vector_length));
        cudaSafeCall(cudaMemcpyVectorDoubleToDevice(cu_spindle_vector, SIM.spindle_vector, vector_length));

        cudaSafeCall(cudaMalloc((void ** )&cu_a0, sizeof(*cu_a0) * vector_length));
        cudaSafeCall(cudaMemcpyVectorDoubleToDevice(cu_a0, SIM.a0, vector_length));

        cudaSafeCall(cudaMalloc((void ** )&cu_b0, sizeof(*cu_b0) * vector_length));
        cudaSafeCall(cudaMemcpyVectorDoubleToDevice(cu_b0, SIM.b0, vector_length));

        cudaSafeCall(cudaMalloc((void ** )&cu_c0, sizeof(*cu_c0) * vector_length));
        cudaSafeCall(cudaMemcpyVectorDoubleToDevice(cu_c0, SIM.c0, vector_length));

        // Unitize polar vector before sending it to the GPU.
        // Optimization do it only once here rather than multiple time per pixel in the GPU.
        double polar_vector_unitized[4];
        cpu_unitize(SIM.polar_vector, polar_vector_unitized);
        cudaSafeCall(cudaMalloc((void ** )&cu_polar_vector, sizeof(*cu_polar_vector) * vector_length));
        cudaSafeCall(cudaMemcpyVectorDoubleToDevice(cu_polar_vector, polar_vector_unitized, vector_length));

        cudaSafeCall(cudaMalloc((void ** )&cu_source_X, sizeof(*cu_source_X) * cu_sources));
        cudaSafeCall(cudaMemcpyVectorDoubleToDevice(cu_source_X, SIM.source_X, cu_sources));

        cudaSafeCall(cudaMalloc((void ** )&cu_source_Y, sizeof(*cu_source_Y) * cu_sources));
        cudaSafeCall(cudaMemcpyVectorDoubleToDevice(cu_source_Y, SIM.source_Y, cu_sources));

        cudaSafeCall(cudaMalloc((void ** )&cu_source_Z, sizeof(*cu_source_Z) * cu_sources));
        cudaSafeCall(cudaMemcpyVectorDoubleToDevice(cu_source_Z, SIM.source_Z, cu_sources));

        cudaSafeCall(cudaMalloc((void ** )&cu_source_I, sizeof(*cu_source_I) * cu_sources));
        cudaSafeCall(cudaMemcpyVectorDoubleToDevice(cu_source_I, SIM.source_I, cu_sources));

        cudaSafeCall(cudaMalloc((void ** )&cu_source_lambda, sizeof(*cu_source_lambda) * cu_sources));
        cudaSafeCall(cudaMemcpyVectorDoubleToDevice(cu_source_lambda, SIM.source_lambda, cu_sources));

        cudaSafeCall(cudaMalloc((void ** )&cu_mosaic_umats, sizeof(*cu_mosaic_umats) * cu_mosaic_domains * 9));
        cudaSafeCall(cudaMemcpyVectorDoubleToDevice(cu_mosaic_umats, SIM.mosaic_umats, cu_mosaic_domains * 9));
  };

  exascale_api::~exascale_api(){
    cudaSafeCall(cudaSetDevice(stash_device_Id));

        cudaSafeCall(cudaFree(cu_beam_vector));
        cudaSafeCall(cudaFree(cu_spindle_vector));
        cudaSafeCall(cudaFree(cu_source_X));
        cudaSafeCall(cudaFree(cu_source_Y));
        cudaSafeCall(cudaFree(cu_source_Z));
        cudaSafeCall(cudaFree(cu_source_I));
        cudaSafeCall(cudaFree(cu_source_lambda));
        cudaSafeCall(cudaFree(cu_a0));
        cudaSafeCall(cudaFree(cu_b0));
        cudaSafeCall(cudaFree(cu_c0));
        cudaSafeCall(cudaFree(cu_mosaic_umats));
        cudaSafeCall(cudaFree(cu_polar_vector));
  }

} // gpu
} // simtbx
