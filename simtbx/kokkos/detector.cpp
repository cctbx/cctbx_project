#include "scitbx/array_family/boost_python/flex_fwd.h"
//#include "cudatbx/cuda_base.cuh"
#include "simtbx/kokkos/detector.h"
#include "scitbx/vec3.h"
#include "scitbx/vec2.h"

#include <unistd.h> // for sleep

#define THREADS_PER_BLOCK_X 128
#define THREADS_PER_BLOCK_Y 1
#define THREADS_PER_BLOCK_TOTAL (THREADS_PER_BLOCK_X * THREADS_PER_BLOCK_Y)

using Kokkos::fence;
using Kokkos::deep_copy;
using Kokkos::create_mirror_view;
using Kokkos::parallel_for;

namespace simtbx { namespace Kokkos {

  template <typename T>
  void
  transfer_vector2kokkos(view_1d_t<T> &dst, const af::shared<T>  &src) {
    if (true) {
      printf("== Transfer %s from %p\n", dst.label().c_str(), (void*) dst.data());
      printf(" - size src|dst: %d|%d\n", src.size(), dst.span() );
    }
    if (dst.span() < src.size()) {
      resize(dst, src.size());
      printf(" - size changed, new size: %d\n", dst.span() );
    }
    auto host_view = create_mirror_view(dst);

    for (int i=0; i<src.size(); ++i) {
      host_view( i ) = src[ i ];
    }
    deep_copy(dst, host_view);
    printf("copied from %p to %p.\n", (void*) host_view.data(), (void*) dst.data());

    print_view(dst);
  }

  packed_metrology::packed_metrology(dxtbx::model::Detector const & arg_detector,
                                   dxtbx::model::Beam const & arg_beam) {

    for (std::size_t panel_id = 0; panel_id < arg_detector.size(); panel_id++){
          // helper code arising from the nanoBragg constructor, with user_beam=True
      typedef scitbx::vec3<double> vec3;

      // DETECTOR properties
      // typically: 1 0 0
      vec3 fdet_vector = arg_detector[panel_id].get_fast_axis();
      fdet_vector = fdet_vector.normalize();

      // typically: 0 -1 0
      vec3 sdet_vector = arg_detector[panel_id].get_slow_axis();
      sdet_vector = sdet_vector.normalize();

      // set orthogonal vector to the detector pixel array
      vec3 odet_vector = fdet_vector.cross(sdet_vector);
      odet_vector = odet_vector.normalize();

      // dxtbx origin is location of outer corner of the first pixel
      vec3 pix0_vector = arg_detector[panel_id].get_origin()/1000.0;

      // what is the point of closest approach between sample and detector?
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
      // set beam centre
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

  vector_double_t
  kokkos_detector::construct_detail(dxtbx::model::Detector const & arg_detector) {
    //1) confirm the size
    SCITBX_ASSERT( m_panel_count == arg_detector.size() );
    SCITBX_ASSERT( m_panel_count >= 1 );

    //2) confirm that array dimensions are similar for each size
    for (int ipanel=1; ipanel < arg_detector.size(); ++ipanel){
      SCITBX_ASSERT( arg_detector[ipanel].get_image_size()[0] == m_slow_dim_size );
      SCITBX_ASSERT( arg_detector[ipanel].get_image_size()[1] == m_fast_dim_size );
    }
    printf(" m_total_pixel_count: %d\n", m_total_pixel_count);
    printf("     m_slow_dim_size: %d\n", m_slow_dim_size);
    printf("     m_fast_dim_size: %d\n", m_fast_dim_size);
    printf("       m_panel_count: %d\n", m_panel_count);

    //3) allocate a cuda array with these dimensions
    // separate accumulator image outside the usual nanoBragg data structure.
    //       1. accumulate contributions from a sequence of source energy channels computed separately
    //       2. represent multiple panels, all same rectangular shape; slowest dimension = n_panels
    vector_double_t view_floatimage( "m_accumulate_floatimage", m_total_pixel_count );
    return view_floatimage;
  };

  kokkos_detector::kokkos_detector(dxtbx::model::Detector const & arg_detector,
                             dxtbx::model::Beam const& arg_beam):
    detector(arg_detector),
    m_panel_count( arg_detector.size() ),
    m_slow_dim_size( arg_detector[0].get_image_size()[0] ),
    m_fast_dim_size( arg_detector[0].get_image_size()[1] ),
    m_total_pixel_count( m_panel_count * m_slow_dim_size * m_fast_dim_size ),
    m_accumulate_floatimage( construct_detail(arg_detector) ),
    cu_active_pixel_list(NULL),
    cu_accumulate_floatimage(NULL),
    metrology(arg_detector, arg_beam) { }

/*  kokkos_detector::kokkos_detector(const simtbx::nanoBragg::nanoBragg& nB):
    metrology(nB),
    cu_active_pixel_list(NULL),
    cu_accumulate_floatimage(NULL){

    //1) determine the size
    m_panel_count = 1;

    //2) confirm that array dimensions are similar for each size
    m_slow_dim_size = nB.spixels;
    m_fast_dim_size = nB.fpixels;
    m_total_pixel_count = m_panel_count * m_slow_dim_size * m_fast_dim_size;

    //3) allocate a cuda array with these dimensions
    // separate accumulator image outside the usual nanoBragg data structure.
    //     1. accumulate contributions from a sequence of source energy channels computed separately
    //     2. represent multiple panels, all same rectangular shape; slowest dimension = n_panels
    cudaSafeCall(cudaMalloc((void ** )&cu_accumulate_floatimage,
                            sizeof(*cu_accumulate_floatimage) * m_total_pixel_count));
    cudaSafeCall(cudaMemset((void *)cu_accumulate_floatimage, 0,
                            sizeof(*cu_accumulate_floatimage) * m_total_pixel_count));
  }
*/
  void
  kokkos_detector::scale_in_place_cuda(const double& factor){
    auto local_accumulate_floatimage = m_accumulate_floatimage;
    parallel_for("scale_in_place", range_policy(0,m_total_pixel_count), KOKKOS_LAMBDA (const int i) {
      local_accumulate_floatimage( i ) = local_accumulate_floatimage( i ) * factor;
    });
  }

  void
  kokkos_detector::write_raw_pixels_cuda(simtbx::nanoBragg::nanoBragg& nB){
    //only implement the monolithic detector case, one panel
    SCITBX_ASSERT(nB.spixels == m_slow_dim_size);
    SCITBX_ASSERT(nB.fpixels == m_fast_dim_size);
    SCITBX_ASSERT(m_panel_count == 1);
    // nB.raw_pixels = af::flex_double(af::flex_grid<>(nB.spixels,nB.fpixels));
    // do not reallocate CPU memory for the data write, as it is not needed
    
    vector_double_t::HostMirror host_floatimage = create_mirror_view(m_accumulate_floatimage);
    deep_copy(host_floatimage, m_accumulate_floatimage);

    printf(" m_total_pixel_count: %d\n", m_total_pixel_count);

    double * double_floatimage = nB.raw_pixels.begin();
    for (int i=0; i<m_total_pixel_count; ++i) {
      double_floatimage[i] = host_floatimage( i ); 
    }
  }

  af::flex_double
  kokkos_detector::get_raw_pixels_cuda(){
    //return the data array for the multipanel detector case
    af::flex_double host_array(af::flex_grid<>(m_panel_count,m_slow_dim_size,m_fast_dim_size), af::init_functor_null<double>());
    double* p_host_array = host_array.begin();
    
    vector_double_t::HostMirror host_floatimage = create_mirror_view(m_accumulate_floatimage);
    deep_copy(host_floatimage, m_accumulate_floatimage);

    for (int i=0; i<m_total_pixel_count; ++i) {
      p_host_array[i] = host_floatimage( i ); 
    }    
    return host_array;
  }

  void
  kokkos_detector::set_active_pixels_on_KOKKOS(af::shared<int> active_pixel_list_value) {
    m_active_pixel_size = active_pixel_list_value.size();
    resize(m_active_pixel_list, m_active_pixel_size);
    vector_int_t::HostMirror host_view = create_mirror_view(m_active_pixel_list);

    int* ptr_active_pixel_list = active_pixel_list.begin();
    for (int i=0; i<m_active_pixel_size; ++i) {
      host_view( i ) = ptr_active_pixel_list[ i ];
    }
    deep_copy(m_active_pixel_list, host_view);

    active_pixel_list = active_pixel_list_value;
    
    // cudaSafeCall(cudaMalloc((void ** )&cu_active_pixel_list, sizeof(*cu_active_pixel_list) * active_pixel_list.size() ));
    // cudaSafeCall(cudaMemcpy(cu_active_pixel_list,
    //                         ptr_active_pixel_list,
    //                         sizeof(*cu_active_pixel_list) * active_pixel_list.size(),
    //                         cudaMemcpyHostToDevice));
  }
/*
  af::shared<double>
  kokkos_detector::get_whitelist_raw_pixels_cuda(af::shared<std::size_t> selection
  ){
    //return the data array for the multipanel detector case, but only for whitelist pixels
    af::shared<double> z(active_pixel_list.size(), af::init_functor_null<double>());
    double* begin = z.begin();
    CUDAREAL * cu_active_pixel_results;
    std::size_t * cu_active_pixel_selection;

    cudaSafeCall(cudaMalloc((void ** )&cu_active_pixel_results, sizeof(*cu_active_pixel_results) * active_pixel_list.size() ));
    cudaSafeCall(cudaMalloc((void ** )&cu_active_pixel_selection, sizeof(*cu_active_pixel_selection) * selection.size() ));
    cudaSafeCall(cudaMemcpy(cu_active_pixel_selection,
                 selection.begin(), sizeof(*cu_active_pixel_selection) * selection.size(),
                 cudaMemcpyHostToDevice));

    int smCount = 84; //deviceProps.multiProcessorCount;
    dim3 threadsPerBlock(THREADS_PER_BLOCK_X, THREADS_PER_BLOCK_Y);
    dim3 numBlocks(smCount * 8, 1);
    int total_pixels = active_pixel_list.size();
    get_active_pixel_selection_CUDAKernel<<<numBlocks, threadsPerBlock>>>(
      cu_active_pixel_results, cu_active_pixel_selection, cu_accumulate_floatimage, total_pixels);

    cudaSafeCall(cudaMemcpy(
      begin,
      cu_active_pixel_results,
      sizeof(*cu_active_pixel_results) * active_pixel_list.size(),
      cudaMemcpyDeviceToHost));
    cudaSafeCall(cudaFree(cu_active_pixel_selection));
    cudaSafeCall(cudaFree(cu_active_pixel_results));
    return z;
  }
*/
  void
  kokkos_detector::each_image_allocate_cuda() {
    resize(m_rangemap, m_total_pixel_count);
    resize(m_omega_reduction, m_total_pixel_count);
    resize(m_max_I_x_reduction, m_total_pixel_count);
    resize(m_max_I_y_reduction, m_total_pixel_count);
    
    resize(m_maskimage, m_total_pixel_count);
    resize(m_floatimage, m_total_pixel_count);

    transfer_vector2kokkos(m_sdet_vector, metrology.sdet);
    transfer_vector2kokkos(m_fdet_vector, metrology.fdet);
    transfer_vector2kokkos(m_odet_vector, metrology.odet);
    transfer_vector2kokkos(m_pix0_vector, metrology.pix0);
    transfer_vector2kokkos(m_distance, metrology.dists);
    transfer_vector2kokkos(m_Xbeam, metrology.Xbeam);
    transfer_vector2kokkos(m_Ybeam, metrology.Ybeam);
    fence();

    metrology.show();
    
    printf(" rangemap size:%d\n", m_rangemap.span());
    printf(" omega_reduction size:%d\n", m_omega_reduction.span());
    printf(" max_I_x_reduction size:%d\n", m_max_I_x_reduction.span());
    printf(" max_I_y_reduction size:%d\n", m_max_I_y_reduction.span());
    printf(" maskimage size:%d\n", m_maskimage.span());
    printf(" floatimage size:%d\n", m_floatimage.span());
    printf(" sdet_vector size:%d\n", m_sdet_vector.span());
    printf(" fdet_vector size:%d\n", m_fdet_vector.span());
    printf(" odet_vector size:%d\n", m_odet_vector.span());
    printf(" pix0_vector size:%d\n", m_pix0_vector.span());
    printf(" distance size:%d\n", m_distance.span());
    printf(" Xbeam size:%d\n", m_Xbeam.span());
    printf(" Ybeam size:%d\n", m_Ybeam.span());
  
    print_view(m_fdet_vector);
    print_view(m_odet_vector, 1, 3);

    printf("DONE.\n");
  }

} // Kokkos
} // simtbx

