#include "scitbx/array_family/boost_python/flex_fwd.h"
//#include "cudatbx/cuda_base.cuh"
#include "simtbx/kokkos/detector.h"
#include "kokkostbx/kokkos_utils.h"
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

auto get_kokkos_vec3 = [](auto&& src) { return vec3(src[0], src[1], src[2]); };

namespace simtbx { namespace Kokkos {

  packed_metrology::packed_metrology(dxtbx::model::Detector const & arg_detector,
                                   dxtbx::model::Beam const & arg_beam) {

    for (std::size_t panel_id = 0; panel_id < arg_detector.size(); panel_id++){
      // helper code arising from the nanoBragg constructor, with user_beam=True

      // DETECTOR properties
      // typically: 1 0 0
      vec3 fdet_vector = get_kokkos_vec3(arg_detector[panel_id].get_fast_axis());
      fdet_vector.normalize();

      // typically: 0 -1 0
      vec3 sdet_vector = get_kokkos_vec3(arg_detector[panel_id].get_slow_axis());
      sdet_vector.normalize();

      // set orthogonal vector to the detector pixel array
      vec3 odet_vector = fdet_vector.cross(sdet_vector);
      odet_vector.normalize();

      // dxtbx origin is location of outer corner of the first pixel
      vec3 pix0_vector = get_kokkos_vec3(arg_detector[panel_id].get_origin()/1000.0);

      // what is the point of closest approach between sample and detector?
      double close_distance = pix0_vector.dot(odet_vector);
      if (close_distance < 0){
        bool verbose = false;
        if(verbose)printf("WARNING: dxtbx model is lefthanded. Inverting odet_vector.\n");
        odet_vector = -1. * odet_vector;
        close_distance = -1 * close_distance;
      }

      sdet.push_back(sdet_vector);
      fdet.push_back(fdet_vector);
      odet.push_back(odet_vector);
      pix0.push_back(pix0_vector);

      // set beam centre
      scitbx::vec2<double> dials_bc=arg_detector[panel_id].get_beam_centre(arg_beam.get_s0());
      dists.push_back(close_distance);
      Xbeam.push_back(dials_bc[0]/1000.0);
      Ybeam.push_back(dials_bc[1]/1000.0);
    }
  };

  packed_metrology::packed_metrology(const simtbx::nanoBragg::nanoBragg& nB){
    // Careful, 4-vectors! [length, x, y, z]
    auto get_kokkos_vec3 = [](auto& src) { return vec3(src[1], src[2], src[3]); };
    sdet.push_back( get_kokkos_vec3(nB.sdet_vector) );
    fdet.push_back( get_kokkos_vec3(nB.fdet_vector) );
    odet.push_back( get_kokkos_vec3(nB.odet_vector) );
    pix0.push_back( get_kokkos_vec3(nB.pix0_vector) );
    dists.push_back(nB.close_distance);
    Xbeam.push_back(nB.Xbeam);
    Ybeam.push_back(nB.Ybeam);
  }

  void
  packed_metrology::show() const {
    for (std::size_t idx_p = 0; idx_p < Xbeam.size(); idx_p++){
      printf(" Panel %3ld\n",idx_p);
      printf(" Panel %3ld sdet %9.6f %9.6f %9.6f fdet %9.6f %9.6f %9.6f\n",
             idx_p,sdet[idx_p][0],sdet[idx_p][1],sdet[idx_p][2],
                   fdet[idx_p][0],fdet[idx_p][1],fdet[idx_p][2]
      );
      printf(" Panel %3ld odet %9.6f %9.6f %9.6f pix0 %9.6f %9.6f %9.6f\n",
             idx_p,odet[idx_p][0],odet[idx_p][1],odet[idx_p][2],
                   pix0[idx_p][0],pix0[idx_p][1],pix0[idx_p][2]
      );
      printf(" Panel %3ld beam %11.8f %11.8f\n",idx_p,Xbeam[idx_p],Ybeam[idx_p]);
    }
  }

  vector_double_t
  kokkos_detector::construct_detail(dxtbx::model::Detector const & arg_detector) {
    //1) confirm the size
    SCITBX_ASSERT( m_panel_count == arg_detector.size() );
    SCITBX_ASSERT( m_panel_count >= 1 );

    //2) confirm that array dimensions are similar for each size
    for (int ipanel=1; ipanel < arg_detector.size(); ++ipanel){
      SCITBX_ASSERT( arg_detector[ipanel].get_image_size()[1] == m_slow_dim_size );
      SCITBX_ASSERT( arg_detector[ipanel].get_image_size()[0] == m_fast_dim_size );
    }
    // printf(" m_total_pixel_count: %d\n", m_total_pixel_count);
    // printf("     m_slow_dim_size: %d\n", m_slow_dim_size);
    // printf("     m_fast_dim_size: %d\n", m_fast_dim_size);
    // printf("       m_panel_count: %d\n", m_panel_count);

    //3) allocate a cuda array with these dimensions
    // separate accumulator image outside the usual nanoBragg data structure.
    //       1. accumulate contributions from a sequence of source energy channels computed separately
    //       2. represent multiple panels, all same rectangular shape; slowest dimension = n_panels
    vector_double_t view_floatimage( "m_accumulate_floatimage", m_total_pixel_count );
    return view_floatimage;
  };

  kokkos_detector::kokkos_detector(int const& arg_device,
                             dxtbx::model::Detector const & arg_detector,
                             dxtbx::model::Beam const& arg_beam):
    h_deviceID(arg_device),
    metrology(arg_detector, arg_beam),
    m_panel_count( arg_detector.size() ),
    m_slow_dim_size( arg_detector[0].get_image_size()[1] ),
    m_fast_dim_size( arg_detector[0].get_image_size()[0] ),
    m_total_pixel_count( m_panel_count * m_slow_dim_size * m_fast_dim_size ),
    m_accumulate_floatimage( construct_detail(arg_detector) ) { }
    // Easy mistake: not realizing that the dxtbx detector model stores (fast,slow) sizes

  kokkos_detector::kokkos_detector(int const& arg_device,
                             const simtbx::nanoBragg::nanoBragg& nB):
    h_deviceID(arg_device),
    metrology(nB),
    m_panel_count(1),
    m_slow_dim_size(nB.spixels),
    m_fast_dim_size(nB.fpixels),
    m_total_pixel_count( m_panel_count * m_slow_dim_size * m_fast_dim_size ),
    m_accumulate_floatimage( vector_double_t( "m_accumulate_floatimage", m_total_pixel_count) ) { }

  void
  kokkos_detector::scale_in_place(const double& factor){
    auto local_accumulate_floatimage = m_accumulate_floatimage;
    parallel_for("scale_in_place", range_policy(0,m_total_pixel_count), KOKKOS_LAMBDA (const int i) {
      local_accumulate_floatimage( i ) = local_accumulate_floatimage( i ) * factor;
    });
  }

  void
  kokkos_detector::write_raw_pixels(simtbx::nanoBragg::nanoBragg& nB) {
    //only implement the monolithic detector case, one panel
    SCITBX_ASSERT(nB.spixels == m_slow_dim_size);
    SCITBX_ASSERT(nB.fpixels == m_fast_dim_size);
    SCITBX_ASSERT(m_panel_count == 1);
    // nB.raw_pixels = af::flex_double(af::flex_grid<>(nB.spixels,nB.fpixels));
    // do not reallocate CPU memory for the data write, as it is not needed

    kokkostbx::transfer_kokkos2flex(nB.raw_pixels, m_accumulate_floatimage);
    // vector_double_t::HostMirror host_floatimage = create_mirror_view(m_accumulate_floatimage);
    // deep_copy(host_floatimage, m_accumulate_floatimage);

    // printf(" m_total_pixel_count: %d\n", m_total_pixel_count);

    // double * double_floatimage = nB.raw_pixels.begin();
    // for (int i=0; i<m_total_pixel_count; ++i) {
    //   double_floatimage[i] = host_floatimage( i );
    // }
  }

  af::flex_double
  kokkos_detector::get_raw_pixels(){
    //return the data array for the multipanel detector case
    af::flex_double output_array(af::flex_grid<>(m_panel_count,m_slow_dim_size,m_fast_dim_size), af::init_functor_null<double>());
    kokkostbx::transfer_kokkos2flex(output_array, m_accumulate_floatimage);

    // vector_double_t::HostMirror host_floatimage = create_mirror_view(m_accumulate_floatimage);
    // deep_copy(host_floatimage, m_accumulate_floatimage);

    // for (int i=0; i<m_total_pixel_count; ++i) {
    //   output_array_ptr[ i ] = host_floatimage( i );
    // }
    return output_array;
  }

  void
  kokkos_detector::set_active_pixels_on_GPU(af::shared<std::size_t> active_pixel_list_value) {
    m_active_pixel_size = active_pixel_list_value.size();
    kokkostbx::transfer_shared2kokkos(m_active_pixel_list, active_pixel_list_value);
    active_pixel_list = active_pixel_list_value;
  }

  af::shared<double>
  kokkos_detector::get_whitelist_raw_pixels(af::shared<std::size_t> selection) {
    //return the data array for the multipanel detector case, but only for whitelist pixels
    vector_size_t active_pixel_selection = vector_size_t("active_pixel_selection", selection.size());
    kokkostbx::transfer_shared2kokkos(active_pixel_selection, selection);

    size_t output_pixel_size = selection.size();
    vector_cudareal_t active_pixel_results = vector_cudareal_t("active_pixel_results", output_pixel_size);

    auto temp = m_accumulate_floatimage;

    parallel_for("get_active_pixel_selection",
                  range_policy(0, output_pixel_size),
                  KOKKOS_LAMBDA (const int i) {
      size_t index = active_pixel_selection( i );
      active_pixel_results( i ) = temp( index );
    });

    af::shared<double> output_array(output_pixel_size, af::init_functor_null<double>());
    kokkostbx::transfer_kokkos2shared(output_array, active_pixel_results);

    SCITBX_ASSERT(output_array.size() == output_pixel_size);
    return output_array;
  }

  void
  kokkos_detector::each_image_allocate() {
    resize(m_rangemap, m_total_pixel_count);
    resize(m_omega_reduction, m_total_pixel_count);
    resize(m_max_I_x_reduction, m_total_pixel_count);
    resize(m_max_I_y_reduction, m_total_pixel_count);

    resize(m_maskimage, m_total_pixel_count);
    resize(m_floatimage, m_total_pixel_count);

    kokkostbx::transfer_shared2kokkos(m_sdet_vector, metrology.sdet);
    kokkostbx::transfer_shared2kokkos(m_fdet_vector, metrology.fdet);
    kokkostbx::transfer_shared2kokkos(m_odet_vector, metrology.odet);
    kokkostbx::transfer_shared2kokkos(m_pix0_vector, metrology.pix0);
    kokkostbx::transfer_shared2kokkos(m_distance, metrology.dists);
    kokkostbx::transfer_shared2kokkos(m_Xbeam, metrology.Xbeam);
    kokkostbx::transfer_shared2kokkos(m_Ybeam, metrology.Ybeam);
    fence();

    // metrology.show();

    // printf(" rangemap size:%d\n", m_rangemap.span());
    // printf(" omega_reduction size:%d\n", m_omega_reduction.span());
    // printf(" max_I_x_reduction size:%d\n", m_max_I_x_reduction.span());
    // printf(" max_I_y_reduction size:%d\n", m_max_I_y_reduction.span());
    // printf(" maskimage size:%d\n", m_maskimage.span());
    // printf(" floatimage size:%d\n", m_floatimage.span());
    // printf(" sdet_vector size:%d\n", m_sdet_vector.span());
    // printf(" fdet_vector size:%d\n", m_fdet_vector.span());
    // printf(" odet_vector size:%d\n", m_odet_vector.span());
    // printf(" pix0_vector size:%d\n", m_pix0_vector.span());
    // printf(" distance size:%d\n", m_distance.span());
    // printf(" Xbeam size:%d\n", m_Xbeam.span());
    // printf(" Ybeam size:%d\n", m_Ybeam.span());

    // print_view(m_fdet_vector);
    // print_view(m_odet_vector, 1, 3);

    // printf("DONE.\n");
  }

} // Kokkos
} // simtbx
