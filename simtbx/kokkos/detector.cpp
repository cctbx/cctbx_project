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

} // Kokkos
} // simtbx
