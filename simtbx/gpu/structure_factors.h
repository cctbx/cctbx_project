#ifndef SIMTBX_GPU_STRUCTURE_FACTORS_H
#define SIMTBX_GPU_STRUCTURE_FACTORS_H
// the *.h file can be included for all compiles regardless of the availability of CUDA

#include <scitbx/array_family/shared.h>
#include <cctbx/miller.h>
#include <simtbx/nanoBragg/nanotypes.h>

using simtbx::nanoBragg::hklParams;
namespace simtbx {
namespace gpu {

namespace af = scitbx::af;
typedef cctbx::miller::index<int> miller_t;
typedef af::shared<miller_t > indices;

class gpu_instance {
  public:
    inline
    gpu_instance(){
      printf("NO OPERATION, NO DEVICE NUMBER, gpu_instance");
    }

    inline
    gpu_instance(int const& t_deviceID) {
      deviceID = t_deviceID;
    }

    inline int get_deviceID() const {
      return deviceID;
    }

  private:
    int deviceID = -1;
};

struct gpu_energy_channels {
  gpu_energy_channels(){printf("NO OPERATION, NO DEVICE NUMBER");};
  gpu_energy_channels(int const&);

  inline
  void structure_factors_to_GPU_direct(
    int, indices pythony_indices, af::shared<double> pythony_amplitudes){
    double F_cell;

    /* load the structure factors */
    SCITBX_ASSERT(pythony_indices.size()>0);
    if(pythony_indices.size()) {
        if(verbose) printf(" noticed pythony_indices.size() = %ld\n",pythony_indices.size());
        /* need to know how much memory to allocate for Fhkl array */
        h_min=k_min=l_min=1e9;
        h_max=k_max=l_max=-1e9;
        miller_t hkl; int h0,k0,l0;
        for (int i_idx=0; i_idx < pythony_indices.size(); ++i_idx){
            hkl = pythony_indices[i_idx];
            h0 = hkl[0]; k0 = hkl[1]; l0 = hkl[2];
            if(h_min > h0) h_min = h0;
            if(k_min > k0) k_min = k0;
            if(l_min > l0) l_min = l0;
            if(h_max < h0) h_max = h0;
            if(k_max < k0) k_max = k0;
            if(l_max < l0) l_max = l0;
        }
        h_range = h_max - h_min + 1;
        k_range = k_max - k_min + 1;
        l_range = l_max - l_min + 1;
        if(verbose) printf("h: %d - %d\n",h_min,h_max);
        if(verbose) printf("k: %d - %d\n",k_min,k_max);
        if(verbose) printf("l: %d - %d\n",l_min,l_max);
        SCITBX_ASSERT (!(h_range < 0 || k_range < 0 || l_range < 0));
    }
    af::shared<double> linear_amplitudes;
    if(pythony_indices.size() && pythony_amplitudes.size()) {
      linear_amplitudes = af::shared<double>( h_range * k_range * l_range, default_F);
      double * raw_ptr = linear_amplitudes.begin();
      miller_t hkl; int h0,k0,l0;
      for (int i_idx=0; i_idx < pythony_indices.size(); ++i_idx) {
        hkl = pythony_indices[i_idx];
        F_cell = pythony_amplitudes[i_idx];
        h0 = hkl[0]; k0 = hkl[1]; l0 = hkl[2];
        raw_ptr[(h0-h_min) * k_range * l_range + (k0-k_min) * l_range + (l0-l_min)] = F_cell;
      }
      raw_ptr[-h_min * k_range * l_range -k_min * l_range -l_min] = F000;
    }
    structure_factors_to_GPU_detail(linear_amplitudes);
  }
  void structure_factors_to_GPU_detail(af::shared<double>);

  inline int get_deviceID(){return h_deviceID;}
  inline int get_nchannels(){return d_channel_Fhkl.size();}
  void free_detail();
  inline ~gpu_energy_channels(){ if (d_channel_Fhkl.size() > 0) {free_detail();} }

    double F000 = 0.0;        // to mark beam center
    double default_F = 0.0;
    int verbose = 0;
    int h_range,k_range,l_range,h_min,h_max,k_min,k_max,l_min,l_max;
    int h_deviceID;

  /* pointers to data on device */
  af::shared<double *> d_channel_Fhkl;
  hklParams * cu_FhklParams;
};
} // gpu
} // simtbx
#endif // SIMTBX_GPU_STRUCTURE_FACTORS_H
