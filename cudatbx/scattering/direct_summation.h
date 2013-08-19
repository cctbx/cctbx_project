#ifndef DIRECT_SUMMATION_H
#define DIRECT_SUMMATION_H

// std::system includes
#include <cmath>
#include <ctime>

#include <complex>

#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include <math.h>
#include <boost/optional.hpp>
#include <cctbx/eltbx/xray_scattering.h>
#include <cctbx/xray/scattering_type_registry.h>
#include <scitbx/constants.h>
#include <scitbx/vec3.h>
#include <scitbx/array_family/shared.h>

// set precision for calculations on GPU
#define fType float

namespace cudatbx {
namespace scattering {

  class direct_summation {

  public:
    direct_summation();
    ~direct_summation();
    void add(const scitbx::af::const_ref<std::string>&,
             const scitbx::af::const_ref<scitbx::vec3<double> >&,
             const scitbx::af::const_ref<double>&,
             const scitbx::af::const_ref<scitbx::vec3<double> >&,
             const scitbx::af::const_ref<double>&,
             const scitbx::af::const_ref<scitbx::vec3<double> >&,
             const cctbx::xray::scattering_type_registry&,
             const bool&);
    scitbx::af::shared<std::complex<double> > get_sum();

    // low level functions for reorganizing and transferring data
    void prepare_saxs(const scitbx::af::const_ref<std::string>&,
                      const scitbx::af::const_ref<scitbx::vec3<double> >&,
                      const scitbx::af::const_ref<double>&,
                      const scitbx::af::const_ref<double>&,
                      const scitbx::af::const_ref<double>&,
                      const scitbx::af::const_ref<double>&,
                      const cctbx::xray::scattering_type_registry&,
                      const bool&);
    void reorganize_xyz(const scitbx::af::const_ref<scitbx::vec3<double> >&);
    void transfer_xyz();
    void clear_xyz();
    void copy_solvent_weights(const scitbx::af::const_ref<double>&);
    void transfer_solvent_weights();
    void clear_solvent_weights();
    void reorganize_hkl(const scitbx::af::const_ref<scitbx::vec3<double> >&);
    void transfer_hkl();
    void clear_hkl();
    void copy_q(const scitbx::af::const_ref<double>&);
    void transfer_q();
    void clear_q();
    void copy_lattice(const scitbx::af::const_ref<double>&,
                      const scitbx::af::const_ref<double>&);
    void transfer_lattice();
    void clear_weights();
    void clear_lattice();
    void reorganize_rotations_translations
      (const scitbx::af::const_ref<double>&,
       const scitbx::af::const_ref<scitbx::vec3<double> >&);
    void transfer_rotations_translations();
    void clear_rotations_translations();
    void convert_scattering_types(const scitbx::af::const_ref<std::string>&,
                                  const cctbx::xray::scattering_type_registry&);
    void transfer_scattering_types();
    void clear_scattering_types();
    void convert_scattering_type_registry(const cctbx::xray::scattering_type_registry&);
    void transfer_scattering_type_registry(const bool&);
    void clear_scattering_type_registry();
    void allocate_amplitudes();
    void reset_amplitudes();
    void clear_amplitudes();
    void allocate_workspace(const int&);
    void clear_workspace();
    void run_kernel();
    void run_saxs_kernel();
    void run_solvent_saxs_kernel();
    void run_collect_solvent_saxs_kernel(const double&,const double&);

  private:
    // functions for reorganizing data
    void reorganize_coordinates
      (const scitbx::af::const_ref<scitbx::vec3<double> >&,
       const scitbx::af::const_ref<double>&);
    void reorganize_h(const scitbx::af::const_ref<scitbx::vec3<double> >&);
    void reorganize_q(const scitbx::af::const_ref<double>&,
                      const scitbx::af::const_ref<double>&,
                      const scitbx::af::const_ref<double>&);
    void reorganize_rt(const scitbx::af::const_ref<double>&,
                       const scitbx::af::const_ref<scitbx::vec3<double> >&);
    void convert_scatterers(const scitbx::af::const_ref<std::string>&,
                            const cctbx::xray::scattering_type_registry&,
                            const bool&);
    void clear_arrays();

    // xyz parameters
    int n_xyz;
    int padded_n_xyz;
    int size_xyz;
    fType * h_xyz, * d_xyz;

    // solvent weight parameters
    fType * h_solvent, * d_solvent;

    // reciprocal space parameters
    int n_h;
    int padded_n_h;
    int size_h;
    fType * h_h, * d_h;

    // rotation & translation parameters
    int n_rt;
    int size_rt;
    fType * h_rt, * d_rt;

    // scaling parameters
    fType * h_weights, * d_weights;

    // scatterer parameters
    int * h_scattering_type, * d_scattering_type;
    int n_types;
    int n_terms;
    int f_size;
    fType * h_a;
    fType * h_b;
    fType * h_c;

    // structure factor parameters
    bool amplitudes_allocated;
    fType * h_real, * h_imag, * d_real, * d_imag;

    // workspace
    bool workspace_allocated;
    int workspace_size;
    fType * d_workspace;
  };

  /* ==========================================================================
   */

}
}
#endif // DIRECT_SUMMATION_H
