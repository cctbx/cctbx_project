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
             const cctbx::xray::scattering_type_registry&);
    scitbx::af::shared<std::complex<double> > get_sum();

  private:
    // functions for reorganizing data
    void reorganize_coordinates
      (const scitbx::af::const_ref<scitbx::vec3<double> >&);
    void copy_solvent_weights(const scitbx::af::const_ref<double>&);
    void reorganize_h(const scitbx::af::const_ref<scitbx::vec3<double> >&);
    void reorganize_rt(const scitbx::af::const_ref<double>&,
                       const scitbx::af::const_ref<scitbx::vec3<double> >&);
    void convert_scatterers(const scitbx::af::const_ref<std::string>&,
                            const cctbx::xray::scattering_type_registry&);

    // xyz parameters
    int n_xyz;
    int padded_n_xyz;
    int size_xyz;
    fType * h_xyz, * d_xyz;

    // solvent weight parameters
    fType * d_solvent;

    // reciprocal space parameters
    int n_h;
    int padded_n_h;
    int size_h;
    fType * h_h, * d_h;

    // rotation & translation parameters
    int n_rt;
    int size_rt;
    fType * h_rt, * d_rt;

    // scatterer parameters
    int * h_scattering_type, * d_scattering_type;
    int n_types;
    int n_terms;
    int f_size;
    fType * h_a, * d_a;
    fType * h_b, * d_b;
    fType * h_c, * d_c;

    // structure factor parameters
    int sf_size;
    fType * sf_real, * sf_imag;
  };

  /* ==========================================================================
   */

}
}
#endif // DIRECT_SUMMATION_H
