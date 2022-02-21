#ifndef SIMTBX_KOKKOS_SIMULATION_H
#define SIMTBX_KOKKOS_SIMULATION_H

#include "scitbx/array_family/shared.h"
#include "simtbx/nanoBragg/nanoBragg.h"
#include "simtbx/nanoBragg/nanotypes.h"
#include "simtbx/kokkos/structure_factors.h"
#include "simtbx/kokkos/detector.h"
#include "kokkos_types.h"

namespace simtbx { namespace Kokkos {

namespace af = scitbx::af;

struct data_holder {
  inline
  data_holder(const simtbx::nanoBragg::nanoBragg& nB) :
    roi_xmin(nB.roi_xmin),
    roi_xmax(nB.roi_xmax),
    roi_ymin(nB.roi_ymin),
    roi_ymax(nB.roi_ymax),
    oversample(nB.oversample),
    point_pixel(nB.point_pixel),
    pixel_size(nB.pixel_size),
    detector_thickstep(nB.detector_thickstep),
    detector_thicksteps(nB.detector_thicksteps),
    detector_thick(nB.detector_thick),
    detector_attnlen(nB.detector_attnlen),
    curved_detector(nB.curved_detector),
    distance(nB.distance),
    close_distance(nB.close_distance),
    dmin(nB.dmin),
    phi0(nB.phi0),
    phistep(nB.phistep),
    phisteps(nB.phisteps),
    sources(nB.sources),
    mosaic_spread(nB.mosaic_spread),
    mosaic_domains(nB.mosaic_domains),
    Na(nB.Na),Nb(nB.Nb),Nc(nB.Nc),
    fluence(nB.fluence),
    spot_scale(nB.spot_scale),
    integral_form(nB.integral_form),
    default_F(nB.default_F),
    interpolate(nB.interpolate),
    nopolar(nB.nopolar),
    polarization(nB.polarization),
    fudge(nB.fudge),
    xtal_shape(nB.xtal_shape),
    V_cell(nB.V_cell),
    stols(nB.stols),
    amorphous_molecules(nB.amorphous_molecules)
  {
    for (int idx = 0 ; idx < sources; ++idx) {
      source_X.push_back(nB.source_X[idx]); source_Y.push_back(nB.source_Y[idx]);
      source_Z.push_back(nB.source_Z[idx]); source_I.push_back(nB.source_I[idx]);
      source_lambda.push_back(nB.source_lambda[idx]);
    }
    for (int idx = 0 ; idx < stols; ++idx) {
      stol_of.push_back(nB.stol_of[idx]); Fbg_of.push_back(nB.Fbg_of[idx]);
    }
    for (int idx = 0 ; idx < 4; ++idx) {
      beam_vector[idx] = nB.beam_vector[idx]; spindle_vector[idx] = nB.spindle_vector[idx];
      a0[idx] = nB.a0[idx]; b0[idx] = nB.b0[idx]; c0[idx] = nB.c0[idx]; polar_vector[idx] = nB.polar_vector[idx];
    }
    for (int idx = 0 ; idx < 9 * mosaic_domains; ++idx) {
      mosaic_umats.push_back ( nB.mosaic_umats[idx] );
    }
  }
    int roi_xmin,roi_xmax,roi_ymin,roi_ymax;
    int oversample;
    bool point_pixel;
    double pixel_size;
    double detector_thickstep;
    int detector_thicksteps;
    double detector_thick;
    double detector_attnlen;
    bool curved_detector;
    double distance;
    double close_distance;
    double dmin;
    double phi0;
    double phistep;
    int phisteps;
    int sources;
    double mosaic_spread;
    int mosaic_domains;
    double Na,Nb,Nc;
    double fluence;
    double spot_scale;
    bool integral_form;
    double default_F;
    int interpolate;
    bool nopolar;
    double polarization; //=0.0;  Kahn "polarization" parameter [-1:1]
    double fudge;
    af::shared<double> source_X,source_Y,source_Z,source_I,source_lambda;
    simtbx::nanoBragg::shapetype xtal_shape;
    double V_cell;
    int stols;
    af::shared<double> stol_of,Fbg_of;
    double amorphous_molecules;
    double beam_vector[4];
    double spindle_vector[4];
    double a0[4],b0[4],c0[4];
    double polar_vector[4];
    af::shared<double> mosaic_umats;
};

struct exascale_api {
  inline
  exascale_api(const simtbx::nanoBragg::nanoBragg& nB) : SIM(nB) { }

  void show();
  void add_energy_channel_from_gpu_amplitudes(int const&,
    simtbx::Kokkos::kokkos_energy_channels &,
    simtbx::Kokkos::kokkos_detector &
  );
  void add_energy_channel_mask_allpanel(int const&,
    simtbx::Kokkos::kokkos_energy_channels &,
    simtbx::Kokkos::kokkos_detector &,
    af::shared<bool>
  );
  void add_energy_channel_mask_allpanel(int const&,
    simtbx::Kokkos::kokkos_energy_channels &,
    simtbx::Kokkos::kokkos_detector &,
    af::shared<int> const
  );
  void add_energy_multichannel_mask_allpanel(af::shared<int> const,
    simtbx::Kokkos::kokkos_energy_channels &,
    simtbx::Kokkos::kokkos_detector &,
    af::shared<int> const
  );

  void add_background(simtbx::Kokkos::kokkos_detector &);
  void allocate();
  //~exascale_api();

  data_holder SIM;

  CUDAREAL m_subpixel_size;
  int m_steps;

  const int m_vector_length = 4;
  vector_cudareal_t m_beam_vector = vector_cudareal_t("m_beam_vector", m_vector_length);
  vector_cudareal_t m_spindle_vector = vector_cudareal_t("m_spindle_vector", m_vector_length);
  vector_cudareal_t m_a0 = vector_cudareal_t("m_a0", m_vector_length);
  vector_cudareal_t m_b0 = vector_cudareal_t("m_b0", m_vector_length);
  vector_cudareal_t m_c0 = vector_cudareal_t("m_c0", m_vector_length);
  vector_cudareal_t m_polar_vector = vector_cudareal_t("m_polar_vector", m_vector_length);
  vector_cudareal_t m_source_X = vector_cudareal_t("m_source_X", 0);
  vector_cudareal_t m_source_Y = vector_cudareal_t("m_source_Y", 0);
  vector_cudareal_t m_source_Z = vector_cudareal_t("m_source_Z", 0);
  vector_cudareal_t m_source_I = vector_cudareal_t("m_source_I", 0);
  vector_cudareal_t m_source_lambda = vector_cudareal_t("m_source_lambda", 0);
  vector_cudareal_t m_mosaic_umats = vector_cudareal_t("m_mosaic_umats", 0);
  CUDAREAL m_water_size = 0;
  CUDAREAL m_water_F = 0;
  CUDAREAL m_water_MW = 0;
};
} // Kokkos
} // simtbx
#endif // SIMTBX_Kokkos_STRUCTURE_FACTORS_H
