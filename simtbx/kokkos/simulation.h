#ifndef SIMTBX_KOKKOS_SIMULATION_H
#define SIMTBX_KOKKOS_SIMULATION_H

#include "scitbx/array_family/shared.h"
#include "simtbx/nanoBragg/nanoBragg.h"
#include "simtbx/kokkos/structure_factors.h"
#include "simtbx/kokkos/detector.h"
#include "kokkos_types.h"

namespace simtbx { namespace Kokkos {

namespace af = scitbx::af;

struct exascale_api {
  inline
  exascale_api(const simtbx::nanoBragg::nanoBragg& nB) : SIM(nB) { }

  void show();
  void add_energy_channel_from_gpu_amplitudes(int const&,
    simtbx::Kokkos::kokkos_energy_channels &,
    simtbx::Kokkos::kokkos_detector &,
    double const&
  );
  void add_energy_channel_mask_allpanel(int const&,
    simtbx::Kokkos::kokkos_energy_channels &,
    simtbx::Kokkos::kokkos_detector &,
    af::shared<bool>
  );
  void add_energy_channel_mask_allpanel(int const&,
    simtbx::Kokkos::kokkos_energy_channels &,
    simtbx::Kokkos::kokkos_detector &,
    af::shared<std::size_t> const
  );
  void add_energy_multichannel_mask_allpanel(af::shared<int> const,
    simtbx::Kokkos::kokkos_energy_channels &,
    simtbx::Kokkos::kokkos_detector &,
    af::shared<std::size_t> const,
    af::shared<double> const
  );

  void add_background(simtbx::Kokkos::kokkos_detector &);
  void allocate();
  //~exascale_api();

  const simtbx::nanoBragg::nanoBragg& SIM;

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
