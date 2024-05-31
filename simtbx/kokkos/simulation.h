#ifndef SIMTBX_KOKKOS_SIMULATION_H
#define SIMTBX_KOKKOS_SIMULATION_H

#include "scitbx/array_family/shared.h"
#include "scitbx/array_family/flex_types.h"
#include "simtbx/nanoBragg/nanoBragg.h"
#include "simtbx/kokkos/structure_factors.h"
#include "simtbx/kokkos/detector.h"
#include "kokkostbx/kokkos_types.h"
#include "kokkostbx/kokkos_vector3.h"
#include "kokkostbx/kokkos_matrix3.h"

using vec3 = kokkostbx::vector3<CUDAREAL>;
using mat3 = kokkostbx::matrix3<CUDAREAL>;
using crystal_orientation_t = Kokkos::View<vec3***, Kokkos::LayoutRight, MemSpace>; // [phisteps, domains, 3]

namespace simtbx { namespace Kokkos {

namespace af = scitbx::af;

struct diffuse_api {
  inline diffuse_api() {};
  inline void show() {};
  bool enable = false;
  mat3 anisoG;
  mat3 anisoU;
  int stencil_size = 1;
  bool symmetrize_diffuse = true;
  int laue_group_num = 12;
  mat3 rotate_principal_axes;
  vec3 a0, b0, c0; // cell basis vectors
};

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

  void add_background(simtbx::Kokkos::kokkos_detector &, int const&);
  af::flex_double add_noise(simtbx::Kokkos::kokkos_detector &);
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
  crystal_orientation_t m_crystal_orientation = crystal_orientation_t("m_crystal_orientation", 0, 0, 3);
  CUDAREAL m_water_size = 0;
  CUDAREAL m_water_F = 0;
  CUDAREAL m_water_MW = 0;
  diffuse_api diffuse;
};
} // Kokkos
} // simtbx
#endif // SIMTBX_Kokkos_STRUCTURE_FACTORS_H
