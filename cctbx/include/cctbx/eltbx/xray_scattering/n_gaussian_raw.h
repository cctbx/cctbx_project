#ifndef CCTBX_ELTBX_XRAY_SCATTERING_N_GAUSSIAN_RAW_H
#define CCTBX_ELTBX_XRAY_SCATTERING_N_GAUSSIAN_RAW_H

namespace cctbx { namespace eltbx { namespace xray_scattering {
namespace n_gaussian { namespace raw {

  struct entry
  {
    const double* max_stols;
    const double* max_relative_errors;
    const double** coeff;
  };

  const char*
  get_tag();

  const char**
  get_labels();

  unsigned int
  get_table_size();

  const entry*
  get_table();

}}}}} // namespace cctbx::eltbx::xray_scattering::n_gaussian::raw

#endif // CCTBX_ELTBX_XRAY_SCATTERING_N_GAUSSIAN_RAW_H
