#ifndef CCTBX_ELTBX_XRAY_SCATTERING_N_GAUSSIAN_RAW_H
#define CCTBX_ELTBX_XRAY_SCATTERING_N_GAUSSIAN_RAW_H

namespace cctbx { namespace eltbx { namespace xray_scattering {
namespace n_gaussian { namespace raw {

  struct entry
  {
    const char* label;
    const double* stols;
    const double** coeff;
  };

}}}}} // namespace cctbx::eltbx::xray_scattering::n_gaussian::raw

#endif // CCTBX_ELTBX_XRAY_SCATTERING_N_GAUSSIAN_RAW_H
