#ifndef CCTBX_ELTBX_XRAY_SCATTERING_N_GAUSSIAN_H
#define CCTBX_ELTBX_XRAY_SCATTERING_N_GAUSSIAN_H

#include <cctbx/eltbx/xray_scattering/gaussian.h>

namespace cctbx { namespace eltbx { namespace xray_scattering {
namespace n_gaussian {

  std::size_t
  table_size();

  std::size_t
  table_index(std::string label);

  //! Entry in table of N-Gaussian approximations.
  class table_entry
  {
    public:
      table_entry() {}

      table_entry(
        std::size_t i_entry,
        std::size_t n_terms) { init(i_entry, n_terms); }

      table_entry(
        std::string const&,
        std::size_t n_terms);

      table_entry(
        std::size_t i_entry,
        double d_min,
        double max_relative_error);

      table_entry(
        std::string const& label,
        double d_min,
        double max_relative_error);

      std::string const&
      label() const { return label_; }

      xray_scattering::gaussian const&
      gaussian() const { return gaussian_; }

      double
      max_stol() const { return max_stol_; }

      double
      d_min() const { return 1/(2*max_stol_); }

      double
      max_relative_error() const { return max_relative_error_; }

    protected:
      void
      init(
        std::size_t i_entry,
        std::size_t n_terms);

      void
      init(
        std::size_t i_entry,
        double d_min,
        double max_relative_error);

      void
      init_core(
        std::size_t i_entry,
        std::size_t n_terms);

      std::string label_;
      xray_scattering::gaussian gaussian_;
      double max_stol_;
      double max_relative_error_;
  };

}}}} // namespace cctbx::eltbx::xray_scattering::n_gaussian

#endif // CCTBX_ELTBX_XRAY_SCATTERING_N_GAUSSIAN_H
