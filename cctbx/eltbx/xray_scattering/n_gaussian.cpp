#include <cctbx/eltbx/xray_scattering/n_gaussian.h>
#include <cctbx/eltbx/xray_scattering/n_gaussian_raw.h>
#include <cctbx/error.h>

namespace cctbx { namespace eltbx { namespace xray_scattering {
namespace n_gaussian {

  std::size_t
  table_size() { return raw::get_table_size(); }

  std::size_t
  table_index(std::string const& label)
  {
    std::size_t i_entry = 0;
    const char** lbl=raw::get_labels();
    for(;*lbl;lbl++,i_entry++) {
      if (label == *lbl) break;
    }
    if (*lbl == 0) {
      throw error("Not in table of N-Gaussian approximations: " + label);
    }
    return i_entry;
  }

  table_entry::table_entry(
    std::string const& label,
    std::size_t n_terms)
  {
    init(table_index(label), n_terms);
  }

  table_entry::table_entry(
    std::size_t i_entry,
    double d_min,
    double max_relative_error)
  {
    init(i_entry, d_min, max_relative_error);
  }

  table_entry::table_entry(
    std::string const& label,
    double d_min,
    double max_relative_error)
  {
    init(table_index(label), d_min, max_relative_error);
  }

  void
  table_entry::init(
    std::size_t i_entry,
    std::size_t n_terms)
  {
    CCTBX_ASSERT(i_entry < table_size());
    CCTBX_ASSERT(n_terms >= 1);
    CCTBX_ASSERT(n_terms <= 6);
    init_core(i_entry, n_terms);
  }

  void
  table_entry::init(
    std::size_t i_entry,
    double d_min,
    double max_relative_error)
  {
    CCTBX_ASSERT(i_entry < table_size());
    if (d_min <= 0) {
      init_core(i_entry, 6);
    }
    else {
      double max_stol = 1/(2*d_min);
      raw::entry const& raw_entry = raw::get_table()[i_entry];
      int i = 5;
      for(;i>=0;i--) {
        if (raw_entry.max_stols[i] >= max_stol) {
          if (   max_relative_error <= 0
              || raw_entry.max_relative_errors[i] <= max_relative_error) {
            init_core(i_entry, 6-i);
            return;
          }
        }
      }
      throw error("No suitable N-Gaussian approximation.");
    }
  }

  void
  table_entry::init_core(
    std::size_t i_entry,
    std::size_t n_terms)
  {
    label_ = raw::get_labels()[i_entry];
    std::size_t i = 6 - n_terms;
    raw::entry const& raw_entry = raw::get_table()[i_entry];
    gaussian_ = xray_scattering::gaussian(af::const_ref<double>(
      raw_entry.coeff[i], 2*n_terms), 0., false),
    max_stol_ = raw_entry.max_stols[i];
    max_relative_error_ = raw_entry.max_relative_errors[i];
  }

}}}} // namespace cctbx::eltbx::xray_scattering::n_gaussian
