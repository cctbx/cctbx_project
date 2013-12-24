#ifndef CCTBX_ELTBX_FP_FDP_H
#define CCTBX_ELTBX_FP_FDP_H

#include <cctbx/eltbx/basic.h>
#include <complex>
#include <ctype.h>

namespace cctbx { namespace eltbx {

  namespace anomalous {

    struct e_fp_fdp
    {
      float e;
      float fp;
      float fdp;
    };

    struct label_z_e_fp_fdp
    {
      const char* label;
      int z;
      const e_fp_fdp* data;
    };

    struct iterator_init_tag {};

  } // namespace anomalous

  //! Indicator for undefined values.
  static const float fp_fdp_undefined = -9999.00;

  //! Helper class for passing f' (f-prime) and f" (f-double-prime).
  class fp_fdp
  {
    public:
      //! Constructor.
      explicit
      fp_fdp(float fp  = fp_fdp_undefined,
             float fdp = fp_fdp_undefined)
      : fp_(fp), fdp_(fdp)
      {}

      //! Tests if f-prime is defined.
      bool is_valid_fp() const { return fp_ != fp_fdp_undefined; }

      //! Tests if f-double-prime is defined.
      bool is_valid_fdp() const { return fdp_ != fp_fdp_undefined; }

      //! Tests if f-prime and f-double-prime are defined.
      bool is_valid() const { return is_valid_fp() && is_valid_fdp(); }

      //! Returns f-prime (f').
      float fp() const { return fp_; }

      //! Returns f-double-prime (f").
      float fdp() const { return fdp_; }

      //! Returns complex(f-prime, f-double-prime).
      /*! Not available in Python.
       */
      std::complex<float> as_complex_float() const
      {
        return std::complex<float>(fp_, fdp_);
      }

      //! Returns complex(f-prime, f-double-prime).
      /*! Python: as_complex()
       */
      std::complex<double> as_complex_double() const
      {
        return std::complex<double>(fp_, fdp_);
      }

    private:
      float fp_;
      float fdp_;
  };

  namespace anomalous {

    template <typename TableType>
    const TableType*
    find_entry(const TableType* tables,
               std::string const& work_label,
               bool exact,
               bool exception_if_no_match)
    {
      // map D to H
      std::string wl = work_label;
      if (wl == "D") wl = "H";
      int m = 0;
      const TableType* m_entry = 0;
      for (const TableType* entry = tables; entry->label; entry++)
      {
        int i = basic::match_labels(wl, entry->label);
        if (i < 0) return entry;
        if (i > m && !isdigit(entry->label[i - 1])) {
          m = i;
          m_entry = entry;
        }
      }
      if (exception_if_no_match && (exact || !m_entry)) {
        throw std::invalid_argument(
          "Unknown scattering type label: " + std::string(work_label));
      }
      return m_entry;
    }

    fp_fdp interpolate(const label_z_e_fp_fdp* m_label_z_e_fp_fdp,
                       double energy);

  } // namespace anomalous

}} // cctbx::eltbx

#endif // CCTBX_ELTBX_FP_FDP_H
