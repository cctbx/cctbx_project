#ifndef CCTBX_ELTBX_SASAKI_H
#define CCTBX_ELTBX_SASAKI_H

#include <cctbx/eltbx/fp_fdp.h>
#include <scitbx/constants.h>

namespace cctbx { namespace eltbx { namespace sasaki {

  namespace detail {

    const long n_raw = 280; // All tables have exactly 280 data points.
    const double first_wide = 0.1; // All tables in fpwide.tbl start at 0.1.
    const double wide_incr = 0.01;
    const double edge_incr = 0.0001;

    struct raw
    {
      float fp;
      float fdp;
    };

    struct info
    {
      char const* label;
      int z;
      raw* wide;
      double first_k;
      raw* k;
      double first_l1;
      raw* l1;
      double first_l2;
      raw* l2;
      double first_l3;
      raw* l3;
    };

  } // namespace detail

  //! Access to Sasaki tables.
  /*! Sasaki tables are available for elements with Z=4-83 and Z=92.
      They are valid in the energy range 4-124 keV and they have a fine
      step size close to the absorption edges (K,L1,L2,L3). The tables
      are therefore suitable for use in connection with anomalous
      diffraction experiments.
      <p>
      Reference: S.Sasaki (1989) Numerical Tables of Anomalous
      Scattering Factors Calculated by the Cromer and Liberman Method,
      KEK Report, 88-14, 1-136<br>
      ftp://pfweis.kek.jp/pub/Sasaki-table/
      <p>
      See also:
        http://www.esrf.fr/computing/scientific/dabax/
   */
  class table
  {
    public:
      //! Default constructor. Calling certain methods may cause crashes!
      table() : info_(0) {}

      //! Searches Sasaki tables for the given scattering factor label.
      /*! If exact == true, the scattering factor label must exactly
          match the tabulated label. However, the lookup is not
          case-sensitive.
          <p>
          If exception_if_no_match == true, an exception is thrown
          if the label cannot be found.<br>
          If exception_if_no_match == false, use is_valid() to test
          if the label was found.<br>
          exception_if_no_match is not available in Python. Use
          try-except instead.
          <p>
          See also: eltbx::basic::strip_label()
       */
      explicit
      table(
        std::string const& label,
        bool exact=false,
        bool exception_if_no_match=true);

      //! Tests if the instance is constructed properly.
      /*! Shorthand for: label() != 0
          <p>
          Not available in Python.
       */
      bool
      is_valid() const { return info_ != 0 && info_->label != 0; }

      //! Returns the scattering factor label.
      const char*
      label() const { return info_->label; }

      //! Returns the atomic number.
      int
      atomic_number() const { return info_->z; }

      //! Computes f-prime (f') and f-double-prime (f") for given energy [eV].
      /*! f-prime and f-double-prime are determined by linear
          interpolation.
          <p>
          See also:
            at_kev(),
            at_angstrom()
       */
      fp_fdp
      at_ev(double energy) const;

      //! Computes f-prime (f') and f-double-prime (f") for given energy [keV].
      /*! See also:
            at_ev()
            at_angstrom()
       */
      fp_fdp
      at_kev(double energy) const { return at_ev(energy * 1000); }

      /*! \brief Computes f-prime (f') and f-double-prime (f") for
          given wavelength [Angstrom].
       */
      /*! See also:
            at_kev(),
            at_ev(),
            scitbx::constants::factor_ev_angstrom
       */
      fp_fdp
      at_angstrom(double wavelength) const
      {
        return at_ev(scitbx::constants::factor_ev_angstrom / wavelength);
      }

    private:
      const detail::info* info_;
      friend class table_iterator;
  };

  /*! \brief Iterator over Sasaki tables.
   */
  class table_iterator
  {
    public:
      //! Initialization of the iterator.
      table_iterator();

      //! Retrieves the next entry from the internal table.
      /*! Use table::is_valid() to detect end-of-iteration.
       */
      table
      next();

    private:
      table current_;
  };

}}} // cctbx::eltbx::sasaki

#endif // CCTBX_ELTBX_SASAKI_H
