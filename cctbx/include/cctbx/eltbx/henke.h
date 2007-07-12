#ifndef CCTBX_ELTBX_HENKE_H
#define CCTBX_ELTBX_HENKE_H

#include <cctbx/eltbx/fp_fdp.h>
#include <scitbx/constants.h>

namespace cctbx { namespace eltbx { namespace henke {

  //! Access to Henke tables.
  /*! Henke tables are available for elements with Z=1-92.
      Each table contains 500+ points on a uniform logarithmic mesh
      from 10 to 30,000 eV with points added 0.1 eV above and below
      "sharp" absorption edges. The atomic scattering factors are
      based upon experimental measurements of the atomic
      photoabsorption cross section. The absorption measurements
      provide values for the imaginary part of the atomic scattering
      factor. The real part is calculated from the absorption
      measurements using the Kramers-Kronig integral relations.
      <p>
      Reference: B. L. Henke, E. M. Gullikson, and J. C. Davis,
      Atomic Data and Nuclear Data Tables Vol. 54 No. 2 (July 1993).<br>
      ftp://grace.lbl.gov/pub/sf/
      <p>
      See also:
        http://www-cxro.lbl.gov/optical_constants/asf.html <br>
        http://www.esrf.fr/computing/scientific/dabax/
   */
  class table
  {
    public:
      //! Default constructor. Calling certain methods may cause crashes!
      table() : label_z_e_fp_fdp_(0) {}

      //! Searches Henke tables for the given scattering factor label.
      /*! If exact == true, the scattering factor label must exactly
          match the tabulated label. However, the search is not
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
      is_valid() const
      {
        return label_z_e_fp_fdp_ != 0 && label_z_e_fp_fdp_->label != 0;
      }

      //! Returns the scattering factor label.
      const char*
      label() const { return label_z_e_fp_fdp_->label; }

      //! Returns the atomic number.
      int
      atomic_number() const { return label_z_e_fp_fdp_->z; }

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
      const anomalous::label_z_e_fp_fdp* label_z_e_fp_fdp_;
      friend class table_iterator;
  };

  /*! \brief Iterator over Henke tables.
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

}}} // cctbx::eltbx::henke

#endif // CCTBX_ELTBX_HENKE_H
