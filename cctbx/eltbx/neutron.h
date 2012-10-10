#ifndef CCTBX_ELTBX_NEUTRON_H
#define CCTBX_ELTBX_NEUTRON_H

#include <string>
#include <complex>
#include <cctbx/eltbx/xray_scattering/gaussian.h>

namespace cctbx { namespace eltbx { namespace neutron {

  namespace detail
  {
    struct raw_record_neutron_news_1992
    {
      const char* label;
      float       bound_coh_scatt_length_real;
      float       bound_coh_scatt_length_imag;
      float       abs_cross_sect; // For 2200 m/s neutrons
    };
  }

  //! Access to neutron bound scattering lengths & cross-sections.
  /*! Reference:<br>
      Neutron News, Vol. 3, No. 3, 1992, pp. 29-37.
      <p>
      http://www.ncnr.nist.gov/resources/n-lengths/list.html
   */
  class neutron_news_1992_table
  {
    public:
      //! Default constructor. Calling certain methods may cause crashes!
      neutron_news_1992_table() : record_(0) {}

      //! Search internal table for the given element label.
      /*! If exact == true, the element label must exactly
          match the tabulated label. However, the lookup is not
          case-sensitive.
          <p>
          See also: eltbx::basic::strip_label()
       */
      explicit
      neutron_news_1992_table(std::string const& label, bool exact=false);

      //! Tests if the instance is constructed properly.
      /*! Shorthand for: label() != 0
          <p>
          Not available in Python.
       */
      bool
      is_valid() const { return record_->label != 0; }

      //! Element label from internal table.
      const char*
      label() const
      {
        return record_->label;
      }

      //! Bound coherent scattering length (fm) as a complex number.
      /*! 1 fm = 1e-15 m
       */
      std::complex<float>
      bound_coh_scatt_length() const
      {
        return std::complex<float>(record_->bound_coh_scatt_length_real,
                                   record_->bound_coh_scatt_length_imag);
      }

      //! Real part of bound coherent scattering length (fm).
      /*! 1 fm = 1e-15 m
          <p>
          Not available in Python.
       */
      float
      bound_coh_scatt_length_real() const
      {
        return record_->bound_coh_scatt_length_real;
      }

      //! Imaginary part of bound coherent scattering length (fm).
      /*! 1 fm = 1e-15 m
          <p>
          Not available in Python.
       */
      float
      bound_coh_scatt_length_imag() const
      {
        return record_->bound_coh_scatt_length_imag;
      }

      //! Absorption cross section (barn) for 2200 m/s neutrons.
      /*! 1 barn = 1e-24 cm^2
       */
      float
      abs_cross_sect() const
      {
        return record_->abs_cross_sect;
      }

      cctbx::eltbx::xray_scattering::gaussian
      fetch() const
      {
          return cctbx::eltbx::xray_scattering::gaussian(record_->bound_coh_scatt_length_real, true);
      }

    private:
      const detail::raw_record_neutron_news_1992* record_;
      friend class neutron_news_1992_table_iterator;
  };

  /*! \brief Iterator over neutron_news_1992_table entries.
   */
  class neutron_news_1992_table_iterator
  {
    public:
      //! Initialization of the iterator.
      neutron_news_1992_table_iterator();

      //! Retrieves the next entry from the internal table.
      /*! Use neutron_news_1992_table_iterator::is_valid() to detect
          end-of-iteration.
       */
      neutron_news_1992_table
      next();

    private:
      neutron_news_1992_table current_;
  };

}}} // cctbx::eltbx::neutron

#endif // CCTBX_ELTBX_NEUTRON_H
