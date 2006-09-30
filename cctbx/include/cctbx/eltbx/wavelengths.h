#ifndef CCTBX_ELTBX_WAVELENGTHS_H
#define CCTBX_ELTBX_WAVELENGTHS_H

#include <scitbx/constants.h>
#include <string>

namespace cctbx { namespace eltbx { namespace wavelengths {

  namespace detail
  {
    struct raw
    {
      const char* label;
      float       value;
    };
  }

  //! Characteristic wavelengths of commonly used x-ray tube target materials.
  class characteristic
  {
    public:
      //! Default constructor. Calling certain methods may cause crashes!
      characteristic() : raw_(0) {}

      //! Retrieves characteristic wavelength for given label.
      /*! The lookup is not case-sensitive.
          <p>
          See also: is_valid()
       */
      characteristic(std::string const& label);

      //! Tests if the label passed to the contructor could be matched.
      /*! Shorthand for: as_angstrom() != 0
       */
      bool
      is_valid() const { return raw_->value != 0; }

      //! Returns label (e.g. Cu or CuA1).
      /*! See also: is_valid()
       */
      const char*
      label() const
      {
        return raw_->label;
      }

      //! Wavelength (Angstrom).
      /*! See also: is_valid()
       */
      float
      as_angstrom() const
      {
        return raw_->value;
      }

      //! Energy [keV].
      /*! See also: is_valid()
       */
      float
      as_kev() const
      {
        if (raw_->value == 0.) return 0.;
        return scitbx::constants::factor_kev_angstrom / raw_->value;
      }

      //! Energy [eV].
      /*! See also: is_valid()
       */
      float
      as_ev() const { return as_kev() * 1000; }

    private:
      const detail::raw* raw_;
      friend class characteristic_iterator;
  };

  //! Iterator over characteristic wavelengths.
  class characteristic_iterator
  {
    public:
      //! Initialization of the iterator.
      characteristic_iterator();

      //! Retrieves the next entry from the internal table.
      /*! Use characteristic::is_valid() to detect end-of-iteration.
       */
      characteristic
      next();

    private:
      characteristic current_;
  };

}}} // cctbx::eltbx::wavelengths

#endif // CCTBX_ELTBX_WAVELENGTHS_H
