// $Id$

#ifndef CCTBX_ELTBX_WAVELENGTHS_H
#define CCTBX_ELTBX_WAVELENGTHS_H

#include <string>
#include <cctbx/constants.h>

namespace eltbx {

  namespace detail {
    struct RawWaveLength {
      const char* Label;
      float       lambda;
    };
  }

  //! Characteristic wavelengths of commonly used x-ray tube target materials.
  class WaveLength {
    public:
      //! Get i'th entry in the internal table.
      /*! Can be used to iterate over all entries.
       */
      WaveLength(int i);
      //! Lookup characteristic wavelength for given label.
      /*! The lookup is not case-sensitive.
       */
      WaveLength(const std::string& Label);
      //! Return label.
      inline const char* Label() const { return m_RawEntry->Label; }
      //! Return wavelength (Angstrom).
      inline float operator()() const { return m_RawEntry->lambda; }
      //! Return energy (keV).
      inline float Energy() const {
        if (m_RawEntry->lambda == 0.) return 0.;
        return cctbx::constants::factor_keV_Angstrom / m_RawEntry->lambda;
      }
    private:
      const detail::RawWaveLength* m_RawEntry;
  };

} // eltbx

#endif // CCTBX_ELTBX_WAVELENGTHS_H
