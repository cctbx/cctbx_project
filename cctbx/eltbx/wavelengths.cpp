// $Id$
/* Copyright (c) 2001 The Regents of the University of California through
   E.O. Lawrence Berkeley National Laboratory, subject to approval by the
   U.S. Department of Energy. See files COPYRIGHT.txt and
   cctbx/LICENSE.txt for further details.

   Revision history:
     2001 May 31: merged from CVS branch sgtbx_type (R.W. Grosse-Kunstleve)
     Apr 2001: SourceForge release (R.W. Grosse-Kunstleve)
 */

#include <ctype.h> // cannot use cctype b/o non-conforming compilers
#include <cctbx/eltbx/basic.h>
#include <cctbx/eltbx/wavelengths.h>

namespace eltbx {
  namespace tables {

    const detail::RawWaveLength RawWaveLengths[] = {
// BEGIN_COMPILED_IN_REFERENCE_DATA
      {"CrA1", 2.28970}, {"CrA2", 2.29361}, {"Cr", 2.2909},
      {"FeA1", 1.93604}, {"FeA2", 1.93998}, {"Fe", 1.9373},
      {"CuA1", 1.54056}, {"CuA2", 1.54439}, {"Cu", 1.5418},
      {"MoA1", 0.70930}, {"MoA2", 0.71359}, {"Mo", 0.7107},
      {"AgA1", 0.55941}, {"AgA2", 0.56380}, {"Ag", 0.5608},
      {0, 0}
// END_COMPILED_IN_REFERENCE_DATA
    };

  } // namespace tables
} // namespace eltbx

namespace eltbx {

  WaveLength::WaveLength(int i)
  {
    static const detail::RawWaveLength null = {"", 0.};
    if (   i < 0
        || i >= (  sizeof   tables::RawWaveLengths
                 / sizeof (*tables::RawWaveLengths)) - 1) {
      m_RawEntry = &null;
    }
    else {
      m_RawEntry = &tables::RawWaveLengths[i];
    }
  }

  WaveLength::WaveLength(const std::string& Label)
  {
    std::string lbl = Label;
    if (lbl.size() > 0) lbl[0] = toupper(lbl[0]);
    if (lbl.size() > 1) lbl[1] = tolower(lbl[1]);
    if (lbl.size() > 2) lbl[2] = toupper(lbl[2]);
    for (m_RawEntry = tables::RawWaveLengths;
         m_RawEntry->Label;
         m_RawEntry++) {
     if (lbl == std::string(m_RawEntry->Label)) return;
    }
    throw error("Unknown label for characteristic wavelength.");
  }

} // namespace eltbx
