/* Copyright (c) 2001-2002 The Regents of the University of California
   through E.O. Lawrence Berkeley National Laboratory, subject to
   approval by the U.S. Department of Energy.
   See files COPYRIGHT.txt and LICENSE.txt for further details.

   Revision history:
     2001 May: merged from CVS branch sgtbx_type (R.W. Grosse-Kunstleve)
     2001 Apr: SourceForge release (R.W. Grosse-Kunstleve)
 */

#include <cctbx/eltbx/wavelengths.h>
#include <cctbx/eltbx/basic.h>
#include <ctype.h> // cannot use cctype b/o non-conforming compilers

namespace cctbx { namespace eltbx { namespace wavelengths {

namespace detail {
namespace {

    const raw table[] = {
// BEGIN_COMPILED_IN_REFERENCE_DATA
      {"CrA1", 2.28970}, {"CrA2", 2.29361}, {"Cr", 2.2909},
      {"FeA1", 1.93604}, {"FeA2", 1.93998}, {"Fe", 1.9373},
      {"CuA1", 1.54056}, {"CuA2", 1.54439}, {"Cu", 1.5418},
      {"MoA1", 0.70930}, {"MoA2", 0.71359}, {"Mo", 0.7107},
      {"AgA1", 0.55941}, {"AgA2", 0.56380}, {"Ag", 0.5608},
      {"", 0}
// END_COMPILED_IN_REFERENCE_DATA
    };

} // namespace <anonymous>
} // namespace detail

  characteristic::characteristic(std::string const& label)
  {
    std::string lbl = label;
    if (lbl.size() > 0) lbl[0] = toupper(lbl[0]);
    if (lbl.size() > 1) lbl[1] = tolower(lbl[1]);
    if (lbl.size() > 2) lbl[2] = toupper(lbl[2]);
    for (raw_ = detail::table; raw_->value; raw_++) {
     if (lbl == std::string(raw_->label)) return;
    }
    throw error("Unknown label for characteristic wavelength.");
  }

  characteristic_iterator::characteristic_iterator()
  :
    current_("CrA1")
  {}

  characteristic
  characteristic_iterator::next()
  {
    characteristic result = current_;
    if (current_.is_valid()) current_.raw_++;
    return result;
  }

}}} // namespace cctbx::eltbx::wavelengths
