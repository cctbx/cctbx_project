#ifndef CCTBX_ELTBX_BASIC_H
#define CCTBX_ELTBX_BASIC_H

#include <cctbx/error.h>

namespace cctbx {
  //! Element Toolbox namespace.
  /*! The eltbx is a collection of tables of various x-ray and
      neutron scattering factors, element names, atomic numbers,
      atomic weights, ionic radii, and characteristic x-ray
      wavelengths. Associated with each table are procedures for
      accessing the tabulated data, e.g. by using interpolation.
   */
  namespace eltbx {
  namespace basic {

  //! Strips element label, scattering factor label, or ion label.
  /*! For internal use only.
      <p>
      Returns leading part of atom label. Ignores leading whitespace.
      Copies input label starting at first non-whitespace. Copying
      stops at whitespace, second digit, or the characters "+" or "-".
      Lower case letters are converted to upper case.<br>
      "na+" is converted to "NA1+", i.e. "1" is implicit.<br>
      "Si4+A" is converted to "SI4+".
      <p>
      If exact == true, copying must stop at a whitespace or
      end-of-string. Otherwise the empty string is returned.
      For example, with exact == true "Si4+A" is converted to "".
   */
  std::string
  strip_label(std::string const& Label, bool exact = false);

  //! Compares an input label with a label from an internal table.
  /*! For internal use only.
      <p>
      Comparison is case-insensitive.<br>
      If the labels match exactly, the return value is the number
      of matching characters multiplied by -1.<br>
      If only the first characters match and the second character of
      tab_label is a letter, the return value is 0.<br>
      Otherwise the return value is the number of matching characters.
   */
  int
  match_labels(std::string const& work_label, const char* tab_label);

}}} // namespace cctbx::eltbx::basic

#endif // CCTBX_ELTBX_BASIC_H
