// $Id$
/* Copyright (c) 2001 The Regents of the University of California through
   E.O. Lawrence Berkeley National Laboratory, subject to approval by the
   U.S. Department of Energy. See files COPYRIGHT.txt and
   cctbx/LICENSE.txt for further details.

   Revision history:
     2001 Jul 02: Merged from CVS branch sgtbx_special_pos (rwgk)
     2001 May 31: merged from CVS branch sgtbx_type (R.W. Grosse-Kunstleve)
     Apr 2001: SourceForge release (R.W. Grosse-Kunstleve)
 */

#ifndef CCTBX_ELTBX_BASIC_H
#define CCTBX_ELTBX_BASIC_H

#include <string>
#include <cctbx/error.h>

//! Element Toolbox namespace.
/*! The eltbx is a collection of tables of various x-ray and
    neutron scattering factors, element names, atomic numbers,
    atomic weights, ionic radii, and characteristic x-ray
    wavelenghts. Associated with each table are procedures for
    accessing the tabluated data, e.g. by using interpolation.
 */
namespace eltbx {

  using cctbx::error;
  using cctbx::error_index;

  //! Strip element label, scattering factor label, or ion label.
  /*! For internal use only.
      <p>
      Returns leading part of atom label. Ignores leading whitespace.
      Copies input label starting at first non-whitespace. Copying
      stops at whitespace, second digit, or the characters "+" or "-".
      Lower case letters are converted to upper case.<br>
      "na+" is converted to "NA1+", i.e. "1" is implicit.<br>
      "Ca+2" is converted to "CA2+", i.e. "+" or "-" are moved
        before the digit.<br>
      "Si4+A" is converted to "SI4+".
      <p>
      If Exact == true, copying must stop at a whitespace or
      end-of-string. Otherwise the empty string is returned.
      For example, with Exact == true "Si4+A" is converted to "".
   */
  std::string StripLabel(const std::string& Label, bool Exact = false);
  //! Compare an input label with a label from an internal table.
  /*! For internal use only.
      <p>
      Comparison is case-insensitive.<br>
      If the labels match exactly, the return value is the number
      of matching characters multiplied by -1.<br>
      If only the first characters match and the second character of
      TabLabel is a letter, the return value is 0.<br>
      Otherwise the return value is the number of matching characters.
   */
  int MatchLabels(const std::string& WorkLabel, const char* TabLabel);

} // namespace eltbx

#endif // CCTBX_ELTBX_BASIC_H
