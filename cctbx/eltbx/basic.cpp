/* Copyright (c) 2001-2002 The Regents of the University of California
   through E.O. Lawrence Berkeley National Laboratory, subject to
   approval by the U.S. Department of Energy.
   See files COPYRIGHT.txt and LICENSE.txt for further details.

   Revision history:
     2001 May: merged from CVS branch sgtbx_type (R.W. Grosse-Kunstleve)
     2001 Dec: SourceForge release (R.W. Grosse-Kunstleve)
 */

#include <cctbx/eltbx/basic.h>
#include <cstddef>
#include <ctype.h>

namespace cctbx { namespace eltbx { namespace basic {

  std::string strip_label(std::string const& label, bool exact)
  {
    std::string result;
    std::string::const_iterator l;
    for (l = label.begin(); l != label.end(); l++) {
      if (!isspace(*l)) break;
    }
    char digit = '\0';
    for (; l != label.end(); l++) {
      if (isspace(*l)) {
        break;
      }
      if (isdigit(*l)) {
        if (digit != '\0') break;
        digit = *l;
      }
      else if (*l == '+' || *l == '-') {
        if (digit == '\0')
            digit = '1';
        result += digit;
        result += *l++;
        break;
      }
      else if (digit != '\0') {
        break;
      }
      else {
        result += toupper(*l);
      }
    }
    if (exact && l != label.end() && !isspace(*l)) {
      return "";
    }
    return result;
  }

  int match_labels(std::string const& work_label, const char* tab_label)
  {
    int i;
    for (i = 0; i < work_label.size() && tab_label[i]; i++)
      if (work_label[i] != toupper(tab_label[i])) break;
    if (i == work_label.size() && tab_label[i] == '\0') return -i;
    if (i == 1 && isalpha(tab_label[1])) return 0;
    return i;
  }

}}} // namespace cctbx::eltbx::basic
