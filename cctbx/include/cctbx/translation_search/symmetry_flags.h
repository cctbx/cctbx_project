/* Copyright (c) 2001-2002 The Regents of the University of California
   through E.O. Lawrence Berkeley National Laboratory, subject to
   approval by the U.S. Department of Energy.
   See files COPYRIGHT.txt and LICENSE.txt for further details.

   Revision history:
     Oct 2002: Modified fragment from phenix/translation_search.h (rwgk)
     Jan 2002: Created (R.W. Grosse-Kunstleve)
 */

#ifndef CCTBX_TRANSLATION_SEARCH_SYMMETRY_FLAGS_H
#define CCTBX_TRANSLATION_SEARCH_SYMMETRY_FLAGS_H

#include <cctbx/maptbx/symmetry_flags.h>

namespace cctbx { namespace translation_search {

  class symmetry_flags : public maptbx::symmetry_flags
  {
    public:
      symmetry_flags() {}

      symmetry_flags(bool is_isotropic_search_model,
                     bool have_f_part)
      :
        maptbx::symmetry_flags(
          is_isotropic_search_model,
          is_isotropic_search_model && !have_f_part,
          !have_f_part)
      {}

      bool
      is_isotropic_search_model() const { return use_space_group_symmetry_; }

      bool
      have_f_part() const { return !use_structure_seminvariants_; }
  };

}} // namespace cctbx::translation_search

#endif // CCTBX_TRANSLATION_SEARCH_SYMMETRY_FLAGS_H
