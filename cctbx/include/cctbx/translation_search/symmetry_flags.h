#ifndef CCTBX_TRANSLATION_SEARCH_SYMMETRY_FLAGS_H
#define CCTBX_TRANSLATION_SEARCH_SYMMETRY_FLAGS_H

#include <cctbx/sgtbx/search_symmetry.h>

namespace cctbx { namespace translation_search {

  class symmetry_flags : public sgtbx::search_symmetry_flags
  {
    public:
      symmetry_flags() {}

      symmetry_flags(bool is_isotropic_search_model,
                     bool have_f_part)
      :
        sgtbx::search_symmetry_flags(
          is_isotropic_search_model,
          0,
          !have_f_part,
          is_isotropic_search_model && !have_f_part,
          false)
      {}

      bool
      is_isotropic_search_model() const { return use_space_group_symmetry_; }

      bool
      have_f_part() const { return !use_seminvariant_; }
  };

}} // namespace cctbx::translation_search

#endif // CCTBX_TRANSLATION_SEARCH_SYMMETRY_FLAGS_H
