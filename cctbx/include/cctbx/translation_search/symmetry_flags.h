#ifndef CCTBX_TRANSLATION_SEARCH_SYMMETRY_FLAGS_H
#define CCTBX_TRANSLATION_SEARCH_SYMMETRY_FLAGS_H

#include <cctbx/sgtbx/search_symmetry.h>

namespace cctbx { namespace translation_search {

  //! Grouping of flags determining the symmetry of the translation function.
  class symmetry_flags : public sgtbx::search_symmetry_flags
  {
    public:
      //! Default constructor. Some data members are not initialized!
      symmetry_flags() {}

      /*! \brief Custom sgtbx::search_symmetry_flags initializer
          for translation searches.
       */
      /*! See implementation.
       */
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

      //! Value as passed to the constructor.
      /*! Equivalent to use_space_group_symmetry().
       */
      bool
      is_isotropic_search_model() const { return use_space_group_symmetry_; }

      //! Value as passed to the constructor.
      /*! Equivalent to !use_seminvariants().
       */
      bool
      have_f_part() const { return !use_seminvariants_; }
  };

}} // namespace cctbx::translation_search

#endif // CCTBX_TRANSLATION_SEARCH_SYMMETRY_FLAGS_H
