/* Copyright (c) 2001-2002 The Regents of the University of California
   through E.O. Lawrence Berkeley National Laboratory, subject to
   approval by the U.S. Department of Energy.
   See files COPYRIGHT.txt and LICENSE.txt for further details.

   Revision history:
     2002 Oct: Fragment from cctbx/maps (rwgk)
     2002 Jan: Created (R.W. Grosse-Kunstleve)
 */

#ifndef CCTBX_MAPTBX_SYMMETRY_FLAGS_H
#define CCTBX_MAPTBX_SYMMETRY_FLAGS_H

#include <cctbx/sgtbx/seminvariant.h>
#include <cctbx/sgtbx/space_group_type.h>

namespace cctbx { namespace maptbx {

  class symmetry_flags
  {
    public:
      symmetry_flags() {}

      explicit
      symmetry_flags(bool use_space_group_symmetry,
                     bool use_normalizer_k2l=false,
                     bool use_structure_seminvariants=false)
      :
        use_space_group_symmetry_(use_space_group_symmetry),
        use_normalizer_k2l_(use_normalizer_k2l),
        use_structure_seminvariants_(use_structure_seminvariants)
      {}

      bool
      use_space_group_symmetry() const { return use_space_group_symmetry_; }

      bool
      use_normalizer_k2l() const { return use_normalizer_k2l_; }

      bool
      use_structure_seminvariants() const
      {
        return use_structure_seminvariants_;
      }

      sgtbx::space_group
      select_sub_space_group(sgtbx::space_group_type const& sg_type) const
      {
        sgtbx::space_group result;
        if (use_space_group_symmetry()) {
          result = sg_type.group();
        }
        else if (use_structure_seminvariants()) {
          for(std::size_t i=1;i<sg_type.group().n_ltr();i++) {
            result.expand_ltr(sg_type.group().ltr(i));
          }
        }
        if (use_normalizer_k2l()) {
          result.expand_smx(
            sg_type.addl_generators_of_euclidean_normalizer(true, false)
              .const_ref());
        }
        return result;
      }

      af::int3
      grid_factors(sgtbx::space_group_type const& sg_type) const
      {
        af::int3 grid_ss(1,1,1);
        if (use_structure_seminvariants()) {
          grid_ss = sgtbx::structure_seminvariant(sg_type.group()).gridding();
        }
        return select_sub_space_group(sg_type).refine_gridding(grid_ss);
      }

      bool
      operator==(symmetry_flags const& rhs) const
      {
        return
             use_space_group_symmetry_ == rhs.use_space_group_symmetry_
          && use_normalizer_k2l_ == rhs.use_normalizer_k2l_
          && use_structure_seminvariants_ == rhs.use_structure_seminvariants_;
      }

      bool
      operator!=(symmetry_flags const& rhs) const { return !((*this) == rhs); }

    protected:
      bool use_space_group_symmetry_;
      bool use_normalizer_k2l_;
      bool use_structure_seminvariants_;
  };

}} // namespace cctbx::maptbx

#endif // CCTBX_MAPTBX_SYMMETRY_FLAGS_H
