#ifndef CCTBX_SGTBX_SEARCH_SYMMETRY_H
#define CCTBX_SGTBX_SEARCH_SYMMETRY_H

#include <cctbx/sgtbx/seminvariant.h>
#include <cctbx/sgtbx/space_group_type.h>

namespace cctbx { namespace sgtbx {

  class search_symmetry_flags
  {
    public:
      search_symmetry_flags() {}

      explicit
      search_symmetry_flags(bool use_space_group_symmetry,
                            bool use_normalizer_k2l=false,
                            bool use_normalizer_l2n=false,
                            bool use_structure_seminvariants=false)
      :
        use_space_group_symmetry_(use_space_group_symmetry),
        use_normalizer_k2l_(use_normalizer_k2l),
        use_normalizer_l2n_(use_normalizer_l2n),
        use_structure_seminvariants_(use_structure_seminvariants)
      {}

      bool
      use_space_group_symmetry() const { return use_space_group_symmetry_; }

      bool
      use_normalizer_k2l() const { return use_normalizer_k2l_; }

      bool
      use_normalizer_l2n() const { return use_normalizer_l2n_; }

      bool
      use_structure_seminvariants() const
      {
        return use_structure_seminvariants_;
      }

      bool
      operator==(search_symmetry_flags const& rhs) const
      {
        return
             use_space_group_symmetry_ == rhs.use_space_group_symmetry_
          && use_normalizer_k2l_ == rhs.use_normalizer_k2l_
          && use_normalizer_l2n_ == rhs.use_normalizer_l2n_
          && use_structure_seminvariants_ == rhs.use_structure_seminvariants_;
      }

      bool
      operator!=(search_symmetry_flags const& rhs) const
      {
        return !((*this) == rhs);
      }

    protected:
      bool use_space_group_symmetry_;
      bool use_normalizer_k2l_;
      bool use_normalizer_l2n_;
      bool use_structure_seminvariants_;
  };

  class search_symmetry
  {
    public:
      search_symmetry() {}

      search_symmetry(
        search_symmetry_flags const& flags,
        space_group_type const& group_type)
      :
        flags_(flags)
      {
        init(group_type);
      }

      search_symmetry(
        search_symmetry_flags const& flags,
        space_group_type const& group_type,
        structure_seminvariant const& seminvariant)
      :
        flags_(flags)
      {
        init(group_type, &seminvariant);
      }

      search_symmetry_flags const&
      flags() const { return flags_; }

      space_group const&
      group() const { return group_; }

      scitbx::af::small<sg_vec3, 3> const&
      continuous_shifts() const { return continuous_shifts_; }

    protected:
      search_symmetry_flags flags_;
      space_group group_;
      scitbx::af::small<sg_vec3, 3> continuous_shifts_;

      void
      init(
        space_group_type const& group_type,
        const structure_seminvariant* seminvariant=0)
      {
        if (flags_.use_space_group_symmetry()) {
          group_ = group_type.group();
        }
        else if (flags_.use_structure_seminvariants()) {
          for(std::size_t i=1;i<group_type.group().n_ltr();i++) {
            group_.expand_ltr(group_type.group().ltr(i));
          }
        }
        if (flags_.use_structure_seminvariants() && seminvariant) {
          af::small<ss_vec_mod, 3> const&
            ss = seminvariant->vectors_and_moduli();
          for(std::size_t i_ss=0;i_ss<ss.size();i_ss++) {
            if (ss[i_ss].m == 0) {
              continuous_shifts_.push_back(ss[i_ss].v);
            }
            else {
              group_.expand_ltr(tr_vec(ss[i_ss].v, ss[i_ss].m)
                .new_denominator(group_.t_den()));
            }
          }
        }
        if (flags_.use_normalizer_k2l() || flags_.use_normalizer_l2n()) {
          group_.expand_smx(
            group_type.addl_generators_of_euclidean_normalizer(
              flags_.use_normalizer_k2l(),
              flags_.use_normalizer_l2n()).const_ref());
        }
      }
  };

}} // namespace cctbx::sgtbx

#endif // CCTBX_SGTBX_SEARCH_SYMMETRY_H
