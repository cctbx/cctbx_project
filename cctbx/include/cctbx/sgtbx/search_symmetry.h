#ifndef CCTBX_SGTBX_SEARCH_SYMMETRY_H
#define CCTBX_SGTBX_SEARCH_SYMMETRY_H

#include <cctbx/sgtbx/seminvariant.h>
#include <cctbx/sgtbx/space_group_type.h>

namespace cctbx { namespace sgtbx {

  //! Grouping of flags determining the symmetry of search spaces.
  /*! See also:
        cctbx::sgtbx::search_symmetry,
        cctbx::translation_search::symmetry_flags,
        cctbx::crystal::close_packing::hexagonal_sampling_generator,
        cctbx::maptbx::grid_tags,
        cctbx.euclidean_model_matching (Python)

      Remark: A "normalizer" (see below) is the term used in
      mathematics and in the International Tables for Crystallography
      to describe what is known colloquially as "Cheshire symmetry."

      Related reference:
        R.W. Grosse-Kunstleve and P.D. Adams (2003),
        On symmetries of substructures,
        Acta Cryst D59, 1974-1977,
        http://journals.iucr.org/d/issues/2003/11/00/ba5051/
   */
  class search_symmetry_flags
  {
    public:
      //! Default constructor. Some data members are not initialized!
      search_symmetry_flags() {}

      //! Definition of the symmetry flags.
      /*! use_space_group_ltr > 1 indicates that the lattice
          translations are to be used even if
          use_space_group_symmetry is false.

          use_space_group_ltr has no effect if use_space_group_symmetry
          is true.

          use_seminvariants implies use_space_group_ltr unless
          use_space_group_ltr < 0.

          use_seminvariants indicates that the "Additional generators"
          in the "Translations" column of Table 15.3.2 in the International
          Tables for Crystallography, Volume A are to be used.

          use_normalizer_k2l indicates that the Additional generator
          in the "Inversion through a centre at" column of
          Table 15.3.2 is to be used.

          use_normalizer_l2n indicates that the Additional generators
          in the "Further generators" column of Table 15.3.2 are to be
          used.
       */
      explicit
      search_symmetry_flags(
        bool use_space_group_symmetry,
        int use_space_group_ltr=0,
        bool use_seminvariants=false,
        bool use_normalizer_k2l=false,
        bool use_normalizer_l2n=false)
      :
        use_space_group_symmetry_(use_space_group_symmetry),
        use_space_group_ltr_(use_space_group_ltr),
        use_seminvariants_(use_seminvariants),
        use_normalizer_k2l_(use_normalizer_k2l),
        use_normalizer_l2n_(use_normalizer_l2n)
      {}

      //! Value as passed to the constructor.
      bool
      use_space_group_symmetry() const { return use_space_group_symmetry_; }

      //! Value as passed to the constructor.
      int
      use_space_group_ltr() const { return use_space_group_ltr_; }

      //! Value as passed to the constructor.
      bool
      use_seminvariants() const { return use_seminvariants_; }

      //! Value as passed to the constructor.
      bool
      use_normalizer_k2l() const { return use_normalizer_k2l_; }

      //! Value as passed to the constructor.
      bool
      use_normalizer_l2n() const { return use_normalizer_l2n_; }

      //! True if all values are exactly identical.
      bool
      operator==(search_symmetry_flags const& rhs) const
      {
        return
             use_space_group_symmetry_ == rhs.use_space_group_symmetry_
          && use_space_group_ltr_ == rhs.use_space_group_ltr_
          && use_seminvariants_ == rhs.use_seminvariants_
          && use_normalizer_k2l_ == rhs.use_normalizer_k2l_
          && use_normalizer_l2n_ == rhs.use_normalizer_l2n_;
      }

      //! Negation of test for equality.
      bool
      operator!=(search_symmetry_flags const& rhs) const
      {
        return !((*this) == rhs);
      }

    protected:
      bool use_space_group_symmetry_;
      int use_space_group_ltr_;
      bool use_seminvariants_;
      bool use_normalizer_k2l_;
      bool use_normalizer_l2n_;
  };

  //! Parameterization of search symmetries.
  class search_symmetry
  {
    public:
      //! Default constructor. Some data members are not initialized!
      search_symmetry() {}

      //! Initialization given search_symmetry_flags and space_group_type.
      search_symmetry(
        search_symmetry_flags const& flags,
        space_group_type const& group_type)
      :
        flags_(flags)
      {
        init(group_type);
      }

      /*! \brief Initialization given search_symmetry_flags,
          space_group_type and structure_seminvariants.
       */
      search_symmetry(
        search_symmetry_flags const& flags,
        space_group_type const& group_type,
        structure_seminvariants const& seminvariant)
      :
        flags_(flags)
      {
        init(group_type, &seminvariant);
      }

      //! Flags as passed to the constructor.
      search_symmetry_flags const&
      flags() const { return flags_; }

      //! Subgroup of the search symmetry.
      /*! Result of performing group multiplication given the
          operations of the original space group and the discrete
          additional generators of the Euclidean normalizer, as
          indicates by flags(). The continuous allowed origin shifts
          of the structure seminvariants are not multiplied into
          subgroup(). See continuous_shifts().
       */
      space_group const&
      subgroup() const { return subgroup_; }

      //! Vectors specifying the continuous allowed origin shifts.
      /*! These are the structure seminvariant vectors with
          modulus zero.
       */
      af::small<scitbx::vec3<int>, 3> const&
      continuous_shifts() const { return continuous_shifts_; }

      //! True only if all continuous_shifts() are principal.
      /*! The three possible principal directions are
          (1,0,0), (0,1,0), (0,0,1).
       */
      bool
      continuous_shifts_are_principal() const
      {
        typedef scitbx::vec3<int> v;
        for(std::size_t i=0;i<continuous_shifts_.size();i++) {
          v const& s = continuous_shifts_[i];
          if (   s != v(1,0,0)
              && s != v(0,1,0)
              && s != v(0,0,1)) {
            return false;
          }
        }
        return true;
      }

      /*! \brief Flags indicating if a given direction is a continuous
          allowed origin shift.
       */
      /*! Union of continuous_shifts() != 0. Intended to
          be used only if continuous_shifts_are_principal(), but
          this is enforced only if assert_principal == true.
       */
      af::tiny<bool, 3>
      continuous_shift_flags(bool assert_principal=true) const
      {
        if (assert_principal) {
          CCTBX_ASSERT(continuous_shifts_are_principal());
        }
        af::tiny<bool, 3> result(false,false,false);
        typedef scitbx::vec3<int> v;
        for(std::size_t i=0;i<continuous_shifts_.size();i++) {
          v const& s = continuous_shifts_[i];
          for(std::size_t j=0;j<3;j++) {
            if (s[j]) result[j] = true;
          }
        }
        return result;
      }

      //! Projection of symmetry operations along continuous shifts.
      /*! An exception is thrown if continuous_shifts_are_principal()
          is false.
       */
      space_group
      projected_subgroup() const
      {
        CCTBX_ASSERT(continuous_shifts_are_principal());
        space_group result;
        for(std::size_t i_smx=1;i_smx<subgroup_.order_z();i_smx++) {
          rt_mx s = subgroup_(i_smx);
          for(std::size_t i_sh=0;i_sh<continuous_shifts_.size();i_sh++) {
            std::size_t i=0;
            for(;i<3;i++) {
              if (continuous_shifts_[i_sh][i] != 0) break;
            }
            for(std::size_t j=0;j<3;j++) {
              if (j != i) s.r().num()(i,j) = 0;
            }
            s.t().num()[i] = 0;
          }
          result.expand_smx(s);
        }
        return result;
      }

    protected:
      search_symmetry_flags flags_;
      space_group subgroup_;
      af::small<scitbx::vec3<int>, 3> continuous_shifts_;

      void
      init(
        space_group_type const& group_type,
        const structure_seminvariants* seminvariant=0)
      {
        if (flags_.use_space_group_symmetry()) {
          subgroup_ = group_type.group();
        }
        else if (   flags_.use_space_group_ltr() > 0
                 || (   flags_.use_space_group_ltr() == 0
                     && flags_.use_seminvariants())) {
          for(std::size_t i=1;i<group_type.group().n_ltr();i++) {
            subgroup_.expand_ltr(group_type.group().ltr(i));
          }
        }
        if (flags_.use_seminvariants()) {
          CCTBX_ASSERT(seminvariant != 0);
          af::small<ss_vec_mod, 3> const&
            ss = seminvariant->vectors_and_moduli();
          for(std::size_t i_ss=0;i_ss<ss.size();i_ss++) {
            if (ss[i_ss].m == 0) {
              continuous_shifts_.push_back(ss[i_ss].v);
            }
            else {
              subgroup_.expand_ltr(tr_vec(ss[i_ss].v, ss[i_ss].m)
                .new_denominator(subgroup_.t_den()));
            }
          }
        }
        if (flags_.use_normalizer_k2l() || flags_.use_normalizer_l2n()) {
          subgroup_.expand_smx(
            group_type.addl_generators_of_euclidean_normalizer(
              flags_.use_normalizer_k2l(),
              flags_.use_normalizer_l2n()).const_ref());
        }
      }
  };

}} // namespace cctbx::sgtbx

#endif // CCTBX_SGTBX_SEARCH_SYMMETRY_H
