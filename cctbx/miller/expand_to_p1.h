#ifndef CCTBX_MILLER_EXPAND_TO_P1_H
#define CCTBX_MILLER_EXPAND_TO_P1_H

#include <cctbx/miller/sym_equiv.h>

namespace cctbx { namespace miller {

  namespace detail {

    struct expand_to_p1_generator
    {
      expand_to_p1_generator() {}

      expand_to_p1_generator(
        sgtbx::space_group const& space_group,
        bool anomalous_flag,
        af::const_ref<index<> > const& indices)
      :
        space_group_(&space_group),
        indices_(indices),
        anomalous_flag_(anomalous_flag),
        start(true)
      {}

      bool
      incr()
      {
        if (!start) goto continue_after_return;
        start = false;
        for(i_index=0;i_index<indices_.size();i_index++) {
          sym_equiv = sym_equiv_indices(*space_group_, indices_[i_index]);
          p1_listing = sym_equiv.p1_listing(anomalous_flag_);
          p1_listing_ref = p1_listing.const_ref();
          for(p1_index=p1_listing_ref.begin();
              p1_index!=p1_listing_ref.end();
              p1_index++) {
            return true;
            continue_after_return:;
          }
        }
        start = true;
        return false;
      }

      const sgtbx::space_group* space_group_;
      af::const_ref<index<> > indices_;
      bool anomalous_flag_;
      bool start;
      std::size_t i_index;
      sym_equiv_indices sym_equiv;
      af::shared<sym_equiv_index> p1_listing;
      af::const_ref<sym_equiv_index> p1_listing_ref;
      sym_equiv_index const* p1_index;
    };

  } // namespace detail

  //! Expands an array of Miller indices to P1 symmetry.
  /*! The symmetry operations are applied to each element
      of the input array.
      <p>
      If anomalous_flag == false, centric indices are treated in a
      special way: Friedel mates are suppressed. If N is the
      number of unique symmetrically equivalent indices for
      a given centric index, only N/2 indices will be generated.
      <p>
      See also: sym_equiv_indices
   */
  inline
  af::shared<index<> >
  expand_to_p1_indices(
    sgtbx::space_group const& space_group,
    bool anomalous_flag,
    af::const_ref<index<> > const& indices)
  {
    af::shared<index<> > result;
    detail::expand_to_p1_generator generator(
      space_group, anomalous_flag, indices);
    while (generator.incr()) {
      result.push_back(generator.p1_index->h());
    }
    return result;
  }

  //! Expands an array of Miller indices to P1 symmetry.
  /*! For details see: expand_to_p1_indices
      <p>
      If build_iselection == true, the iselection array keeps
      track of integer indices into the original array.
   */
  struct expand_to_p1_iselection
  {
    expand_to_p1_iselection() {}

    expand_to_p1_iselection(
      sgtbx::space_group const& space_group,
      bool anomalous_flag,
      af::const_ref<index<> > const& indices_,
      bool build_iselection)
    {
      std::size_t n_reserve = indices_.size() * space_group.order_p();
      indices.reserve(n_reserve);
      if (build_iselection) {
        iselection.reserve(n_reserve);
      }
      detail::expand_to_p1_generator generator(
        space_group, anomalous_flag, indices_);
      while (generator.incr()) {
        indices.push_back(generator.p1_index->h());
        if (build_iselection) {
          iselection.push_back(generator.i_index);
        }
      }
    }

    af::shared<index<> > indices;
    af::shared<std::size_t> iselection;
  };

  /*! \brief Expands an array of Miller indices and associated
      complex data to P1 symmetry.
   */
  /*! See also: expand_to_p1_indices
   */
  template <typename FloatType>
  struct expand_to_p1_complex
  {
    expand_to_p1_complex() {}

    expand_to_p1_complex(
      sgtbx::space_group const& space_group,
      bool anomalous_flag,
      af::const_ref<index<> > const& indices_,
      af::const_ref<std::complex<FloatType> > const& data_)
    {
      CCTBX_ASSERT(data_.size() == indices_.size());
      detail::expand_to_p1_generator generator(
        space_group, anomalous_flag, indices_);
      while (generator.incr()) {
        indices.push_back(generator.p1_index->h());
        data.push_back(generator.p1_index->complex_eq(
          data_[generator.i_index]));
      }
    }

    af::shared<index<> > indices;
    af::shared<std::complex<FloatType> > data;
  };

  /*! \brief Expands an array of Miller indices and associated
      Hendrickson-Lattman coefficients to P1 symmetry.
   */
  /*! See also: expand_to_p1_indices
   */
  template <typename FloatType>
  struct expand_to_p1_hendrickson_lattman
  {
    expand_to_p1_hendrickson_lattman() {}

    expand_to_p1_hendrickson_lattman(
      sgtbx::space_group const& space_group,
      bool anomalous_flag,
      af::const_ref<index<> > const& indices_,
      af::const_ref<hendrickson_lattman<FloatType> > const& data_)
    {
      CCTBX_ASSERT(data_.size() == indices_.size());
      detail::expand_to_p1_generator generator(
        space_group, anomalous_flag, indices_);
      while (generator.incr()) {
        indices.push_back(generator.p1_index->h());
        data.push_back(generator.p1_index->hendrickson_lattman_eq(
          data_[generator.i_index]));
      }
    }

    af::shared<index<> > indices;
    af::shared<hendrickson_lattman<FloatType> > data;
  };

  /*! \brief Expands an array of Miller indices and associated
      phases to P1 symmetry.
   */
  /*! See also: expand_to_p1_indices
   */
  template <typename FloatType>
  struct expand_to_p1_phases
  {
    expand_to_p1_phases() {}

    expand_to_p1_phases(
      sgtbx::space_group const& space_group,
      bool anomalous_flag,
      af::const_ref<index<> > const& indices_,
      af::const_ref<FloatType> const& data_,
      bool deg)
    {
      CCTBX_ASSERT(data_.size() == indices_.size());
      detail::expand_to_p1_generator generator(
        space_group, anomalous_flag, indices_);
      while (generator.incr()) {
        indices.push_back(generator.p1_index->h());
        data.push_back(generator.p1_index->phase_eq(
          data_[generator.i_index], deg));
      }
    }

    af::shared<index<> > indices;
    af::shared<FloatType> data;
  };

}} // namespace cctbx::miller

#endif // CCTBX_MILLER_EXPAND_TO_P1_H
