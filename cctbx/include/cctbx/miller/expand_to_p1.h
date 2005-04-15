#ifndef CCTBX_MILLER_EXPAND_TO_P1_H
#define CCTBX_MILLER_EXPAND_TO_P1_H

#include <cctbx/miller/sym_equiv.h>

namespace cctbx { namespace miller {

  //! XXX OBSOLETE
  /*! Expands an array of Miller indices to P1 symmetry.
      <p>
      The symmetry operations are applied to each element
      of the input array.
      <p>
      If anomalous_flag == false, centric indices are treated in a
      special way: Friedel mates are suppressed. If N is the
      number of unique symmetrically equivalent indices for
      a given centric index, only N/2 indices will be generated.
      <p>
      See also: sym_equiv_indices
   */
  template <typename FloatType = double>
  class expand_to_p1
  {
    public:
      //! Default constructor. All output arrays have zero length.
      /*! Not available in Python.
       */
      expand_to_p1() {}

      //! Expands a list of Miller indices.
      expand_to_p1(
        sgtbx::space_group const& space_group,
        bool anomalous_flag,
        af::const_ref<index<> > const& indices)
      {
        for(std::size_t i_indices=0;i_indices<indices.size();i_indices++) {
          sym_equiv_indices eq(space_group, indices[i_indices]);
          af::shared<sym_equiv_index> p1_list = eq.p1_listing(anomalous_flag);
          for(sym_equiv_index const* p1=p1_list.begin();
                                     p1!=p1_list.end();
                                     p1++) {
            indices_.push_back(p1->h());
          }
        }
      }

      //! Indices in P1.
      af::shared<index<> > const&
      indices() const { return indices_; }

    protected:
      af::shared<index<> > indices_;
  };

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

  /*! \brief Expands an array of Miller indices and associated
      scalar data (e.g. bool, int, double) to P1 symmetry.
   */
  /*! See also: expand_to_p1_indices
   */
  template <typename ScalarType>
  struct expand_to_p1_scalar
  {
    expand_to_p1_scalar() {}

    expand_to_p1_scalar(
      sgtbx::space_group const& space_group,
      bool anomalous_flag,
      af::const_ref<index<> > const& indices_,
      af::const_ref<ScalarType> const& data_)
    {
      CCTBX_ASSERT(data_.size() == indices_.size());
      detail::expand_to_p1_generator generator(
        space_group, anomalous_flag, indices_);
      while (generator.incr()) {
        indices.push_back(generator.p1_index->h());
        data.push_back(data_[generator.i_index]);
      }
    }

    af::shared<index<> > indices;
    af::shared<ScalarType> data;
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
      observations (data, sigmas) to P1 symmetry.
   */
  /*! See also: expand_to_p1_indices
   */
  template <typename FloatType>
  struct expand_to_p1_obs
  {
    expand_to_p1_obs() {}

    expand_to_p1_obs(
      sgtbx::space_group const& space_group,
      bool anomalous_flag,
      af::const_ref<index<> > const& indices_,
      af::const_ref<FloatType> const& data_,
      af::const_ref<FloatType> const& sigmas_)
    {
      CCTBX_ASSERT(data_.size() == indices_.size());
      CCTBX_ASSERT(sigmas_.size() == indices_.size());
      detail::expand_to_p1_generator generator(
        space_group, anomalous_flag, indices_);
      while (generator.incr()) {
        indices.push_back(generator.p1_index->h());
        data.push_back(data_[generator.i_index]);
        sigmas.push_back(sigmas_[generator.i_index]);
      }
    }

    af::shared<index<> > indices;
    af::shared<FloatType> data;
    af::shared<FloatType> sigmas;
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
