/* Copyright (c) 2001-2002 The Regents of the University of California
   through E.O. Lawrence Berkeley National Laboratory, subject to
   approval by the U.S. Department of Energy.
   See files COPYRIGHT.txt and LICENSE.txt for further details.

   Revision history:
     2002 Oct: Refactored (rwgk)
     2002 Jul: Created from fragments of cctbx/sgtbx/miller_asu.h (rwgk)
 */

#ifndef CCTBX_MILLER_EXPAND_TO_P1_H
#define CCTBX_MILLER_EXPAND_TO_P1_H

#include <cctbx/miller/sym_equiv.h>

namespace cctbx { namespace miller {

  //! Expands an array of Miller indices and associated data to P1 symmetry.
  /*! The symmetry operations are applied to each element
      of the input array in. The unique symmetrically
      equivalent indices and associated data are available through
      member functions. Unused output arrays have length zero.
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
        af::const_ref<index<> > const& indices);

      //! Expands a list of Miller indices and associated amplitudes.
      expand_to_p1(
        sgtbx::space_group const& space_group,
        bool anomalous_flag,
        af::const_ref<index<> > const& indices,
        af::const_ref<FloatType> const& amplitudes);

      //! Expands a list of Miller indices and associated phases.
      expand_to_p1(
        sgtbx::space_group const& space_group,
        bool anomalous_flag,
        af::const_ref<index<> > const& indices,
        af::const_ref<FloatType> const& phases, bool phase_degrees);

      /*! \brief Expands a list of Miller indices and associated amplitudes
          and phases.
       */
      expand_to_p1(
        sgtbx::space_group const& space_group,
        bool anomalous_flag,
        af::const_ref<index<> > const& indices,
        af::const_ref<FloatType> const& amplitudes,
        af::const_ref<FloatType> const& phases, bool phase_degrees=false);

      //! \brief Expands a list of structure factors.
      expand_to_p1(
        sgtbx::space_group const& space_group,
        bool anomalous_flag,
        af::const_ref<index<> > const& indices,
        af::const_ref<std::complex<FloatType> > const& structure_factors);

      //! Indices in P1.
      af::shared<index<> > const&
      indices() const { return indices_; }

      //! Amplitudes in P1.
      af::shared<FloatType> const&
      amplitudes() const { return amplitudes_; }

      //! Phases in P1.
      af::shared<FloatType> const&
      phases() const { return phases_; }

      //! Structure factors in P1.
      af::shared<std::complex<FloatType> > const&
      structure_factors() const { return structure_factors_; }

    protected:
      af::shared<index<> > indices_;
      af::shared<FloatType> amplitudes_;
      af::shared<FloatType> phases_;
      af::shared<std::complex<FloatType> > structure_factors_;

      void
      work(
        sgtbx::space_group const& space_group,
        bool anomalous_flag,
        af::const_ref<index<> > const& indices,
        af::const_ref<FloatType> const& amplitudes,
        af::const_ref<FloatType> const& phases,
        bool phase_degrees,
        af::const_ref<std::complex<FloatType> > const& structure_factors);
  };

  template <typename FloatType>
  expand_to_p1<FloatType>::expand_to_p1(
    sgtbx::space_group const& space_group,
    bool anomalous_flag,
    af::const_ref<index<> > const& indices)
  {
    work(space_group, anomalous_flag, indices,
         af::const_ref<FloatType>(0,0),
         af::const_ref<FloatType>(0,0), false,
         af::const_ref<std::complex<FloatType> >(0,0));
  }

  template <typename FloatType>
  expand_to_p1<FloatType>::expand_to_p1(
    sgtbx::space_group const& space_group,
    bool anomalous_flag,
    af::const_ref<index<> > const& indices,
    af::const_ref<FloatType> const& amplitudes)
  {
    work(space_group, anomalous_flag, indices,
         amplitudes,
         af::const_ref<FloatType>(0,0), false,
         af::const_ref<std::complex<FloatType> >(0,0));
  }

  template <typename FloatType>
  expand_to_p1<FloatType>::expand_to_p1(
    sgtbx::space_group const& space_group,
    bool anomalous_flag,
    af::const_ref<index<> > const& indices,
    af::const_ref<FloatType> const& phases, bool phase_degrees)
  {
    work(space_group, anomalous_flag, indices,
         af::const_ref<FloatType>(0,0),
         phases, phase_degrees,
         af::const_ref<std::complex<FloatType> >(0,0));
  }

  template <typename FloatType>
  expand_to_p1<FloatType>::expand_to_p1(
    sgtbx::space_group const& space_group,
    bool anomalous_flag,
    af::const_ref<index<> > const& indices,
    af::const_ref<FloatType> const& amplitudes,
    af::const_ref<FloatType> const& phases,
    bool phase_degrees)
  {
    work(space_group, anomalous_flag, indices,
         amplitudes, phases, phase_degrees,
         af::const_ref<std::complex<FloatType> >(0,0));
  }

  template <typename FloatType>
  expand_to_p1<FloatType>::expand_to_p1(
    sgtbx::space_group const& space_group,
    bool anomalous_flag,
    af::const_ref<index<> > const& indices,
    af::const_ref<std::complex<FloatType> > const& structure_factors)
  {
    work(space_group, anomalous_flag, indices,
         af::const_ref<FloatType>(0,0),
         af::const_ref<FloatType>(0,0), false,
         structure_factors);
  }

  template <typename FloatType>
  void
  expand_to_p1<FloatType>::work(
    sgtbx::space_group const& space_group,
    bool anomalous_flag,
    af::const_ref<index<> > const& indices,
    af::const_ref<FloatType> const& amplitudes,
    af::const_ref<FloatType> const& phases, bool phase_degrees,
    af::const_ref<std::complex<FloatType> > const& structure_factors)
  {
    CCTBX_ASSERT(   amplitudes.size() == indices.size()
                 || amplitudes.size() == 0);
    CCTBX_ASSERT(   phases.size() == indices.size()
                 || phases.size() == 0);
    CCTBX_ASSERT(   structure_factors.size() == indices.size()
                 || structure_factors.size() == 0);
    for(std::size_t i_indices=0;i_indices<indices.size();i_indices++) {
      sym_equiv_indices eq(space_group, indices[i_indices]);
      af::shared<sym_equiv_index> p1_list = eq.p1_listing(anomalous_flag);
      for(sym_equiv_index const* p1=p1_list.begin();p1!=p1_list.end();p1++) {
        indices_.push_back(p1->h());
        if (amplitudes.size()) {
          amplitudes_.push_back(amplitudes[i_indices]);
        }
        if (phases.size()) {
          phases_.push_back(
            p1->phase_eq(phases[i_indices], phase_degrees));
        }
        if (structure_factors.size()) {
          structure_factors_.push_back(
            p1->complex_eq(structure_factors[i_indices]));
        }
      }
    }
  }

}} // namespace cctbx::miller

#endif // CCTBX_MILLER_EXPAND_TO_P1_H
