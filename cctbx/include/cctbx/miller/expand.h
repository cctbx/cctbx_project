// $Id$
/* Copyright (c) 2001 The Regents of the University of California through
   E.O. Lawrence Berkeley National Laboratory, subject to approval by the
   U.S. Department of Energy. See files COPYRIGHT.txt and
   cctbx/LICENSE.txt for further details.

   Revision history:
     2002 Jul: Created from fragments of cctbx/sgtbx/miller_asu.h (rwgk)
 */

#ifndef CCTBX_MILLER_EXPAND_H
#define CCTBX_MILLER_EXPAND_H

#include <cctbx/sgtbx/miller_asu.h>

namespace cctbx { namespace miller {

  //! Expand an array of Miller indices to P1 symmetry.
  /*! The symmetry operations are applied to each element
      of the input array in. The unique symmetrically
      equivalent indices are appended to the output array out.
      <p>
      If friedel_flag == true, centric indices are treated in a
      special way: Friedel mates are suppressed. If N is the
      number of unique symmetrically equivalent indices for
      a given centric index, only N/2 indices will be generated.
      <p>
      See also: class SymEquivMillerIndices
   */
  template <typename MillerIndexArrayType>
  void
  expand_to_p1(
    sgtbx::SpaceGroup const& SgOps,
    bool friedel_flag,
    MillerIndexArrayType const& h_in,
    MillerIndexArrayType& h_out)
  {
    for(std::size_t i_in = 0; i_in < h_in.size(); i_in++) {
      sgtbx::SymEquivMillerIndices
      h_seq = SgOps.getEquivMillerIndices(h_in[i_in]);
      af::shared<miller::SymEquivIndex>
      p1_listing = h_seq.p1_listing(friedel_flag);
      for (int i_eq = 0; i_eq < p1_listing.size(); i_eq++) {
        h_out.push_back(p1_listing[i_eq].H());
      }
    }
  }

  //! XXX
  template <typename MillerIndexArrayType,
            typename AmplitudeArrayType,
            typename PhaseArrayType>
  void
  expand_to_p1(
    sgtbx::SpaceGroup const& SgOps,
    bool friedel_flag,
    MillerIndexArrayType const& h_in,
    AmplitudeArrayType const& ampl_in,
    PhaseArrayType const& phase_in,
    MillerIndexArrayType& h_out,
    AmplitudeArrayType& ampl_out,
    PhaseArrayType& phase_out,
    bool phase_degrees = false)
  {
    cctbx_assert(h_in.size() == ampl_in.size() || ampl_in.size() == 0);
    cctbx_assert(h_in.size() == phase_in.size() || phase_in.size() == 0);
    for(std::size_t i_in = 0; i_in < h_in.size(); i_in++) {
      sgtbx::SymEquivMillerIndices
      h_seq = SgOps.getEquivMillerIndices(h_in[i_in]);
      af::shared<miller::SymEquivIndex>
      p1_listing = h_seq.p1_listing(friedel_flag);
      for (int i_eq = 0; i_eq < p1_listing.size(); i_eq++) {
        h_out.push_back(p1_listing[i_eq].H());
        if (ampl_in.size()) {
          ampl_out.push_back(ampl_in[i_in]);
        }
        if (phase_in.size()) {
          phase_out.push_back(
            p1_listing[i_eq].phase_eq(phase_in[i_in], phase_degrees));
        }
      }
    }
  }

}} // namespace cctbx::miller

#endif // CCTBX_MILLER_EXPAND_H
