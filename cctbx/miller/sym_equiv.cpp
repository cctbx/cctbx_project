/* Copyright (c) 2001-2002 The Regents of the University of California
   through E.O. Lawrence Berkeley National Laboratory, subject to
   approval by the U.S. Department of Energy.
   See files COPYRIGHT.txt and LICENSE.txt for further details.

   Revision history:
     2002 Oct: Refactored parts of miller/miller_lib.cpp (rwgk)
 */

#include <cctbx/miller/sym_equiv.h>
#include <cctbx/sgtbx/miller_ops.h>
#include <cctbx/sgtbx/reciprocal_space_reference_asu.h>

namespace cctbx { namespace miller {

  sym_equiv_indices::sym_equiv_indices(
    sgtbx::space_group const& space_group,
    index<> const& h_in)
  :
    t_den_(space_group.t_den()),
    order_p_(space_group.order_p()),
    ht_restriction_(-1)
  {
    using namespace sgtbx;
    for(std::size_t i_inv=0;i_inv<space_group.f_inv();i_inv++) {
      for(std::size_t i_smx=0;i_smx<space_group.n_smx();i_smx++) {
        rt_mx s = space_group(0, i_inv, i_smx);
        index<> hr = h_in * s.r();
        bool found = false;
        for(std::size_t i=0;i<indices_.size();i++) {
          if (indices_[i].hr() == hr) {
            found = true;
            break;
          }
        }
        if (!found) {
          add(sym_equiv_index(hr, ht_mod_1(h_in, s.t()), t_den_, false));
        }
      }
    }
    CCTBX_ASSERT((space_group.n_smx() * space_group.f_inv()) % indices_.size()
                 == 0);
    CCTBX_ASSERT(!is_centric() || indices_.size() % 2 == 0);
  }

  void sym_equiv_indices::add(sym_equiv_index const& eq)
  {
    indices_.push_back(eq);
    if (indices_.size()) {
      if (eq.hr() == -indices_[0].hr()) {
        CCTBX_ASSERT(ht_restriction_ < 0 || ht_restriction_ == eq.ht());
        ht_restriction_ = eq.ht();
      }
    }
  }

  sym_equiv_index
  sym_equiv_indices::operator()(
    std::size_t i_mate,
    std::size_t i_indices) const
  {
    if (i_mate >= f_mates(false) || i_indices >= indices_.size()) {
      throw error_index();
    }
    return indices_[i_indices].mate(i_mate);
  }

  sym_equiv_indices::index_mate_indices_decomposition
  sym_equiv_indices::decompose_index_mate_indices(std::size_t i) const
  {
    // i = i_mate * indices_.size() + i_indices
    if (i >= multiplicity(false)) {
      throw error_index();
    }
    return index_mate_indices_decomposition(i / indices_.size(),
                                            i % indices_.size());
  }

  sym_equiv_index
  sym_equiv_indices::operator()(std::size_t i) const
  {
    index_mate_indices_decomposition d = decompose_index_mate_indices(i);
    return operator()(d.i_mate, d.i_indices);
  }

  af::shared<sym_equiv_index>
  sym_equiv_indices::p1_listing(bool anomalous_flag) const
  {
    af::shared<sym_equiv_index> result;
    if (anomalous_flag) {
      result.assign(indices_.begin(), indices_.end());
    }
    else {
      if (is_centric()) result.reserve(indices_.size() / 2);
      else              result.reserve(indices_.size());
      for(std::size_t i=0;i<multiplicity(false);i++) {
        sym_equiv_index h_eq = operator()(i);
        if (sgtbx::reciprocal_space::is_in_reference_asu_1b(h_eq.h())) {
          result.push_back(h_eq);
        }
      }
      CCTBX_ASSERT(result.size() == result.capacity());
    }
    return result;
  }

}} // namespace cctbx::miller
