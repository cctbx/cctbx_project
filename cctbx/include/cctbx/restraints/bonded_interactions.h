#ifndef CCTBX_RESTRAINTS_BONDED_INTERACTIONS_H
#define CCTBX_RESTRAINTS_BONDED_INTERACTIONS_H

#include <cctbx/import_scitbx_af.h>
#include <scitbx/array_family/ref.h>
#include <set>

namespace cctbx { namespace restraints {

  class bonded_interactions
  {
    protected:
      typedef std::set<unsigned>::const_iterator set_const_iter;

    public:
      bonded_interactions() {}

      bonded_interactions(
        af::const_ref<std::set<unsigned> > const& bond_sets,
        unsigned i_seq_0,
        bool eliminate_redundant_interactions=false)
      :
        i_seq_0_(i_seq_0),
        eliminate_redundant_interactions_(eliminate_redundant_interactions)
      {
        find_0_1(bond_sets);
        for(set_const_iter
            p=interactions_1_2_.begin(); p!=interactions_1_2_.end(); p++) {
          find_1_2(bond_sets, *p);
        }
        for(set_const_iter
            p=interactions_1_3_.begin(); p!=interactions_1_3_.end(); p++) {
          find_2_3(bond_sets, *p);
        }
        if (eliminate_redundant_interactions_) {
          do_eliminate_redundant_interactions(interactions_1_2_);
          do_eliminate_redundant_interactions(interactions_1_3_);
          do_eliminate_redundant_interactions(interactions_1_4_);
        }
      }

      unsigned
      i_seq_0() const { return i_seq_0_; }

      bool
      eliminate_redundant_interactions() const
      {
        return eliminate_redundant_interactions_;
      }

      std::set<unsigned> const&
      interactions_1_2() const { return interactions_1_2_; }

      std::set<unsigned> const&
      interactions_1_3() const { return interactions_1_3_; }

      std::set<unsigned> const&
      interactions_1_4() const { return interactions_1_4_; }

      int
      interaction_type_of(unsigned j_seq)
      {
        if (interactions_1_2_.find(j_seq) != interactions_1_2_.end()) return 2;
        if (interactions_1_3_.find(j_seq) != interactions_1_3_.end()) return 3;
        if (interactions_1_4_.find(j_seq) != interactions_1_4_.end()) return 4;
        return 0;
      }

    protected:
      bool eliminate_redundant_interactions_;
      unsigned i_seq_0_;
      std::set<unsigned> interactions_1_2_;
      std::set<unsigned> interactions_1_3_;
      std::set<unsigned> interactions_1_4_;

      void
      do_eliminate_redundant_interactions(std::set<unsigned>& interactions)
      {
        interactions.erase(
          interactions.begin(),
          interactions.upper_bound(i_seq_0_));
      }

      void
      find_0_1(
        af::const_ref<std::set<unsigned> > const& bond_sets)
      {
        std::set<unsigned> const& bonds = bond_sets[i_seq_0_];
        for(set_const_iter p=bonds.begin(); p!=bonds.end(); p++) {
          unsigned i_seq_1 = *p;
          if (i_seq_1 != i_seq_0_) {
            interactions_1_2_.insert(i_seq_1);
          }
        }
      }

      void
      find_1_2(
        af::const_ref<std::set<unsigned> > const& bond_sets,
        unsigned i_seq_1)
      {
        std::set<unsigned> const& bonds = bond_sets[i_seq_1];
        for(set_const_iter p=bonds.begin(); p!=bonds.end(); p++) {
          unsigned i_seq_2 = *p;
          if (   i_seq_2 != i_seq_0_
              && interactions_1_2_.find(i_seq_2) == interactions_1_2_.end()) {
            interactions_1_3_.insert(i_seq_2);
          }
        }
      }

      void
      find_2_3(
        af::const_ref<std::set<unsigned> > const& bond_sets,
        unsigned i_seq_2)
      {
        std::set<unsigned> const& bonds = bond_sets[i_seq_2];
        for(set_const_iter p=bonds.begin(); p!=bonds.end(); p++) {
          unsigned i_seq_3 = *p;
          if (   i_seq_3 != i_seq_0_
              && interactions_1_2_.find(i_seq_3) == interactions_1_2_.end()
              && interactions_1_3_.find(i_seq_3) == interactions_1_3_.end()) {
            interactions_1_4_.insert(i_seq_3);
          }
        }
      }
  };

}} // namespace cctbx::restraints

#endif // CCTBX_RESTRAINTS_BONDED_INTERACTIONS_H
