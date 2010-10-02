#ifndef CCTBX_SGTBX_TR_GROUP_H
#define CCTBX_SGTBX_TR_GROUP_H

#include <cctbx/sgtbx/change_of_basis_op.h>
#include <vector>

namespace cctbx { namespace sgtbx {

  class tr_group
  {
    friend class space_group;

    public:
      explicit
      tr_group(int t_den=sg_t_den)
      {
        elems_.push_back(tr_vec(t_den));
      }

      void reset(int t_den=sg_t_den)
      {
        elems_.clear();
        elems_.push_back(tr_vec(t_den));
      }

      std::size_t size() const { return elems_.size(); }

      std::vector<tr_vec> const& elems() const { return elems_; }

      std::vector<tr_vec>&       elems()       { return elems_; }

      tr_vec const& operator[](std::size_t i) const { return elems_[i]; }

      tr_vec&       operator[](std::size_t i)       { return elems_[i]; }

      int t_den() const { return elems_[0].den(); }

      bool
      contains(tr_vec const& t) const;

      bool operator==(tr_group const& rhs)
      {
        return this->elems_ == rhs.elems_;
      }

      bool operator!=(tr_group const& rhs)
      {
        return !((*this) == rhs);
      }

      bool expand(tr_vec const& new_t);

      tr_group change_basis(change_of_basis_op const& cb_op) const;

      char conventional_centring_type_symbol() const;

      change_of_basis_op
      conventional_z2p_op(int r_den=cb_r_den, int t_den=cb_t_den) const;

      tr_vec tidy(tr_vec const& t) const;

      void
      find_best_equiv_in_place(
        vec3_rat& t) const;

    private:
      std::vector<tr_vec> elems_;
      bool add(tr_vec const& new_t);
  };

}} // namespace cctbx::sgtbx

#endif // CCTBX_SGTBX_TR_GROUP_H
