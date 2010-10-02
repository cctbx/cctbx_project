#include <cctbx/sgtbx/tr_group.h>
#include <cctbx/sgtbx/lattice_tr.h>
#include <cctbx/sgtbx/utils.h>

namespace cctbx { namespace sgtbx {

  bool
  tr_group::contains(tr_vec const& t) const
  {
    if (!t.is_valid()) return false;
    return (   std::find(elems_.begin(), elems_.end(), t.mod_positive())
            != elems_.end());
  }

  bool tr_group::add(tr_vec const& new_t)
  {
    if (!new_t.is_valid()) return false;
    tr_vec t = new_t.mod_positive();
    if (std::find(elems_.begin(), elems_.end(), t) != elems_.end()) {
      return false;
    }
    CCTBX_ASSERT(t.den() == elems_[0].den());
    elems_.push_back(t);
    return true;
  }

  bool tr_group::expand(tr_vec const& new_t)
  {
    tr_vec const *p_new_t = &new_t;
    tr_vec trial_t; // XXX avoid the pointer stuff: trial_t = new_t;
    int old_size = size();
    int i = old_size;
    int j = 1;
    for (;;) {
      add(*p_new_t);
      if (j > i) {
        i++;
        j = 1;
      }
      if (i == size()) break;
      trial_t = elems_[j] + elems_[i];
      p_new_t = &trial_t;
      j++;
    }
    return size() != old_size;
  }

  tr_group tr_group::change_basis(change_of_basis_op const& cb_op) const
  {
    int t_den = elems_[0].den();
    tr_group result(t_den);
    for(std::size_t i=0;i<3;i++) {
      tr_vec bv(t_den);
      bv[i] = t_den;
      result.expand(cb_op(bv, 1));
    }
    for(std::size_t i=1;i<size();i++) {
      result.expand(cb_op(elems_[i], 1));
    }
    return result;
  }

  char
  tr_group::conventional_centring_type_symbol() const
  {
    using namespace lattice_tr::conventional_centring_types;
    for (const table_entry* e = table(); e->symbol != '\0'; e++) {
      if (e->n_translations == size()) {
        af::small<bool, 4> match_flags(size(), false);
        std::size_t n_matches = 0;
        for(std::size_t i=0;i<size();i++) {
          for(std::size_t j=0;j<size();j++) {
            if (!match_flags[j]
                && e->t[i].new_denominator(t_den()) == elems_[j]) {
              match_flags[j] = true;
              n_matches++;
              break;
            }
          }
        }
        if (n_matches == size()) return e->symbol;
      }
    }
    return '\0';
  }

  change_of_basis_op
  tr_group::conventional_z2p_op(int r_den, int t_den) const
  {
    char z_symbol = conventional_centring_type_symbol();
    rot_mx const&
      z2p_mx = lattice_tr::conventional_z2p_matrices::get(z_symbol);
    if (!z2p_mx.is_valid()) return change_of_basis_op(0, 0);
    return change_of_basis_op(
      rt_mx(z2p_mx.new_denominator(r_den), tr_vec(t_den)));
  }

  tr_vec
  tr_group::tidy(tr_vec const& t) const
  {
    int lcm_den = boost::lcm(elems_[0].den(), t.den());
    int f_den = lcm_den / elems_[0].den();
    tr_vec t_lcm = t.scale(lcm_den / t.den());
    tr_vec t_best = t_lcm.mod_short();
    for(std::size_t i=1;i<size();i++) {
      tr_vec t_trial = (t_lcm + elems_[i].scale(f_den)).mod_short();
      if (utils::cmp_i_vec(3)(t_trial.num().begin(), t_best.num().begin())) {
        t_best = t_trial;
      }
    }
    return t_best.new_denominator(t.den()).mod_positive();
  }

  void
  tr_group::find_best_equiv_in_place(
    vec3_rat& t) const
  {
    vec3_rat best;
    vec3_rat curr;
    for(std::size_t i=0;i<size();i++) {
      tr_vec const& lt = elems_[i];
      bool update_best = false;
      for(unsigned j=0;j<3;j++) {
        rat tt(lt.num()[j], lt.den());
        tt += t[j];
        int d = tt.denominator();
        tt = rat(scitbx::math::mod_positive(tt.numerator(), d), d);
        if (i == 0) {
          best[j] = tt;
        }
        else {
          if (!update_best) {
            if (best[j] < tt) {
              break;
            }
            update_best = (tt < best[j]);
          }
          curr[j] = tt;
        }
      }
      if (update_best) {
        best = curr;
      }
    }
    t = best;
  }

}} // namespace cctbx::sgtbx
