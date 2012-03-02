#include <cctbx/sgtbx/space_group.h>
#include <cctbx/sgtbx/lattice_tr.h>
#include <cctbx/sgtbx/utils.h>
#include <cctbx/sgtbx/reference_settings.h>
#include <cctbx/sgtbx/rot_mx_info.h>

namespace cctbx { namespace sgtbx {

  void
  space_group::reset(int t_den)
  {
    n_lsl_ = 1;
    n_ssl_ = 1;
    f_inv_ = 1;
    ltr_.reset(t_den);
    inv_t_ = tr_vec(0);
    smx_.clear();
    smx_.push_back(rt_mx(1, t_den));
    is_tidy_ = false;
  }

  void
  space_group::add_inv(tr_vec const& new_inv_t)
  {
    if (is_centric()) {
      if (ltr_.add(inv_t_ - new_inv_t)) is_tidy_ = false;
      return;
    }
    inv_t_ = new_inv_t.mod_positive();
    f_inv_ = 2;
    is_tidy_ = false;
    if (!no_expand_) {
      for(std::size_t i=1;i<n_smx();i++) {
        if (ltr_.add(      smx_[i].r() * inv_t_
                     + 2 * smx_[i].t() - inv_t_)) {
          is_tidy_ = false;
        }
      }
    }
  }

  void
  space_group::add_smx(rt_mx const& new_smx)
  {
    rot_mx minus_r = -new_smx.r();
    for(std::size_t i=0;i<n_smx();i++) {
      if (smx_[i].r() == new_smx.r()) {
        if (ltr_.add(smx_[i].t() - new_smx.t())) is_tidy_ = false;
        return;
      }
      if (smx_[i].r() == minus_r) {
        add_inv(smx_[i].t() + new_smx.t());
        return;
      }
    }
    int d = new_smx.r().num().determinant();
    if (n_smx() >= smx_.capacity() || (d != -1 && d != 1))
      throw error_non_crystallographic_rotation_matrix_encountered();
    CCTBX_ASSERT(new_smx.t().den() == smx_[0].t().den());
    smx_.push_back(new_smx.mod_positive());
    rt_mx const& s = smx_.back();
    if (!no_expand_ && is_centric()) {
      ltr_.add(      s.r() * inv_t_
               + 2 * s.t() - inv_t_);
    }
    is_tidy_ = false;
  }

  space_group&
  space_group::expand_ltr(tr_vec const& new_ltr)
  {
    if (no_expand_) {
      if (ltr_.add(new_ltr)) is_tidy_ = false;
      return *this;
    }
    for (std::size_t i_smx=n_ssl_;i_smx<n_smx();i_smx++) {
      for (std::size_t i_ltr=1;i_ltr<n_lsl_;i_ltr++) {
        if (ltr_.add(smx_[i_smx].r() * ltr_[i_ltr])) is_tidy_ = false;
      }
    }
    n_ssl_ = n_smx();
    tr_vec trial_ltr = new_ltr;
    std::size_t i = n_lsl_;
    std::size_t j = 1;
    for (;;)
    {
      if (ltr_.add(trial_ltr)) is_tidy_ = false;
      for (std::size_t i_smx=1;i_smx<n_smx();i_smx++) {
        for (std::size_t i_ltr=n_lsl_;i_ltr<ltr_.size(); i_ltr++) {
          if (ltr_.add(smx_[i_smx].r() * ltr_[i_ltr])) is_tidy_ = false;
        }
      }
      n_lsl_ = ltr_.size();
      if (j > i) {
        i++;
        j = 1;
      }
      if (i == ltr_.size()) break;
      trial_ltr = ltr_[j] + ltr_[i];
      j++;
    }
    return *this;
  }

  space_group&
  space_group::expand_inv(tr_vec const& new_inv_t)
  {
    add_inv(new_inv_t);
    expand_ltr(tr_vec(0));
    return *this;
  }

  space_group&
  space_group::expand_smx(rt_mx const& new_smx)
  {
    if (new_smx.r().den() != 1) {
      throw error(
        "sgtbx::space_group::expand_smx():"
        " rotation-part denominator must be 1 (implementation limitation).");
    }
    if (new_smx.t().den() != smx_[0].t().den()) {
      throw error(
        "sgtbx::space_group::expand_smx():"
        " incompatible translation-part denominator.");
    }
    if (no_expand_) {
      add_smx(new_smx);
      return *this;
    }
    const rt_mx *pnew_smx = &new_smx; // XXX avoid pointer stuff
    rt_mx trial_smx;
    std::size_t i = n_smx();
    std::size_t j = 1;
    for (;;) {
      add_smx(*pnew_smx);
      if (j > i) {
        i++;
        j = 1;
      }
      if (i == n_smx()) break;
      trial_smx = smx_[j] * smx_[i];
      pnew_smx = &trial_smx;
      j++;
    }
    expand_ltr(tr_vec(0));
    return *this;
  }

  std::size_t
  space_group::expand_conventional_centring_type(char symbol)
  {
    using namespace lattice_tr::conventional_centring_types;
    const table_entry *e = get(symbol);
    if (e == 0) {
      throw error("Illegal symbol for centring type of cell.");
    }
    for(std::size_t i=0;i<e->n_translations;i++) {
      expand_ltr(e->t[i].new_denominator(t_den()));
    }
    return e->n_translations;
  }

  rt_mx
  space_group::operator()(
    std::size_t i_ltr,
    std::size_t i_inv,
    std::size_t i_smx) const
  {
    if (   i_ltr >= ltr_.size()
        || i_inv >= f_inv_
        || i_smx >= n_smx()) {
      throw error_index();
    }
    if (i_inv == 0) return smx_[i_smx] + ltr_[i_ltr];
    return -smx_[i_smx] + inv_t_ + ltr_[i_ltr];
  }

  rt_mx
  space_group::operator()(std::size_t i_op) const
  {
    // i_op = ((i_ltr * f_inv()) + i_inv) * n_smx() + i_smx
    if (i_op >= order_z()) {
      throw error_index();
    }
    std::size_t i_smx = i_op % n_smx();
    std::size_t i_inv = (i_op / n_smx()) % f_inv_;
    std::size_t i_ltr = i_op / (f_inv_ * n_smx());
    return operator()(i_ltr, i_inv, i_smx);
  }

  bool
  space_group::contains(rt_mx const& smx) const
  {
    CCTBX_ASSERT(smx.r().den() == r_den());
    CCTBX_ASSERT(smx.t().den() == t_den());
    rot_mx const& r = smx.r();
    for(std::size_t i=0;i<n_smx();i++) {
      if (smx_[i].r() == r) {
        return ltr_.contains(smx_[i].t() - smx.t());
      }
    }
    if (!is_centric()) return false;
    rot_mx minus_r = -smx.r();
    for(std::size_t i=0;i<n_smx();i++) {
      if (smx_[i].r() == minus_r) {
        return ltr_.contains(smx_[i].t() + smx.t() - inv_t_);
      }
    }
    return false;
  }

  bool
  space_group::operator==(space_group const& rhs) const
  {
    CCTBX_ASSERT(r_den() == rhs.r_den());
    CCTBX_ASSERT(t_den() == rhs.t_den());
    if (n_ltr() != rhs.n_ltr()) return false;
    if (f_inv() != rhs.f_inv()) return false;
    if (n_smx() != rhs.n_smx()) return false;
    space_group tidy_lhs = *this; tidy_lhs.make_tidy();
    space_group tidy_rhs = rhs; tidy_rhs.make_tidy();
    if (tidy_lhs.inv_t_ != tidy_rhs.inv_t_) return false;
    if (tidy_lhs.ltr_ != tidy_rhs.ltr_) return false;
    for(std::size_t i=0;i<tidy_lhs.n_smx();i++) {
      if (tidy_lhs.smx(i) != tidy_rhs.smx(i)) return false;
    }
    return true;
  }

  space_group
  space_group::build_derived_group(bool discard_z, bool add_inv) const
  {
    space_group result(false, t_den());
    if (!discard_z) {
      for(std::size_t i=0;i<n_ltr();i++) {
        result.expand_ltr(ltr_[i]);
      }
    }
    if (is_centric() || add_inv) {
      result.expand_inv(tr_vec(t_den()));
    }
    for(std::size_t i=0;i<n_smx();i++) {
      result.expand_smx(rt_mx(smx_[i].r(), t_den()));
    }
    return result;
  }

  space_group
  space_group::build_derived_acentric_group() const
  {
    if (!is_centric()) return *this;
    space_group tidy_sg = *this;
    tidy_sg.make_tidy();
    space_group result(false, tidy_sg.t_den());
    for(std::size_t i=0;i<tidy_sg.n_ltr();i++) {
      result.expand_ltr(tidy_sg.ltr_[i]);
    }
    for(std::size_t i=0;i<tidy_sg.n_smx();i++) {
      result.expand_smx(tidy_sg.smx_[i]);
    }
    return result;
  }

  space_group_symbols
  space_group::match_tabulated_settings() const
  {
    matrix_group::code point_group = point_group_type();
    space_group tidy_sg = *this;
    tidy_sg.make_tidy();
    space_group_symbol_iterator symbol_iter;
    space_group_symbols symbols;
    for (;;) {
      symbols = symbol_iter.next();
      if (symbols.number() == 0) return symbols;
      if (point_group
          != reference_settings::matrix_group_code_table(symbols.number())
               .point_group_type())
        continue;
      try {
        space_group tab_sg(symbols.hall(), true, false, false, t_den());
        if (tab_sg == tidy_sg) return symbols;
      }
      catch (error const&) {
        throw CCTBX_INTERNAL_ERROR();
      }
    }
  }

  namespace {

    class cmp_tr_vec
    {
      public:
        cmp_tr_vec() : civ_(3) {}

        bool operator()(tr_vec const& a, tr_vec const& b)
        {
          return civ_(a.num().begin(), b.num().begin());
        }

      private:
        utils::cmp_i_vec civ_;
    };

    bool first_is_shorter(sg_vec3 const& a, sg_vec3 const& b)
    {
      using scitbx::fn::absolute;
      for(std::size_t i=0;i<3;i++) {
        if (a[i]) {
          if (absolute(a[i]) > absolute(b[i])) return false;
          return true;
        }
      }
      return true;
    }

    af::shared<tr_vec>
    build_list_tot_tr(tr_group const& ltr, int t_den)
    {
      af::shared<tr_vec> tlt;
      for (std::size_t i_ltr=1;i_ltr<ltr.size();i_ltr++) {
        sg_vec3 n_utr(1,1,1);
        for(std::size_t i=0;i<3;i++) if (ltr[i_ltr][i]) n_utr[i] = 2;
        sg_vec3 unit_tr;
        for(unit_tr[0]=0;unit_tr[0]<n_utr[0];unit_tr[0]++)
        for(unit_tr[1]=0;unit_tr[1]<n_utr[1];unit_tr[1]++)
        for(unit_tr[2]=0;unit_tr[2]<n_utr[2];unit_tr[2]++)
        {
          tr_vec v = ltr[i_ltr]
                   - tr_vec(unit_tr, 1).new_denominator(ltr[0].den());
          v = v.new_denominator(t_den);
          std::size_t i = 0;
          for (;i<tlt.size();i++) {
            if (tlt[i].num().cross(v.num()) == 0) {
              if (!first_is_shorter(tlt[i].num(), v.num())) {
                tlt[i] = v;
              }
              break;
            }
          }
          if (i == tlt.size()) tlt.push_back(v);
        }
      }
      std::sort(tlt.begin(), tlt.end(), cmp_tr_vec());
      for(std::size_t i=0;i<3;i++) {
        tr_vec v(t_den);
        v[i] = t_den;
        tlt.push_back(v);
      }
      return tlt;
    }

  } // namespace <anonymous>

  change_of_basis_op
  space_group::construct_z2p_op(int r_den, int t_den) const
  {
    change_of_basis_op result;
    space_group primitive_sg(false, this->t_den());
    const int r_den_3 = r_den * r_den * r_den;
    af::shared<tr_vec> tlt = build_list_tot_tr(ltr_, r_den);
    std::size_t i_tlt[3];
    sg_mat3 basis;
    for (i_tlt[0] =           0; i_tlt[0] < tlt.size() - 2; i_tlt[0]++) {
      basis.set_column(0, tlt[i_tlt[0]].num());
    for (i_tlt[1] = i_tlt[0] + 1; i_tlt[1] < tlt.size() - 1; i_tlt[1]++) {
      basis.set_column(1, tlt[i_tlt[1]].num());
    for (i_tlt[2] = i_tlt[1] + 1; i_tlt[2] < tlt.size();     i_tlt[2]++) {
      basis.set_column(2, tlt[i_tlt[2]].num());
      int f = basis.determinant() * n_ltr();
      if (f == r_den_3 || -f == r_den_3) {
        if (f < 0) for(std::size_t i=0;i<3;i++) basis[i * 3] *= -1;
        try {
          result = change_of_basis_op(
            rt_mx(rot_mx(basis, r_den), t_den)).inverse();
          primitive_sg = change_basis(result);
        }
        catch (error const&) {
          continue;
        }
        if (primitive_sg.n_ltr() == 1) {
          CCTBX_ASSERT(
            result.c().r().num().determinant() == n_ltr() * r_den_3);
          return result;
        }
      }
    }}}
    throw CCTBX_INTERNAL_ERROR();
  }

  change_of_basis_op
  space_group::z2p_op(int r_den, int t_den) const
  {
    change_of_basis_op cb_op = ltr_.conventional_z2p_op(r_den, t_den);
    if (cb_op.is_valid()) return cb_op;
    return construct_z2p_op(r_den, t_den);
  }

  std::map<int, int>
  space_group::count_rotation_part_types() const
  {
    std::map<int, int> result;
    for(std::size_t i=0;i<n_smx();i++) {
      result[smx_[i].r().type()]++;
    }
    return result;
  }

  matrix_group::code
  space_group::point_group_type() const
  {
    using namespace matrix_group;

    std::map<int, int> counter = count_rotation_part_types();

    if      (counter[-3] + counter[3] == 8) {
      if      (n_smx() == 12) {
        if (!is_centric()) return code_23;
        else               return code_m3b;
      }
      else if (n_smx() == 24) {
        if (!is_centric()) {
          if (counter[ 4] == 6) return code_432;
          if (counter[-4] == 6) return code_4b3m;
        }
        else return code_m3bm;
      }
    }
    else if (counter[-6] + counter[6] == 2) {
      if      (n_smx() ==  6) {
        if (!is_centric()) {
          if (counter[ 6] == 2) return code_6;
          if (counter[-6] == 2) return code_6b;
        }
        else return code_6_m;
      }
      else if (n_smx() == 12) {
        if (!is_centric()) {
          if (counter[ 6] == 2) {
            if (counter[ 2] == 7) return code_622;
            if (counter[-2] == 6) return code_6mm;
          }
          else if (counter[-6] == 2) return code_6bm2;
        }
        else return code_6_mmm;
      }
    }
    else if (counter[-3] + counter[3] == 2) {
      if      (n_smx() ==  3) {
        if (!is_centric()) return code_3;
        else return code_3b;
      }
      else if (n_smx() ==  6) {
        if (!is_centric()) {
          if (counter[ 2] == 3) return code_32;
          if (counter[-2] == 3) return code_3m;
        }
        else return code_3bm;
      }
    }
    else if (counter[-4] + counter[4] == 2) {
      if (n_smx() ==  4) {
        if (!is_centric()) {
          if (counter[ 4] == 2) return code_4;
          if (counter[-4] == 2) return code_4b;
        }
        else return code_4_m;
      }
      else if (n_smx() ==  8) {
        if (!is_centric()) {
          if (counter[ 4] == 2) {
            if (counter[ 2] == 5) return code_422;
            if (counter[-2] == 4) return code_4mm;
          }
          else if (counter[-4] == 2) return code_4bm2;
        }
        else return code_4_mmm;
      }
    }
    else if (counter[-2] + counter[2] == 3) {
      if (!is_centric()) {
        if (counter[ 2] == 3) return code_222;
        if (counter[-2] == 2) return code_mm2;
      }
      else return code_mmm;
    }
    else if (counter[-2] + counter[2] == 1) {
      if (!is_centric()) {
        if (counter[ 2] == 1) return code_2;
        if (counter[-2] == 1) return code_m;
      }
      else return code_2_m;
    }
    else if (n_smx() == 1) {
      if (!is_centric()) return code_1;
      else return code_1b;
    }

    throw CCTBX_INTERNAL_ERROR();
  }

  namespace {

    struct cmp_ltr
    {
      bool operator()(tr_vec const& a, tr_vec const& b)
      {
        for(std::size_t i=0;i<3;i++) {
          if (a[i] < b[i]) return true;
          if (a[i] > b[i]) return false;
        }
        return false;
      }
    };

    struct cmp_smx
    {
      bool operator()(rt_mx const& a, rt_mx const& b)
      {
        using scitbx::fn::absolute;
        using utils::cmp_i_vec;
        rot_mx_info ri_a(a.r());
        rot_mx_info ri_b(b.r());
        if (absolute(ri_a.type()) > absolute(ri_b.type())) return true;
        if (absolute(ri_a.type()) < absolute(ri_b.type())) return false;
        if (ri_a.type() > ri_b.type()) return true;
        if (ri_a.type() < ri_b.type()) return false;
        if (cmp_i_vec(3)(ri_a.ev().begin(), ri_b.ev().begin())) return true;
        if (cmp_i_vec(3)(ri_b.ev().begin(), ri_a.ev().begin())) return false;
        if (ri_a.sense() > ri_b.sense()) return true;
        if (ri_a.sense() < ri_b.sense()) return false;
        if (a.t().den() == b.t().den()) {
          if (cmp_i_vec(3)(
            a.t().num().begin(), b.t().num().begin())) return true;
          if (cmp_i_vec(3)(
            b.t().num().begin(), a.t().num().begin())) return false;
        }
        else {
          throw std::runtime_error("Not implemented.");
        }
        for(std::size_t i=0;i<9;i++) {
          if (a.r()[i] < b.r()[i]) return true;
          if (a.r()[i] > b.r()[i]) return false;
        }
        if (a.t().den() == b.t().den()) {
          for(std::size_t i=0;i<3;i++) {
            if (a.t()[i] < b.t()[i]) return true;
            if (a.t()[i] > b.t()[i]) return false;
          }
        }
        else {
          throw std::runtime_error("Not implemented.");
        }
        return false;
      }
    };

  } // namespace <anonymous>

  bool
  rt_mx::operator<(rt_mx const& rhs) const
  {
    return cmp_smx()(*this, rhs);
  }

  bool
  rt_mx::operator>(rt_mx const& rhs) const
  {
    return cmp_smx()(rhs, *this);
  }

  space_group&
  space_group::make_tidy()
  {
    if (is_tidy_) return *this;
    if (is_centric()) {
      inv_t_ = inv_t(true);
      for (std::size_t i=1;i<n_smx();i++) {
        if (smx_[i].r().num().determinant() < 0) {
          smx_[i] = smx_[i].pre_multiply_inv_t(inv_t_);
        }
      }
    }
    for (std::size_t i=1;i<n_smx();i++) {
      smx_[i] = rt_mx(smx_[i].r(), ltr_.tidy(smx_[i].t()));
    }
    if (n_ltr() > 2)
      std::sort(ltr_.elems_.begin() + 1, ltr_.elems_.end(), cmp_ltr());
    if (n_smx() > 2)
      std::sort(smx_.begin() + 1, smx_.begin() + n_smx(), cmp_smx());
    is_tidy_ = true;
    return *this;
  }

  space_group
  space_group::change_basis(const change_of_basis_op& cb_op) const
  {
    space_group result(no_expand_, t_den());
    result.ltr_ = ltr_.change_basis(cb_op);
    if (is_centric()) {
      result.expand_inv(cb_op(inv_t_, -1));
    }
    for (std::size_t i=1;i<n_smx();i++) {
      result.expand_smx(cb_op(smx_[i]));
    }
    return result;
  }

  change_of_basis_op
  space_group::change_of_origin_realising_origin_centricity() const
  {
    if(!is_centric() || is_origin_centric()) {
      return change_of_basis_op();
    }
    return change_of_basis_op(rt_mx(
      tr_vec(inv_t().num(), -inv_t().den()*2).new_denominator(cb_t_den),
      cb_r_den));
  }

  bool
  space_group::is_chiral() const
  {
    if (is_centric()) return false;
    for (std::size_t i=1;i<n_smx();i++) {
      if (smx_[i].r().type() < 0) return false;
    }
    return true;
  }

  af::shared<rt_mx>
  space_group::all_ops(int mod, bool cancel) const
  {
    af::shared<rt_mx> result((af::reserve(order_z())));
    for(std::size_t i=0;i<order_z();i++) {
      rt_mx s = (*this)(i);
      if (cancel) s = s.cancel();
      if      (mod > 0) s.mod_positive_in_place();
      else if (mod < 0) s.mod_short_in_place();
      result.push_back(s);
    }
    return result;
  }

  af::shared<rt_mx>
  space_group::unique(rt_mx const& special_op) const
  {
    if (special_op.is_unit_mx()) return all_ops(1, true);
    af::shared<rt_mx> result;
    for(std::size_t i=0;i<order_z();i++) {
      rt_mx s = (*this)(i).multiply(special_op).mod_positive();
      if (std::find(result.begin(), result.end(), s) == result.end()) {
        result.push_back(s);
      }
    }
    return result;
  }

  af::shared<double>
  space_group::nearest_valid_phases(
    af::const_ref<miller::index<> > const& miller_indices,
    af::const_ref<double> const& phases,
    bool deg) const
  {
    CCTBX_ASSERT(phases.size() == miller_indices.size());
    af::shared<double> result((af::reserve(miller_indices.size())));
    for(std::size_t i=0;i<miller_indices.size();i++) {
      miller::index<>
#if !defined(CCTBX_SGTBX_PHASE_INFO_APPLE_LLVM2335_WORKAROUND)
      const&
#endif
      h = miller_indices[i];
      result.push_back(
        phase_restriction(h).nearest_valid_phase(phases[i], deg));
    }
    return result;
  }

  int
  space_group::multiplicity(
    vec3_rat const& site) const
  {
    int result = 0;
    vec3_rat pivot = site;
    ltr_.find_best_equiv_in_place(pivot);
    for(unsigned i=0;i<n_smx();i++) {
      vec3_rat sym_site = smx_[i] * site;
      ltr_.find_best_equiv_in_place(sym_site);
      if (sym_site == pivot) result++;
      if (is_centric()) {
        for(unsigned j=0;j<3;j++) {
          sym_site[j] *= -1;
          sym_site[j] += rat(inv_t_.num()[j], inv_t_.den());
        }
        ltr_.find_best_equiv_in_place(sym_site);
        if (sym_site == pivot) result++;
      }
    }
    CCTBX_ASSERT(order_z() % result == 0);
    return order_z() / result;
  }

}} // namespace cctbx::sgtbx
