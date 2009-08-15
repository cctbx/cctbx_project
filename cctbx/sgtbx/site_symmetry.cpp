#include <cctbx/sgtbx/site_symmetry_table.h>

namespace cctbx { namespace sgtbx {

  rt_point_group::rt_point_group(sgtbx::space_group const& sg,
                                 rt_mx const& projection)
  :
    is_valid_(true)
  {
    reset(sg(0));
    if (projection.is_unit_mx()) return;
    rt_mx pcm = projection.cancel().mod_positive();
    for(std::size_t i_op=0;i_op<sg.order_z();i_op++) {
      rt_mx s = sg(i_op);
      rt_mx sp = s.multiply(projection);
      if (sp.mod_positive() == pcm) {
        tr_vec unit_t = projection.t().minus(sp.t());
        CCTBX_ASSERT(unit_t.mod_positive().num().is_zero());
        expand(rt_mx(s.r(), s.t() + unit_t.new_denominator(s.t().den())));
      }
    }
  }

  void rt_point_group::reset(rt_mx const& s)
  {
    matrices_.clear();
    matrices_.push_back(s);
  }

  void rt_point_group::add(rt_mx const& s)
  {
    for (const rt_mx* i=matrices_.begin(); i != matrices_.end(); i++) {
      if (i->r() == s.r()) {
        if (i->t() != s.t()) is_valid_ = false;
        return;
      }
    }
    matrices_.push_back(s);
  }

  void rt_point_group::expand(rt_mx const& s)
  {
    rt_mx trial_s = s;
    std::size_t i = matrices_.size();
    std::size_t j = 1;
    for (;;) {
      add(trial_s);
      if (!is_valid_) return;
      if (j > i) {
        i++;
        j = 1;
      }
      if (i == matrices_.size()) break;
      trial_s = matrices_[j] * matrices_[i];
      j++;
    }
  }

  bool rt_point_group::try_expand(rt_mx const& s)
  {
    std::size_t size_before_expand = matrices_.size();
    expand(s);
    if (!is_valid_) {
      matrices_.resize(size_before_expand);
      is_valid_ = true;
      return false;
    }
    return true;
  }

  rt_mx rt_point_group::accumulate() const
  {
    rt_mx result = matrices_[0];
    for(std::size_t i=1;i<matrices_.size();i++) result += matrices_[i];
    result.pseudo_divide(matrices_.size());
    return result.cancel();
  }

  sgtbx::space_group rt_point_group::space_group() const
  {
    sgtbx::space_group sg(false, matrices_[0].t().den());
    sg.expand_smx(matrices_.const_ref());
    CCTBX_ASSERT(sg.n_ltr() == 1);
    return sg;
  }

  matrix_group::code rt_point_group::type() const
  {
    return space_group().point_group_type();
  }

  namespace {

    struct close_mate
    {
      close_mate(rt_mx const& s, double cart_delta_sq)
      : s_(s), cart_delta_sq_(cart_delta_sq)
      {}

      bool
      operator<(close_mate const& other) const
      {
        if (cart_delta_sq_ < other.cart_delta_sq_) return true;
        return false;
      }

      rt_mx s_;
      double cart_delta_sq_;
    };

  } // namespace <anonymous>

  void
  site_symmetry::
  build_special_op()
  {
    int t_den = space_group_.t_den();
    shortest_distance_sq_ = unit_cell_.longest_vector_sq();
    std::vector<close_mate> close_mates;
    close_mates.push_back(close_mate(space_group_.smx(0), 0.));
    for(std::size_t i_smx=0;i_smx<space_group_.order_p();i_smx++) {
      rt_mx s = space_group_(i_smx);
      rot_mx cum_r = s.r().accumulate();
      fractional<> mate0 = s * exact_site_;
      for (std::size_t i_ltr=0;i_ltr<space_group_.n_ltr();i_ltr++) {
        fractional<> mate = mate0 + space_group_.ltr(i_ltr).as_double();
        fractional<> delta0 = mate - exact_site_;
        delta0 = delta0.mod_short();
        tr_vec u_shifts(
          fractional<>((exact_site_ + delta0) - mate).unit_shifts(), 1);
        u_shifts = u_shifts.scale(t_den);
        tr_vec st = s.t() + space_group_.ltr(i_ltr) + u_shifts;
        bool special = false;
        double new_shortest_distance_sq = shortest_distance_sq_;
        sg_vec3& u_num = u_shifts.num();
        for (u_num[0] = -t_den; u_num[0] <= t_den; u_num[0] += t_den)
        for (u_num[1] = -t_den; u_num[1] <= t_den; u_num[1] += t_den)
        for (u_num[2] = -t_den; u_num[2] <= t_den; u_num[2] += t_den) {
          fractional<> delta = delta0 + u_shifts.as_double();
          double cart_delta_sq = unit_cell_.length_sq(delta);
          scitbx::math::update_max(new_shortest_distance_sq, cart_delta_sq);
          if (   min_distance_sym_equiv_sq_ != 0
              && cart_delta_sq <= min_distance_sym_equiv_sq_) {
            tr_vec stu = st + u_shifts;
            tr_vec intrinsic_part = cum_r * stu;
            if (intrinsic_part.num().is_zero()) {
              close_mates.push_back(
                close_mate(rt_mx(s.r(), stu), cart_delta_sq));
              special = true;
            }
          }
        }
        if (!special) shortest_distance_sq_ = new_shortest_distance_sq;
      }
    }
    if (close_mates.size() > 1) {
      std::sort(close_mates.begin() + 1, close_mates.end());
    }
    point_group_.reset(close_mates[0].s_);
    for(std::size_t i=1;i<close_mates.size();i++) {
      if (!point_group_.try_expand(close_mates[i].s_)) {
        scitbx::math::update_min(
          shortest_distance_sq_, close_mates[i].cart_delta_sq_);
      }
    }
    CCTBX_ASSERT(space_group_.order_z() % point_group_.matrices().size() == 0);
    multiplicity_ = space_group_.order_z() / point_group_.matrices().size();
    special_op_ = point_group_.accumulate();
    exact_site_ = special_op_ * original_site_;
  }

  site_symmetry::
  site_symmetry(uctbx::unit_cell const& unit_cell,
                sgtbx::space_group const& space_group,
                fractional<> const& original_site,
                double min_distance_sym_equiv,
                bool assert_min_distance_sym_equiv)
  :
    site_symmetry_ops(0, 1, 1),
    unit_cell_(unit_cell),
    space_group_(space_group),
    original_site_(original_site),
    min_distance_sym_equiv_sq_(scitbx::fn::pow2(min_distance_sym_equiv)),
    shortest_distance_sq_(-1.),
    exact_site_(original_site)
  {
    for (;;) {
      rt_mx last_special_op = special_op_;
      build_special_op();
      if (special_op_ == last_special_op) break;
    }
    if (assert_min_distance_sym_equiv && !check_min_distance_sym_equiv()) {
      throw error("site_symmetry: min_distance_sym_equiv too large.");
    }
    matrices_ = point_group_.matrices();
  }

  af::shared<rt_mx>
  site_symmetry::unique_ops()
  {
    af::shared<rt_mx> result = space_group_.unique(special_op_);
    CCTBX_ASSERT(result.size() == multiplicity_);
    return result;
  }

  void
  site_symmetry_table::
  process(
    std::size_t insert_at_index,
    site_symmetry_ops const& site_symmetry_ops_)
  {
    CCTBX_ASSERT(indices_const_ref_.end() == indices_.end());
    CCTBX_ASSERT(table_const_ref_.end() == table_.end());
    CCTBX_ASSERT(insert_at_index <= indices_const_ref_.size());
    if (table_const_ref_.size() == 0) {
      table_.push_back(site_symmetry_ops_.make_point_group_1());
      table_const_ref_ = table_.const_ref();
    }
    std::size_t i_tab = 0;
    bool is_pg_1 = site_symmetry_ops_.is_point_group_1();
    if (!is_pg_1) {
      for(i_tab=1;i_tab<table_const_ref_.size();i_tab++) {
        if (table_const_ref_[i_tab] == site_symmetry_ops_) {
          break;
        }
      }
      if (i_tab == table_const_ref_.size()) {
        table_.push_back(site_symmetry_ops_);
        table_const_ref_ = table_.const_ref();
      }
    }
    if (insert_at_index == indices_const_ref_.size()) {
      indices_.push_back(i_tab);
      if (!is_pg_1) special_position_indices_.push_back(insert_at_index);
    }
    else {
      indices_.insert(indices_.begin()+insert_at_index, i_tab);
      std::size_t n = special_position_indices_.size();
      if (!is_pg_1) special_position_indices_.resize(n+1);
      std::size_t* si = special_position_indices_.begin();
      std::size_t i = n;
      while (i != 0) {
        i--;
        if (si[i] < insert_at_index) break;
        if (is_pg_1) si[i]++;
        else         si[i+1] = si[i] + 1;
      }
      if (!is_pg_1) si[i] = insert_at_index;
    }
    indices_const_ref_ = indices_.const_ref();
  }

}} // namespace cctbx::sgtbx
