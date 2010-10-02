#include <cctbx/sgtbx/wyckoff.h>
#include <cctbx/sgtbx/reference_settings.h>
#include <cctbx/sgtbx/smith_normal_form.h>

namespace cctbx { namespace sgtbx { namespace wyckoff {

  namespace {
    static const char letter_table[] = "abcdefghijklmnopqrstuvwxyz@";
  }

  position::position(wyckoff::table const* table,
                     int multiplicity,
                     char letter,
                     rt_mx const& special_op)
  : table_(table),
    multiplicity_(multiplicity),
    letter_(letter),
    special_op_(special_op)
  {}

  matrix_group::code
  position::point_group_type() const
  {
    return rt_point_group(table_->space_group_type().group(), special_op_)
             .type();
  }

  af::shared<rt_mx>
  position::unique_ops(const sgtbx::space_group& space_group)
  {
    af::shared<rt_mx> result = space_group.unique(special_op_);
    CCTBX_ASSERT(result.size() == multiplicity_);
    return result;
  }

  table::table(sgtbx::space_group_type const& space_group_type)
  :
    space_group_type_(space_group_type)
  {
    using reference_settings::wyckoff::general_position_multiplicities;
    using reference_settings::wyckoff::raw_table;
    using reference_settings::wyckoff::raw_tables;
    int sg_number = space_group_type_.number();
    CCTBX_ASSERT(1 <= sg_number && sg_number <= 230);
    rot_mx const& cb_r = space_group_type_.cb_op().c().r();
    rat factor_mult(cb_r.num().determinant(), scitbx::fn::pow3(cb_r.den()));
    rat mult = general_position_multiplicities(sg_number) * factor_mult;
    CCTBX_ASSERT(mult.denominator() == 1);
    raw_table const& raw_tab = raw_tables(sg_number);
    CCTBX_ASSERT(raw_tab.n < 27);
    char letter = letter_table[raw_tab.n];
    positions_.push_back(wyckoff::position(
      this, mult.numerator(), letter, rt_mx(1, 1)));
    change_of_basis_op cb_op_inv = space_group_type_.cb_op().inverse();
    rt_mx reference_special_op;
    for(int i=0;i<raw_tab.n;i++) {
      try {
        reference_special_op = rt_mx(raw_tab.op[i].xyz, "", 6, 24);
      }
      catch (error const&) {
        throw CCTBX_INTERNAL_ERROR();
      }
      mult = raw_tab.op[i].m * factor_mult;
      CCTBX_ASSERT(mult.denominator() == 1);
      letter = letter_table[raw_tab.n - 1 - i];
      rt_mx special_op = cb_op_inv.apply(reference_special_op);
      positions_.push_back(wyckoff::position(
        this, mult.numerator(), letter, special_op));
    }
  }

  std::size_t
  table::lookup_index(char letter) const
  {
    for (std::size_t i=0; letter_table[i]; i++) {
      if (letter_table[i] == letter) {
        if (i < size()) return size() - 1 - i;
        throw error("Wyckoff letter out of range.");
      }
    }
    throw error("Not a Wyckoff letter.");
  }

  wyckoff::mapping
  table::mapping(
    uctbx::unit_cell const& unit_cell,
    fractional<> const& x,
    double special_position_radius) const
  {
    space_group const& sg = space_group_type_.group();
    fractional<> norm_x = x.mod_short();
    tr_vec norm_u(fractional<>(norm_x - x).unit_shifts(), 1);
    for (std::size_t i_pos=size()-1;i_pos>0;i_pos--) {
      wyckoff::mapping result;
      double shortest_distance_sq = unit_cell.longest_vector_sq();
      for(std::size_t i_op = 0;i_op<sg.order_z();i_op++) {
        rt_mx sym_op = sg(i_op).mod_short();
        fractional<> sx = sym_op * norm_x;
        tr_vec u(1);
        sg_vec3& u_n = u.num();
        for (u_n[0] = -1; u_n[0] <= 1; u_n[0]++)
        for (u_n[1] = -1; u_n[1] <= 1; u_n[1]++)
        for (u_n[2] = -1; u_n[2] <= 1; u_n[2]++) {
          fractional<> usx = sx + u.as_double();
          fractional<> exact_usx = positions_[i_pos].special_op() * usx;
          double dist_sq = unit_cell.distance_sq(exact_usx, usx);
          if (shortest_distance_sq > dist_sq) {
            shortest_distance_sq = dist_sq;
            result = wyckoff::mapping(
              unit_cell,
              x,
              positions_[i_pos],
              rt_mx(sym_op.r(),
                    sym_op.t().plus(u).plus(sym_op.r().multiply(norm_u))));
          }
        }
      }
      if (shortest_distance_sq <= scitbx::fn::pow2(special_position_radius)) {
        return result;
      }
    }
    return wyckoff::mapping(unit_cell, x, positions_[0], rt_mx(1, 1));
  }

  namespace {

    tr_vec solve_in_z(rot_mx m, tr_vec b)
    {
      int f = boost::lcm(m.den(), b.den());
      m = m.scale(f / m.den());
      b = b.scale(f / b.den());
      rot_mx p(1);
      rot_mx q(1);
      af::ref<int, af::mat_grid> m_ref(m.num().begin(), 3, 3);
      af::ref<int, af::mat_grid> p_ref(p.num().begin(), 3, 3);
      af::ref<int, af::mat_grid> q_ref(q.num().begin(), 3, 3);
      smith_normal_form(m_ref, p_ref, q_ref);
      CCTBX_ASSERT(m_ref.is_square());
      std::size_t nd = m_ref.n_rows();
      CCTBX_ASSERT(nd <= 3);
      tr_vec pb = p * b;
      for(std::size_t i=nd;i<3;i++) {
        if (pb[i] != 0) return tr_vec(0);
      }
      tr_vec xp(1);
      for(std::size_t i=0;i<nd;i++) {
        int d = m_ref(i,i);
        if (pb[i] % d != 0) return tr_vec(0);
        xp[i] = pb[i] / d;
      }
      return q * xp;
    }

  } // namespace <anonymous>

  wyckoff::mapping
  table::mapping(sgtbx::site_symmetry const& site_symmetry) const
  {
    CCTBX_ASSERT(positions_.size() > 0);
    CCTBX_ASSERT(site_symmetry.space_group() == space_group_type_.group());
    const wyckoff::position* pos = positions_.begin();
    if (pos->multiplicity() == site_symmetry.multiplicity()) {
      return wyckoff::mapping(
        site_symmetry.unit_cell(),
        site_symmetry.original_site(),
        *pos,
        rt_mx(1, 1));
    }
    rot_mx const& r_s = site_symmetry.special_op().r();
    tr_vec const& t_s = site_symmetry.special_op().t();
    space_group const& sg = site_symmetry.space_group();
    for(;pos!=positions_.end();pos++) {
      if (pos->multiplicity() == site_symmetry.multiplicity()) {
        rot_mx const& r_w = pos->special_op().r();
        tr_vec const& t_w = pos->special_op().t();
        for(std::size_t i_op=0;i_op<sg.order_p();i_op++) {
          rt_mx sym_op = sg(i_op);
          rot_mx const& r = sym_op.r();
          if (r_w.multiply(r) == r.multiply(r_s)) {
            for(std::size_t i_ltr=0;i_ltr<sg.n_ltr();i_ltr++) {
              tr_vec t = sym_op.t() + sg.ltr(i_ltr);
              rot_mx r_w_m_i = r_w.minus_unit_mx().cancel();
              tr_vec b = r.multiply(t_s).plus(t)
                          .minus(r_w.multiply(t)).minus(t_w);
              tr_vec u = solve_in_z(r_w_m_i, b);
              if (u.is_valid()) {
                return wyckoff::mapping(
                  site_symmetry.unit_cell(),
                  site_symmetry.original_site(),
                  *pos,
                  rt_mx(r.cancel(),
                  t.plus(u)));
              }
            }
          }
        }
      }
    }
    throw error("Cannot determine mapping to Wyckoff position.");
  }

}}} // namespace cctbx::sgtbx::wyckoff
