#include <cctbx/sgtbx/space_group_type.h>
#include <cctbx/sgtbx/select_generators.h>
#include <cctbx/sgtbx/reference_settings.h>
#include <cctbx/sgtbx/rot_mx_info.h>
#include <cctbx/sgtbx/smith_normal_form.h>
#include <cctbx/sgtbx/row_echelon_solve.h>
#include <cctbx/sgtbx/symbols.h>
#include <cctbx/math/loop_n_from_m.h>
#include <cctbx/sgtbx/utils.h>

namespace cctbx { namespace sgtbx {

  namespace construct_cb_op_r {
  namespace {

    matrix_group::code
    get_matrix_group_type(space_group const& std_sg,
                          matrix_group::code const& point_group)
    {
      using namespace matrix_group;

      bool two_fold = false;
      bool mirror = false;

      if      (   point_group == code_4bm2
               || point_group == code_6bm2) {
        two_fold = true;
      }
      else if (std_sg.n_ltr() == 1) {
        if      (point_group == code_32)  two_fold = true;
        else if (point_group == code_3m)  mirror  = true;
        else if (point_group == code_3bm) two_fold = mirror = true;
      }

      if (!(two_fold || mirror)) return point_group;

      for(std::size_t i_smx=1;i_smx<std_sg.n_smx();i_smx++) {
        int r_type = std_sg.smx(i_smx).r().type();
        CCTBX_ASSERT(r_type != 0);
        if (   (r_type ==  2 && two_fold)
            || (r_type == -2 && mirror)) {
          rot_mx_info ri(std_sg.smx(i_smx).r());
          const sg_vec3 ev_100(1, 0, 0);
          if (ri.ev() == ev_100) {
            if (point_group == code_4bm2) return code_4b2m;
            if (point_group == code_32)   return code_321;
            if (point_group == code_3m)   return code_3m1;
            if (point_group == code_3bm)  return code_3bm1;
            if (point_group == code_6bm2) return code_6b2m;
          }
        }
      }

      if (point_group == code_4bm2) return code_4bm2;
      if (point_group == code_32)   return code_312;
      if (point_group == code_3m)   return code_31m;
      if (point_group == code_3bm)  return code_3b1m;
      if (point_group == code_6bm2) return code_6bm2;

      throw CCTBX_INTERNAL_ERROR();
    }

    struct basis_builder
    {
      basis_builder(matrix_group::code const& laue_group);

      void pick_eigen_vectors(space_group const& raw_sg);

      rot_mx uni_axial_basis();

      rot_mx two_axial_basis();

      const crystal_system::code xtal_system;
      std::size_t n_wtd;
      int ord[3];
      rot_mx r[3];
      rot_mx_info ri[3];
    };

    basis_builder::basis_builder(matrix_group::code const& laue_group)
    : xtal_system(laue_group.crystal_system())
    {
      using namespace matrix_group;

      if      (laue_group == code_1b) {
        n_wtd = 0;
      }
      else if (laue_group == code_2_m) {
        n_wtd = 1; ord[0] = 2;
      }
      else if (laue_group == code_mmm) {
        n_wtd = 3; ord[0] = 2; ord[1] = 2; ord[2] = 2;
      }
      else if (laue_group == code_4_m) {
        n_wtd = 1; ord[0] = 4;
      }
      else if (laue_group == code_4_mmm) {
        n_wtd = 2; ord[0] = 4; ord[1] = 2;
      }
      else if (laue_group == code_3b) {
        n_wtd = 1; ord[0] = 3;
      }
      else if (laue_group == code_3bm) {
        n_wtd = 2; ord[0] = 3; ord[1] = 2;
      }
      else if (laue_group == code_6_m) {
        n_wtd = 1; ord[0] = 3;
      }
      else if (laue_group == code_6_mmm) {
        n_wtd = 2; ord[0] = 3; ord[1] = 2;
      }
      else if (laue_group == code_m3b) {
        n_wtd = 2; ord[0] = 3; ord[1] = 2;
      }
      else if (laue_group == code_m3bm) {
        n_wtd = 2; ord[0] = 3; ord[1] = 4;
      }
      else {
        throw CCTBX_INTERNAL_ERROR();
      }
    }

    void basis_builder::pick_eigen_vectors(space_group const& raw_sg)
    {
      using scitbx::fn::absolute;
      for(std::size_t i_wtd=0;i_wtd<n_wtd;i_wtd++) r[i_wtd] = rot_mx(0);
      std::size_t n_found = 0;
      for (;;) {
        bool restart = false;
        for (std::size_t i_smx=0;i_smx<raw_sg.n_smx();i_smx++) {
          int r_type = raw_sg.smx(i_smx).r().type();
          CCTBX_ASSERT(r_type != 0);
          for (std::size_t i_wtd=0;i_wtd<n_wtd;i_wtd++) {
            if (   !r[i_wtd].is_valid()
                && (r_type == ord[i_wtd] || -r_type == ord[i_wtd])) {
              ri[i_wtd] = raw_sg.smx(i_smx).r().info();
              if (ri[i_wtd].sense() >= 0) {
                bool use_this_smx = true;
                for (std::size_t j_wtd=0;j_wtd<n_wtd;j_wtd++) {
                  if (   r[j_wtd].is_valid()
                      && ri[i_wtd].ev() == ri[j_wtd].ev()) {
                    if (absolute(ord[j_wtd]) >= absolute(r_type)) {
                      use_this_smx = false;
                    }
                    else {
                      r[j_wtd] = rot_mx(0);
                      n_found--;
                      restart = true;
                    }
                    break;
                  }
                }
                if (use_this_smx) {
                  ord[i_wtd] = r_type;
                  r[i_wtd] = raw_sg.smx(i_smx).r();
                      n_found++;
                  if (n_found == n_wtd)
                    return;
                }
              }
              break;
            }
          }
        }
        if (!restart) break;
      }
      throw CCTBX_INTERNAL_ERROR();
    }

    rot_mx basis_builder::uni_axial_basis()
    {
      using scitbx::fn::absolute;
      if (ord[0] < 0) {
        ord[0] *= -1;
        r[0] *= -1;
        ri[0] = r[0].info();
      }
      rot_mx cum_mx = r[0].accumulate(ord[0]);
      af::ref<int, af::mat_grid> re_mx(cum_mx.num().begin(), 3, 3);
      CCTBX_ASSERT(scitbx::matrix::row_echelon::form(re_mx) == 1);
      af::tiny<sg_vec3, 4> sol = row_echelon::solve::homog_rank_1(re_mx);
      std::size_t n_ix = 1; if (ord[0] == 2) n_ix++;
      rot_mx trial_basis;
      trial_basis.num().set_column(2, ri[0].ev());
      int min_det = 0;
      rot_mx result;
      for (math::loop_n_from_m<2> ix(4, n_ix); !ix.over(); ix.incr()) {
        for(std::size_t i=0;i<ix.n();i++) {
          trial_basis.num().set_column(i, sol[ix[i]]);
        }
        if (ix.n() == 1) trial_basis.num().set_column(1, r[0] * sol[ix[0]]);
        int det = trial_basis.num().determinant();
        if (det != 0
            && (min_det == 0 || absolute(min_det) > absolute(det))) {
          min_det = det;
          result = trial_basis;
        }
      }
      CCTBX_ASSERT(min_det != 0);
      if (min_det < 0) result.num().swap_columns(0, 1);
      return result;
    }

    rot_mx basis_builder::two_axial_basis()
    {
      rot_mx result;
      result.num().set_column(0, ri[1].ev());
      if (xtal_system == crystal_system::cubic) {
        for(std::size_t i=0;i<2;i++) {
          result.num().set_column(i + 1, r[0] * result.num().get_column(i));
        }
        if (result.num().determinant() < 0) result.num().swap_columns(0, 1);
      }
      else {
        result.num().set_column(2, ri[0].ev());
        if (n_wtd == 3) {
          result.num().set_column(1, ri[2].ev());
          if (result.num().determinant() < 0) result.num().swap_columns(0, 1);
        }
        else {
          if (ord[0] > 0) {
            result.num().set_column(1,  r[0] * result.num().get_column(0));
          }
          else {
            result.num().set_column(1, -r[0] * result.num().get_column(0));
          }
          CCTBX_ASSERT(result.num().determinant() > 0);
        }
      }
      return result;
    }

    rot_mx std_basis(space_group const& raw_sg, matrix_group::code const& mgc)
    {
      basis_builder bb(mgc.laue_group_type());
      if (bb.n_wtd == 0) return rot_mx();
      bb.pick_eigen_vectors(raw_sg);
      if (bb.n_wtd == 1) return bb.uni_axial_basis();
      return bb.two_axial_basis();
    }

    bool is_m3b_100_glide(space_group const& work_sg)
    {
      using scitbx::fn::absolute;
      CCTBX_ASSERT(work_sg.is_centric());
      for (std::size_t i_smx=1;i_smx<work_sg.n_smx();i_smx++) {
        rot_mx const& r = work_sg.smx(i_smx).r();
        int r_type = r.type();
        CCTBX_ASSERT(r_type != 0);
        if (absolute(r_type) == 2) {
          rot_mx_info ri(r);
          const sg_vec3 ev_100(1, 0, 0);
          if (ri.ev() == ev_100) {
            std::size_t i_inv = 0; if (r_type == 2) i_inv = 1;
            rt_mx smx = work_sg(0, i_inv, i_smx);
            tr_vec wi = smx.t_intrinsic_part();
            if (wi[2] % wi.den() != 0) return true;
            return false;
          }
        }
      }
      throw CCTBX_INTERNAL_ERROR();
    }

    rot_mx get_adj_rmx(matrix_group::code const& laue_group, char Z)
    {
      using namespace matrix_group;

      if      (    laue_group == code_2_m) {
        // monoclinic unique c -> unique b
        return rot_mx(  0,  1,  0,
                        0,  0,  1,
                        1,  0,  0,  1);
      }
      if ((laue_group == code_4_m || laue_group == code_4_mmm) && Z == 'C') {
        // C -> P
        return rot_mx(  1,  1,  0,
                        1, -1,  0,
                        0,  0, -1,  1);
      }
      if ((laue_group == code_4_m || laue_group == code_4_mmm) && Z == 'F') {
        // F -> I
        return rot_mx(  1,  1,  0,
                       -1,  1,  0,
                        0,  0,  1,  1);
      }
      if ((laue_group == code_3b || laue_group == code_3bm) && Z == 'Q') {
        // reverse -> obverse
        return rot_mx( -1,  0,  0,
                        0, -1,  0,
                        0,  0,  1,  1);
      }
      if ((laue_group == code_3bm || laue_group == code_6_mmm) && Z == 'H') {
        // H -> P
        return rot_mx(  1,  1,  0,
                       -1,  2,  0,
                        0,  0,  1,  1);
      }
      return rot_mx(0);
    }

  } // namespace <anonymous>
  } // namespace construct_cb_op_r

  namespace construct_cb_op_t {
  namespace {

    bool solve_inhomog_mod_z(int *m, std::size_t nr, std::size_t nc,
                             int *b, int den,
                             int *x)
    {
      const std::size_t maxr = 9;
      const std::size_t maxc = 6;

      CCTBX_ASSERT(nr <= maxr);
      CCTBX_ASSERT(nc <= maxc);

      int p[maxr * maxr];
      int q[maxc * maxc];
      af::ref<int, af::mat_grid> m_ref(m, nr, nc);
      af::ref<int, af::mat_grid> p_ref(p, nr, nr);
      af::ref<int, af::mat_grid> q_ref(q, nc, nc);
      smith_normal_form(m_ref, p_ref, q_ref);
      CCTBX_ASSERT(m_ref.is_square());
      std::size_t nd = m_ref.n_rows();
      CCTBX_ASSERT(nd <= nc);

      int pb[maxr];
      scitbx::matrix::multiply(p, b, nr, nr, 1, pb);
      for(std::size_t i=nd;i<nr;i++) {
        if (pb[i] % den != 0) return false;
      }

      if (x) {
        int xp[maxc];
        for(std::size_t i=0;i<nc;i++) {
          xp[i] = 0;
          int d = m[i * nd + i];
          if (d) {
            CCTBX_ASSERT(pb[i] % d == 0);
            xp[i] = pb[i] / d;
          }
        }
        scitbx::matrix::multiply(q, xp, nc, nc, 1, x);
      }

      return true;
    }

    tr_vec find_origin_shift(
      select_generators::standard const& tab_generators,
      select_generators::standard const& tst_generators,
      int t_den)
    {
      /*    (I|K)(R|T)(I|-K)=(R|S)
         => S=-RK+T+K=-(R-I)K+T
         => S=-(R-I)K+T
         => (R-I)K=T-S
         => (R-I)^-1(T-S)=K
       */

      rot_mx rmi[3];
      tr_vec delta_t[3];
      for(std::size_t i=0;i<tab_generators.n_gen;i++) {
        CCTBX_ASSERT(   tst_generators.p_gen[i].r()
                     == tab_generators.p_gen[i].r());
        rmi[i] = tst_generators.p_gen[i].r() - rot_mx();
        delta_t[i] = (  tst_generators.p_gen[i].t()
                      - tab_generators.p_gen[i].t()).new_denominator(t_den);
      }
      CCTBX_ASSERT(   tst_generators.p_inv_t.is_valid()
                   == tab_generators.p_inv_t.is_valid());
      std::size_t n_gen = tab_generators.n_gen;
      if (tst_generators.p_inv_t.is_valid()) {
        rmi[n_gen] = rot_mx(1, -2);
        delta_t[n_gen] = (  tst_generators.p_inv_t
                          - tab_generators.p_inv_t).new_denominator(t_den);
        n_gen++;
      }
      std::size_t nr_snf = n_gen * 3;
      int snf[9 * 3], v[3 * 3];
      for (std::size_t iGen = 0; iGen < n_gen; iGen++) {
        for(std::size_t i=0;i<9;i++) snf[iGen * 9 + i] = rmi[iGen][i];
        for(std::size_t i=0;i<3;i++) v[iGen * 3 + i] = delta_t[iGen][i];
      }
      sg_vec3 x;
      if (solve_inhomog_mod_z(snf, nr_snf, 3, v, t_den, x.begin())) {
        tr_vec cb_t = tst_generators.z2p_op.c_inv().r() * tr_vec(x, t_den);
        return cb_t.new_denominator(t_den);
      }
      return tr_vec(0);
    }

    change_of_basis_op
    match_generators(
      int r_den,
      int t_den,
      space_group const& work_sg,
      matrix_group::code const& point_group,
      char tab_z,
      select_generators::standard const& tab_generators)
    {
      if (tab_generators.n_gen == 0 && !tab_generators.z_inv_t.is_valid()) {
        return change_of_basis_op(r_den, t_den); // space group P 1
      }
      rot_mx r_2fold(0);
      rot_mx r_3fold(0);
      if      (point_group.crystal_system() == crystal_system::monoclinic) {
        r_2fold = rot_mx( // 2 [101]
          0,  0,  1,
          0, -1,  0,
          1,  0,  0,  1);
        r_3fold = rot_mx( // 3 [010]
         -1,  0,  1,
          0,  1,  0,
         -1,  0,  0,  1);
      }
      else if (point_group.crystal_system() == crystal_system::orthorhombic) {
        r_2fold = rot_mx( // 2 [110]
          0,  1,  0,
          1,  0,  0,
          0,  0, -1,  1);
        r_3fold = rot_mx( // 3 [111]
          0,  0,  1,
          1,  0,  0,
          0,  1,  0,  1);
      }
      if (r_2fold.is_valid())
      {
        change_of_basis_op
          cb_op_2fold(rt_mx(r_2fold.new_denominator(r_den), t_den));
        change_of_basis_op
          cb_op_3fold(rt_mx(r_3fold.new_denominator(r_den), t_den));
        change_of_basis_op cb_op(r_den, t_den);
        for(std::size_t i2fold=0;i2fold<2;i2fold++) {
          if (i2fold) cb_op = cb_op_2fold;
          for (std::size_t i3fold=0;i3fold<3;i3fold++) {
            if (i3fold) cb_op.update(cb_op_3fold);
            space_group tst_sg = work_sg.change_basis(cb_op);
            char tst_z = tst_sg.ltr().conventional_centring_type_symbol();
            CCTBX_ASSERT(tst_z != '\0' && tst_z != 'Q');
            if (tst_z != tab_z)
              continue;
            select_generators::standard tst_generators(
              tst_sg, r_den, t_den, point_group);
            CCTBX_ASSERT(tst_generators.n_gen == tab_generators.n_gen);
            if (    tab_generators.n_gen != 2
                ||  (      tab_generators.z_gen[0].r()[8]
                        == tst_generators.z_gen[0].r()[8]
                     &&    tab_generators.z_gen[1].r()[0]
                        == tst_generators.z_gen[1].r()[0])) {
              tst_generators.set_primitive();
              tr_vec cb_t = find_origin_shift(tab_generators,
                                              tst_generators, t_den);
              if (cb_t.is_valid()) {
                cb_op.update(cb_t);
                return cb_op;
              }
            }
          }
        }
      }
      else {
        select_generators::standard tst_generators(
          work_sg, r_den, t_den, point_group);
        CCTBX_ASSERT(tst_generators.n_gen == tab_generators.n_gen);
        tst_generators.set_primitive();
        tr_vec cb_t = find_origin_shift(tab_generators, tst_generators, t_den);
        if (cb_t.is_valid()) {
          return change_of_basis_op(rt_mx(cb_t, r_den));
        }
      }
      return change_of_basis_op(0, 0);
    }

    void
    tidy_cb_op_t(
      select_generators::standard const& target_generators,
      space_group const& given_sg,
      matrix_group::code const& point_group,
      change_of_basis_op& trial_cb_op)
    {
      // set translation parts to zero
      trial_cb_op = change_of_basis_op(
        rt_mx(trial_cb_op.c().r(), trial_cb_op.c().t().den()),
        rt_mx(trial_cb_op.c_inv().r(), trial_cb_op.c_inv().t().den()));

      // done if space group is P 1
      if (given_sg.n_smx() == 1 && !given_sg.is_centric()) return;

      space_group transformed_sg = given_sg.change_basis(trial_cb_op);
      select_generators::standard transformed_generators(
        transformed_sg,
        trial_cb_op.c().r().den(),
        trial_cb_op.c().t().den(),
        point_group);
      transformed_generators.set_primitive();
      tr_vec cb_t = find_origin_shift(
        target_generators, transformed_generators, trial_cb_op.c().t().den());
      CCTBX_ASSERT(cb_t.is_valid());
      trial_cb_op.update(cb_t.mod_short());
    }

    class cmp_change_of_basis_mx
    {
      public:
        cmp_change_of_basis_mx() : cmp_r_(9), cmp_t_(3) {}

        bool operator()(rt_mx const& a, rt_mx const& b) const
        {
          using std::size_t;
          using scitbx::fn::absolute;

          rot_mx const& ar = a.r(); tr_vec const& at = a.t();
          rot_mx const& br = b.r(); tr_vec const& bt = b.t();

          bool ba = ar.is_unit_mx();
          bool bb = br.is_unit_mx();
          if ( ba && !bb) return true;
          if (!ba &&  bb) return false;

          ba = ar.num().is_diagonal();
          bb = br.num().is_diagonal();
          if ( ba && !bb) return true;
          if (!ba &&  bb) return false;

          ba = at.num().is_zero();
          bb = bt.num().is_zero();
          if ( ba && !bb) return true;
          if (!ba &&  bb) return false;

          size_t na = 0; for(size_t i=0;i<9;i++) if (ar[i] == 0) na++;
          size_t nb = 0; for(size_t i=0;i<9;i++) if (br[i] == 0) nb++;
          if (na > nb) return true;
          if (na < nb) return false;

          na = 0; for(size_t i=0;i<9;i++) if (absolute(ar[i])==ar.den()) na++;
          nb = 0; for(size_t i=0;i<9;i++) if (absolute(br[i])==br.den()) nb++;
          if (na > nb) return true;
          if (na < nb) return false;

          na = 0; for(size_t i=0;i<9;i++) if (ar[i] > 0) na++;
          nb = 0; for(size_t i=0;i<9;i++) if (br[i] > 0) nb++;
          if (na > nb) return true;
          if (na < nb) return false;

          if (cmp_t_(at.num().begin(), bt.num().begin())) return true;
          if (cmp_t_(bt.num().begin(), at.num().begin())) return false;

          return cmp_r_(br.num().begin(), ar.num().begin());
        }
      private:
        const utils::cmp_i_vec cmp_r_;
        const utils::cmp_i_vec cmp_t_;
    };

    /* Uses affine normalizer operations to find the "best"
       change-of-basis matrix.
       For a given space group representation the resulting
       change-of-basis matrix should always be identical,
       independent of the order of symmetry operations
       and indepenent of the raw_cb_op.
     */
    change_of_basis_op
    find_best_cb_op(
      space_group const& given_sg,
      matrix_group::code const& point_group,
      int sg_number,
      change_of_basis_op const& reference_cb_op,
      space_group const& target_sg,
      select_generators::standard const& target_generators,
      change_of_basis_op const& raw_cb_op)
    {
      af::shared<rt_mx>
        addl_g = reference_settings::normalizer::get_addl_generators(
          sg_number, true, true, true);
      change_of_basis_op cb_op = reference_cb_op.inverse();
      space_group norm_sg = target_sg;
      std::size_t old_order_p = norm_sg.order_p();
      for(std::size_t i=0;i<addl_g.size();i++) {
        rt_mx s = cb_op(addl_g[i].new_denominators(norm_sg.smx()[0]));
        norm_sg.expand_smx(s);
        std::size_t new_order_p = norm_sg.order_p();
        CCTBX_ASSERT(   old_order_p * s.r().order()
                     == new_order_p);
        old_order_p = new_order_p;
      }
      change_of_basis_op best_cb_op = raw_cb_op;
      bool cmp_inv = false;
      if (  best_cb_op.c().r().num().determinant()
          < best_cb_op.c_inv().r().num().determinant()) {
        cmp_inv = true;
      }
      cmp_change_of_basis_mx cmp_cb_mx;
      for(std::size_t i_inv=0;i_inv<norm_sg.f_inv();i_inv++)
      for(std::size_t i_smx=0;i_smx<norm_sg.n_smx();i_smx++)
      {
        rt_mx s = norm_sg(0, i_inv, i_smx);
        if (s.r().num().determinant() < 0) continue;
        change_of_basis_op norm_cb_op(s.new_denominators(raw_cb_op.c()));
        change_of_basis_op trial_cb_op = norm_cb_op * raw_cb_op;
        if (sg_number < 3 || sg_number > 15) {
          tidy_cb_op_t(target_generators, given_sg, point_group, trial_cb_op);
          if (!cmp_cb_mx(best_cb_op.select(cmp_inv),
                         trial_cb_op.select(cmp_inv))) {
            best_cb_op = trial_cb_op;
          }
        }
        else {
          int r_den = raw_cb_op.c().r().den();
          int t_den = raw_cb_op.c().t().den();
          int r00, r22;
          reference_settings::normalizer::get_monoclinic_affine_trial_ranges(
            trial_cb_op.c().r(), r00, r22);
          sg_mat3 r, r_inv;
          r.fill(0);
          r_inv.fill(0);
#define loop(i, rxx) \
          for (r[i] = -rxx * r_den; r[i] <= rxx * r_den; r[i] += r_den)
          loop(0, r00)
          loop(2, r22)
          loop(6, r00)
          loop(8, r22) {
#undef loop
            if (!reference_settings::normalizer
                  ::check_monoclinic_affine_restrictions(
                     sg_number, rot_mx(r, r_den))) {
              continue;
            }
            r[4] = r[0] * r[8] - r[2] * r[6];
            if (   r[4] != -r_den * r_den
                && r[4] !=  r_den * r_den) continue;
            r[4] /= r_den;
            int f = r[4] / r_den;
            r_inv[0] =  f * r[8];
            r_inv[2] = -f * r[2];
            r_inv[4] =      r[4];
            r_inv[6] = -f * r[6];
            r_inv[8] =  f * r[0];
            change_of_basis_op r_trial_cb_op = change_of_basis_op(
                rt_mx(rot_mx(r, r_den), t_den),
                rt_mx(rot_mx(r_inv, r_den), t_den)) * trial_cb_op;
            tidy_cb_op_t(
              target_generators, given_sg, point_group, r_trial_cb_op);
            if (!cmp_cb_mx(best_cb_op.select(cmp_inv),
                           r_trial_cb_op.select(cmp_inv))) {
              best_cb_op = r_trial_cb_op;
            }
          }
        }
      }
      return best_cb_op;
    }

  } // namespace <anonymous>
  } // namespace construct_cb_op_t

  space_group_type::space_group_type(
    space_group const& group,
    bool tidy_cb_op,
    int r_den, int t_den)
  :
    group_(group),
    number_(0),
    cb_op_(r_den, t_den),
    cb_op_is_tidy_(tidy_cb_op)
  {
    matrix_group::code point_group = group_.point_group_type();
    matrix_group::code laue_group = point_group.laue_group_type();
    crystal_system::code xtal_system = point_group.crystal_system();

    space_group z_point_group_sg = group_.build_derived_group(false, false);
    space_group work_sg = z_point_group_sg;
    char centring_type = '\0';
    std::size_t run_away_counter = 0;

    do // align the point group operations
    {
      CCTBX_ASSERT(run_away_counter++ < 10);

      change_of_basis_op addl_cb_op = work_sg.z2p_op(r_den, t_den);
      cb_op_.update(addl_cb_op);
      work_sg = z_point_group_sg.change_basis(cb_op_);
      CCTBX_ASSERT(work_sg.n_ltr() == 1);

      rot_mx std_basis = construct_cb_op_r::std_basis(work_sg, laue_group);
      std_basis = std_basis.new_denominator(r_den);
      addl_cb_op = change_of_basis_op(rt_mx(std_basis, t_den)).inverse();
      cb_op_.update(addl_cb_op);
      work_sg = z_point_group_sg.change_basis(cb_op_);
      centring_type = work_sg.ltr().conventional_centring_type_symbol();

      rot_mx adj_rmx = construct_cb_op_r::get_adj_rmx(
        laue_group, centring_type);
      if (adj_rmx.is_valid()) {
        adj_rmx = adj_rmx.new_denominator(r_den);
        addl_cb_op = change_of_basis_op(rt_mx(adj_rmx, t_den));
        cb_op_.update(addl_cb_op);
        work_sg = z_point_group_sg.change_basis(cb_op_);
        centring_type = work_sg.ltr().conventional_centring_type_symbol();
      }
    }
    while (centring_type == '\0');
    CCTBX_ASSERT(centring_type != 'Q');

    // transform original symmetry operations
    work_sg = group_.change_basis(cb_op_);

    if (   point_group == matrix_group::code_m3b
        && centring_type == 'P'
        && construct_cb_op_r::is_m3b_100_glide(work_sg)) {
      // rotate by 90 degree (4 [0 0 1])
      rot_mx r_4_001(  0, -1,  0,
                       1,  0,  0,
                       0,  0,  1,  1);
      change_of_basis_op addl_cb_op(
        rt_mx(r_4_001.new_denominator(r_den), t_den));
      cb_op_.update(addl_cb_op);
      work_sg = group_.change_basis(cb_op_);
    }

    matrix_group::code
      mx_group = construct_cb_op_r::get_matrix_group_type(
        work_sg, point_group);
    CCTBX_ASSERT(point_group == mx_group.point_group_type());

    bool match_centring_type_symbol = (
             xtal_system != crystal_system::monoclinic
      && (   xtal_system != crystal_system::orthorhombic
          || (centring_type == 'I' || centring_type == 'F')));

    for (number_=1;number_<=230;number_++)
    {
      const char*
        hall_symbol = reference_settings::hall_symbol_table(number_);

      if (match_centring_type_symbol && centring_type != hall_symbol[1])
        continue;

      if ((centring_type == 'P') != (hall_symbol[1] == 'P'))
        continue;

      if (reference_settings::matrix_group_code_table(number_) != mx_group)
        continue;

      space_group tab_sg;
      try {
        tab_sg = space_group(hall_symbol, true, false, false, work_sg.t_den());
      }
      catch (error const&) {
        throw CCTBX_INTERNAL_ERROR();
      }

      if (tab_sg.n_ltr() != work_sg.n_ltr())
        continue;

      select_generators::standard tab_generators(
        tab_sg, r_den, t_den, point_group);
      tab_generators.set_primitive();
      change_of_basis_op addl_cb_op = construct_cb_op_t::match_generators(
        r_den, t_den, work_sg, point_group, hall_symbol[1], tab_generators);
      if (addl_cb_op.is_valid()) {
        cb_op_.update(addl_cb_op);
        if (tidy_cb_op) {
          cb_op_ = construct_cb_op_t::find_best_cb_op(
            group_, point_group, number_,
            change_of_basis_op(r_den, t_den), tab_sg, tab_generators, cb_op_);
        }
        return;
      }
    }
    throw CCTBX_INTERNAL_ERROR();
  }

  std::string
  space_group_type::hall_symbol(bool tidy_cb_op) const
  {
    if (tidy_cb_op) {
      if (hall_symbol_tidy_true_.size()) return hall_symbol_tidy_true_;
    }
    else {
      if (hall_symbol_tidy_false_.size()) return hall_symbol_tidy_false_;
    }
    std::string hall_symbol(reference_settings::hall_symbol_table(number_));
    parse_string ps(hall_symbol);
    space_group target_sg(false, group_.t_den());
    change_of_basis_op reference_cb_op = cb_op_.identity_op();
    target_sg.parse_hall_symbol_cb_op(ps, reference_cb_op, true);
    change_of_basis_op cb_op = reference_cb_op;
    if (cb_op.is_valid()) {
      cb_op = cb_op.inverse() * cb_op_;
    }
    else {
      cb_op = cb_op_;
    }
    if (tidy_cb_op) {
      matrix_group::code point_group = target_sg.point_group_type();
      select_generators::standard target_generators(
        target_sg,
        cb_op.c().r().den(),
        cb_op.c().t().den(),
        point_group);
      target_generators.set_primitive();
      cb_op = construct_cb_op_t::find_best_cb_op(
        group_, point_group, number_,
        reference_cb_op, target_sg, target_generators, cb_op);
    }

    std::string::size_type par = hall_symbol.find(" (");
    if (par != std::string::npos) {
      hall_symbol.resize(par);
    }

    if (!cb_op.is_identity_op()) {
      hall_symbol += " (" + cb_op.c_inv().mod_short().as_xyz() + ")";
    }

    if (tidy_cb_op) {
      hall_symbol_tidy_true_ = hall_symbol;
    }
    else {
      hall_symbol_tidy_false_ = hall_symbol;
    }

    return hall_symbol;
  }

  namespace {

    std::string
    uhm_reference_symbol(std::size_t space_group_number)
    {
      std::string result(
        reference_settings::hermann_mauguin_symbol_table(space_group_number));
      std::size_t i = result.size();
      if (i > 1) {
        i--;
        std::size_t j = i--;
        if (result[i] == ':' && result[j] == 'h') result[j] = 'H';
      }
      return result;
    }

    std::string
    uhm_change_of_basis_symbol(change_of_basis_op const& cb_op)
    {
      return " (" + cb_op.inverse().mod_short().symbol() + ")";
    }

  }

  std::string
  space_group_type::universal_hermann_mauguin_symbol(bool tidy_cb_op) const
  {
    if (tidy_cb_op) {
      if (uhm_symbol_tidy_true_.size() == 0) {
        uhm_symbol_tidy_true_ = uhm_reference_symbol(number_);
        if (!cb_op_.is_identity_op()) {
          if (cb_op_is_tidy_) {
            uhm_symbol_tidy_true_ += uhm_change_of_basis_symbol(cb_op_);
          }
          else {
            space_group tab_sg(
              reference_settings::hall_symbol_table(number_),
              true, false, false, group_.t_den());
            matrix_group::code point_group = tab_sg.point_group_type();
            select_generators::standard target_generators(
              tab_sg,
              cb_op_.c().r().den(),
              cb_op_.c().t().den(),
              point_group);
            target_generators.set_primitive();
            change_of_basis_op cb_op = construct_cb_op_t::find_best_cb_op(
              group_, point_group, number_,
              cb_op_.identity_op(), tab_sg, target_generators, cb_op_);
            if (!cb_op.is_identity_op()) {
              uhm_symbol_tidy_true_ += uhm_change_of_basis_symbol(cb_op);
            }
          }
        }
      }
      return uhm_symbol_tidy_true_;
    }
    else {
      if (uhm_symbol_tidy_false_.size() == 0) {
        uhm_symbol_tidy_false_ = uhm_reference_symbol(number_);
        if (!cb_op_.is_identity_op()) {
          uhm_symbol_tidy_false_ += uhm_change_of_basis_symbol(cb_op_);
        }
      }
      return uhm_symbol_tidy_false_;
    }
  }

  std::string
  space_group_type::lookup_symbol(
    bool ad_hoc_1992) const
  {
    // simlar table in symbols.cpp, ad_hoc_1992_symbol_as_a1983_symbol()
    static char const* adh_a38_pairs[] = {
      // No. 39
      "A e m 2", "A b m 2",
      "B m e 2", "B m a 2",
      "B 2 e m", "B 2 c m",
      "C 2 m e", "C 2 m b",
      "C m 2 e", "C m 2 a",
      "A e 2 m", "A c 2 m",
      // No. 41
      "A e a 2", "A b a 2",
      "B b e 2", "B b a 2",
      "B 2 e b", "B 2 c b",
      "C 2 c e", "C 2 c b",
      "C c 2 e", "C c 2 a",
      "A e 2 a", "A c 2 a",
      // No. 64
      "C m c e", "C m c a",
      "C c m e", "C c m b",
      "A e m a", "A b m a",
      "A e a m", "A c a m",
      "B b e m", "B b c m",
      "B m e b", "B m a b",
      // No. 67
      "C m m e", "C m m a", // ambiguous: Cmmb
      "A e m m", "A b m m", // ambiguous: Acmm
      "B m e m", "B m c m", // ambiguous: Bmam
      // No. 68
      "C c c e", "C c c a", // ambiguous: Cccb
      "A e a a", "A b a a", // ambiguous: Acaa
      "B b e b", "B b c b"  // ambiguous: Bbab
    };
    if (!lookup_symbol_.size()) {
      space_group_symbols symbols = group_.match_tabulated_settings();
      if (symbols.number() != 0) {
        lookup_symbol_ = symbols.universal_hermann_mauguin();
      }
      else {
        lookup_symbol_ = universal_hermann_mauguin_symbol();
      }
    }
    if (   ad_hoc_1992
        && lookup_symbol_.size() >= 7
        && (lookup_symbol_.size() == 7 || lookup_symbol_[7] == ' ')) {
      for(int i=0;i<48;i+=2) {
        if (std::strncmp(lookup_symbol_.data(), adh_a38_pairs[i+1], 7) == 0) {
          std::string result = lookup_symbol_;
          for(int j=0;j<7;j++) {
            result[j] = adh_a38_pairs[i][j];
          }
          return result;
        }
      }
    }
    return lookup_symbol_;
  }

  space_group_type
  space_group::type() const
  {
    return space_group_type(*this);
  }

  space_group_type::space_group_type(
    std::string const& symbol,
    std::string const& table_id,
    bool tidy_cb_op)
  {
    *this = space_group_type(
      space_group(space_group_symbols(symbol, table_id)),
      tidy_cb_op);
  }

  af::shared<rt_mx>
  space_group_type::addl_generators_of_euclidean_normalizer(
    bool flag_k2l,
    bool flag_l2n) const
  {
    af::shared<rt_mx>
      result = reference_settings::normalizer::get_addl_generators(
        number_, false, flag_k2l, flag_l2n);
    change_of_basis_op cb_op_inv = cb_op_.inverse();
    for(std::size_t i=0;i<result.size();i++) {
      result[i] = cb_op_inv(result[i]);
    }
    return result;
  }

  bool
  space_group_type::is_enantiomorphic() const
  {
    if (group_.is_centric()) return false;
    af::shared<rt_mx>
      addl_g = reference_settings::normalizer::get_addl_generators(
        number_, false, true, false);
    if (addl_g.size() == 1) return false;
    CCTBX_ASSERT(addl_g.size() == 0);
    return true;
  }

  bool
  space_group_type::is_symmorphic() const
  {
    static const char tab[] = {
      0,1,1,1,0,1,1,0,1,0,1,0,1,0,0,0,1,0,0,0,0,1,1,1,0,1,0,0,0,0,0,0,0,0,0,1,
      0,0,1,0,0,0,1,0,1,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,1,0,1,
      0,0,0,1,0,0,0,1,0,1,1,1,0,0,0,1,0,1,0,0,0,0,0,0,0,1,0,1,0,0,0,0,0,0,0,1,
      0,0,0,1,0,0,0,1,0,0,0,1,0,1,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,1,
      0,0,1,1,1,1,1,0,0,0,0,1,1,1,0,0,1,0,1,0,1,0,1,0,1,0,0,0,0,0,1,1,0,1,0,0,
      0,0,0,1,0,0,0,1,0,1,0,1,0,0,0,1,1,1,0,0,1,0,1,0,1,0,0,1,0,1,0,1,0,0,0,1,
      1,1,0,0,0,1,0,0,0,1,0,0,0,1,0};
    return (tab[number_] != tab[0]);
  }

  change_of_basis_op
  space_group_type::change_of_hand_op() const
  {
    if (group_.is_centric()) return change_of_basis_op(1, group_.t_den());
    af::shared<rt_mx>
      addl_g = addl_generators_of_euclidean_normalizer(true, false);
    if (addl_g.size() == 1) return change_of_basis_op(addl_g[0]);
    CCTBX_ASSERT(addl_g.size() == 0);
    return change_of_basis_op(
      cb_op_.inverse()(rt_mx(rot_mx(1, -1), group_.t_den())));
  }

}} // namespace cctbx::sgtbx
