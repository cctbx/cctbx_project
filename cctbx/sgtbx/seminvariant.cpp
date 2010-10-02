#include <cctbx/sgtbx/seminvariant.h>
#include <cctbx/sgtbx/select_generators.h>
#include <cctbx/sgtbx/smith_normal_form.h>
#include <cctbx/sgtbx/row_echelon_solve.h>
#include <cctbx/sgtbx/utils.h>
#include <cctbx/math/loop_n_from_m.h>
#include <scitbx/array_family/loops.h>

namespace cctbx { namespace sgtbx {

  namespace {

    inline void
    copy(rot_mx const& source, int* target, std::size_t n)
    {
      for(std::size_t i=0;i<n;i++) target[i] = source[i];
    }

    af::tiny<int, 3 * 3*3>
    construct_gen_rmi(select_generators::any const& gen, bool primitive)
    {
      af::tiny<int, 3 * 3*3> result;
      for(std::size_t i=0;i<gen.n_gen;i++)
      {
        rot_mx const& r = gen.z_gen[i].r();
        if (!primitive) {
          copy(r.minus_unit_mx(), &result[i * 3*3], 3*3);
        }
        else {
          copy(gen.z2p_op(r).minus_unit_mx(), &result[i * 3*3], 3*3);
        }
      }
      if (gen.z_inv_t.is_valid()) {
        copy(rot_mx(1, -2), &result[gen.n_gen * 3*3], 3*3);
      }
      return result;
    }

    class cmp_o_len_sq
    {
      public:
        cmp_o_len_sq() : cmp_tr(3) {}

        bool operator()(af::int3 const& a, af::int3 const& b) const
        {
          using std::size_t;
          int o_len_sq_a=0; for(size_t i=0;i<3;i++) o_len_sq_a += a[i] * a[i];
          int o_len_sq_b=0; for(size_t i=0;i<3;i++) o_len_sq_b += b[i] * b[i];
          if (o_len_sq_a < o_len_sq_b) return true;
          if (o_len_sq_a > o_len_sq_b) return false;
          return cmp_tr(a.begin(), b.begin());
        }
      private:
        const utils::cmp_i_vec cmp_tr;
    };

    af::small<ss_vec_mod, 3>
    get_cont_null_space(select_generators::any const& gen)
    {
      af::small<ss_vec_mod, 3> result;
      af::tiny<int, 3 * 3 * 3> gen_rmi = construct_gen_rmi(gen, false);
      af::ref<int, af::mat_grid> re_mx(gen_rmi.begin(), gen.n_all() * 3, 3);
      CCTBX_ASSERT(scitbx::matrix::row_echelon::form(re_mx) <= 3);
      scitbx::matrix::row_echelon::independent<int> indep(re_mx);
      if (indep.indices.size() != 2) {
        for (std::size_t i_indep=0;i_indep<indep.indices.size();i_indep++) {
          ss_vec_mod vm;
          if (gen.n_all() == 0) vm.v.fill(0);
          vm.v[indep.indices[i_indep]] = 1;
          int* n_a = 0;
          CCTBX_ASSERT(scitbx::matrix::row_echelon::back_substitution_int(
            re_mx, n_a, vm.v.begin()) > 0);
          vm.m = 0;
          result.push_back(vm);
        }
      }
      else {
        af::tiny<sg_vec3, 4>
          sol = row_echelon::solve::homog_rank_1(re_mx, indep);
        std::sort(sol.begin(), sol.end(), cmp_o_len_sq());
        for(std::size_t i_indep=0;i_indep<2;i_indep++) {
          ss_vec_mod vm;
          vm.v = sol[i_indep];
          vm.m = 0;
          result.push_back(vm);
        }
      }
      return result;
    }

    typedef af::small<tr_vec, 8> tr_vec_8;

    void update_best_z(tr_vec_8 const& orig_z_f,
                       tr_vec const& shift,
                       tr_vec_8& best_z_f,
                       tr_vec_8& best_z_c)
    {
      for (std::size_t i_dl=1;i_dl<orig_z_f.size();i_dl++) {
        tr_vec z_f = (orig_z_f[i_dl] + shift).mod_positive();
        tr_vec z_c = z_f.cancel();
        for(std::size_t i=0;i<3;i++) {
          if (z_f[i]) {
            if (   cmp_o_len_sq()(z_c.num(), best_z_c[i_dl].num())
                || (   z_c.num() == best_z_c[i_dl].num()
                    && z_c.den() <  best_z_c[i_dl].den())) {
              best_z_f[i_dl] = z_f;
              best_z_c[i_dl] = z_c;
            }
            break;
          }
        }
      }
    }

    void best_vectors(space_group const& sg,
                      af::small<ss_vec_mod, 3> const& continuous_vm,
                      tr_vec_8& discr_z)
    {
      if (sg.n_ltr() == 1 && continuous_vm.size() == 0) return;
      int lt_den = 1;
      for(std::size_t i_dl=1;i_dl<discr_z.size();i_dl++) {
        int den = discr_z[i_dl].den();
        for(std::size_t i=0;i<3;i++) {
          int g = scitbx::math::gcd_int(discr_z[i_dl][i], den);
          lt_den = boost::lcm(lt_den, den / g);
        }
      }
      for(std::size_t i_ltr=1;i_ltr<sg.n_ltr();i_ltr++) {
        int den = sg.ltr(i_ltr).den();
        for(std::size_t i=0;i<3;i++) {
          int g = scitbx::math::gcd_int(sg.ltr(i_ltr)[i], den);
          lt_den = boost::lcm(lt_den, den / g);
        }
      }
      int f_grd = 1;
      for(std::size_t i_vm = 0; i_vm < continuous_vm.size(); i_vm++) {
        for(std::size_t i=0;i<3;i++) {
          if (continuous_vm[i_vm].v[i]) {
            f_grd = boost::lcm(f_grd, continuous_vm[i_vm].v[i]);
          }
        }
      }
      lt_den *= f_grd;
      CCTBX_ASSERT(lt_den > 0);
      if (lt_den > 6) lt_den = boost::lcm(lt_den,  6);
      else            lt_den = boost::lcm(lt_den, 12);
      tr_vec_8 orig_z_f;
      tr_vec_8 best_z_f;
      tr_vec_8 best_z_c;
      for(std::size_t i_dl=0;i_dl<discr_z.size();i_dl++) {
        orig_z_f.push_back(discr_z[i_dl].new_denominator(lt_den));
        best_z_f.push_back(orig_z_f[i_dl]);
        best_z_c.push_back(best_z_f[i_dl].cancel());
      }
      CCTBX_ASSERT(continuous_vm.size() < 3);
      af::small<int, 2> loop_end(continuous_vm.size());
      std::fill(loop_end.begin(), loop_end.end(), lt_den);
      tr_vec ltr[3];
      for(std::size_t i_ltr=0;i_ltr<sg.n_ltr();i_ltr++) {
        ltr[0] = sg.ltr(i_ltr).new_denominator(lt_den);
        af::nested_loop<af::small<int, 2> > loop(loop_end);
        af::small<int, 2> const& f = loop();
        do {
          for(std::size_t i_vm=0;i_vm<continuous_vm.size();i_vm++) {
            ltr[i_vm + 1] = ltr[i_vm]
                          + f[i_vm] * tr_vec(continuous_vm[i_vm].v, lt_den);
          }
          update_best_z(orig_z_f, ltr[continuous_vm.size()],
                        best_z_f, best_z_c);
          loop.incr();
        }
        while (loop.over() == 0);
      }
      for(std::size_t i_dl=1;i_dl<discr_z.size();i_dl++) {
        discr_z[i_dl] = best_z_f[i_dl].new_denominator(discr_z[i_dl].den());
      }
    }

    af::small<tr_vec, 3>
    select_discrete_generators(tr_vec_8 const& discr_p,
                               tr_vec_8 const& discr_z)
    {
      if (discr_p.size() == 1) return af::small<tr_vec, 3>();
      for(std::size_t n_d_vm=1;
          n_d_vm<=discr_p.size()-1 && n_d_vm<=3;
          n_d_vm++)
      {
        for(math::loop_n_from_m<3> ix(discr_p.size() - 1, n_d_vm);
            !ix.over();
            ix.incr())
        {
          tr_group discr_group(discr_p[0].den());
          for(std::size_t i_ix=0;i_ix<ix.n();i_ix++) {
            discr_group.expand(discr_p[ix[i_ix] + 1]);
          }
          if (discr_group.size() == discr_p.size()) {
            af::small<tr_vec, 3> result;
            for(std::size_t i_ix=0;i_ix<ix.n();i_ix++) {
              result.push_back(discr_z[ix[i_ix] + 1]);
            }
            return result;
          }
          CCTBX_ASSERT(discr_group.size() < discr_p.size());
        }
      }
      throw CCTBX_INTERNAL_ERROR();
    }

    class cmp_discr_z
    {
      public:
        cmp_discr_z() : cmp_tr(3) {}

        bool
        operator()(tr_vec const& a, tr_vec const& b) const
        {
          return cmp_tr(a.num().begin(), b.num().begin());
        }

      private:
        const utils::cmp_i_vec cmp_tr;
    };

    class cmp_ss_vec_mod
    {
      public:
        cmp_ss_vec_mod() : cmp_tr(3) {}

        bool
        operator()(ss_vec_mod const& a, ss_vec_mod const& b) const
        {
          return cmp_tr(a.v.begin(), b.v.begin());
        }

      private:
        const utils::cmp_i_vec cmp_tr;
    };

  } // namespace <anonymous>

  structure_seminvariants::structure_seminvariants(space_group const& sg)
  {
    select_generators::any gen(sg, cb_r_den, cb_t_den);
    vec_mod_ = get_cont_null_space(gen);
    if (vec_mod_.size() == 3) return; // space group P1
    af::tiny<int, 3 * 3 * 3> snf = construct_gen_rmi(gen, true);
    sg_mat3 q;
    af::ref<int, af::mat_grid> m_ref(snf.begin(), gen.n_all() * 3, 3);
    af::ref<int, af::mat_grid> p_ref(0, 0, 0);
    af::ref<int, af::mat_grid> q_ref(q.begin(), 3, 3);
    smith_normal_form(m_ref, p_ref, q_ref);
    CCTBX_ASSERT(m_ref.is_square());
    std::size_t nd = m_ref.n_rows();
    CCTBX_ASSERT(nd <= 3);
    int d_t_den = 1;
    for(std::size_t id=0;id<nd;id++) {
      d_t_den = boost::lcm(d_t_den, snf[(nd + 1) * id]);
    }
    tr_group discr_group(d_t_den);
    for(std::size_t id=0;id<nd;id++) {
      int d = snf[(nd + 1) * id];
      for(int f=1;f<d;f++) {
        sg_vec3 xp(0,0,0);
        xp[id] = f * d_t_den / d;
        tr_vec x(q * xp, d_t_den);
        discr_group.expand(x);
      }
    }
    tr_vec_8 discr_z;
    for(std::size_t i_dl=0;i_dl<discr_group.size();i_dl++) {
      discr_z.push_back(
        (gen.z2p_op.c_inv().r() * discr_group[i_dl]).mod_positive());
    }
    best_vectors(sg, vec_mod_, discr_z);
    std::sort(discr_z.begin(), discr_z.end(), cmp_discr_z());
    tr_vec_8 discr_p;
    for(std::size_t i_dl=0;i_dl<discr_z.size();i_dl++) {
      discr_p.push_back(
        (gen.z2p_op.c().r() * discr_z[i_dl]).new_denominator(d_t_den));
    }
    af::small<tr_vec, 3>
      discr_gen = select_discrete_generators(discr_p, discr_z);
    for(std::size_t i_g=0;i_g<discr_gen.size();i_g++) {
      CCTBX_ASSERT(vec_mod_.size() < 3);
      tr_vec g = discr_gen[i_g].cancel();
      ss_vec_mod vm;
      vm.v = g.num();
      vm.m = g.den();
      vec_mod_.push_back(vm);
    }
    std::sort(vec_mod_.begin(), vec_mod_.end(), cmp_ss_vec_mod());
  }

  bool
  structure_seminvariants::is_ss(miller::index<> const& h) const
  {
    for(std::size_t i=0;i<vec_mod_.size();i++) {
      int u = vec_mod_[i].v * h;
      if (vec_mod_[i].m) {
        if (u % vec_mod_[i].m) return false;
      }
      else {
        if (u) return false;
      }
    }
    return true;
  }

  af::small<int, 3>
  structure_seminvariants::apply_mod(miller::index<> const& h) const
  {
    af::small<int, 3> result;
    for(std::size_t i=0;i<vec_mod_.size();i++) {
      result.push_back(vec_mod_[i].v * h);
      if (vec_mod_[i].m) result[i] %= vec_mod_[i].m;
    }
    return result;
  }

}} // namespace cctbx::sgtbx
