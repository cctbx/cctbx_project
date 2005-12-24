#include <cctbx/sgtbx/find_affine.h>
#include <scitbx/matrix/row_echelon.h>
#include <scitbx/math/unimodular_generator.h>
#include <scitbx/array_family/loops.h>

namespace cctbx { namespace sgtbx {

  find_affine::find_affine(
    space_group const& group,
    int range,
    bool use_p1_algorithm)
  {
    if (use_p1_algorithm || group.n_smx() == 1) {
      p1_algorithm(group, range);
    }
    else {
      sg_algorithm(group, range);
    }
  }

  void
  find_affine::p1_algorithm(space_group const& group, int range)
  {
    space_group tidy_group(group);
    tidy_group.make_tidy();
    scitbx::math::unimodular_generator<int> unimodular_generator(range);
    while (!unimodular_generator.at_end()) {
      change_of_basis_op cb_op((rt_mx(rot_mx(unimodular_generator.next()))));
      if (   tidy_group.n_smx() == 1
          || tidy_group.change_basis(cb_op) == tidy_group) {
        cb_mx_.push_back(cb_op.c());
      }
    }
  }

  namespace {

    void
    setup_affine_row_reduced_echelon_form(
      rot_mx const& r,
      std::vector<int>& m)
    {
      CCTBX_ASSERT(r.den() == 1);
      sg_mat3 const& r_num = r.num();
      // setup system of linear equations for coefficients of the matrix c,
      // given the condition c*r-r*c==0
      for(std::size_t i=0;i<3;i++) {
        for(std::size_t k=0;k<3;k++) {
          int* row = (&*m.begin())+m.size();
          m.resize(m.size()+9, 0);
          for(std::size_t j=0;j<3;j++) {
            if (i == j && j == k) continue;
            row[i*3+j] += r_num[j*3+k];
            row[j*3+k] -= r_num[i*3+j];
          }
        }
      }
    }

  }

  void
  find_affine::sg_algorithm(space_group const& group, int range)
  {
    space_group tidy_group(group);
    tidy_group.make_tidy();
    std::size_t n_rows = group.order_z() * 9;
    std::vector<int> m; m.reserve(n_rows * 9);
    for(int i=0;i<group.order_z();i++) {
      setup_affine_row_reduced_echelon_form(group(i).r(), m);
    }
    scitbx::mat_ref<int> m_ref(&*m.begin(), n_rows, 9);
    std::size_t r = scitbx::matrix::row_echelon::form(m_ref);
    scitbx::matrix::row_echelon::independent<int, 9> indep(m_ref);
    typedef af::nested_loop<af::small<int, 9> > loop_t;
    af::small<int, 9> loop_begin(indep.indices.size(), -range);
    af::small<int, 9> loop_end(indep.indices.size(), range+1);
    for(loop_t loop(loop_begin, loop_end);!loop.over();loop.incr()) {
      rot_mx c(1, 0);
      for(std::size_t i=0;i<indep.indices.size();i++) {
        c[indep.indices[i]] = loop()[i];
      }
      int den = scitbx::matrix::row_echelon::back_substitution_int(
        m_ref,
        static_cast<int*>(0),
        c.num().begin());
      CCTBX_ASSERT(den > 0);
      if (c.num().determinant() == 1) {
        std::size_t i;
        for(i=0;i<9;i++) {
          if (c[i] < -range) break;
          if (c[i] > range) break;
        }
        if (i == 9) {
          change_of_basis_op cb_op((rt_mx(c)));
          if (   tidy_group.n_smx() == 1
              || tidy_group.change_basis(cb_op) == tidy_group) {
            cb_mx_.push_back(cb_op.c());
          }
        }
      }
    }
  }

}} // namespace cctbx::sgtbx
