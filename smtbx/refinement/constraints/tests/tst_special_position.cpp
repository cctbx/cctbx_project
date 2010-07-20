#include <smtbx/refinement/constraints/special_position.h>
#include <smtbx/import_cctbx.h>
#include <scitbx/sparse/approx_equal.h>
#include <scitbx/sparse/io.h>
#include <iostream>

using namespace smtbx;
using namespace smtbx::refinement::constraints;
typedef crystallographic_parameter::scatterer_type sc_t;

namespace smtbx { namespace refinement { namespace constraints {

  #define site_approx_equal(s1, s2) \
    s1.ref().all_approx_equal(s2.const_ref(), 1e-15)

  void exercise_special_position(sgtbx::space_group const &sg) {
    uctbx::unit_cell uc(af::double6(10, 10, 10, 90, 90, 90));
    af::shared<sc_t> sc(3);
    sc_t *sc0 = &sc[0], *sc1 = &sc[1], *sc2 = &sc[2];
    sc0->site = frac_t(0.1 , 0.01, 0.01 );
    sc0->flags.set_grad_site(true);
    sc1->site = frac_t(0.01, 0.1 , 0.01 );
    sc1->flags.set_grad_site(true);
    sc2->site = frac_t(0.01, 0.01, 0.1  );
    sc2->flags.set_grad_site(true);
    sgtbx::site_symmetry site_symm_0(uc, sg, sc0->site),
                         site_symm_1(uc, sg, sc1->site),
                         site_symm_2(uc, sg, sc2->site);
    special_position_site *s[] = {
      new special_position_site(site_symm_0, sc0),
      new special_position_site(site_symm_1, sc1),
      new special_position_site(site_symm_2, sc2)
    };

    reparametrisation reparam(uc, boost::make_iterator_range(s, s+3));
    reparam.linearise();
    reparam.store();

    SMTBX_ASSERT(site_approx_equal(sc0->site, frac_t(0.1, 0. , 0.  )));
    SMTBX_ASSERT(site_approx_equal(sc1->site, frac_t(0. , 0.1, 0.  )));
    SMTBX_ASSERT(site_approx_equal(sc2->site, frac_t(0. , 0. , 0.1 )));

    scitbx::sparse::approx_equal<double> sparse_approx_equal(1.e-15);
    sparse_matrix_type jt = reparam.jacobian_transpose;
    sparse_matrix_type jt0(3, 12);

    for (int k=0; k<3; ++k) {
      for (int l=0; l<s[k]->independent_params().size(); ++l) {
        int idx = s[k]->independent_params().index() + l;
        jt0(idx, idx) = 1.;
        for (int i=0; i<3; ++i) {
          jt0(idx, s[k]->index() + i) = i == k ? 1. : 0.;
        }
      }
    }
    SMTBX_ASSERT(sparse_approx_equal(jt, jt0));
  }

}}}

int main() {
  using namespace smtbx::refinement::constraints;
  using cctbx::sgtbx::space_group;
  exercise_special_position(space_group("P 2 2"));
  std::cout << "OK\n";
}
