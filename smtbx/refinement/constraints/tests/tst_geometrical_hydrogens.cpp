#include <smtbx/refinement/constraints/geometrical_hydrogens.h>

#include <scitbx/math/approx_equal.h>
#include <scitbx/sparse/approx_equal.h>
#include <scitbx/sparse/io.h>

using namespace smtbx;
using namespace smtbx::refinement::constraints;

typedef crystallographic_parameter::scatterer_type sc_t;


void exercise_ch3() {
  uctbx::unit_cell uc(af::double6(1, 2, 3, 90, 90, 90));

  boost::shared_ptr<sc_t> c1(new sc_t()), c0(new sc_t());
  af::small<boost::shared_ptr<sc_t>, 3> h;
  h.push_back(boost::shared_ptr<sc_t>( new sc_t()));
  h.push_back(boost::shared_ptr<sc_t>( new sc_t()));
  h.push_back(boost::shared_ptr<sc_t>( new sc_t()));

  c1->site = frac_t(-1., -1., -1.);
  c0->site = frac_t( 0.,  0.,  0.);
  c0->flags.set_grad_site(true);
  independent_site_parameter *is_c1 = new independent_site_parameter(c1),
                             *is_c0 = new independent_site_parameter(c0);
  independent_scalar_parameter
    *azimuth = new independent_scalar_parameter(constants::pi/3.),
    *length  = new independent_scalar_parameter(0.9);

  scitbx::vec3<double> e_zero_azimuth(1., 0., 0.);
  terminal_tetrahedral_xhn_sites
  *ch3 = new terminal_tetrahedral_xhn_sites(is_c0, is_c1, azimuth, length,
                                            e_zero_azimuth, h);

  reparametrisation reparam(uc, boost::make_iterator_range(&ch3, &ch3 + 1));
  reparam.linearise();
  reparam.store();

  scitbx::math::approx_equal_absolutely<double> scalar_approx_equal(1e-15);

  // Check geometry
  af::tiny<cart_t, 3> ch;
  for (int k=0; k<3; ++k) ch[k] = uc.orthogonalize(h[k]->site - c0->site);
  cart_t cc = uc.orthogonalize(c1->site - c0->site);
  SMTBX_ASSERT(scalar_approx_equal(ch[0].length(), length->value));
  SMTBX_ASSERT(scalar_approx_equal(ch[1].length(), length->value));
  SMTBX_ASSERT(scalar_approx_equal(ch[2].length(), length->value));
  SMTBX_ASSERT(scalar_approx_equal(ch[0].angle(ch[1]),
                                   constants::tetrahedral_angle));
  SMTBX_ASSERT(scalar_approx_equal(ch[1].angle(ch[2]),
                                   constants::tetrahedral_angle));
  SMTBX_ASSERT(scalar_approx_equal(ch[2].angle(ch[0]),
                                   constants::tetrahedral_angle));
  SMTBX_ASSERT(scalar_approx_equal(ch[0].angle(cc),
                                   constants::tetrahedral_angle));
  SMTBX_ASSERT(scalar_approx_equal(ch[1].angle(cc),
                                   constants::tetrahedral_angle));
  SMTBX_ASSERT(scalar_approx_equal(ch[2].angle(cc),
                                   constants::tetrahedral_angle));

  // Check gradients
  sparse_matrix_type jt = reparam.jacobian_transpose.clone();
  sparse_matrix_type jt0(5, 17);
  double eps = 1.e-6;
  af::tiny<frac_t, 3> xm_h, xp_h;
  scitbx::sparse::approx_equal<double> sparse_approx_equal(eps);

  for (int i=0; i<3; ++i) {
    jt0(is_c0->index() + i, is_c0->index() + i) = 1.;
  }
  jt0(azimuth->index(), azimuth->index()) = 1.;
  jt0(length->index(), length->index()) = 1.;

  for (int k=0; k<3; ++k) for (int i=0; i<3; ++i) {
    jt0(is_c0->index() + i, ch3->index() + 3*k + i) = 1.;
  }

  azimuth->value += eps;
  reparam.linearise();
  reparam.store();
  for (int k=0; k<3; ++k) xp_h[k] = h[k]->site;

  azimuth->value -= 2*eps;
  reparam.linearise();
  reparam.store();
  for (int k=0; k<3; ++k) xm_h[k] = h[k]->site;

  for (int k=0; k<3; ++k) for (int i=0; i<3; ++i) {
    jt0(azimuth->index(),
        ch3->index() + 3*k + i) = (xp_h[k][i] - xm_h[k][i])/(2*eps);
  }

  length->value += eps;
  reparam.linearise();
  reparam.store();
  for (int k=0; k<3; ++k) xp_h[k] = h[k]->site;

  length->value -= 2*eps;
  reparam.linearise();
  reparam.store();
  for (int k=0; k<3; ++k) xm_h[k] = h[k]->site;

  for (int k=0; k<3; ++k) for (int i=0; i<3; ++i) {
    jt0(length->index(),
        ch3->index() + 3*k + i) = (xp_h[k][i] - xm_h[k][i])/(2*eps);
  }

  SMTBX_ASSERT(sparse_approx_equal(jt, jt0));
}

int main() {
  exercise_ch3();
  std::cout << "OK\n";
  return 0;
}
