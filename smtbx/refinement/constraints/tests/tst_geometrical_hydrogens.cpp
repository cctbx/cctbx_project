#include <smtbx/refinement/constraints/geometrical_hydrogens.h>

#include <scitbx/math/r3_rotation.h>
#include <scitbx/math/approx_equal.h>
#include <scitbx/sparse/approx_equal.h>
#include <scitbx/sparse/io.h>

using namespace smtbx;
using namespace smtbx::refinement::constraints;

typedef crystallographic_parameter::scatterer_type sc_t;


void exercise_ch3() {
  uctbx::unit_cell uc(af::double6(1, 2, 3, 90, 90, 90));
  af::shared<sc_t> sc(5);
  sc_t *c0 = &sc[0], *c1 = &sc[1];

  af::small<sc_t *, 3> h;
  h.push_back(&sc[2]);
  h.push_back(&sc[3]);
  h.push_back(&sc[4]);

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

void exercise_secondary_ch2() {
  uctbx::unit_cell uc(af::double6(1, 2, 3, 90, 90, 90));

  af::shared<sc_t> sc(5);
  sc_t *c0 = &sc[0], *c1 = &sc[1], *c2 = &sc[2], *h0 = &sc[3], *h1 = &sc[4];
  c0->site = frac_t( 0.,  0.,  0.);
  c0->flags.set_grad_site(true);
  double delta = 0.1;
  double c = std::cos(constants::tetrahedral_angle/2 + delta),
         s = std::sin(constants::tetrahedral_angle/2 + delta);
  cart_t u = cart_t(1, 2, -1).normalize(), v = u.ortho();
  c1->site = uc.fractionalize(-s*u + c*v);
  c2->site = uc.fractionalize( s*u + c*v);
  independent_site_parameter *is_c0 = new independent_site_parameter(c0),
                             *is_c1 = new independent_site_parameter(c1),
                             *is_c2 = new independent_site_parameter(c2);
  independent_scalar_parameter
  *length = new independent_scalar_parameter(0.9, false);
  angle_starting_tetrahedral *h_c_h = new angle_starting_tetrahedral();
  secondary_ch2_sites
  *ch2 = new secondary_ch2_sites(is_c0, is_c1, is_c2, length, h_c_h,
                                 h0, h1);

  reparametrisation reparam(uc, boost::make_iterator_range(&ch2, &ch2 + 1));
  reparam.linearise();
  reparam.store();

  scitbx::math::approx_equal_absolutely<double> scalar_approx_equal(1e-15);

  // Check geometry
  cart_t ch0  = uc.orthogonalize(h0->site - c0->site),
         ch1  = uc.orthogonalize(h1->site - c0->site),
         c0c1 = uc.orthogonalize(c1->site - c0->site),
         c0c2 = uc.orthogonalize(c2->site - c0->site);

  SMTBX_ASSERT(scalar_approx_equal(ch0.angle(ch1),
                                   constants::tetrahedral_angle));
  SMTBX_ASSERT(scalar_approx_equal(ch0.angle(c0c1), ch0.angle(c0c2)));
  SMTBX_ASSERT(scalar_approx_equal(ch1.angle(c0c1), ch1.angle(c0c2)));

  // Check gradient
  sparse_matrix_type jt = reparam.jacobian_transpose.clone();
  sparse_matrix_type jt0(4, 17);
  double eps = 1.e-6;
  frac_t xm_h0, xm_h1, xp_h0, xp_h1;
  scitbx::sparse::approx_equal<double> sparse_approx_equal(eps);

  for (int i=0; i<3; ++i) {
    jt0(is_c0->index() + i, is_c0->index() + i) = 1.;
  }
  jt0(h_c_h->index(), h_c_h->index()) = 1.;

  for (int i=0; i<3; ++i) {
    jt0(is_c0->index() + i, ch2->index() + i) = 1.;
    jt0(is_c0->index() + i, ch2->index() + 3 + i) = 1.;
  }

  h_c_h->value += eps;
  reparam.linearise();
  reparam.store();
  xp_h0 = h0->site;
  xp_h1 = h1->site;

  h_c_h->value -= 2*eps;
  reparam.linearise();
  reparam.store();
  xm_h0 = h0->site;
  xm_h1 = h1->site;

  for (int i=0; i<3; ++i) {
    jt0(h_c_h->index(), ch2->index() + i)     = (xp_h0[i] - xm_h0[i])/(2*eps);
    jt0(h_c_h->index(), ch2->index() + 3 + i) = (xp_h1[i] - xm_h1[i])/(2*eps);
  }

  SMTBX_ASSERT(sparse_approx_equal(jt, jt0));
}

void exercise_tertiary_ch() {
  using constants::pi;
  uctbx::unit_cell uc(af::double6(1, 2, 3, 90, 90, 90));

  af::shared<sc_t> sc(5);
  sc_t *c = &sc[0], *x = &sc[1], *y = &sc[2], *z = &sc[3], *h = &sc[4];

  c->site = frac_t( 0.,  0.,  0.);
  c->flags.set_grad_site(true);
  cart_t u = cart_t(2, -1, 1).normalize(), v = u.ortho(), w = u.cross(v);
  // classic embedding of a tetrahedron inside a cube
  cart_t delta = 0.1*cart_t(1, -2, 3).normalize();
  x->site = uc.fractionalize(-u -v +w + delta);
  y->site = uc.fractionalize(-u +v -w + delta);
  z->site = uc.fractionalize(+u -v -w + delta);

  independent_site_parameter *is_c = new independent_site_parameter(c),
                             *is_x = new independent_site_parameter(x),
                             *is_y = new independent_site_parameter(y),
                             *is_z = new independent_site_parameter(z);
  independent_scalar_parameter
  *length = new independent_scalar_parameter(0.9, false);
  tertiary_ch_site
  *tert_ch = new tertiary_ch_site(is_c, is_x, is_y, is_z, length, h);

  reparametrisation reparam(uc,
                            boost::make_iterator_range(&tert_ch, &tert_ch + 1));
  reparam.linearise();
  reparam.store();

  // Test geometry
  scitbx::math::approx_equal_absolutely<double> scalar_approx_equal(1e-15);
  cart_t ch = uc.orthogonalize(h->site - c->site),
         cx = uc.orthogonalize(x->site - c->site),
         cy = uc.orthogonalize(y->site - c->site),
         cz = uc.orthogonalize(z->site - c->site);
  SMTBX_ASSERT(scalar_approx_equal(ch.angle(cx), ch.angle(cy)));
  SMTBX_ASSERT(scalar_approx_equal(ch.angle(cy), ch.angle(cz)));
  SMTBX_ASSERT(scalar_approx_equal(ch.angle(cz), ch.angle(cx)));
  SMTBX_ASSERT(scalar_approx_equal(ch.length(), length->value));

  // Jacobian
  scitbx::sparse::approx_equal<double> sparse_approx_equal(1e-15);
  sparse_matrix_type jt = reparam.jacobian_transpose.clone();
  sparse_matrix_type jt0(3, 16);
  for (int i=0; i<3; ++i) {
    jt0(is_c->index() + i, is_c->index() + i) = 1.;
  }
  for (int i=0; i<3; ++i) {
    jt0(is_c->index() + i, tert_ch->index() + i) = 1.;
  }
  SMTBX_ASSERT(sparse_approx_equal(jt, jt0));
}

void exercise_aromatic_ch() {
  using constants::pi;
  uctbx::unit_cell uc(af::double6(1, 2, 3, 90, 90, 90));
  af::shared<sc_t> sc(4);
  sc_t *x = &sc[0], *y = &sc[1], *z = &sc[2], *h = &sc[3];
  x->site = frac_t( 0.,  0.,  0.);
  x->flags.set_grad_site(true);
  cart_t xy = cart_t(-2, -1, 1).normalize();
  cart_t v = cart_t(1, 1, 1).normalize();
  cart_t n = xy.cross(v);
  scitbx::mat3<double>
  r = scitbx::math::r3_rotation::axis_and_angle_as_matrix(n, 3*pi/4);
  cart_t xz = r*xy;
  y->site = uc.fractionalize(xy);
  z->site = uc.fractionalize(xz);

  independent_site_parameter *is_x = new independent_site_parameter(x),
                             *is_y = new independent_site_parameter(y),
                             *is_z = new independent_site_parameter(z);
  independent_scalar_parameter
  *length = new independent_scalar_parameter(1.1);
  secondary_planar_xh_site
  *arom_ch = new secondary_planar_xh_site(is_x, is_y, is_z, length, h);
  reparametrisation reparam(uc,
                            boost::make_iterator_range(&arom_ch,
                                                       &arom_ch + 1));
  reparam.linearise();
  reparam.store();

  // Check geometry
  scitbx::math::approx_equal_absolutely<double> scalar_approx_equal(1e-15);
  cart_t xh = uc.orthogonalize(h->site - x->site);
  SMTBX_ASSERT(scalar_approx_equal(xh.angle(xy), xh.angle(xy)));
  SMTBX_ASSERT(scalar_approx_equal(xh*xy.cross(xz), 0));

  // Check Jacobian
  scitbx::sparse::approx_equal<double> sparse_approx_equal(1e-15);
  sparse_matrix_type jt = reparam.jacobian_transpose.clone();
  sparse_matrix_type jt0(4, 13);
  for (int i=0; i<3; ++i) {
    jt0(is_x->index() + i, is_x->index() + i) = 1.;
  }
  jt0(length->index(), length->index()) = 1.;
  for (int i=0; i<3; ++i) {
    jt0(is_x->index() + i, arom_ch->index() + i) = 1.;
  }
  frac_t u_h = uc.fractionalize(xh.normalize());
  for (int i=0; i<3; ++i) {
    jt0(length->index(), arom_ch->index() + i) = u_h[i];
  }
  SMTBX_ASSERT(sparse_approx_equal(jt, jt0));
}

void exercise_terminal_xh2() {
  using constants::pi;
  uctbx::unit_cell uc(af::double6(1, 2, 3, 90, 90, 90));
  af::shared<sc_t> sc(5);
  sc_t *x = &sc[0], *y = &sc[1], *z = &sc[2], *h0 = &sc[3], *h1 = &sc[4];
  x->site = frac_t( 0.,  0.,  0.);
  x->flags.set_grad_site(true);
  y->site = frac_t( 1., 1., 1.);
  z->site = frac_t( 2., 0., 2.);

  independent_site_parameter *is_x = new independent_site_parameter(x),
                             *is_y = new independent_site_parameter(y),
                             *is_z = new independent_site_parameter(z);
  independent_scalar_parameter
  *length = new independent_scalar_parameter(1.1, false);
  terminal_planar_xh2_sites
  *ch2 = new terminal_planar_xh2_sites(is_x, is_y, is_z, length,
                                       h0, h1);
  reparametrisation reparam(uc,
                            boost::make_iterator_range(&ch2, &ch2 + 1));
  reparam.linearise();
  reparam.store();

  // Check geometry
  scitbx::math::approx_equal_absolutely<double> scalar_approx_equal(5.e-15);
    // Need that slightly larger tolerance for gcc 4.4.0 on 64-bit Linux
  cart_t xh0 = uc.orthogonalize(h0->site - x->site),
         xh1 = uc.orthogonalize(h1->site - x->site),
         xy  = uc.orthogonalize(y->site - x->site),
         yz  = uc.orthogonalize(z->site - y->site);
  SMTBX_ASSERT(scalar_approx_equal(xh0*xy.cross(yz), 0));
  SMTBX_ASSERT(scalar_approx_equal(xh1*xy.cross(yz), 0));
  SMTBX_ASSERT(scalar_approx_equal(xh0.angle(xy), 2*pi/3));
  SMTBX_ASSERT(scalar_approx_equal(xh1.angle(xy), 2*pi/3));
  SMTBX_ASSERT(scalar_approx_equal(xh0.length(), length->value));
  SMTBX_ASSERT(scalar_approx_equal(xh1.length(), length->value));

  // Jacobian
  scitbx::sparse::approx_equal<double> sparse_approx_equal(1e-15);
  sparse_matrix_type jt = reparam.jacobian_transpose.clone();
  sparse_matrix_type jt0(3, 16);
  for (int i=0; i<3; ++i) {
    jt0(is_x->index() + i, is_x->index() + i) = 1.;
  }
  for (int k=0; k<2; ++k) for (int i=0; i<3; ++i) {
    jt0(is_x->index() + i, ch2->index() + 3*k + i) = 1.;
  }
  SMTBX_ASSERT(sparse_approx_equal(jt, jt0));
}

void exercise_acetylenic_ch() {
  using constants::pi;
  uctbx::unit_cell uc(af::double6(1, 2, 3, 90, 90, 90));

  af::shared<sc_t> sc(3);
  sc_t *c = &sc[0], *x = &sc[1], *h = &sc[2];
  c->site = frac_t( 0.,  0.,  0.);
  c->flags.set_grad_site(true);
  x->site = frac_t( 1., 1., 1.);

  independent_site_parameter *is_c = new independent_site_parameter(c),
                             *is_x = new independent_site_parameter(x);
  independent_scalar_parameter *length = new independent_scalar_parameter(1.1);
  terminal_linear_ch_site
  *term_ch = new terminal_linear_ch_site(is_c, is_x, length, h);

  reparametrisation reparam(uc,
                            boost::make_iterator_range(&term_ch, &term_ch + 1));
  reparam.linearise();
  reparam.store();

  // Check geometry
  scitbx::math::approx_equal_absolutely<double> scalar_approx_equal(1e-15);
  cart_t ch = uc.orthogonalize(h->site - c->site),
         cx = uc.orthogonalize(x->site - c->site);
  SMTBX_ASSERT(scalar_approx_equal(ch.length(), length->value));
  SMTBX_ASSERT(scalar_approx_equal(ch*cx, -ch.length()*cx.length()))
              (ch*cx + ch.length()*cx.length());

  // Jacobian
  scitbx::sparse::approx_equal<double> sparse_approx_equal(1e-15);
  sparse_matrix_type jt = reparam.jacobian_transpose.clone();
  sparse_matrix_type jt0(4, 10);
  for (int i=0; i<3; ++i) {
    jt0(is_c->index() + i, is_c->index() + i) = 1.;
  }
  jt0(length->index(), length->index()) = 1.;
  for (int i=0; i<3; ++i) {
    jt0(is_c->index() + i, term_ch->index() + i) = 1.;
  }
  frac_t u_h = uc.fractionalize(ch.normalize());
  for (int i=0; i<3; ++i) {
    jt0(length->index(), term_ch->index() + i) = u_h[i];
  }
  SMTBX_ASSERT(sparse_approx_equal(jt, jt0));
}

int main() {
  exercise_ch3();
  exercise_secondary_ch2();
  exercise_tertiary_ch();
  exercise_aromatic_ch();
  exercise_terminal_xh2();
  exercise_acetylenic_ch();
  std::cout << "OK\n";
  return 0;
}
