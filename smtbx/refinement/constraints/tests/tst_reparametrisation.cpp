#include <smtbx/refinement/constraints/reparametrisation.h>

#include <scitbx/array_family/ref_reductions.h>
#include <scitbx/sparse/io.h>
#include <scitbx/sparse/approx_equal.h>
#include <iostream>

using namespace smtbx;
using namespace smtbx::refinement::constraints;

class dependent_site_1 : public site_parameter
{
public:
  dependent_site_1(scatterer_type *scatterer,
                   site_parameter *site_1,
                   site_parameter *site_2)
    : site_parameter(scatterer, 2)
  {
    set_arguments(site_1, site_2);
  }

  virtual void linearise(uctbx::unit_cell const &unit_cell,
                         sparse_matrix_type *jt)
  {
    site_parameter *s1 = (site_parameter *)argument(0),
                   *s2 = (site_parameter *)argument(1);
    value = s1->value + s2->value;
    if (!jt) return;
    std::size_t j_s1=s1->index(), j_s2=s2->index(), j=index();
    for (std::size_t k=0; k<3; ++k) {
      jt->col(j+k) = jt->col(j_s1+k) + jt->col(j_s2+k);
    }
  }
};


class dependent_site_2 : public site_parameter
{
public:
  dependent_site_2(scatterer_type *scatterer,
                   site_parameter *site_1,
                   site_parameter *site_2)
  : site_parameter(scatterer, 2)
  {
    set_arguments(site_1, site_2);
  }

  virtual void linearise(uctbx::unit_cell const &unit_cell,
                         sparse_matrix_type *jt)
  {
    site_parameter *s1 = (site_parameter *)argument(0),
                   *s2 = (site_parameter *)argument(1);
    value = s1->value - s2->value;
    if (!jt) return;
    std::size_t j_s1=s1->index(), j_s2=s2->index(), j=index();
    for (std::size_t k=0; k<3; ++k) {
      jt->col(j+k) = jt->col(j_s1+k) - jt->col(j_s2+k);
    }
  }
};


class dependent_site_3 : public site_parameter
{
public:
  dependent_site_3(scatterer_type *scatterer,
                   site_parameter *site,
                   independent_scalar_parameter *x)
    : site_parameter(scatterer, 2)
  {
    set_arguments(site, x);
  }

  virtual void linearise(uctbx::unit_cell const &unit_cell,
                         sparse_matrix_type *jt)
  {
    site_parameter *s = (site_parameter *)argument(0);
    independent_scalar_parameter *x = (independent_scalar_parameter *)argument(1);
    value = x->value * s->value;
    if (!jt) return;
    std::size_t j_s=s->index(), j_x=x->index(), j=index();
    for (std::size_t k=0; k<3; ++k) {
      jt->col(j+k) = x->value*jt->col(j_s+k) + jt->col(j_x)*s->value[k];
    }
  }
};

#define site_approx_equal(s1, s2) \
  s1.ref().all_approx_equal(s2.const_ref(), 1e-15)


class test_case
{
public:
  typedef crystallographic_parameter::scatterer_type sc_t;

  uctbx::unit_cell uc;
  af::shared<sc_t> sc;
  independent_site_parameter *is1, *is2, *is3, *is4;

  site_parameter *s1, *s2, *s3, *s4;
  site_parameter *s5, *s6, *s7, *s8;

  independent_scalar_parameter *s11;
  site_parameter *s10, *s9;


  test_case()
    : uc(af::double6(1, 2, 3, 90, 90, 90)),
      sc(11)
  {
    sc[1].site = frac_t( 0.1,  0.2,  0.3);
    sc[1].flags.set_grad_site(true);
    sc[2].site = frac_t(-0.1, -0.2, -0.3);
    sc[2].flags.set_grad_site(true);
    sc[3].site = frac_t( 0.4,  0.3,  0.2);
    sc[4].site = frac_t(-0.5, -0.6, -0.7);

    sc[5].site = frac_t(1, 2, 3);
    sc[6].site = frac_t(4, 5, 6);
    sc[7].site = frac_t(-1, -2, -3);
    sc[8].site = frac_t(-4, -5, -6);
    sc[9].site = frac_t(10, 11, 12);
    sc[10].site = frac_t(20, 21, 22);

    is1 = new independent_site_parameter(&sc[1]);
    is2 = new independent_site_parameter(&sc[2]);
    is3 = new independent_site_parameter(&sc[3]);
    is4 = new independent_site_parameter(&sc[4]);
    s1 = is1; s2 = is2; s3 = is3; s4 = is4;

    s5 = new dependent_site_1(&sc[5], is1, is2);
    s6 = new dependent_site_1(&sc[6], is3, is4);
    s7 = new dependent_site_2(&sc[7], s5, s2);
    s8 = new dependent_site_2(&sc[8], s7, s6);

    s11 = new independent_scalar_parameter(2.);
    s10 = new dependent_site_3(&sc[10], s6, s11);
    s9 = new dependent_site_1(&sc[9], s10, s3);
  }

  void check_parameter_status() {
    SMTBX_ASSERT(!s1->is_root());
    SMTBX_ASSERT(!s2->is_root());
    SMTBX_ASSERT(!s3->is_root());
    SMTBX_ASSERT(!s4->is_root());
    SMTBX_ASSERT(!s5->is_root());
    SMTBX_ASSERT(!s6->is_root());
    SMTBX_ASSERT(!s7->is_root());
    SMTBX_ASSERT(s8->is_root());
    SMTBX_ASSERT(s9->is_root());
    SMTBX_ASSERT(!s10->is_root());
    SMTBX_ASSERT(!s11->is_root());

    SMTBX_ASSERT(s8->is_variable());
    SMTBX_ASSERT(s7->is_variable());
    SMTBX_ASSERT(s5->is_variable());
    SMTBX_ASSERT(s1->is_variable());
    SMTBX_ASSERT(s2->is_variable());
    SMTBX_ASSERT(!s6->is_variable());
    SMTBX_ASSERT(!s3->is_variable());
    SMTBX_ASSERT(!s4->is_variable());
  }
};


class test_case_1 : public test_case
{
public:
  void run() {
    std::cout << "sizeof(parameter) = " << sizeof(parameter) << "\n";

    //*** check root and variability status are correctly assigned ***
    std::vector<crystallographic_parameter *> cryst;
    cryst.push_back(s8);
    cryst.push_back(s9);
    reparametrisation
    *reparam = new reparametrisation(uc, boost::make_iterator_range(cryst));
    reparametrisation::iterator p = reparam->parameters().begin();

    // Expected from DFS
    SMTBX_ASSERT(*p++ == s1);
    SMTBX_ASSERT(*p++ == s2);
    SMTBX_ASSERT(*p++ == s5);
    SMTBX_ASSERT(*p++ == s7);
    SMTBX_ASSERT(*p++ == s3);
    SMTBX_ASSERT(*p++ == s4);
    SMTBX_ASSERT(*p++ == s6);
    SMTBX_ASSERT(*p++ == s8);
    SMTBX_ASSERT(*p++ == s11);
    SMTBX_ASSERT(*p++ == s10);
    SMTBX_ASSERT(*p++ == s9);

    check_parameter_status();

    //*** Independent site variability ***
    SMTBX_ASSERT(s1->is_variable());
    SMTBX_ASSERT(s2->is_variable());
    SMTBX_ASSERT(!s3->is_variable());
    SMTBX_ASSERT(!s4->is_variable());

    is3->set_variable(true);
    SMTBX_ASSERT(sc[3].flags.grad_site());
    sc[4].flags.set_grad_site(true);
    SMTBX_ASSERT(s4->is_variable());

    //*** Shifts ***
    sc[3].flags.set_grad_site(false);
    sc[4].flags.set_grad_site(false);
    double s11_old = s11->value;
    af::shared<double> s(7, 1.);
    is1->evaluate(uc); is2->evaluate(uc); is3->evaluate(uc); is4->evaluate(uc);
    reparam->apply_shifts(s.const_ref());
    SMTBX_ASSERT(site_approx_equal(s1->value, (sc[1].site + frac_t(1., 1., 1.))));
    SMTBX_ASSERT(site_approx_equal(s2->value, (sc[2].site + frac_t(1., 1., 1.))));
    SMTBX_ASSERT(site_approx_equal(s3->value, (sc[3].site)));
    SMTBX_ASSERT(site_approx_equal(s4->value, (sc[4].site)));
    SMTBX_ASSERT(std::abs(s11->value - (s11_old + 1.)) < 1.e-15);

    // Clean-up
    delete reparam;
  }
};


class test_case_2 : public test_case
{
public:
  reparametrisation *reparam;

  void check_linearisation(bool s3_s4_are_variable) {
    reparam->linearise();
    SMTBX_ASSERT(site_approx_equal(s8->value, (frac_t(0.2, 0.5, 0.8))));
    SMTBX_ASSERT(site_approx_equal(s9->value, (frac_t(0.2, -0.3, -0.8))));

    reparam->store();
    SMTBX_ASSERT(site_approx_equal(sc[1].site, (frac_t(0.1, 0.2, 0.3))));
    SMTBX_ASSERT(site_approx_equal(sc[2].site, (frac_t(-0.1, -0.2, -0.3))));
    SMTBX_ASSERT(site_approx_equal(sc[3].site, (frac_t(0.4,  0.3,  0.2))));
    SMTBX_ASSERT(site_approx_equal(sc[4].site, (frac_t(-0.5, -0.6, -0.7))));
    SMTBX_ASSERT(site_approx_equal(sc[5].site, (frac_t(0, 0, 0))));
    SMTBX_ASSERT(site_approx_equal(sc[6].site, (frac_t(-0.1, -0.3, -0.5))));
    SMTBX_ASSERT(site_approx_equal(sc[7].site, (frac_t(0.1, 0.2, 0.3))));
    SMTBX_ASSERT(site_approx_equal(sc[8].site, (frac_t(0.2, 0.5, 0.8))));

    /* [ jt(i,j) ]_ij = [ dx_j / dx_i ]_ij
     where i runs through independent parameters, and
     where j runs through all parameters
     */
    sparse_matrix_type &jt = reparam->jacobian_transpose;

    // Construct reference Jacobian transpose
    int m = s3_s4_are_variable ? 4*3 + 1 : 2*3 + 1;
    int n = 10*3 + 1;
    sparse_matrix_type jt0(m, n);

    // columns corresponding to independent parameters
    for (int k=0; k<3; ++k) {
      jt0(s1->index() + k, s1->index() + k) = 1.;
      jt0(s2->index() + k, s2->index() + k) = 1.;
      if (s3_s4_are_variable) {
        jt0(s3->index() + k, s3->index() + k) = 1.;
        jt0(s4->index() + k, s4->index() + k) = 1.;
      }
      jt0(s11->index()   , s11->index())    = 1.;
    }

    // s6 = s3 + s4
    if (s3_s4_are_variable) {
      for (int k=0; k<3; ++k) {
        jt0(s1->index() + k, s6->index() + k) = 0.;
        jt0(s2->index() + k, s6->index() + k) = 0.;
        jt0(s3->index() + k, s6->index() + k) = 1.;
        jt0(s4->index() + k, s6->index() + k) = 1.;
        jt0(s11->index()   , s6->index())     = 0.;
      }
    }

    // s10 = s11*s6 = s11*s3 + s11*s4
    for (int k=0; k<3; ++k) {
      jt0(s1->index() + k, s10->index() + k) = 0.;
      jt0(s2->index() + k, s10->index() + k) = 0.;
      if (s3_s4_are_variable) {
        jt0(s3->index() + k, s10->index() + k) = s11->value;
        jt0(s4->index() + k, s10->index() + k) = s11->value;
      }
      jt0(s11->index()   , s10->index() + k) = s6->value[k];
    }

    // s9 = s10 + s3 = (s11 + 1)*s3 + s11*s4
    for (int k=0; k<3; ++k) {
      jt0(s1->index() + k, s9->index() + k) = 0.;
      jt0(s2->index() + k, s9->index() + k) = 0.;
      if (s3_s4_are_variable) {
        jt0(s3->index() + k, s9->index() + k) = s11->value + 1.;
        jt0(s4->index() + k, s9->index() + k) = s11->value;
      }
      jt0(s11->index()   , s9->index() + k) = s6->value[k];
    }

    // s5 = s1 + s2
    for (int k=0; k<3; ++k) {
      jt0(s1->index() + k, s5->index() + k) = 1.;
      jt0(s2->index() + k, s5->index() + k) = 1.;
      if (s3_s4_are_variable) {
        jt0(s3->index() + k, s5->index() + k) = 0.;
        jt0(s4->index() + k, s5->index() + k) = 0.;
      }
      jt0(s11->index()   , s5->index())     = 0.;
    }

    // s7 = s5 - s2 = s1 (cancellation of s2)
    for (int k=0; k<3; ++k) {
      jt0(s1->index() + k, s7->index() + k) = 1.;
      jt0(s2->index() + k, s7->index() + k) = 0.;
      if (s3_s4_are_variable) {
        jt0(s3->index() + k, s7->index() + k) = 0.;
        jt0(s4->index() + k, s7->index() + k) = 0.;
      }
      jt0(s11->index()   , s7->index())     = 0.;
    }

    // s8 = s7 - s6 = s1 - s3 - s4
    for (int k=0; k<3; ++k) {
      jt0(s1->index() + k, s8->index() + k) =  1.;
      jt0(s2->index() + k, s8->index() + k) =  0.;
      if (s3_s4_are_variable) {
        jt0(s3->index() + k, s8->index() + k) = -1.;
        jt0(s4->index() + k, s8->index() + k) = -1.;
      }
      jt0(s11->index()   , s8->index())     =  0.;
    }

    scitbx::sparse::approx_equal<double> sparse_approx_equal(1.e-15);
    #if 0
      std::cout << scitbx::sparse::dense_display(jt) << std::endl;
      std::cout << scitbx::sparse::dense_display(jt0) << std::endl;
    #endif
    SMTBX_ASSERT(sparse_approx_equal(jt, jt0));
  }

  void run() {
    // Construct
    std::vector<crystallographic_parameter *> cryst;
    cryst.push_back(s6);
    cryst.push_back(s1);
    cryst.push_back(s2);
    cryst.push_back(s8);
    cryst.push_back(s4);
    cryst.push_back(s7);
    cryst.push_back(s9);

    reparam = new reparametrisation(uc, boost::make_iterator_range(cryst));

    // Test
    check_linearisation(/*s3 and s4 are variable:*/ false);

    sc[3].flags.set_grad_site(true);
    sc[4].flags.set_grad_site(true);
    reparam->analyse_variability();
    check_linearisation(/*s3 and s4 are variable:*/ true);

    // Clean-up
    delete reparam;
  }

};



/// ** Don't forget to run with valgrind to check leaks **
int main() {
  test_case_1 t1;
  t1.run();
  test_case_2 t2;
  t2.run();
  std::cout << "OK\n";
  return 0;
}
