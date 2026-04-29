#include <cctbx/symmetry_search/translation_refinement.h>

#include <scitbx/array_family/boost_python/shared_flat_conversions.h>

#include <boost/python/module.hpp>
#include <boost/python/class.hpp>
#include <boost/python/return_value_policy.hpp>
#include <boost/python/return_by_value.hpp>

namespace boost {
  template <>
  struct has_trivial_destructor<scitbx::af::tiny<std::complex<double>, 3> > {
    static const bool value = true;
  };
}

namespace cctbx { namespace symmetry_search { namespace boost_python {

  template <typename FloatType>
  struct symmetrised_shifted_structure_factors_wrapper
  {
    typedef symmetrised_shifted_structure_factors<FloatType> wt;
    typedef typename wt::real_type real_type;
    typedef typename wt::complex_type complex_type;
    typedef typename wt::complex_grad_type complex_grad_type;
    typedef typename wt::vector_type vector_type;

    static void wrap(char const *name) {
      using namespace boost::python;
      return_value_policy<return_by_value> rbv;
      af::boost_python::flat_shared_conversions<complex_grad_type>();
      class_<wt>(name, no_init)
        .def(init<sgtbx::space_group const &,
                  af::const_ref<miller::index<> > const &,
                  miller::f_calc_map<real_type> const &,
                  vector_type const &,
                  bool>
             ((arg("space_group"), arg("indices"), arg("f_c"),
               arg("x"), arg("compute_gradient")=false)))
        .add_property("f_x", make_getter(&wt::f_x, rbv))
        .add_property("grad_f_x", make_getter(&wt::grad_f_x, rbv))
        ;
    }
  };

  template <typename FloatType>
  struct ls_with_scale_and_bias_wrapper
  {
    typedef ls_with_scale_and_bias<FloatType> wt;
    typedef typename wt::real_type real_type;
    typedef typename wt::complex_type complex_type;
    typedef typename wt::complex_grad_type complex_grad_type;

    static void wrap(char const *name) {
      using namespace boost::python;
      return_value_policy<return_by_value> rbv;
      class_<wt>(name, no_init)
        .def(init<af::const_ref<complex_type> const &,
                  af::const_ref<complex_grad_type> const &,
                  af::const_ref<real_type> const &,
                  af::const_ref<real_type> const &>
             ((arg("f_x"), arg("grad_f_x"), arg("f_o_sq"), arg("weight"))))
        .def_readonly("value", &wt::q)
        .def_readonly("correlation", &wt::c)
        .add_property("gradient", make_getter(&wt::grad_q, rbv))
        .def_readonly("scale", &wt::lambda)
        .def_readonly("bias", &wt::mu)
        ;
    }
  };


  void init_module() {
    using namespace boost::python;
    symmetrised_shifted_structure_factors_wrapper<double>::wrap(
      "symmetrised_shifted_structure_factors");
    ls_with_scale_and_bias_wrapper<double>::wrap("ls_with_scale_and_bias");
  }

}}} // boost_python

BOOST_PYTHON_MODULE(cctbx_symmetry_search_ext)
{
  cctbx::symmetry_search::boost_python::init_module();
}
