#include <cctbx/symmetry_search/translation_refinement.h>

#include <boost/python/module.hpp>
#include <boost/python/class.hpp>
#include <boost/python/return_value_policy.hpp>
#include <boost/python/return_by_value.hpp>


namespace cctbx { namespace symmetry_search { namespace boost_python {
  
  template <typename FloatType>
  struct goodness_of_symmetry_wrapper
  {
    typedef goodness_of_symmetry<FloatType> wt;
    typedef typename wt::real_type real_type;
    typedef typename wt::vector_type vector_type;

    static void wrap(char const *name) {
      using namespace boost::python;
      return_value_policy<return_by_value> rbv;
      class_<wt>(name, no_init)
        .def(init<sgtbx::space_group const &,
                  af::const_ref<miller::index<> > const &,
                  af::const_ref<real_type> const &,
                  miller::f_calc_map<real_type> &,
                  vector_type const &>
             ((arg("space_group"), arg("indices"), arg("f_o"), arg("f_c"),
               arg("x"))))
        .def_readonly("value", &wt::q)
        .add_property("gradient", make_getter(&wt::dq, rbv))
        .def_readonly("scale", &wt::lambda)
        .def_readonly("bias", &wt::mu)
        ;
    }
  };


  void init_module() {
    using namespace boost::python;
    goodness_of_symmetry_wrapper<double>::wrap("goodness_of_symmetry");
  }
  
}}} // boost_python

BOOST_PYTHON_MODULE(cctbx_symmetry_search_ext)
{
  cctbx::symmetry_search::boost_python::init_module();
}
