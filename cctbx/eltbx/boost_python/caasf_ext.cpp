#include <boost/python/module.hpp>
#include <boost/python/class.hpp>
#include <boost/python/return_value_policy.hpp>
#include <boost/python/copy_const_reference.hpp>
#include <cctbx/eltbx/caasf.h>
#include <scitbx/boost_python/iterator_wrappers.h>

namespace cctbx { namespace eltbx { namespace caasf { namespace boost_python {

namespace {

  template <std::size_t N>
  struct base_wrappers
  {
    typedef base<N> w_t;

    static af::small<float, custom::max_n_ab>
    a(w_t const& o)
    {
      af::small<float, custom::max_n_ab> result;
      for(std::size_t i=0;i<N;i++) result.push_back(o.a(i));
      return result;
    }

    static af::small<float, custom::max_n_ab>
    b(w_t const& o)
    {
      af::small<float, custom::max_n_ab> result;
      for(std::size_t i=0;i<N;i++) result.push_back(o.b(i));
      return result;
    }

    static void
    wrap(const char* python_name)
    {
      using namespace boost::python;
      class_<w_t>(python_name, no_init)
        .def("table", &w_t::table)
        .def("label", &w_t::label)
        .def("a", a)
        .def("b", b)
        .def("c", &w_t::c)
        .def("at_stol_sq", &w_t::at_stol_sq)
        .def("at_stol", &w_t::at_stol)
        .def("at_d_star_sq", &w_t::at_d_star_sq)
      ;
    }
  };

  struct it1992_wrappers
  {
    typedef it1992 w_t;

    static void
    wrap()
    {
      using namespace boost::python;
      base_wrappers<4>::wrap("base_4");
      class_<w_t, bases<base<4> > >("it1992", no_init)
        .def(init<std::string const&, optional<bool> >())
      ;
    }
  };

  struct wk1995_wrappers
  {
    typedef wk1995 w_t;

    static void
    wrap()
    {
      using namespace boost::python;
      base_wrappers<5>::wrap("base_5");
      class_<w_t, bases<base<5> > >("wk1995", no_init)
        .def(init<std::string const&, optional<bool> >())
      ;
    }
  };

  struct custom_wrappers
  {
    typedef custom w_t;

    static void
    wrap()
    {
      using namespace boost::python;
      typedef return_value_policy<copy_const_reference> ccr;
      class_<w_t>("custom", no_init)
        .def(init<std::string const&,
                  af::small<float, w_t::max_n_ab> const&,
                  af::small<float, w_t::max_n_ab> const&,
                  float>())
        .def("label", &w_t::label, ccr())
        .def("n_ab", &w_t::n_ab)
        .def("a", (af::small<float, w_t::max_n_ab> const&(w_t::*)()const)
          &w_t::a, ccr())
        .def("b", (af::small<float, w_t::max_n_ab> const&(w_t::*)()const)
          &w_t::b, ccr())
        .def("c", &w_t::c)
        .def("at_stol_sq", &w_t::at_stol_sq)
        .def("at_stol", &w_t::at_stol)
        .def("at_d_star_sq", &w_t::at_d_star_sq)
      ;
    }
  };

  void init_module()
  {
    using namespace boost::python;

    it1992_wrappers::wrap();
    scitbx::boost_python::iterator_wrappers<
      it1992, it1992_iterator>::wrap("it1992_iterator");

    wk1995_wrappers::wrap();
    scitbx::boost_python::iterator_wrappers<
      wk1995, wk1995_iterator>::wrap("wk1995_iterator");

    custom_wrappers::wrap();
  }

} // namespace <anonymous>
}}}} // namespace cctbx::eltbx::caasf::boost_python

BOOST_PYTHON_MODULE(caasf_ext)
{
  cctbx::eltbx::caasf::boost_python::init_module();
}
