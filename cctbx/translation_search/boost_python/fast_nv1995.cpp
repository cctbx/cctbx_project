#include <cctbx/boost_python/flex_fwd.h>

#include <cctbx/translation_search/fast_nv1995.h>
#include <boost/python/class.hpp>
#include <boost/python/args.hpp>

namespace cctbx { namespace translation_search { namespace boost_python {

namespace {

  struct fast_nv1995_wrappers
  {
    typedef fast_nv1995<> w_t;

    static void
    wrap()
    {
      using namespace boost::python;
      typedef boost::python::arg arg_; // gcc 2.96 workaround
      class_<w_t>("fast_nv1995", no_init)
        .def(init<af::int3 const&,
                  sgtbx::space_group const&,
                  bool,
                  af::const_ref<miller::index<> > const&,
                  af::const_ref<double> const&,
                  af::const_ref<std::complex<double> > const&,
                  af::const_ref<miller::index<> > const&,
                  af::const_ref<std::complex<double> > >(
          (arg_("gridding"),
           arg_("space_group"),
           arg_("anomalous_flag"),
           arg_("miller_indices_f_obs"),
           arg_("f_obs"),
           arg_("f_part"),
           arg_("miller_indices_p1_f_calc"),
           arg_("p1_f_calc"))))
        .def("target_map", &w_t::target_map)
      ;
    }
  };

} // namespace <anoymous>

  void wrap_fast_nv1995()
  {
    fast_nv1995_wrappers::wrap();
  }

}}} // namespace cctbx::translation_search::boost_python
