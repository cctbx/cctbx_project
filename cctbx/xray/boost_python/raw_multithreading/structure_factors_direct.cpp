#include <cctbx/boost_python/flex_fwd.h>

#include <cctbx/math/cos_sin_table.h>
#include <cctbx/xray/raw_multithreading/structure_factors_direct.h>
#include <boost/python/class.hpp>
#include <boost/python/return_value_policy.hpp>
#include <boost/python/copy_const_reference.hpp>

namespace cctbx { namespace xray { namespace structure_factors {
  namespace boost_python {

    namespace {

      struct raw_multithreaded_direct_wrapper
      {
        typedef raw_multithreaded_direct<> w_t;
        typedef w_t::scatterer_type scatterer_type;
        typedef w_t::float_type float_type;

        static void
        wrap()
        {
          using namespace boost::python;
          typedef return_value_policy<copy_const_reference> ccr;
          class_<w_t>("structure_factors_raw_multithreaded_direct", no_init)
          .def(init<math::cos_sin_table<float_type> const&,
                    uctbx::unit_cell const&,
                    sgtbx::space_group const&,
                    af::const_ref<miller::index<> > const&,
                    af::const_ref<scatterer_type> const&,
                    scattering_type_registry const&,
                    unsigned>((arg("cos_sin_table"),
                               arg("unit_cell"),
                               arg("space_group"),
                               arg("miller_indices"),
                               arg("scatterers"),
                               arg("scattering_type_registry"),
                               arg("n_threads"))))
          .def(init<uctbx::unit_cell const&,
                    sgtbx::space_group const&,
                    af::const_ref<miller::index<> > const&,
                    af::const_ref<scatterer_type> const&,
                    scattering_type_registry const&,
                    unsigned>((arg("unit_cell"),
                               arg("space_group"),
                               arg("miller_indices"),
                               arg("scatterers"),
                               arg("scattering_type_registry"),
                               arg("n_threads"))))
          .def("f_calc", &w_t::f_calc, ccr())
          ;
        }
      };

    } // namespace <anoymous>

  }} // namespace structure_factors::boost_python

  namespace boost_python {

    void wrap_structure_factors_raw_multithreaded_direct()
    {
      structure_factors::boost_python::raw_multithreaded_direct_wrapper::wrap();
    }

  }}} // namespace cctbx::xray::boost_python
