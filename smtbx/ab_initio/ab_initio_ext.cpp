#include <smtbx/boost_python/flex_fwd.h>

#include <boost/python/module.hpp>
#include <boost/python/class.hpp>
#include <boost/python/def.hpp>
#include <boost/python/args.hpp>

#include <scitbx/array_family/accessors/c_grid_padded.h>
#include <scitbx/array_family/accessors/flex_grid.h>

#include <smtbx/ab_initio/density_modification.h>

namespace smtbx { namespace ab_initio { namespace boost_python {

  namespace af = scitbx::af;

  template<class FloatType, class AccessorType>
  struct density_modification_wrapper
  {
    typedef void (*f_t)(af::ref<FloatType, AccessorType>, FloatType);
    static void wrap() {
      using namespace boost::python;
      using namespace density_modification;
      def("flip_charges_in_place",
          static_cast<f_t>(flip_charges_in_place));
      def("low_density_elimination_in_place_tanaka_et_al_2001",
          static_cast<f_t>(low_density_elimination_in_place_tanaka_et_al_2001));
    }
  };

  void init_module() {
    density_modification_wrapper<double, af::c_grid_padded<3> >::wrap();
    density_modification_wrapper<double, af::flex_grid<> >::wrap();
  }

}}} // namespace smtbx::ab_initio::boost_python


BOOST_PYTHON_MODULE(smtbx_ab_initio_ext)
{
  smtbx::ab_initio::boost_python::init_module();
}
