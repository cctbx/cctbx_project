#include <boost/python/module.hpp>
#include <boost/python/class.hpp>
#include <boost/python/def.hpp>
#include <boost/python/args.hpp>

#include <scitbx/array_family/accessors/c_grid_padded.h>
#include <scitbx/array_family/accessors/flex_grid.h>
#include <smtbx/ab_initio/charge_flipping.h>

namespace smtbx { namespace ab_initio { namespace boost_python {                            

  template<class FloatType, class AccessorType>
  struct charge_flipping_wrapper 
  {
    typedef void (*f_t)(af::ref<FloatType, AccessorType>, FloatType);
    static void wrap() {
      using namespace boost::python;
      def("flip_charges_in_place", static_cast<f_t>(flip_charges_in_place));
    }
  };
  
  namespace {
    
    void init_module() {
      charge_flipping_wrapper<double, af::c_grid_padded<3> >::wrap();
      charge_flipping_wrapper<double, af::flex_grid<> >::wrap();
    }
    
  } // namespace anonymous
}}} // namespace smtbx::ab_initio::boost_python


BOOST_PYTHON_MODULE(smtbx_ab_initio_ext)
{
  smtbx::ab_initio::boost_python::init_module();
}
