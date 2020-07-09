#include <cctbx/boost_python/flex_fwd.h>

#include <boost/python/module.hpp>
#include <boost/python/class.hpp>
#include <boost/python/def.hpp>
#include <boost/python/args.hpp>
#include <mmtbx/bulk_solvent/mosaic.h>
#include <mmtbx/f_model/f_model.h>
#include <mmtbx/bulk_solvent/bulk_solvent.h>

#include <boost/python/return_value_policy.hpp>
#include <boost/python/return_by_value.hpp>
#include <boost/python/list.hpp>

namespace mmtbx { namespace mosaic {
namespace {

  void init_module()
  {
    using namespace boost::python;
    using boost::python::arg;
    typedef return_value_policy<return_by_value> rbv;


   def("alg4",
      (af::shared<double>(*)
        (boost::python::list const& F,
         af::const_ref<std::complex<double> > const& f_obs,
         double                const& k
         //af::const_ref<bool>   const& selection,
         //af::shared<double>           data
         )) alg4);
   ;


  }

} // namespace <anonymous>
}} // namespace mmtbx::mosaic

BOOST_PYTHON_MODULE(mmtbx_mosaic_ext)
{
  mmtbx::mosaic::init_module();
}
