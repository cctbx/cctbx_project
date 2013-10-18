#include <cctbx/boost_python/flex_fwd.h>

#include <boost/python.hpp>
#include <rstbx/indexing_api/indexing_api.h>


using namespace boost::python;
using namespace rstbx::indexing_api;

BOOST_PYTHON_MODULE(rstbx_indexing_api_ext)
{

   typedef return_value_policy<return_by_value> rbv;
   typedef default_call_policies dcp;

   def("cpp_absence_test",cpp_absence_test);

   class_<dps_extended, bases<rstbx::dps_core> >("dps_extended",init< >())
     .def("getData",&dps_extended::getData)
     .def("setData",&dps_extended::setData)
     .def("refine_direction",&dps_extended::refine_direction,
          (arg("candidate"),arg("current_grid"),
           arg("target_grid")))
   ;

}
