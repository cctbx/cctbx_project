#include <boost/python.hpp>
#include <rstbx/bpcx/detector_model/sensor.h>

namespace rstbx { namespace detector_model {
    
    struct detector_model_wrappers
    {
        static void wrap()
        {
            using namespace boost::python;
            
            class_<sensor>("sensor", init<const scitbx::vec3<double>&,
                                          const scitbx::vec3<double>&,
                                          const scitbx::vec3<double>&,
                                          const scitbx::vec2<double>&,
                                          const scitbx::vec2<double>&>())
                .add_property("distance", &sensor::get_distance)
                .add_property("origin", &sensor::get_origin)
                .add_property("normal", &sensor::get_normal)
                .add_property("dir1", &sensor::get_dir1)
                .add_property("dir2", &sensor::get_dir2)
                .add_property("lim1", &sensor::get_lim1)
                .add_property("lim2", &sensor::get_lim2)
                .add_property("D", &sensor::get_D)
                .add_property("d", &sensor::get_d)
            ;
        }
    };
    
    void init_module()
    {
        using namespace boost::python;
        detector_model_wrappers::wrap();
    }
    
}} //namespace rstbx::detector_model

BOOST_PYTHON_MODULE(rstbx_bpcx_detector_model_ext)
{
    rstbx::detector_model::init_module();
}
