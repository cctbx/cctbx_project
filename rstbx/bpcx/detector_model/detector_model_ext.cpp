#include <boost/python.hpp>
#include <rstbx/bpcx/detector_model/sensor.h>
#include <boost/python/overloads.hpp>

namespace rstbx { namespace detector_model {

    using namespace boost::python;

    struct detector_model_wrappers
    {
        typedef scitbx::vec3<double> v3;
        typedef scitbx::vec2<double> v2;

        static void wrap()
        {
            class_<sensor>("sensor", init<const v3&,
                                          const v3&,
                                          const v3&,
                                          const v2&,
                                          const v2&>())
                .add_property("distance", &sensor::get_distance)
                .add_property("origin", &sensor::get_origin)
                .add_property("normal", &sensor::get_normal)
                .add_property("dir1", &sensor::get_dir1)
                .add_property("dir2", &sensor::get_dir2)
                .add_property("lim1", &sensor::get_lim1)
                .add_property("lim2", &sensor::get_lim2)
                .add_property("D", &sensor::get_D)
                .add_property("d", &sensor::get_d)
                .def("set_frame", &sensor::set_frame)
                .def("set_origin", &sensor::set_origin);
        }
    };

    void init_module()
    {
        detector_model_wrappers::wrap();
    }

}} //namespace rstbx::detector_model

BOOST_PYTHON_MODULE(rstbx_bpcx_detector_model_ext)
{
    rstbx::detector_model::init_module();
}
