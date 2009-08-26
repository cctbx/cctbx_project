#include <scitbx/array_family/boost_python/flex_fwd.h>
#include <boost/python.hpp>
#include <spotfinder/core_toolbox/hough.h>


using namespace boost::python;

BOOST_PYTHON_MODULE(spotfinder_hough_ext)
{

    typedef return_value_policy<return_by_value> rbv;

    class_<spotfinder::hough>("hough",init<>())
      .def("importData",&spotfinder::hough::importData)
      .def("exportData",&spotfinder::hough::exportData)
      .def("setGeometry",&spotfinder::hough::setGeometry)
      .def("getGeometry",&spotfinder::hough::getGeometry)
      .def("gaussianBlur",&spotfinder::hough::gaussianBlur)
      .def("sobelEdge",&spotfinder::hough::sobelEdge)
      .def("cannyEdge",&spotfinder::hough::cannyEdge)
      .def("findEllipse",&spotfinder::hough::findEllipse)
      .def("getDistance",&spotfinder::hough::getDistance)
      .def("getRings",&spotfinder::hough::getRings)
      ;
}
