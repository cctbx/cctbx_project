#include <cctbx/miller/f_calc_map.h>

#include <boost/python/class.hpp>

namespace cctbx { namespace miller { namespace boost_python {

  template <typename FloatType>
  struct f_calc_map_wrapper
  {
    typedef f_calc_map<FloatType> wt;

    static void wrap(char const *name) {
      using namespace boost::python;
      class_<wt>(name, no_init)
        .def(init<af::const_ref<miller::index<> > const&,
                  af::const_ref<typename wt::complex_type> const &,
                  bool>
             ((arg("indices"),
               arg("f_calc"),
               arg("anomalous_flag"))))
        .add_property("anomalous_flag", &wt::anomalous_flag)
        .def("import", &wt::import)
        .def("__getitem__", &wt::operator[])
        ;
    }
  };

  void wrap_f_calc_map() {
    f_calc_map_wrapper<double>::wrap("f_calc_map");
  }


}}}
