#include <boost/python/class.hpp>
#include <boost/python/def.hpp>
#include <cctbx/sgtbx/lattice_symmetry.h>

namespace cctbx { namespace sgtbx { namespace boost_python {

  void wrap_lattice_symmetry()
  {
    typedef lattice_symmetry::group_search gs;
    using namespace boost::python;
    class_<gs>("lattice_symmetry_group_search", no_init)
      .def(init<int>())
      .def("n_potential_axes", &gs::n_potential_axes)
      .def("__call__", &gs::operator())
    ;

    def("lattice_symmetry_find_max_delta", lattice_symmetry::find_max_delta);
  }

}}} // namespace cctbx::sgtbx::boost_python
