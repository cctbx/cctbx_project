#include <boost/python/class.hpp>
#include <boost/python/def.hpp>
#include <boost/python/args.hpp>
#include <cctbx/sgtbx/lattice_symmetry.h>

namespace cctbx { namespace sgtbx { namespace boost_python {

  uc_mat3
  N_fold_operator_from_axis_direction(
    uc_vec3 const& cart, int const& type, int const& sense){
#define N_fold lattice_symmetry::N_fold_operator_from_axis_direction
    switch (type){
      case 1: return N_fold<1>(cart,sense);
      case 2: return N_fold<2>(cart,sense);
      case 3: return N_fold<3>(cart,sense);
      case 4: return N_fold<4>(cart,sense);
      case 6: return N_fold<6>(cart,sense);
      default: throw error("Rotation type must be one of 1,2,3,4,6");
    }
  };

  void wrap_lattice_symmetry()
  {
    typedef lattice_symmetry::group_search gs;
    using namespace boost::python;
    def("N_fold_operator_from_axis_direction",
      N_fold_operator_from_axis_direction);

    class_<gs>("lattice_symmetry_group_search", no_init)
      .def(init<int>((arg_("modulus"))))
      .def("n_potential_axes", &gs::n_potential_axes)
      .def("__call__", &gs::operator(), (
        arg_("niggli_cell"),
        arg_("max_delta"),
        arg_("only_test_generators")))
    ;

    def("lattice_symmetry_find_max_delta", lattice_symmetry::find_max_delta);
  }

}}} // namespace cctbx::sgtbx::boost_python
