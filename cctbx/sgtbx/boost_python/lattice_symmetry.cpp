#include <boost/python/class.hpp>
#include <boost/python/def.hpp>
#include <boost/python/args.hpp>
#include <cctbx/sgtbx/lattice_symmetry.h>

namespace cctbx { namespace sgtbx { namespace boost_python {

  uc_mat3
  n_fold_operator_from_axis_direction(
    uc_vec3 const& cart, int const& type, int const& sense){
#define N_FOLD lattice_symmetry::n_fold_operator_from_axis_direction
    switch (type){
      case 1: return N_FOLD<1>(cart,sense);
      case 2: return N_FOLD<2>(cart,sense);
      case 3: return N_FOLD<3>(cart,sense);
      case 4: return N_FOLD<4>(cart,sense);
      case 6: return N_FOLD<6>(cart,sense);
      default: throw error("Rotation type must be one of 1,2,3,4,6");
    }
  };

  void wrap_lattice_symmetry()
  {
    using namespace boost::python;
    def("n_fold_operator_from_axis_direction",
      n_fold_operator_from_axis_direction);

    typedef lattice_symmetry::evaluated_axis_t eat;
    typedef return_value_policy<return_by_value> rbv;
    class_<eat>("evaluated_axis_t",no_init)
      .def("u",&eat::get_u)
      .def("h",&eat::get_h)
      .add_property("t",make_getter(&eat::t, rbv()))
      .add_property("tau",make_getter(&eat::tau, rbv()))
      .def_readonly("delta",&eat::delta)
    ;

    typedef lattice_symmetry::group_search gs;
    class_<gs>("lattice_symmetry_group_search", no_init)
      .def(init<int>((arg_("modulus"))))
      .def("n_potential_axes", &gs::n_potential_axes)
      .def("__call__", &gs::operator(), (
        arg_("niggli_cell"),
        arg_("max_delta"),
        arg_("only_test_generators"),
        arg_("introspection")
        ))
      .add_property("candidates",make_getter(&gs::candidates, rbv()))
    ;

    def("lattice_symmetry_find_max_delta", lattice_symmetry::find_max_delta);
  }

}}} // namespace cctbx::sgtbx::boost_python
