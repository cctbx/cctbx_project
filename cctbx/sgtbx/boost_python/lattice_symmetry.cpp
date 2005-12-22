#include <boost/python/class.hpp>
#include <boost/python/def.hpp>
#include <boost/python/overloads.hpp>
#include <boost/python/args.hpp>
#include <cctbx/sgtbx/lattice_symmetry.h>

namespace cctbx { namespace sgtbx { namespace boost_python {

  BOOST_PYTHON_FUNCTION_OVERLOADS(
    group_search_fast_overloads, lattice_symmetry::group_search_fast, 1, 3);

  void wrap_lattice_symmetry()
  {
    using namespace boost::python;
    def("lattice_symmetry_n_fold_operator_from_axis_direction",
      lattice_symmetry::n_fold_operator_from_axis_direction, (
        arg_("cart"), arg_("n"), arg_("sense")=1));

    typedef return_value_policy<return_by_value> rbv;

    typedef lattice_symmetry::potential_axis_t pat;
    class_<pat>("lattice_symmetry_potential_axis_t", no_init)
      .add_property("u", make_getter(&pat::u, rbv()))
      .add_property("h", make_getter(&pat::h, rbv()))
      .def_readonly("abs_uh", &pat::abs_uh)
    ;

    typedef lattice_symmetry::group_search gs;
    class_<gs>("lattice_symmetry_group_search", no_init)
      .def(init<int>((arg_("modulus"))))
      .def("n_potential_axes", &gs::n_potential_axes)
      .def("__call__", &gs::operator(), (
        arg_("niggli_cell"),
        arg_("max_delta"),
        arg_("enforce_max_delta_for_generated_two_folds")))
    ;

    def("lattice_symmetry_find_max_delta", lattice_symmetry::find_max_delta);

    def("lattice_symmetry_group_search_fast",
      lattice_symmetry::group_search_fast,
        group_search_fast_overloads((
          arg_("niggli_cell"),
          arg_("max_delta")=3.,
          arg_("enforce_max_delta_for_generated_two_folds")=true)));
  }

}}} // namespace cctbx::sgtbx::boost_python
