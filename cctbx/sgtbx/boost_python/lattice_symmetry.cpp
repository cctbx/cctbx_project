#include <boost/python/class.hpp>
#include <boost/python/def.hpp>
#include <boost/python/overloads.hpp>
#include <boost/python/args.hpp>
#include <cctbx/sgtbx/lattice_symmetry.h>

namespace cctbx { namespace sgtbx { namespace boost_python {

  BOOST_PYTHON_FUNCTION_OVERLOADS(
    lattice_symmetry_group_overloads, lattice_symmetry::group, 1, 3)

  void wrap_lattice_symmetry()
  {
    using namespace boost::python;
    def("lattice_symmetry_n_fold_operator_from_axis_direction",
      lattice_symmetry::n_fold_operator_from_axis_direction, (
        arg_("cart"), arg_("n"), arg_("sense")=1));

    def("lattice_symmetry_find_max_delta", lattice_symmetry::find_max_delta, (
      arg_("reduced_cell"), arg_("space_group")));

    def("lattice_symmetry_group",
      lattice_symmetry::group,
        lattice_symmetry_group_overloads((
          arg_("reduced_cell"),
          arg_("max_delta")=3.,
          arg_("enforce_max_delta_for_generated_two_folds")=true)));
  }

}}} // namespace cctbx::sgtbx::boost_python
