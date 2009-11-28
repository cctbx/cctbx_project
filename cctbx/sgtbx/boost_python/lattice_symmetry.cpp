#include <boost/python/def.hpp>
#include <boost/python/args.hpp>
#include <cctbx/sgtbx/lattice_symmetry.h>

namespace cctbx { namespace sgtbx { namespace boost_python {

  void wrap_lattice_symmetry()
  {
    using namespace boost::python;

    def("lattice_symmetry_find_max_delta", lattice_symmetry::find_max_delta, (
      arg("reduced_cell"), arg("space_group")));

    def("lattice_symmetry_group", lattice_symmetry::group, (
      arg("reduced_cell"),
      arg("max_delta")=3.,
      arg("enforce_max_delta_for_generated_two_folds")=true));
  }

}}} // namespace cctbx::sgtbx::boost_python
