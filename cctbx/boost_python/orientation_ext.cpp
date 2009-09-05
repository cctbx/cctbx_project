#include <cctbx/boost_python/flex_fwd.h>

#include <boost/python/tuple.hpp>
#include <boost/python/module.hpp>
#include <boost/python/def.hpp>
#include <boost/python/class.hpp>
#include <boost/python/args.hpp>
#include <boost/python/return_value_policy.hpp>
#include <boost/python/return_internal_reference.hpp>
#include <cctbx/crystal_orientation.h>

namespace cctbx { namespace boost_python { namespace {

  struct crystal_orientation_wrappers
  {
    typedef crystal_orientation w_t;

    static void
    wrap()
    {
      using namespace boost::python;
      class_<w_t>("crystal_orientation", init<w_t>())
        .enable_pickling()
        .def(init<cctbx::oc_mat3 const&, bool const&>())
        .def("unit_cell", &w_t::unit_cell)
        .def("unit_cell_inverse", &w_t::unit_cell_inverse)
        .def("direct_matrix", &w_t::direct_matrix)
        .def("reciprocal_matrix", &w_t::reciprocal_matrix)
        .def("change_basis",
          (w_t (w_t::*)(
            sgtbx::change_of_basis_op const&) const)
              &w_t::change_basis)
        .def("change_basis",
          (w_t (w_t::*)(oc_mat3 const&) const)
            &w_t::change_basis)
        .def("rotate_thru",
          (w_t (w_t::*)(oc_vec3 const&, double const&) const)
            &w_t::rotate_thru, (
              arg("unit_axis"), arg("angle")))
        .def("direct_mean_square_difference",
          &w_t::direct_mean_square_difference, (
            arg("other")))
        .def("difference_Z_score",
          &w_t::difference_Z_score, (
            arg("other")))
        .def("best_similarity_transformation",
          &w_t::best_similarity_transformation, (
            arg("other"),
            arg("fractional_length_tolerance"),
            arg("unimodular_generator_range")))
      ;
    }
  };

  void init_module()
  {
    crystal_orientation_wrappers::wrap();
  }

}}} // namespace cctbx::boost_python::orientation_ext

BOOST_PYTHON_MODULE(cctbx_orientation_ext)
{
  cctbx::boost_python::init_module();
}
