#include <cctbx/boost_python/flex_fwd.h>

#include <boost/python/tuple.hpp>
#include <boost/python/module.hpp>
#include <boost/python/def.hpp>
#include <boost/python/class.hpp>
#include <boost/python/args.hpp>
#include <boost/python/overloads.hpp>
#include <boost/python/return_value_policy.hpp>
#include <boost/python/return_internal_reference.hpp>
#include <cctbx/crystal_orientation.h>

namespace cctbx { namespace boost_python { namespace {

  struct crystal_orientation_wrappers: boost::python::pickle_suite {
    BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(
      rotate_thru_overloads, rotate_thru, 2, 2)

    static void
    wrap()
    {
      using namespace boost::python;
      class_<crystal_orientation>("crystal_orientation",
                                  init<crystal_orientation>())
        .enable_pickling()
        .def(init<cctbx::oc_mat3 const&, bool const&>())
        .def("unit_cell",&crystal_orientation::unit_cell)
        .def("unit_cell_inverse",&crystal_orientation::unit_cell_inverse)
        .def("direct_matrix",&crystal_orientation::direct_matrix)
        .def("reciprocal_matrix",&crystal_orientation::reciprocal_matrix)
        .def("change_basis",
          (void (crystal_orientation::*)(cctbx::sgtbx::change_of_basis_op const&))
          &crystal_orientation::change_basis)
        .def("change_basis",
          (void (crystal_orientation::*)(cctbx::oc_mat3 const&))
          &crystal_orientation::change_basis)
        .def("rotate_thru",
          (cctbx::crystal_orientation
          (cctbx::crystal_orientation::*)(cctbx::oc_vec3 const&,double const&) const) 0,
          rotate_thru_overloads((
            arg_("unit_axis"), arg_("angle"))))
        .def("direct_mean_square_difference",
          &crystal_orientation::direct_mean_square_difference,
            (arg_("other")))
        .def("best_similarity_transformation",
          &crystal_orientation::best_similarity_transformation,
            (arg_("other"), arg_("unimodular_generator_range")))
      ;
    }
  };

  void init_module()
  {
    using namespace boost::python;
    crystal_orientation_wrappers::wrap();
  }

}}} // namespace cctbx::boost_python::orientation_ext

BOOST_PYTHON_MODULE(cctbx_orientation_ext)
{
  cctbx::boost_python::init_module();
}
