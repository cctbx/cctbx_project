#include <cctbx/boost_python/flex_fwd.h>

#include <boost/python/module.hpp>
#include <boost/python/class.hpp>
#include <boost/python/def.hpp>
#include <boost/python/dict.hpp>
#include <boost/python/list.hpp>
#include <scitbx/array_family/flex_types.h>
#include <scitbx/array_family/shared.h>
#include <scitbx/math/mean_and_variance.h>
#include <scitbx/array_family/boost_python/shared_wrapper.h>
#include <scitbx/examples/bevington/prototype_core.h>
#include <Eigen/Sparse>
#include <boost/python/return_internal_reference.hpp>

using namespace boost::python;
namespace scitbx{
namespace example{
namespace boost_python { namespace {

  void
  scaling_init_module() {
    using namespace boost::python;

    class_<linear_ls_eigen_wrapper>("linear_ls_eigen_wrapper", no_init)
      .def(init<int>(arg("n_parameters")))
      .def("n_parameters", &linear_ls_eigen_wrapper::n_parameters)
      .def("reset", &linear_ls_eigen_wrapper::reset)
      .def("right_hand_side", &linear_ls_eigen_wrapper::right_hand_side)
      .def("solve", &linear_ls_eigen_wrapper::solve)
      .def("solved", &linear_ls_eigen_wrapper::solved)
      .def("normal_matrix", &linear_ls_eigen_wrapper::normal_matrix)
      .def("solution", &linear_ls_eigen_wrapper::solution)
    ;

    typedef non_linear_ls_eigen_wrapper nllsew;
    boost::python::return_internal_reference<> rir;
    class_<nllsew,
           bases<scitbx::lstbx::normal_equations::non_linear_ls<double> > >(
           "non_linear_ls_eigen_wrapper", no_init)
      .def(init<int const&>(arg("n_parameters")))
      .def("reset", &nllsew::reset)
      .def("step_equations",&nllsew::step_equations, rir)
      .def("add_constant_to_diagonal",&nllsew::add_constant_to_diagonal)
      .def("show_eigen_summary",&nllsew::show_eigen_summary)
      .def("get_normal_matrix_diagonal",&nllsew::get_normal_matrix_diagonal)
      .def("get_normal_matrix",&nllsew::get_normal_matrix)
      .def("get_cholesky_lower", &nllsew::get_cholesky_lower)
      .def("get_cholesky_diagonal", &nllsew::get_cholesky_diagonal)
      .def("get_eigen_permutation_ordering", &nllsew::get_eigen_permutation_ordering)
      .def("add_equations", &nllsew::add_equations)
      .def("solved", &nllsew::solved)
      .def("get_normal_matrix_ncols", &nllsew::get_normal_matrix_ncols)
      .def("n_parameters", &nllsew::n_parameters)
      .def("get_normal_matrix_nnonZeros", &nllsew::get_normal_matrix_nnonZeros)
      .def("get_lower_cholesky_nnonZeros", &nllsew::get_lower_cholesky_nnonZeros)
    ;

    typedef bevington_silver silver;
    class_<silver>( "bevington_silver")
      .def(init<>())
      .def("set_cpp_data", &silver::set_cpp_data)
      .def("fvec_callable", &silver::fvec_callable)
      .def("gvec_callable", &silver::gvec_callable)
      .def("curvatures", &silver::curvatures)
      .def("functional", &silver::functional)
      .def("get_xobs", &silver::get_xobs)
      .def("get_wobs", &silver::get_wobs)
    ;

    typedef dense_base_class dbc;
    class_<dbc,bases<silver, scitbx::lstbx::normal_equations::non_linear_ls<double> > >( "dense_base_class", no_init)
      .def(init<int>(arg("n_parameters")))
      .def("access_cpp_build_up_directly_dense", &dbc::access_cpp_build_up_directly_dense,
        (arg("objective_only"),arg("current_values")))
    ;

    typedef eigen_base_class bev;
    class_<bev,bases<silver, nllsew> >( "eigen_base_class", no_init)
      .def(init<int>(arg("n_parameters")))
      .def("access_cpp_build_up_directly_eigen_eqn", &bev::access_cpp_build_up_directly_eigen_eqn,
        (arg("objective_only"),arg("current_values")))
    ;
  }

}
}}} // namespace xfel::boost_python::<anonymous>

BOOST_PYTHON_MODULE(scitbx_examples_bevington_ext)
{
  scitbx::example::boost_python::scaling_init_module();

}
