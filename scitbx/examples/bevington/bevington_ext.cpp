#include <cctbx/boost_python/flex_fwd.h>
#include <boost/python/enum.hpp>
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
#define _OPENMP 1
#include <omp.h>

#ifdef _STRUMPACK_
#include "StrumpackSparseSolver.hpp"
#include <scitbx/examples/bevington/prototype_core_strumpack.h>
#endif

using namespace boost::python;
namespace scitbx{
namespace example{
namespace boost_python { namespace {

  void
  scaling_init_module() {
    using namespace boost::python;

    typedef return_value_policy<return_by_value> rbv;
    typedef default_call_policies dcp;

    enum_<linear_ls_eigen_wrapper::SolverAlgo>("linsolver_backend")
      .value("ldlt", linear_ls_eigen_wrapper::LDLT)
      .value("llt", linear_ls_eigen_wrapper::LLT)
      .value("cg", linear_ls_eigen_wrapper::CG)
      .value("bicgstab", linear_ls_eigen_wrapper::BICGSTAB)
      .export_values()
    ;

    class_<linear_ls_eigen_wrapper>("linear_ls_eigen_wrapper", no_init)
      .def(init<int>(arg("n_parameters")))
      .def("n_parameters", &linear_ls_eigen_wrapper::n_parameters)
      .def("reset", &linear_ls_eigen_wrapper::reset)
      .def("right_hand_side", &linear_ls_eigen_wrapper::right_hand_side)
      .def("solve", &linear_ls_eigen_wrapper::solve, ( arg("solverIdx")=linear_ls_eigen_wrapper::LDLT ) )
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

  //#############################

#ifdef _STRUMPACK_

    enum_<strumpack::KrylovSolver>("strumpack_algo")
      .value("automatic", strumpack::KrylovSolver::AUTO)
      .value("direct", strumpack::KrylovSolver::DIRECT)
      .value("refine", strumpack::KrylovSolver::REFINE)
      .value("prec_gmres", strumpack::KrylovSolver::PREC_GMRES)
      .value("gmres", strumpack::KrylovSolver::GMRES)
      .value("prec_bicgstab", strumpack::KrylovSolver::PREC_BICGSTAB)
      .value("bicgstab", strumpack::KrylovSolver::BICGSTAB)
      .export_values()
    ;

    enum_<strumpack::ReorderingStrategy>("strumpack_reorder")
      .value("natural", strumpack::ReorderingStrategy::NATURAL)
      .value("metis", strumpack::ReorderingStrategy::METIS)
      .value("parmetis", strumpack::ReorderingStrategy::PARMETIS)
      .value("scotch", strumpack::ReorderingStrategy::SCOTCH)
      .value("ptscotch", strumpack::ReorderingStrategy::PTSCOTCH)
      .value("rcm", strumpack::ReorderingStrategy::RCM)
      .value("geometric", strumpack::ReorderingStrategy::GEOMETRIC)
      .export_values()
    ;

    class_<linear_ls_strumpack_wrapper>("linear_ls_strumpack_wrapper", no_init)
      .def(init<int>(arg("n_parameters")))
      .def("n_parameters", &linear_ls_strumpack_wrapper::n_parameters)
      .def("reset", &linear_ls_strumpack_wrapper::reset)
      .def("right_hand_side", &linear_ls_strumpack_wrapper::right_hand_side)
      //.def("solve", &linear_ls_strumpack_wrapper::solve)
      .def("solve", &linear_ls_strumpack_wrapper::solve, 
        ( 
          arg("solverAlgo")=strumpack::KrylovSolver::AUTO, 
          arg("solverReordering")=strumpack::ReorderingStrategy::NATURAL,
          arg("verbose")=false,
          arg("enableHSS")=false 
        ) 
      )
      .def("solved", &linear_ls_strumpack_wrapper::solved)
      .def("normal_matrix", &linear_ls_strumpack_wrapper::normal_matrix)
      .def("solution", &linear_ls_strumpack_wrapper::solution)
    ;
    typedef non_linear_ls_strumpack_wrapper nllssw;
    class_<nllssw,
           bases<scitbx::lstbx::normal_equations::non_linear_ls<double> > >(
           "non_linear_ls_strumpack_wrapper", no_init)
      .def(init<int const&>(arg("n_parameters")))
      .def("reset", &nllssw::reset)
      .def("step_equations",&nllssw::step_equations, rir)
      .def("add_constant_to_diagonal",&nllssw::add_constant_to_diagonal)
      .def("show_eigen_summary",&nllssw::show_eigen_summary)
      .def("get_normal_matrix_diagonal",&nllssw::get_normal_matrix_diagonal)
      .def("get_normal_matrix",&nllssw::get_normal_matrix)
      .def("get_cholesky_lower", &nllssw::get_cholesky_lower)
      .def("get_cholesky_diagonal", &nllssw::get_cholesky_diagonal)
      .def("add_equations", &nllssw::add_equations)
      .def("solved", &nllssw::solved)
      .def("get_normal_matrix_ncols", &nllssw::get_normal_matrix_ncols)
      .def("n_parameters", &nllssw::n_parameters)
      .def("get_normal_matrix_nnonZeros", &nllssw::get_normal_matrix_nnonZeros)
      .def("get_lower_cholesky_nnonZeros", &nllssw::get_lower_cholesky_nnonZeros)
    ;

    typedef strumpack_base_class bevs;
    class_<bevs,bases<silver, nllssw> >( "strumpack_base_class", no_init)
      .def(init<int>(arg("n_parameters")))
      .def("access_cpp_build_up_directly_strumpack_eqn", &bevs::access_cpp_build_up_directly_strumpack_eqn,
        (arg("objective_only"),arg("current_values")))
    ;

#endif

  }

}
}}} // namespace xfel::boost_python::<anonymous>

BOOST_PYTHON_MODULE(scitbx_examples_bevington_ext)
{
  scitbx::example::boost_python::scaling_init_module();

}
