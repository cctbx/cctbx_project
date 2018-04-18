#include <iostream>
#include "StrumpackSparseSolver.hpp"

#include "StrumpackSparseSolverMPIDist.hpp"
#include <getopt.h>
#include <mpi4py/mpi4py.h>

#include <boost/python/enum.hpp>
#include <boost/python.hpp>
#include <boost/python/def.hpp>
#include <scitbx/array_family/versa.h>
#include <scitbx/array_family/shared.h>

using namespace boost::python;
using namespace strumpack;
typedef double numType;
typedef int intType;

namespace sparse_solver {
  struct strumpack_mpi_dist_solver{

    /**
      Constructor for the STRUMPACK distributed MPI Ax=b sparse solver strumpack_mpi_dist_solver.
      This object will be callable from python upon importing the appropriate extension module using
      'from scitbx.examples.bevington import sparse_solver as ss'
      Result vector x_res will be initialised by constructor to the same size as n_rows to contain solution.

      @param[in]  n_rows_local Rank-local row count of the A matrix.
      @param[in]  n_cols Column count of the A matrix.
      @param[in]  A_row_offset_local Rank-local row offset indices for CSR sparse format matrix.
      @param[in]  A_col_idx_local Rank-local column indices for CSR sparse format matrix.
      @param[in]  A_values_local Rank-local A matrix values in CSR sparse format matrix.
      @param[in]  b_local Rank-local right hand side vector of linear system.
      @return Description of returned value.
    */
    strumpack_mpi_dist_solver(
              const int n_rows_local, const int n_cols,
              boost::python::object py_comm,
              const scitbx::af::shared<intType> A_row_offset_local,
              const scitbx::af::shared<intType> A_col_idx_local,
              const scitbx::af::shared<numType> A_values_local,
              const scitbx::af::shared<numType> b_local,
              const scitbx::af::shared<intType> len_rows_local,
              const int reorderEnum,
              const int solverEnum) : x_strum_local(n_rows_local, 0.){

      int thread_level, rank, size;
/*      //Pass MPI communicator from mpi4py to C++
      PyObject* py_obj = py_comm.ptr();
std::cout <<"CPP Here 1a"<< py_obj <<"  " << py_comm.ptr() << std::endl;
      MPI_Comm *comm_p = PyMPIComm_Get(py_comm.ptr());
std::cout <<"CPP Here 1b" << std::endl;
      //Create a duplicate communicator for extension work; may not be necessary
*/    MPI_Comm newcomm;
      MPI_Comm_dup(MPI_COMM_WORLD, &newcomm);

      MPI_Comm_size(newcomm, &size);
      MPI_Comm_rank(newcomm, &rank);

      //Check the initialisation of the MPI environment is using MPI_THREAD_FUNNELED
      int flag;
      MPI_Initialized(&flag);
      if (!flag){
        std::cout << "Library not initialised." << std::endl;
        MPI_Init_thread(NULL,NULL, MPI_THREAD_FUNNELED, &thread_level);
        //atexit(library_onexit);
      }

      {
        std::size_t nnz_local = A_values_local.size();
        std::size_t gridSize_local = n_rows_local * n_cols;

        //Taking raw pointer to flex data types
        const intType* r_ptr = A_row_offset_local.begin();
        const intType* c_ptr = A_col_idx_local.begin();
        const numType* val_ptr = A_values_local.begin();

        StrumpackSparseSolverMPIDist <numType,intType> spss(newcomm,true);

        spss.set_distributed_csr_matrix(n_rows_local, r_ptr, c_ptr, val_ptr, len_rows_local.begin());

        spss.options().set_reordering_method( ReorderingStrategy(reorderEnum) );
        spss.options().set_Krylov_solver( KrylovSolver(solverEnum) );

        spss.reorder();
        spss.solve(b_local.begin(), x_strum_local.begin());
      }
      MPI_Barrier(newcomm);
      MPI_Comm_free( &newcomm );
    }
    scitbx::af::shared<numType> x_strum_local; //Resulting x will be stored here
  }; //end of struct sstrumpack_mpi_dist_solver

namespace boost_python{
namespace {
  void export_strumpack_mpi_dist_solver() {
    typedef return_value_policy<return_by_value> rbv;

    enum_<KrylovSolver>("strumpack_mpi_dist_solver")
     .value("auto", KrylovSolver::AUTO)
     .value("direct", KrylovSolver::DIRECT)
     .value("refine", KrylovSolver::REFINE)
     .value("prec_gmres", KrylovSolver::PREC_GMRES)
     .value("gmres", KrylovSolver::GMRES)
     .value("prec_bicgstab", KrylovSolver::PREC_BICGSTAB)
     .value("bicgstab", KrylovSolver::BICGSTAB)
     .export_values()
    ;
    enum_<ReorderingStrategy>("strumpack_mpi_dist_solver")
     .value("natural", ReorderingStrategy::NATURAL)
     .value("metis", ReorderingStrategy::METIS)
     .value("parmetis", ReorderingStrategy::PARMETIS)
     .value("scotch", ReorderingStrategy::SCOTCH)
     .value("ptscotch", ReorderingStrategy::PTSCOTCH)
     .value("rcm", ReorderingStrategy::RCM)
     .value("geometric", ReorderingStrategy::GEOMETRIC)
     .export_values()
    ;

    class_<sparse_solver::strumpack_mpi_dist_solver> ("strumpack_mpi_dist_solver", no_init)
      .def( init< int, int,
                  boost::python::object,
                  scitbx::af::shared<intType>,
                  scitbx::af::shared<intType>,
                  scitbx::af::shared<numType>,
                  scitbx::af::shared<numType>,
                  scitbx::af::shared<intType>,
                  int, int
                >
      ((
        boost::python::arg("n_rows_local"),
        boost::python::arg("n_cols"),
        boost::python::arg("mpi_comm"),
        boost::python::arg("A_row_offset_local"),
        boost::python::arg("A_col_idx_local"),
        boost::python::arg("A_values_local"),
        boost::python::arg("b_local"),
        boost::python::arg("len_rows_local"),
        boost::python::arg("reorderEnum")=ReorderingStrategy::SCOTCH,
        boost::python::arg("solverEnum")=KrylovSolver::AUTO
      ))
      )
      .add_property("x",make_getter(&sparse_solver::strumpack_mpi_dist_solver::x_strum_local, rbv()))
    ;
  }
}//anon
}//boost_python
}//sparse_solver
BOOST_PYTHON_MODULE(scitbx_examples_strumpack_mpi_dist_solver_ext)
{
   sparse_solver::boost_python::export_strumpack_mpi_dist_solver();
}
