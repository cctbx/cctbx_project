#include <iostream>
#include <boost/python.hpp>
#include <mpi4py/mpi4py.h>
#include <mpi.h>
#include "StrumpackSparseSolverMPIDist.hpp"
//#include "StrumpackSparseSolver.hpp"
#include <boost/python/def.hpp>
#include <scitbx/array_family/versa.h>
#include <scitbx/array_family/shared.h>
#include <Eigen/Sparse>
#include <Eigen/Core>

using namespace boost::python;
using namespace strumpack;
typedef double numType;
typedef int intType;

namespace sparse_solver {
  struct strumpack_solver_mpi{

    /**
      Constructor for the STRUMPACK Ax=b sparse solver strumpack_solver.
      This object will be callable from python upon importing the appropriate extension module using
      'from scitbx.examples.bevington import sparse_solver as ss'
      Result vector x_res will be initialised by constructor to the same size as n_rows to contain solution.

      @param[in]  n_rows Row count of the A matrix.
      @param[in]  n_cols Column count of the A matrix.
      @param[in]  A_row_offset Row offset indices for CSR sparse format matrix.
      @param[in]  A_col_idx Column indices for CSR sparse format matrix.
      @param[in]  A_values A matrix values in CSR sparse format matrix.
      @param[in]  b Right hand side vector of linear system.
      @return Description of returned value.
    */
    strumpack_solver_mpi( boost::python::object py_comm, const int n_rows, const int n_cols,
              const scitbx::af::shared<intType> A_row_offset,
              const scitbx::af::shared<intType> A_col_idx,
              const scitbx::af::shared<numType> A_values,
              const scitbx::af::shared<numType> b,
              const scitbx::af::shared<intType> len_rows_local) : x_strum(n_rows, 0.){//int argc, char* argv[])

      int thread_level, rank, size;

      //Pass MPI communicator from mpi4py to C++
      PyObject* py_obj = py_comm.ptr();
      MPI_Comm *comm_p = PyMPIComm_Get(py_obj);
      //Create a duplicate communicator for extension work
      MPI_Comm newcomm;
//      MPI_Comm_dup(*comm_p, &newcomm);
      MPI_Comm_dup(MPI_COMM_WORLD, &newcomm);

      MPI_Comm_size(newcomm, &size);
      MPI_Comm_rank(newcomm, &rank);

      std::cout << "CPP::Rank " << rank << " at Pos3a" << std::endl;
      int flag;
      std::cout << "SIZE=" << size << " RANK=" << rank <<std::endl;
      MPI_Initialized(&flag);
      if (!flag){
        std::cout << "Library not initialised." << std::endl;
        MPI_Init_thread(NULL,NULL, MPI_THREAD_FUNNELED, &thread_level);
        //atexit(library_onexit);
      }

      std::cout << "CPP::Rank " << rank << " at Pos3b" << std::endl;

      //if (rank==0){
     { 
        std::size_t nnz = A_values.size();
        std::size_t gridSize = n_rows * n_cols;

        //Taking raw pointer to flex data types
        const intType* r_ptr = A_row_offset.begin();
        const intType* c_ptr = A_col_idx.begin();
        const numType* val_ptr = A_values.begin();
        numType* xx = new numType[n_rows];

        StrumpackSparseSolverMPIDist <double,int> spss(newcomm,false);

        //StrumpackSparseSolver<numType,intType> spss(MPI_COMM_WORLD);
        std::cout << "CPP::Rank " << rank << " at Pos3c" << std::endl;

        //Set the sparse matrix values using the CSR formatted data
        /*for (int ii=0; ii < nnz; ++ii){
          std::cout << "R=" << rank << " A[" << ii+len_rows_local[rank] <<"]=" << val_ptr[ii] << std::endl;
        }*/
 
        spss.set_distributed_csr_matrix(n_rows, r_ptr, c_ptr, val_ptr, len_rows_local.begin());
        std::cout <<  "CPP::Rank " << rank << " at Pos3d" << std::endl;
        //spss.options().set_mc64job(0);
        //spss.options().set_reordering_method(ReorderingStrategy::METIS);

        std::cout <<  "CPP::Rank " << rank << " at Pos3e" << std::endl;
        //spss.reorder();

        std::cout <<  "CPP::Rank " << rank << " at Pos3f" << std::endl;
        spss.factor();

        std::cout <<  "CPP::Rank " << rank << " at Pos3g" << std::endl;

        spss.solve(b.begin(), x_strum.begin());
        //spss.solve(b.begin(), xx);
        std::cout << "CPP::Rank " << rank << " at Pos3h" << std::endl;
      }
      std::cout << "CPP::Rank " << rank << " at Pos3i" << std::endl;
      MPI_Barrier(newcomm);
      std::cout << "CPP::Rank " << rank << " at Pos3j" << std::endl;
      MPI_Comm_free( &newcomm );
      std::cout << "CPP::Rank " << rank << " at Pos3k" << std::endl;
      strumpack::scalapack::Cblacs_exit(1); //
      std::cout << "CPP::Rank " << rank << " at Pos3l" << std::endl;
      //MPI_Finalize();
      std::cout << "CPP::Rank " << rank << " at Pos3m" << std::endl;
    }
    scitbx::af::shared<numType> x_strum; //Resulting x will be stored here
  }; //end of struct sparse_solver

  struct eigen_solver{
    //Eigen solver comparison
    eigen_solver(const int n_rows, const int n_cols,
              const scitbx::af::shared<intType> A_row_offset,
              const scitbx::af::shared<intType> A_col_idx,
              const scitbx::af::shared<numType> A_values,
              const scitbx::af::shared<numType> b) : x_eig(n_rows, 0.){

      Eigen::SparseMatrix<double, Eigen::RowMajor> spMat(n_rows,n_cols);
            std::vector<Eigen::Triplet<double,int> > tList;

      Eigen::VectorXd b_internal(n_cols);
      for (int i = 0; i<n_cols; ++i){
         b_internal[i] = *(b.begin()+i);
      }

      int c_l, r_i;
      for( int row = 0; row < n_cols; ++row ){
        r_i = *(A_row_offset.begin() + row); //Row value at index i
        c_l = ( *( A_row_offset.begin() + 1 + row ) - r_i ); //Column length between the given row offsets
        for( int col = 0; col < c_l; ++col ){
           tList.push_back(
             Eigen::Triplet<double,int>( row, *(A_col_idx.begin() + col + r_i), *(A_values.begin() + col + r_i) )
           );
        }
      }
      spMat.setFromTriplets( tList.begin(), tList.end() ); //Form the Eigen matrix from the given (i,j,value) sparse format
      spMat.makeCompressed();

      Eigen::BiCGSTAB<Eigen::SparseMatrix<double,Eigen::RowMajor> > solver;
      solver.compute(spMat);

      //Eigen::SimplicialLDLT<Eigen::SparseMatrix<double> > chol(spMat);
      //Eigen::VectorXd x = chol.solve(b_internal);
      Eigen::VectorXd x = solver.solve(b_internal);
      double* solnptr = x_eig.begin();
      for (int i = 0; i<n_rows; ++i){
        *solnptr++ = x[i];
      }
    }
    scitbx::af::shared<numType> x_eig; //Resulting x will be stored here
  }; //end of struct sparse_solver

namespace boost_python{
namespace {
  void export_strumpack_solver_mpi() {
    typedef return_value_policy<return_by_value> rbv;
    class_<sparse_solver::strumpack_solver_mpi> ("strumpack_solver_mpi", no_init)
      .def( init< boost::python::object, int, int,
                 scitbx::af::shared<intType>,
                 scitbx::af::shared<intType>,
                 scitbx::af::shared<numType>,
                 scitbx::af::shared<numType>,
                 scitbx::af::shared<intType>
                >
      ((
        boost::python::arg("comm"),
        boost::python::arg("n_rows"),
        boost::python::arg("n_cols"),
        boost::python::arg("A_row_offset"),
        boost::python::arg("A_col_idx"),
        boost::python::arg("A_values"),
        boost::python::arg("b"),
        boost::python::arg("len_rows_local")
      ))
      )
      .add_property("x",make_getter(&sparse_solver::strumpack_solver_mpi::x_strum, rbv()))
    ;


    class_<sparse_solver::eigen_solver> ("eigen_solver", no_init)
      .def( init< int, int,
                 scitbx::af::shared<intType>,
                 scitbx::af::shared<intType>,
                 scitbx::af::shared<numType>,
                 scitbx::af::shared<numType>
                >
      ((
        boost::python::arg("n_rows"),
        boost::python::arg("n_cols"),
        boost::python::arg("A_row_offset"),
        boost::python::arg("A_col_idx"),
        boost::python::arg("A_values"),
        boost::python::arg("b")
      ))
      )
      .add_property("x",make_getter(&sparse_solver::eigen_solver::x_eig, rbv()))
    ;

  }
}//anon
}//boost_python
}//sparse_solver
BOOST_PYTHON_MODULE(scitbx_examples_strumpack_solver_mpi_ext)
{
   if (import_mpi4py() < 0) return; /* Python 2.X */
   sparse_solver::boost_python::export_strumpack_solver_mpi();
}
