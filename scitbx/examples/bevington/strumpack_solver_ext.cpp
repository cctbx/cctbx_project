/*
 */
#include<Eigen/SparseLU>
#include <iostream>
#include "StrumpackSparseSolver.hpp"
#include <boost/python/enum.hpp>
#include <boost/python.hpp>
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
  struct strumpack_solver{

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
    strumpack_solver(
              const int n_rows, const int n_cols,
              const scitbx::af::shared<intType> A_row_offset,
              const scitbx::af::shared<intType> A_col_idx,
              const scitbx::af::shared<numType> A_values,
              const scitbx::af::shared<numType> b,
              const int reorderEnum,//ReorderingStrategy reorderEnum,
              const int solverEnum //KrylovSolver solverEnum
           ) : x_strum(n_rows, 0.){
      StrumpackSparseSolver<numType,intType> spss(false,true);

      spss.options().set_reordering_method( ReorderingStrategy(reorderEnum) );
      spss.options().set_Krylov_solver( KrylovSolver(solverEnum) );

      std::size_t nnz = A_values.size();
      std::size_t gridSize = n_rows * n_cols;

      //Taking raw pointer to flex data types
      const intType* r_ptr = A_row_offset.begin();
      const intType* c_ptr = A_col_idx.begin();
      const numType* val_ptr = A_values.begin();

      //Set the sparse matrix values using the CSR formatted data
      spss.set_csr_matrix(n_rows, r_ptr, c_ptr, val_ptr, true);
      spss.reorder( n_rows, n_cols );

      //Solve Ax=b
      spss.solve(b.begin(), x_strum.begin());
    }
    scitbx::af::shared<numType> x_strum; //Resulting x will be stored here
  }; //end of struct sparse_solver

  struct eigen_solver{
    //Eigen solver comparison
    eigen_solver(const int solver_Eigen, const int n_rows, const int n_cols,
              const scitbx::af::shared<intType> A_row_offset,
              const scitbx::af::shared<intType> A_col_idx,
              const scitbx::af::shared<numType> A_values,
              const scitbx::af::shared<numType> b) : x_eig(n_rows, 0.){

      //Eigen::SparseMatrix<double> spMat(n_rows,n_cols);
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
           //std::cout << "ROW=" << row << " COL=" <<  *(A_col_idx.begin() + col + r_i) << " VAL=" << *(A_values.begin() + col + r_i) << std::endl;
           tList.push_back(
             Eigen::Triplet<double,int>( row, *(A_col_idx.begin() + col + r_i), *(A_values.begin() + col + r_i) )
           );
        }
      }
      spMat.setFromTriplets( tList.begin(), tList.end() ); //Form the Eigen matrix from the given (i,j,value) sparse format
      spMat.makeCompressed();

      Eigen::VectorXd x(n_rows);

      if(solver_Eigen == 0){
        Eigen::SimplicialLLT<Eigen::SparseMatrix<double> > chol(spMat.transpose());
        x = chol.solve(b_internal);
      }
      else if(solver_Eigen == 1){
        Eigen::SimplicialLDLT<Eigen::SparseMatrix<double> > chol(spMat.transpose());
        x = chol.solve(b_internal);
      }
      else if(solver_Eigen == 2){
        Eigen::SparseMatrix<double> eigen_normal_matrix_full = spMat.selfadjointView<Eigen::Upper>();
        Eigen::BiCGSTAB<Eigen::SparseMatrix<double, Eigen::RowMajor> > solver(eigen_normal_matrix_full);
        x = solver.solve(b_internal);
      }
      else {
        Eigen::SparseMatrix<double> eigen_normal_matrix_full = spMat.selfadjointView<Eigen::Upper>();
        Eigen::ConjugateGradient<Eigen::SparseMatrix<double>, Eigen::Lower|Eigen::Upper> solver(eigen_normal_matrix_full);
        x = solver.solve(b_internal);
      }
      double* solnptr = x_eig.begin();
      for (int i = 0; i<n_rows; ++i){
        *solnptr++ = x[i];
      }
    }
    scitbx::af::shared<numType> x_eig; //Resulting x will be stored here
  }; //end of struct sparse_solver

namespace boost_python{
namespace {
  void export_strumpack_solver() {
    typedef return_value_policy<return_by_value> rbv;

    enum_<KrylovSolver>("KrylovSolver")
      .value("auto", KrylovSolver::AUTO)
      .value("direct", KrylovSolver::DIRECT)
      .value("refine", KrylovSolver::REFINE)
      .value("prec_gmres", KrylovSolver::PREC_GMRES)
      .value("gmres", KrylovSolver::GMRES)
      .value("prec_bicgstab", KrylovSolver::PREC_BICGSTAB)
      .value("bicgstab", KrylovSolver::BICGSTAB)
      .export_values()
    ;

    enum_<ReorderingStrategy>("ReorderingStrategy")
      .value("natural", ReorderingStrategy::NATURAL)
      .value("metis", ReorderingStrategy::METIS)
      .value("parmetis", ReorderingStrategy::PARMETIS)
      .value("scotch", ReorderingStrategy::SCOTCH)
      .value("ptscotch", ReorderingStrategy::PTSCOTCH)
      .value("rcm", ReorderingStrategy::RCM)
      .value("geometric", ReorderingStrategy::GEOMETRIC)
      .export_values()
    ;

    class_<sparse_solver::strumpack_solver> ("strumpack_solver", no_init)
      .def( init< int, int,
                  scitbx::af::shared<intType>,
                  scitbx::af::shared<intType>,
                  scitbx::af::shared<numType>,
                  scitbx::af::shared<numType>,
                  int,
                  int
                >
      ((
        boost::python::arg("n_rows"),
        boost::python::arg("n_cols"),
        boost::python::arg("A_row_offset"),
        boost::python::arg("A_col_idx"),
        boost::python::arg("A_values"),
        boost::python::arg("b"),
        boost::python::arg("reorderEnum")=ReorderingStrategy::SCOTCH,
        boost::python::arg("solverEnum")=KrylovSolver::AUTO
      ))
      )
      .add_property("x",make_getter(&sparse_solver::strumpack_solver::x_strum, rbv()))
    ;

    class_<sparse_solver::eigen_solver> ("eigen_solver", no_init)
      .def( init< int, int, int,
                 scitbx::af::shared<intType>,
                 scitbx::af::shared<intType>,
                 scitbx::af::shared<numType>,
                 scitbx::af::shared<numType>
                >
      ((
        boost::python::arg("solver_type"),
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
BOOST_PYTHON_MODULE(scitbx_examples_strumpack_solver_ext)
{
   sparse_solver::boost_python::export_strumpack_solver();
}
