#ifndef SCITBX_EXAMPLES_BEVINGTON_H
#define SCITBX_EXAMPLES_BEVINGTON_H

#include <scitbx/array_family/shared.h>
#include <boost/python/tuple.hpp>
#include <scitbx/lstbx/normal_equations.h>
#include <vector>
#include <Eigen/Sparse>

using std::size_t;

namespace scitbx{
namespace example{

class linear_ls_eigen_wrapper
  {

  public:
    typedef Eigen::SparseMatrix<double> sparse_matrix_t;

    /// Construct a least-squares problem with the given number of unknowns.
    linear_ls_eigen_wrapper(int n_parameters)
      : solved_(false),
        eigen_normal_matrix(n_parameters,n_parameters),
        scitbx_normal_matrix(0),
        right_hand_side_(n_parameters),
        solution_(n_parameters)
    {}

    /// Number of unknown parameters
    int n_parameters() const { return right_hand_side_.size(); }


    /// Reset the state to construction time, i.e. no equations accumulated
    void reset() {
      solved_ = false;
      eigen_normal_matrix = sparse_matrix_t(n_parameters(),n_parameters());
      std::fill(scitbx_normal_matrix.begin(), scitbx_normal_matrix.end(), double(0));
      std::fill(right_hand_side_.begin(), right_hand_side_.end(), double(0));
      std::fill(solution_.begin(), solution_.end(), double(0));
    }

    /// Only available if the equations have not been solved yet
    scitbx::af::shared<double>  right_hand_side() const {
      SCITBX_ASSERT(!solved());
      return right_hand_side_;
    }

    void solve() {
      int N = n_parameters();
      //SCITBX_EXAMINE(eigen_normal_matrix.nonZeros());
      int matsize = eigen_normal_matrix.cols() * (eigen_normal_matrix.cols()+1)/2;
      //printf("NM edge %10d, upper triangle size %10d, non-zeros %10d, %6.2f%%\n",
      //  n_parameters(),
      //  matsize,
      //  eigen_normal_matrix.nonZeros(),
      //  eigen_normal_matrix.nonZeros()/double(matsize));

      Eigen::SimplicialCholesky<sparse_matrix_t> chol(eigen_normal_matrix.transpose());
      // XXX pack the right hand side in a eigen vector type
      Eigen::VectorXd b(n_parameters());
      double* rhsptr = right_hand_side_.begin();
      for (int i = 0; i<N; ++i){
        b[i] = *rhsptr++;
      }

      Eigen::VectorXd x = chol.solve(b);

      //  XXX put the solution in the solution_ variable
      double* solnptr = solution_.begin();
      for (int i = 0; i<N; ++i){
        *solnptr++ = x[i];
      }
      solved_ = true;
    }

    // Only available if the equations have not been solved yet
    scitbx::af::versa<double, scitbx::af::packed_u_accessor> normal_matrix() const {
      SCITBX_ASSERT(!solved());
      int N = n_parameters();

      scitbx::af::versa<double, scitbx::af::packed_u_accessor> result(n_parameters());
      double* ptr = result.begin();
      //SCITBX_EXAMINE(result.size());
      //SCITBX_EXAMINE(n_parameters());
      //SCITBX_EXAMINE(eigen_normal_matrix.nonZeros());
      // loop only thru non-zero elements to populate the result array.
      for (int k=0; k<eigen_normal_matrix.outerSize(); ++k) { // column major, so outer (slow) means loop over column
        for (sparse_matrix_t::InnerIterator it(eigen_normal_matrix,k); it; ++it) {
          int irow = it.row();   // row index
          int icol = it.col();   // col index (here it is equal to k)
          std::size_t offset_slow = N * irow - ( irow * (irow - 1) ) / 2;
          std::size_t offset_fast = icol - irow;
          ptr[ offset_slow + offset_fast ] = it.value();
        }
      }
      return result;
    }

    scitbx::af::shared<double> get_cholesky_diagonal() const {
      SCITBX_ASSERT (!solved_);
      int N = n_parameters();
      scitbx::af::shared<double> diagonal = scitbx::af::shared<double>(N);
      Eigen::SimplicialLDLT<sparse_matrix_t> chol(eigen_normal_matrix.transpose());
      Eigen::VectorXd diagonal_eigen = chol.vectorD();

      for (int k=0; k<N; ++k) { // column major, so outer (slow) means loop over column
        diagonal[k] = diagonal_eigen[k]; // copy over the diagonal into persistent array_family type
      }
      return diagonal;
    }

    scitbx::af::shared<double> get_cholesky_lower() const {
      SCITBX_ASSERT (!solved_);
      int N = n_parameters();
      scitbx::af::versa<double, scitbx::af::packed_l_accessor> triangular_result(n_parameters());
      Eigen::SimplicialLDLT<sparse_matrix_t> chol(eigen_normal_matrix.transpose());
      sparse_matrix_t lower = chol.matrixL();

      double* ptr = triangular_result.begin();
      for (int k=0; k<lower.outerSize(); ++k) { // column major, so outer (slow) means loop over column
        for (sparse_matrix_t::InnerIterator it(lower,k); it; ++it) {
          int irow = it.row();   // row index
          int icol = it.col();   // col index (here it is equal to k)
          std::size_t offset_slow = ( irow * irow +  irow )  / 2;
          std::size_t offset_fast = icol;
          ptr[ offset_slow + offset_fast ] = it.value();
        }
      }
      return triangular_result;
    }

    scitbx::af::shared<int> get_eigen_permutation_ordering() const {
      SCITBX_ASSERT (!solved_);
      int N = n_parameters();
      scitbx::af::shared<int> one_D_result(n_parameters());
      Eigen::SimplicialLDLT<sparse_matrix_t> chol(eigen_normal_matrix.transpose());
      for (int k=0; k<N; ++k){
         one_D_result[k] = chol.permutationP().indices()[k];
      }
      return one_D_result;
    }

    bool solved() const {
      return solved_;
    }

    /// Only available after the equations have been solved
    scitbx::af::shared<double> solution() const {
      SCITBX_ASSERT(solved());
      return solution_;
    }

  public:
    bool solved_;
    sparse_matrix_t eigen_normal_matrix;
    scitbx::af::ref_owning_versa<double, scitbx::af::packed_u_accessor>
                                                                 scitbx_normal_matrix;
    scitbx::af::shared<double> right_hand_side_;
    scitbx::af::shared<double> solution_;
  };


class non_linear_ls_eigen_wrapper:  public scitbx::lstbx::normal_equations::non_linear_ls<double> {
  public:
    typedef Eigen::Triplet<double> triplet_t;

    non_linear_ls_eigen_wrapper(int const& n_parameters):
      scitbx::lstbx::normal_equations::non_linear_ls<double>(0),
      eigen_wrapper(n_parameters)
      {}

    inline void
    reset(){
      scitbx::lstbx::normal_equations::non_linear_ls<double>::reset();
      eigen_wrapper.reset();
      tripletList = scitbx::af::shared<triplet_t>();
    }

    linear_ls_eigen_wrapper& step_equations() {
      return eigen_wrapper;
    }

    /// Respecting encapsulation, add constant to all diagonal elements.
    void add_constant_to_diagonal(double const& increment) {
      // loop only thru non-zero elements to update the eigen array.
      // column major, so outer (slow) means loop over column
      for (int k=0; k<eigen_wrapper.eigen_normal_matrix.outerSize(); ++k) { // column major, so outer (slow) means loop over column
        for (linear_ls_eigen_wrapper::sparse_matrix_t::InnerIterator it(eigen_wrapper.eigen_normal_matrix,k);it;++it) {
          int irow = it.row();   // row index
          int icol = it.col();   // col index (here it is equal to k)
          if (irow!=icol){continue;}
          it.valueRef() = it.value() + increment;
        }
      }
    }

    /// get normal matrix
    scitbx::af::versa<double, scitbx::af::packed_u_accessor> get_normal_matrix() const {
      SCITBX_ASSERT(!eigen_wrapper.solved());
      return eigen_wrapper.normal_matrix();
    }

    /// get diagonal elements of the normal matrix
    scitbx::af::shared<double> get_normal_matrix_diagonal() const {
      SCITBX_ASSERT(!eigen_wrapper.solved());
      int N = eigen_wrapper.n_parameters();

      scitbx::af::shared<double> result(N, scitbx::af::init_functor_null<double>());
      double* ptr = result.begin();
      //SCITBX_EXAMINE(result.size());
      //SCITBX_EXAMINE(N);
      //SCITBX_EXAMINE(eigen_wrapper.eigen_normal_matrix.nonZeros());
      for (int k=0; k<eigen_wrapper.eigen_normal_matrix.outerSize(); ++k) { // column major, so outer (slow) means loop over column
        for (linear_ls_eigen_wrapper::sparse_matrix_t::InnerIterator it(eigen_wrapper.eigen_normal_matrix,k);it;++it) {
          int irow = it.row();   // row index
          int icol = it.col();   // col index (here it is equal to k)
          if (irow!=icol){continue;}
          ptr[ k ] = it.value();
        }
      }
      return result;
    }

    /// Add the equation \f$ A_{i.} x = b_i \f$ with the given weight
    inline
    void add_equation_eigen(double b_i,
                             scitbx::af::const_ref<std::size_t> const &row_idx,
                             scitbx::af::const_ref<double> const &row_data,
                             double w)
    {
      int ndata = row_idx.size();
      for (int i=0; i<ndata; i++)  {
        std::size_t idx_i = row_idx[i];
        eigen_wrapper.right_hand_side_[idx_i] += w * row_data[i] * b_i;
        for (int j=i; j<ndata; ++j) {
          //push this term into the stack, later to be added to normal matrix
          tripletList.push_back( triplet_t(idx_i, row_idx[j], w * row_data[i] * row_data[j]) );
        }
      }
    }

    inline
    void set_from_triplets(){
      eigen_wrapper.eigen_normal_matrix.setFromTriplets(tripletList.begin(), tripletList.end());
    }

    inline void wipe_triplets(){
      //critical to release this memory
      tripletList = scitbx::af::shared<triplet_t>();
    }
    scitbx::af::shared<double> get_cholesky_lower(){
      return eigen_wrapper.get_cholesky_lower();
    }
    scitbx::af::shared<double> get_cholesky_diagonal(){
      return eigen_wrapper.get_cholesky_diagonal();
    }
    scitbx::af::shared<int> get_eigen_permutation_ordering(){
      return eigen_wrapper.get_eigen_permutation_ordering();
    }
  public: /* data */
    linear_ls_eigen_wrapper eigen_wrapper;

  private:
    scitbx::af::shared<triplet_t> tripletList;
};

class bevington_silver {
  // encode the functional form of the silver decay target, also gradients and curvatures
  public:
    typedef scitbx::af::shared<double> vecd;
    bevington_silver(){} // explicitly instantiate the class so it is available for Python derivation

    void set_cpp_data(vecd x, vecd y, vecd w) // pass reference, avoid big copy
    {
       x_obs = x;
       y_obs = y;
       w_obs = w;
    }

    vecd get_xobs(){ return x_obs; }
    vecd get_wobs(){ return w_obs; }

    vecd
    fvec_callable(vecd current_values) {
      vecd y_diff = vecd(x_obs.size());
      for (int i = 0; i < x_obs.size(); ++i){
        double y_calc=current_values[0] +
                      current_values[1] * std::exp( - x_obs[i] / current_values[3]) +
                      current_values[2] * std::exp( - x_obs[i] / current_values[4]);
        y_diff[i] = y_obs[i] - y_calc;
      }
      return y_diff;
    }

    vecd
    gvec_callable(vecd current_values) {
      vecd y_diff = fvec_callable(current_values);
      vecd gvec = vecd(current_values.size());
      cvec = vecd(current_values.size());
      for (int ix = 0; ix < x_obs.size(); ++ix) {
        double prefactor = -2. * w_obs[ix] * y_diff[ix];
        gvec[0] += prefactor;
        gvec[1] += prefactor * std::exp(-x_obs[ix]/current_values[3]);
        gvec[2] += prefactor * std::exp(-x_obs[ix]/current_values[4]);
        gvec[3] += prefactor * current_values[1] * std::exp(-x_obs[ix]/current_values[3]) *
                               (x_obs[ix]/(current_values[3]*current_values[3]));
        gvec[4] += prefactor * current_values[2] * std::exp(-x_obs[ix]/current_values[4]) *
                               (x_obs[ix] * std::pow(current_values[4],-2));
        cvec[0] += 2 * w_obs[ix];
        cvec[1] += 2 * w_obs[ix] * std::pow(std::exp(-x_obs[ix]/current_values[3]),2);
        cvec[2] += 2 * w_obs[ix] * std::pow(std::exp(-x_obs[ix]/current_values[4]),2);
        cvec[3] += 2 * w_obs[ix] * std::pow(current_values[1] * std::exp(-x_obs[ix]/current_values[3]) *
                                   (x_obs[ix]/(current_values[3]*current_values[3])),2);
        cvec[4] += 2 * w_obs[ix] * std::pow(current_values[2] * std::exp(-x_obs[ix]/current_values[4]) *
                                   (x_obs[ix] * std::pow(current_values[4],-2)),2);
      }
      return gvec;
    }

    double functional(vecd current_values) {
      double result = 0;
      vecd fvec = fvec_callable(current_values);
      for (int i = 0; i < fvec.size(); ++i) {
        result += fvec[i]*fvec[i]*w_obs[i];
      }
      return result;
    }

    vecd curvatures() const{
      return cvec;
    }

    vecd x_obs, y_obs, w_obs, cvec;
};

class dense_base_class: public bevington_silver, public scitbx::lstbx::normal_equations::non_linear_ls<double> {
  public:
    dense_base_class(int n_parameters):
      scitbx::lstbx::normal_equations::non_linear_ls<double>(n_parameters)
      {
      }

    void access_cpp_build_up_directly_dense(bool objective_only, scitbx::af::shared<double> current_values) {

        vecd residuals = fvec_callable(current_values);
        if (objective_only){
          add_residuals(residuals.const_ref(), w_obs.const_ref());
          return;
        }

        // add one of the normal equations per each observation
        for (int ix = 0; ix < x_obs.size(); ++ix) {

          scitbx::af::shared<double> jacobian_one_row_data;

          jacobian_one_row_data.push_back( 1. );

          jacobian_one_row_data.push_back( std::exp( -x_obs[ix]/ current_values[3]) );

          jacobian_one_row_data.push_back( std::exp( -x_obs[ix]/ current_values[4]) );

          jacobian_one_row_data.push_back( current_values[1] * std::exp( -x_obs[ix]/ current_values[3]) *
                                           ( x_obs[ix] / (current_values[3]*current_values[3]) ));

          jacobian_one_row_data.push_back( current_values[2] * std::exp( -x_obs[ix]/ current_values[4]) *
                                           ( x_obs[ix] / (current_values[4]*current_values[4]) ));

          add_equation(-residuals[ix], jacobian_one_row_data.const_ref(), w_obs[ix]);
        }

    }
};

class eigen_base_class: public bevington_silver, public non_linear_ls_eigen_wrapper {
  public:
    eigen_base_class(int n_parameters):
      non_linear_ls_eigen_wrapper(n_parameters)
      {
      }

    void access_cpp_build_up_directly_eigen_eqn(bool objective_only, scitbx::af::shared<double> current_values) {

        vecd residuals = fvec_callable(current_values);
        if (objective_only){
          add_residuals(residuals.const_ref(), w_obs.const_ref());
          return;
        }

        // add one of the normal equations per each observation
        for (int ix = 0; ix < x_obs.size(); ++ix) {

          scitbx::af::shared<std::size_t> jacobian_one_row_indices;
          scitbx::af::shared<double> jacobian_one_row_data;

          jacobian_one_row_indices.push_back( 0 );
          jacobian_one_row_data.push_back( 1. );

          jacobian_one_row_indices.push_back( 1 );
          jacobian_one_row_data.push_back( std::exp( -x_obs[ix]/ current_values[3]) );

          jacobian_one_row_indices.push_back( 2 );
          jacobian_one_row_data.push_back( std::exp( -x_obs[ix]/ current_values[4]) );

          jacobian_one_row_indices.push_back( 3 );
          jacobian_one_row_data.push_back( current_values[1] * std::exp( -x_obs[ix]/ current_values[3]) *
                                           ( x_obs[ix] / (current_values[3]*current_values[3]) ));

          jacobian_one_row_indices.push_back( 4 );
          jacobian_one_row_data.push_back( current_values[2] * std::exp( -x_obs[ix]/ current_values[4]) *
                                           ( x_obs[ix] / (current_values[4]*current_values[4]) ));

          //add_equation(residuals[ix], jacobian_one_row.const_ref(), weights[ix]);
          add_residual(-residuals[ix], w_obs[ix]);
          add_equation_eigen(residuals[ix], jacobian_one_row_indices.const_ref(), jacobian_one_row_data.const_ref(), w_obs[ix]);
        }
        set_from_triplets();
        wipe_triplets();
    }
};


}}

#endif // SCITBX_EXAMPLES_BEVINGTON_H
