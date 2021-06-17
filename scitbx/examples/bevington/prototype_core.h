#ifndef SCITBX_EXAMPLES_BEVINGTON_H
#define SCITBX_EXAMPLES_BEVINGTON_H

#include <scitbx/array_family/shared.h>
#include <boost/python/tuple.hpp>
#include <scitbx/lstbx/normal_equations.h>
#include <scitbx/sparse/matrix.h>
#include <vector>
#include <Eigen/Sparse>

using std::size_t;

namespace scitbx{
namespace example{

/* Enforce a standard workflow from triplets to normal matrix to Cholesky factor:

Any time:
  n_parameters()
  reset()
  solved()
  formed_normal_matrix()

Upon instantiation, !solved(), !formed_normal_matrix():
  non_linear::add_equations()
  non_linear::add_equation_eigen()
  non_linear::form_normal_matrix()--> set_from_triplets() + wipe_triplets()

After formation of eigen_normal_matrix, !solved():
  solve()
  normal_matrix()
  show_eigen_summary()
  get_cholesky_diagonal()
  get_cholesky_lower()
  get_eigen_permutation_ordering()
  non_linear::add_constant_to_diagonal()
  non_linear::get_normal_matrix()
  non_linear::get_normal_matrix_diagonal()

After solution:
  solution()

Note: non_linear_ls_eigen_wrapper *has a* linear_ls_eigen_wrapper (not *is_a*)
*/

class linear_ls_eigen_wrapper
  {

  public:
    typedef Eigen::SparseMatrix<double> sparse_matrix_t;

    /// Construct a least-squares problem with the given number of unknowns.
    linear_ls_eigen_wrapper(int n_parameters)
      : solved_(false),
        formed_normal_matrix_(false),
        eigen_normal_matrix(n_parameters,n_parameters),
        scitbx_normal_matrix(0),
        right_hand_side_(n_parameters),
        solution_(n_parameters)
    {}

    /// Number of unknown parameters
    long n_parameters() const { return right_hand_side_.size(); }


    /// Reset the state to construction time, i.e. no equations accumulated
    void reset() {
      solved_ = false;
      formed_normal_matrix_ = false;
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
      SCITBX_ASSERT(formed_normal_matrix());
      int N = n_parameters();
      Eigen::SimplicialLDLT<sparse_matrix_t> chol(eigen_normal_matrix.transpose());
      // XXX pack the right hand side in a eigen vector type
      Eigen::VectorXd b(n_parameters());
      double* rhsptr = right_hand_side_.begin();
      for (int i = 0; i<N; ++i){
        b[i] = *rhsptr++;
      }

      Eigen::VectorXd x = chol.solve(b);

      //Try to record state of lower Cholesky factor without incurring cost of computing # non-Zeros again
      sparse_matrix_t lower = chol.matrixL();
      last_computed_matrixL_nonZeros_ = lower.nonZeros();
      // XXX Not sure if this adds measurable time, if so it should be refactored to a separate function

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
      SCITBX_ASSERT(formed_normal_matrix());
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

    void show_eigen_summary() const {
      SCITBX_ASSERT(formed_normal_matrix());
      long matsize = long(eigen_normal_matrix.cols()) * (eigen_normal_matrix.cols()+1)/2;
      printf("Number of parameters      %12ld\n",n_parameters());
      printf("Normal matrix square size %12ld\n",long(eigen_normal_matrix.cols()) * eigen_normal_matrix.cols());
      printf("Upper triangle size       %12ld\n",matsize);
      printf("Normal matrix non-zeros   %12ld, %6.2f%%\n",
              long(eigen_normal_matrix.nonZeros()),
              100. * long(eigen_normal_matrix.nonZeros())/double(matsize));
      Eigen::SimplicialLDLT<sparse_matrix_t> chol(eigen_normal_matrix.transpose());
      sparse_matrix_t lower = chol.matrixL();
      printf("Cholesky factor non-zeros %12ld, %6.2f%%\n",
              long(lower.nonZeros()),
              100. * long(lower.nonZeros())/double(matsize));
    }

    scitbx::af::shared<double> get_cholesky_diagonal() const {
      SCITBX_ASSERT (!solved_);
      SCITBX_ASSERT(formed_normal_matrix());
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
      SCITBX_ASSERT(formed_normal_matrix());
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
      SCITBX_ASSERT(formed_normal_matrix());
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

    bool formed_normal_matrix() const {
      return formed_normal_matrix_;
    }

    /// Only available after the equations have been solved
    scitbx::af::shared<double> solution() const {
      SCITBX_ASSERT(solved());
      return solution_;
    }

  public:
    bool solved_;
    bool formed_normal_matrix_;
    long last_computed_matrixL_nonZeros_;
    sparse_matrix_t eigen_normal_matrix;
    scitbx::af::ref_owning_versa<double, scitbx::af::packed_u_accessor>
                                                                 scitbx_normal_matrix;
    scitbx::af::shared<double> right_hand_side_;
    scitbx::af::shared<double> solution_;
  };


class non_linear_ls_eigen_wrapper:  public scitbx::lstbx::normal_equations::non_linear_ls<double> {
  public:
    typedef Eigen::Triplet<double> triplet_t;
    typedef std::vector<triplet_t> triplet_list_t;

    non_linear_ls_eigen_wrapper(int const& n_parameters):
      scitbx::lstbx::normal_equations::non_linear_ls<double>(0),
      eigen_wrapper(n_parameters)
      {}

    inline void
    reset(){
      scitbx::lstbx::normal_equations::non_linear_ls<double>::reset();
      eigen_wrapper.reset();
      tripletList = triplet_list_t();
    }

    linear_ls_eigen_wrapper& step_equations() {
      return eigen_wrapper;
    }

    /// Respecting encapsulation, add constant to all diagonal elements.
    void add_constant_to_diagonal(double const& increment) {
      form_normal_matrix();
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
    scitbx::af::versa<double, scitbx::af::packed_u_accessor> get_normal_matrix() {
      SCITBX_ASSERT(!eigen_wrapper.solved());
      form_normal_matrix();
      return eigen_wrapper.normal_matrix();
    }

    /// get diagonal elements of the normal matrix
    scitbx::af::shared<double> get_normal_matrix_diagonal() {
      SCITBX_ASSERT(!eigen_wrapper.solved());
      form_normal_matrix();
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

    /// Add equations A x = b given a CCTBX-sparse matrix jacobian
    /** w[i] weights the i-th equation, i.e. the row \f$ A_{i.} \f$.
        The right_hand_side is negated, see lstbx/normal_equations.h
        Function intended for use within SparseLevMar algorithm of DIALS refinement.
     */
    inline
    void add_equations(af::const_ref<scalar_t> const &r,
                       sparse::matrix<scalar_t> const &jacobian,
                       af::const_ref<scalar_t> const &w)
    {
      typedef sparse::matrix<scalar_t>::row_iterator row_iterator;
      SCITBX_ASSERT(!eigen_wrapper.formed_normal_matrix());
      SCITBX_ASSERT(   r.size() == jacobian.n_rows()
                    && (!w.size() || r.size() == w.size()))
                   (r.size())(jacobian.n_rows())(w.size());
      SCITBX_ASSERT(jacobian.n_cols() == eigen_wrapper.n_parameters())
                   (jacobian.n_cols())(eigen_wrapper.n_parameters());
      add_residuals(r, w);
      sparse::matrix<scalar_t> a = jacobian.transpose();
      for (int ieqn = 0; ieqn < a.n_cols(); ++ieqn){
        row_iterator iend = a.col(ieqn).end();
        for (row_iterator i = a.col(ieqn).begin(); i != iend; ++i)  {
          std::size_t idx_i = i.index();
          /* Need to skip over zero values, since the sparse::matrix Jacobian does
             not guarantee that it contains only non-zeroes.  For true sparse
             problems it creates a MemoryError to fill up the tripletList with
             > 10^9 entries.
          */
          if ( (*i)!=0. ){
            eigen_wrapper.right_hand_side_[idx_i] -= w[ieqn] * (*i) * r[ieqn];

            for (row_iterator j=i; j != iend; ++j) {
              if ( (*j)!=0. ) {
                tripletList.push_back( triplet_t(idx_i, j.index(), w[ieqn] * (*i) * (*j)) );
              }
            }
          }
        }
      }
    }

    /// Add the equation \f$ A_{i.} x = b_i \f$ with the given weight
    inline
    void add_equation_eigen(double b_i,
                             scitbx::af::const_ref<std::size_t> const &row_idx,
                             scitbx::af::const_ref<double> const &row_data,
                             double w)
    {
      SCITBX_ASSERT(!eigen_wrapper.formed_normal_matrix());
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

    inline void form_normal_matrix() {
      if (!eigen_wrapper.formed_normal_matrix()){
        set_from_triplets();
        wipe_triplets();
      }
    }

    inline
    void set_from_triplets(){
      SCITBX_ASSERT(!eigen_wrapper.formed_normal_matrix());
      eigen_wrapper.eigen_normal_matrix.setFromTriplets(tripletList.begin(), tripletList.end());
      eigen_wrapper.formed_normal_matrix_=true;
    }

    inline void wipe_triplets(){
      SCITBX_ASSERT(eigen_wrapper.formed_normal_matrix());
      //critical to release this memory
      tripletList = triplet_list_t();
    }
    void show_eigen_summary(){
      form_normal_matrix();
      eigen_wrapper.show_eigen_summary();
    }
    scitbx::af::shared<double> get_cholesky_lower(){
      form_normal_matrix();
      return eigen_wrapper.get_cholesky_lower();
    }
    scitbx::af::shared<double> get_cholesky_diagonal(){
      form_normal_matrix();
      return eigen_wrapper.get_cholesky_diagonal();
    }
    scitbx::af::shared<int> get_eigen_permutation_ordering(){
      form_normal_matrix();
      return eigen_wrapper.get_eigen_permutation_ordering();
    }
    bool solved() const{
      return eigen_wrapper.solved();
    }
    long get_normal_matrix_ncols() const{
      return long(eigen_wrapper.eigen_normal_matrix.cols());
    }
    long n_parameters() const{
      return long(eigen_wrapper.n_parameters());
    }
    long get_normal_matrix_nnonZeros() const{
      return long(eigen_wrapper.eigen_normal_matrix.nonZeros());
    }
    long get_lower_cholesky_nnonZeros() const{
      return long(eigen_wrapper.last_computed_matrixL_nonZeros_);
    }

  public: /* data */
    linear_ls_eigen_wrapper eigen_wrapper;

  private:
    triplet_list_t tripletList;
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
    }
};


}}

#endif // SCITBX_EXAMPLES_BEVINGTON_H
