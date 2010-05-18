#ifndef SCITBX_SPARSE_RANDOM_H
#define SCITBX_SPARSE_RANDOM_H

#include<scitbx/sparse/matrix.h>
#include<scitbx/random.h>

namespace scitbx { namespace sparse {

/// A sequence of random sparse matrices
template <typename T>
class random_matrix_generator
{
  private:
    random::mersenne_twister gen;
    std::size_t non_zeroes_;

  public:
    typedef typename matrix<T>::index_type index_type;
    typedef typename matrix<T>::value_type value_type;

    /// A matrix with the given number of rows and columns
    /// and random elements within the given bounds and the given ratio
    /// of non-zero elements wrt the total number of elements.
    matrix<T> operator()(index_type rows,
                         index_type cols,
                         value_type lower_bound,
                         value_type upper_bound,
                         value_type sparsity)
    {
      matrix<T> result(rows, cols);
      non_zeroes_ = static_cast<std::size_t>(rows*cols*sparsity);
      std::size_t nz = non_zeroes_;
      while (nz) {
        typename matrix<T>::index_type i = gen.random_size_t() % rows;
        typename matrix<T>::index_type j = gen.random_size_t() % cols;
        if (result.is_structural_zero(i,j)) {
          result(i,j) = lower_bound
                        + gen.random_double()*(upper_bound - lower_bound);
          --nz;
        }
      }
      return result;
    }

    /// Number of non_zero entries in the last matrix returned by operator()
    std::size_t non_zeroes() { return non_zeroes_; }
};

}} // scitbx::sparse

#endif // GUARD
