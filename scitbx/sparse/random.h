#ifndef SCITBX_SPARSE_RANDOM_H
#define SCITBX_SPARSE_RANDOM_H

#include <scitbx/sparse/matrix.h>
#include <scitbx/random/variate_generator.h>
#include <scitbx/math/utils.h>
#include <boost/random/uniform_int.hpp>
#include <algorithm>

namespace scitbx { namespace sparse {

namespace details {

  // Generator of non-zero locations in sparse vectors and matrices
  template <class IndexType>
  struct random_non_zero_locations
  {
    typedef IndexType index_type;
    typedef math::float_int_conversions<double, index_type> convert_to;

    // For vectors, mere indices
    // For matrices, compounded indices made from (i,j)
    af::shared<index_type> indices;

    // Number of non-zero locations
    index_type nz;

    // For vector, its dimension
    // For matrices, rows x columns
    index_type range;

    index_type range_with_leeway(index_type r) {
      return std::min(index_type(r*1.1), index_type(r + 10));
    }

    // Construct a generator of nz-tuple of indices in the range [0, r)
    random_non_zero_locations(index_type nz, index_type r)
      : nz(nz), range(r),
        indices(range_with_leeway(r),
                af::init_functor_null<index_type>())
    {
      SCITBX_ASSERT(0 < nz && nz < range)(nz)(range);
    }

    // Alternative constructor with a density d: s = d*r
    random_non_zero_locations(double density, index_type r)
      : nz(convert_to::nearest_integer(density*r)), range(r),
        indices(range_with_leeway(r),
                af::init_functor_null<index_type>())
    {
      SCITBX_ASSERT(0 < density && density < 1)(density);
    }

    /* Generate about s and never more than s unique indices,
     the exact number being stored in member n_uniques
     */
    template <class Engine>
    af::const_ref<index_type> generate(Engine &eng) {
      typedef boost::uniform_int<index_type> indices_dist_t;
      typedef random::variate_generator<Engine &, indices_dist_t>
              indices_variate_t;
      indices_dist_t indices_dist(0, range-1);
      indices_variate_t indices_variate(eng, indices_dist);
      std::generate(indices.begin(), indices.end(), indices_variate);
      std::sort(indices.begin(), indices.end());
      index_type unique_nz = std::unique(indices.begin(), indices.end())
                              - indices.begin();
      if (unique_nz > nz) unique_nz = nz;
      return af::const_ref<index_type>(indices.begin(), unique_nz);
    }
  };


}


/// Random sparse matrix distribution.
/** This abides to the Distribution concept of boost::random.
 */
template <typename FloatType, class ElementDistribution>
class matrix_distribution
{
public:
  typedef typename matrix<FloatType>::index_type index_type;
  typedef typename matrix<FloatType>::value_type value_type;

  /// boost::random Distribution requirement
  typedef matrix<FloatType> result_type;

  /// boost::random Distribution requirement
  typedef typename ElementDistribution::input_type input_type;

  /// Construct a distribution of matrices with the given dimensions
  /// and density.
  /** The density is the proportion of non-zero elements.
   */
  matrix_distribution(index_type n_rows, index_type n_cols,
                      double density,
                      ElementDistribution &d)
    : n_rows_(n_rows), n_cols_(n_cols),
      random_non_zero_locations(density, n_rows*n_cols),
      element_distribution(d)
  {}

  /// Construct a distribution of matrices with the given dimensions
  /// and number of structural zeroes.
  matrix_distribution(index_type n_rows, index_type n_cols,
                      index_type non_zeroes,
                      ElementDistribution &d)
    : n_rows_(n_rows), n_cols_(n_cols),
      random_non_zero_locations(non_zeroes, n_rows*n_cols),
      element_distribution(d)
  {}

  /// boost::random Distribution requirement
  template <class Engine>
  result_type operator()(Engine &eng) {
    index_type m = n_rows(), n = n_cols();
    af::const_ref<index_type> indices = random_non_zero_locations.generate(eng);
    result_type result(m, n);
    for (index_type k=0; k < indices.size(); ++k) {
      index_type j = indices[k]/m, i = indices[k] % m;
      result(i,j) = element_distribution(eng);
    }
    return result;
  }

  /// boost::random Distribution requirement
  void reset() { element_distribution.reset(); }

  /// Number of rows of the generated matrices
  index_type n_rows() { return n_rows_; }

  /// Number of rows of the generated matrices
  index_type n_cols() { return n_cols_; }

  /// Exact number of structural zeroes.
  index_type non_zeroes() { return random_non_zero_locations.nz; }

private:
  index_type n_rows_, n_cols_;
  details::random_non_zero_locations<index_type> random_non_zero_locations;
  ElementDistribution element_distribution;
};


/// Random sparse vector distribution.
/** This abides to the Distribution concept of boost::random.
 */
template <typename FloatType, class ElementDistribution>
class vector_distribution
{
public:
  typedef typename vector<FloatType>::index_type index_type;
  typedef typename vector<FloatType>::value_type value_type;

  /// boost::random Distribution requirement
  typedef vector<FloatType> result_type;

  /// boost::random Distribution requirement
  typedef typename ElementDistribution::input_type input_type;

  /// Construct a distribution of matrices with the given dimensions
  /// and density.
  /** The density is the proportion of non-zero elements.
   */
  vector_distribution(index_type size,
                      double density,
                      ElementDistribution &d)
  : size_(size),
    random_non_zero_locations(density, size),
    element_distribution(d)
  {}

  /// Construct a distribution of matrices with the given dimensions
  /// and number of structural zeroes.
  vector_distribution(index_type size,
                      index_type non_zeroes,
                      ElementDistribution &d)
  : size_(size),
    random_non_zero_locations(non_zeroes, size),
    element_distribution(d)
  {}

  /// boost::random Distribution requirement
  template <class Engine>
  result_type operator()(Engine &eng) {
    af::const_ref<index_type> indices = random_non_zero_locations.generate(eng);
    result_type result(size());
    for (index_type k=0; k < indices.size(); ++k) {
      result[ indices[k] ] = element_distribution(eng);
    }
    return result;
  }

  /// boost::random Distribution requirement
  void reset() { element_distribution.reset(); }

  /// Number of rows of the generated matrices
  index_type size() { return size_; }

  /// Exact number of structural zeroes.
  index_type non_zeroes() { return random_non_zero_locations.nz; }

private:
  index_type size_;
  details::random_non_zero_locations<index_type> random_non_zero_locations;
  ElementDistribution element_distribution;
};

}} // scitbx::sparse

#endif // GUARD
