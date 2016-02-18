#pragma once

#include <scitbx/array_family/tiny_plain.h>
#include <scitbx/array_family/accessors/f_grid.h>

#include <boost/static_assert.hpp>

namespace scitbx { namespace af {

  enum upper_lower_tag { upper, lower };

  /// Accessor for the rectangular full packed format (RFPF)
  /** This is a data format for storing triangular, symmetric, and Hermitian
      matrices. It utilises the same amount of memory as the packed format
      provided by af::packed_u_accessor and af::packed_l_accessor but it
      enables much more efficient algorithm by making use of BLAS level 3.

      Ref: Gustavson, F. G., Wasniewski, J., Dongarra, J. J., and Langou, J.
      Rectangular full packed format for cholesky’s algorithm:
      Factorization, solution, and inversion.
      ACM Trans. Math. Softw. 37, 2 (Apr. 2010), 18:1–18:21

      and: file dtfttp.f from LAPACK
  */
  template<upper_lower_tag uplo>
  class rectangular_full_packed_accessor
  {
  private:
    BOOST_STATIC_ASSERT_MSG(uplo == upper,
                            "Lower diagonal rectangular full packed "
                            "not implemented");
  public:
    typedef unsigned index_value_type;
    typedef af::tiny_plain<unsigned, 2> index_type;

    index_value_type n;

    /// Construct an accessor to a n x n matrix
    explicit
    rectangular_full_packed_accessor(unsigned n=0)
    : n(n),
      ld(n&1 ? n : n+1), tc(n/2), o(ld - tc)
    {};

    /// The size of the storage
    std::size_t size_1d() const { return n*(n+1)/2; }

    /// The index in the storage array for element (i,j) of the matrix
    /** Precondition: i <= j
        Not enforced for the sake of efficiency
      */
    std::size_t operator()(index_value_type i, index_value_type j) const {
      return j < tc ? o + j + i*ld : i + (j-tc)*ld;
    }

    /// The index in the storage array for element (i[0], i[1]) of the matrix
    std::size_t operator()(index_type const &i) const {
      return (*this)(i[0], i[1]);
    }

  private:
    // leading dimension,
    // size of lower triangle,
    // pos of upper left corner of lower triangle
    // C.f. comments in dtfttp.f and diagrams therein:
    // * for N=6: the lower triangle is 00 to 22
    // * for N=5: it is 00 to 11
    index_value_type ld, tc, o;
  };

}}
