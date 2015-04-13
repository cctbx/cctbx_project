/// Eigenvalues and eigenvectors computation

#ifndef SCITBX_MATRIX_EIGENSYSTEM_H
#define SCITBX_MATRIX_EIGENSYSTEM_H

#include <scitbx/array_family/accessors/mat_grid.h>
#include <scitbx/array_family/versa.h>
#include <scitbx/array_family/shared.h>
#include <scitbx/array_family/accessors/c_grid.h>
#include <scitbx/array_family/ref_matrix.h>
#include <scitbx/sym_mat3.h>
#include <scitbx/sym_mat2.h>
#include <boost/scoped_array.hpp>

namespace scitbx { namespace matrix { namespace eigensystem {

  namespace detail {

    template<typename FloatType>
    FloatType
    real_symmetric_given_lower_triangle(
      FloatType* a, // size of memory pointed to by a must be n*(n+1)/2
      std::size_t n,
      FloatType* eigenvectors,
      FloatType* eigenvalues,
      FloatType relative_epsilon,
      FloatType absolute_epsilon)
    {
      SCITBX_ASSERT(relative_epsilon >= 0);
      SCITBX_ASSERT(absolute_epsilon >= 0);
      if (n == 0) return 0;
      // The matrix that will hold the results is initially = I.
      std::size_t i;
      for (i=0; i< (n*n); i++) {
        *(eigenvectors+i) = 0.0;
      }
      for (i=0; i< (n*n); i += (n+1)) {
        *(eigenvectors+i) = 1.0;
      }
      // Setup variables
      std::size_t il, ilq, ilr, im, imq, imr, ind, iq;
      std::size_t j, k, km, l, ll, lm, lq, m, mm, mq;
      FloatType am, anorm, anrmx, cosx, cosx2, sincs, sinx, sinx2, thr, x, y;
      // Initial and final norms (anorm & anrmx).
      anorm=0.0;
      iq=0;
      for (i=0; i<n; i++) for (j=0; j<=i; j++) {
        if (j!=i) anorm+=a[iq]*a[iq];
        ++iq;
      }
      anorm=std::sqrt(2.0*anorm);
      anrmx=relative_epsilon*anorm/n;
      if (anrmx < absolute_epsilon) anrmx = absolute_epsilon;
      if (anorm>0.0) {
        // Compute threshold and initialise flag.
        thr=anorm;
        while (thr>anrmx) { // Compare threshold with final norm
          thr/=n;
          ind=1;
          while (ind) {
            ind=0;
            l=0;
            while (l != n-1) { // Test for l beyond penultimate column
              lq=l*(l+1)/2;
              ll=l+lq;
              m=l+1;
              ilq=n*l;
              while (m != n) { // Test for m beyond last column
                // Compute sin & cos.
                mq=m*(m+1)/2;
                lm=l+mq;
                if (a[lm]*a[lm]>thr*thr) {
                  ind=1;
                  mm=m+mq;
                  x=0.5*(a[ll]-a[mm]);
                  FloatType denominator=std::sqrt(a[lm]*a[lm]+x*x);
                  SCITBX_ASSERT(denominator != 0);
                  y=-a[lm]/denominator;
                  if (x<0.0) y=-y;
                  sinx=y/std::sqrt(2.0*(1.0+(std::sqrt(1.0-y*y))));
                  sinx2=sinx*sinx;
                  cosx=std::sqrt(1.0-sinx2);
                  cosx2=cosx*cosx;
                  sincs=sinx*cosx;
                  // Rotate l & m columns.
                  imq=n*m;
                  for (i=0; i<n; i++) {
                    iq=i*(i+1)/2;
                    if (i!=l && i!=m) {
                      if (i<m) im=i+mq;
                      else     im=m+iq;
                      if (i<l) il=i+lq;
                      else     il=l+iq;
                      x=a[il]*cosx-a[im]*sinx;
                      a[im]=a[il]*sinx+a[im]*cosx;
                      a[il]=x;
                    }
                    ilr=ilq+i;
                    imr=imq+i;
                    x = (*(eigenvectors+ilr)*cosx)
                      - (*(eigenvectors+imr)*sinx);
                    *(eigenvectors+imr) = (*(eigenvectors+ilr)*sinx)
                                        + (*(eigenvectors+imr)*cosx);
                    *(eigenvectors+ilr) = x;
                  }
                  x=2.0*a[lm]*sincs;
                  y=a[ll]*cosx2+a[mm]*sinx2-x;
                  x=a[ll]*sinx2+a[mm]*cosx2+x;
                  a[lm]=(a[ll]-a[mm])*sincs+a[lm]*(cosx2-sinx2);
                  a[ll]=y;
                  a[mm]=x;
                }
                m++;
              }
              l++;
            }
          }
        }
      }
      // Sort eigenvalues & eigenvectors in order of descending eigenvalue.
      k=0;
      for (i=0; i<n-1; i++) {
        im=i;
        km=k;
        am=a[k];
        l=0;
        for (j=0; j<n; j++) {
          if (j>i && a[l]>am) {
            im=j;
            km=l;
            am=a[l];
          }
          l+=j+2;
        }
        if (im!=i) {
          a[km]=a[k];
          a[k]=am;
          l=n*i;
          m=n*im;
          for (j=0; j<n; j++) {
            am=*(eigenvectors+l);
            *(eigenvectors+(l++)) = *(eigenvectors+m);
            *(eigenvectors+(m++)) = am;
          }
        }
        k+=i+2;
      }
      // place sorted eigenvalues into the matrix_vector structure
      for (j=0, k=0; j<n; j++) {
        eigenvalues[j]=a[k];
        k+=j+2;
      }
      return anrmx;
    }

    /* Based on code from BTL3: btl/btl_matrix_algorithms.h
       Original comments:

       [Description="Eigenvalues and eigenvectors of a symmetric matrix"]
       [Restrictions="the input matrix must be real, square and
        symmetric, the output vector of eigen values will be nrows in size
        and the output matrix of eigenvectors will be nrows*nrows in size."]

       Routine originally from IBM SSP manual (see p165) Ian Tickle April 1992,
       (modified by  David Moss February 1993 and Mark Williams November 1998).

       n - number of rows in input matrix
       a - an array of size n*(n+1)/2 containing lower triangle of the original
           n*n matrix in the order:

                  1      2      3    ...
           1    a[0]
           2    a[1]   a[2]
           3    a[3]   a[4]   a[5]   ...

           NOTE a is used as working space and is overwritten.
           Eigenvalues are written into the diagonal elements of a
           i.e.  a[0]  a[2]  a[5]  for a 3*3 matrix.

       r - Resultant matrix of eigenvectors stored columnwise in the same
           order as eigenvalues, initially set equal to identity matrix.
    */
    template<typename FloatType>
    FloatType
    real_symmetric_given_full_matrix(
      const FloatType* first,
      const FloatType* /*last*/,
      std::size_t n,
      FloatType* eigenvectors,
      FloatType* eigenvalues,
      FloatType relative_epsilon,
      FloatType absolute_epsilon)
    {
      boost::scoped_array<FloatType> a(new FloatType[n * (n+1) / 2]);
      FloatType* trng = a.get();
      // Copy lower triangle of the input matrix to a numeric vector.
      for (std::size_t row = 1; row <= n; row++) {
        for (const FloatType* in=(first+(n*(row-1)));
             in!=(first+(n*(row-1)+row));
             in++, trng++) {
          *trng = *in;
        }
      }
      return real_symmetric_given_lower_triangle(
        a.get(), n, eigenvectors, eigenvalues,
        relative_epsilon, absolute_epsilon);
    }

  } // namespace detail

  //! Group of associated eigenvectors and eigenvalues.
  /**
    The Cyclic Jacobi algorith is used (Algorithm 8.4.3 in the Golub
    and Van Loan, 3rd edition).

    It is not generally competitive with the symmetric QR method (c.f. Golub and
    Van Loan, section 8.3), which can be found in LAPACK, but this should not be
    very dramatic for small matrices. The main focus of the modern interest in
    the Jacobi method is that it can be heavily parallelised, but the
    implementation in this class does not take advantage of that.

    The algorithm could be tuned to achieve the same precision with less
    iterations by using the element-wise convergence criteria described in
    section 8.4.5, instead of the norm-wise one used here.
  */
  template <typename FloatType = double>
  class real_symmetric
  {
    public:
      //! Default constructor.
      real_symmetric() {}

      /*! \brief Determines the eigenvectors and eigenvalues of the
          real-symmetric square matrix.
       */
      real_symmetric(
        af::const_ref<FloatType, af::mat_grid> const& m,
        FloatType relative_epsilon=1.e-10,
        FloatType absolute_epsilon=0)
      {
        initialize(m, relative_epsilon, absolute_epsilon);
      }

      /*! \brief Determines the eigenvectors and eigenvalues of the
          real-symmetric square matrix.
       */
      real_symmetric(
        scitbx::sym_mat3<FloatType> const& m,
        FloatType relative_epsilon=1.e-10,
        FloatType absolute_epsilon=0);

      real_symmetric(
        scitbx::sym_mat2<FloatType> const& m,
        FloatType relative_epsilon=1.e-10,
        FloatType absolute_epsilon=0);

      //! Threshold used in the cyclic Jacobi algorithm.
      FloatType
      min_abs_pivot() const { return min_abs_pivot_; }

      //! The list of eigenvectors.
      af::versa<FloatType, af::c_grid<2> >
      vectors() const { return vectors_; }

      //! The eigenvalues.
      af::shared<FloatType>
      values() const { return values_; }

      //! Generalized inverse, using min_abs_pivot() in eigenvalue filtering.
      af::shared<FloatType>
      generalized_inverse_as_packed_u() const;

    private:
      FloatType min_abs_pivot_;
      af::versa<FloatType, af::c_grid<2> > vectors_;
      af::shared<FloatType> values_;

      void
      initialize(
        af::const_ref<FloatType, af::mat_grid> const& m,
        FloatType relative_epsilon,
        FloatType absolute_epsilon);
  };

  template <typename FloatType>
  void
  real_symmetric<FloatType>::initialize(
    af::const_ref<FloatType, af::mat_grid> const& m,
    FloatType relative_epsilon,
    FloatType absolute_epsilon)
  {
    SCITBX_ASSERT(m.is_square());
    vectors_.resize(af::c_grid<2>(m.n_rows(), m.n_rows()));
    values_.resize(m.n_rows());
    min_abs_pivot_ = detail::real_symmetric_given_full_matrix(
      m.begin(),
      m.end(),
      m.n_rows(),
      vectors_.begin(),
      values_.begin(),
      relative_epsilon,
      absolute_epsilon);
  }

  template <typename FloatType>
  real_symmetric<FloatType>::real_symmetric(
    scitbx::sym_mat3<FloatType> const& m,
    FloatType relative_epsilon,
    FloatType absolute_epsilon)
  {
    FloatType a[6];
    a[0] = m(0,0);
    a[1] = m(1,0);
    a[2] = m(1,1);
    a[3] = m(2,0);
    a[4] = m(2,1);
    a[5] = m(2,2);
    vectors_.resize(af::c_grid<2>(3, 3));
    values_.resize(3);
    min_abs_pivot_ = detail::real_symmetric_given_lower_triangle(
      a,
      std::size_t(3),
      vectors_.begin(),
      values_.begin(),
      relative_epsilon,
      absolute_epsilon);
  }

  template <typename FloatType>
  real_symmetric<FloatType>::real_symmetric(
    scitbx::sym_mat2<FloatType> const& m,
    FloatType relative_epsilon,
    FloatType absolute_epsilon)
  {
    FloatType a[3];
    a[0] = m(0,0);
    a[1] = m(1,0);
    a[2] = m(1,1);
    vectors_.resize(af::c_grid<2>(2, 2));
    values_.resize(2);
    min_abs_pivot_ = detail::real_symmetric_given_lower_triangle(
      a,
      std::size_t(2),
      vectors_.begin(),
      values_.begin(),
      relative_epsilon,
      absolute_epsilon);
  }

  template <typename FloatType>
  af::shared<FloatType>
  real_symmetric<FloatType>::generalized_inverse_as_packed_u() const
  {
    unsigned n = values_.size();
    af::shared<FloatType>
      result(n*(n+1)/2, af::init_functor_null<FloatType>());
    boost::scoped_array<FloatType> diagonal_elements(new FloatType[n]);
    const FloatType* v = values_.begin();
    for(unsigned i=0;i<n;i++) {
      FloatType const vi = v[i];
      if ((-min_abs_pivot_ < vi && vi < min_abs_pivot_) || vi == 0) {
        diagonal_elements[i] = 0;
      }
      else {
        diagonal_elements[i] = 1 / vi;
      }
    }
    matrix::transpose_multiply_diagonal_multiply_as_packed_u(
      vectors_.begin(), diagonal_elements.get(), n, result.begin());
    return result;
  }

}}} // namespace scitbx::matrix::eigensystem


#endif // GUARD
