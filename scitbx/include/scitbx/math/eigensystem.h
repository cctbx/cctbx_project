/* Copyright (c) 2001-2002 The Regents of the University of California
   through E.O. Lawrence Berkeley National Laboratory, subject to
   approval by the U.S. Department of Energy.
   See files COPYRIGHT.txt and LICENSE.txt for further details.

   The real_symmetric implementation is based on code from
   btl/btl_matrix_algorithms.h. Included here under
   the cctbx license with kind permission by the authors
   (Ian Tickle, David Moss, Mark Williams).

   Revision history:
     2003 Mar: Created (R.W. Grosse-Kunstleve)
 */

#ifndef SCITBX_MATH_EIGENSYSTEM_H
#define SCITBX_MATH_EIGENSYSTEM_H

#include <scitbx/mat_ref.h>
#include <scitbx/array_family/versa.h>
#include <scitbx/array_family/shared.h>
#include <scitbx/array_family/accessors/c_grid.h>
#include <vector>

namespace scitbx { namespace math { namespace eigensystem {

  namespace detail {

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
    void
    real_symmetric(
      const FloatType* first,
      const FloatType* last,
      std::size_t n,
      FloatType* eigenvectors,
      FloatType* eigenvalues,
      FloatType epsilon)
    {
      std::vector<FloatType> a(n * (n+1) / 2);
      FloatType* trng = &*a.begin();
      // Copy lower triangle of the input matrix to a numeric vector.
      for (std::size_t row = 1; row <= n; row++) {
        for (const FloatType* in=(first+(n*(row-1)));
             in!=(first+(n*(row-1)+row));
             in++, trng++) {
          *trng = *in;
        }
      }
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
      if (anorm>0.0) {
        anorm=std::sqrt(2.0*anorm);
        anrmx=epsilon*anorm/n;
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
    }

  } // namespace detail

  //! Group of associated eigenvectors and eigenvalues.
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
        mat_const_ref<FloatType> const& m,
        FloatType epsilon=1.e-10);

      /*! \brief Determines the eigenvectors and eigenvalues of the
          real-symmetric square matrix.
       */
      real_symmetric(
        af::const_ref<FloatType, af::c_grid<2> > const& m,
        FloatType epsilon=1.e-10);

      //! The list of eigenvectors.
      af::versa<FloatType, af::c_grid<2> >
      vectors() const { return vectors_; }

      //! The eigenvalues.
      af::shared<FloatType>
      values() const { return values_; }

    private:
      af::versa<FloatType, af::c_grid<2> > vectors_;
      af::shared<FloatType> values_;

      void
      initialize(mat_const_ref<FloatType> const& m, FloatType epsilon);
  };

  template <typename FloatType>
  real_symmetric<FloatType>::real_symmetric(
    mat_const_ref<FloatType> const& m,
    FloatType epsilon)
  {
    initialize(m, epsilon);
  }

  template <typename FloatType>
  real_symmetric<FloatType>::real_symmetric(
    af::const_ref<FloatType, af::c_grid<2> > const& m,
    FloatType epsilon)
  {
    initialize(mat_const_ref<FloatType>(
      m.begin(),
      m.accessor()[0],
      m.accessor()[1]),
      epsilon);
  }

  template <typename FloatType>
  void
  real_symmetric<FloatType>::initialize(
    mat_const_ref<FloatType> const& m,
    FloatType epsilon)
  {
    SCITBX_ASSERT(m.is_square());
    vectors_.resize(af::c_grid<2>(m.n_rows(), m.n_rows()));
    values_.resize(m.n_rows());
    detail::real_symmetric(
      m.begin(),
      m.end(),
      m.n_rows(),
      vectors_.begin(),
      values_.begin(),
      epsilon);
  }

}}} // namespace scitbx::math::eigensystem

#endif // SCITBX_MATH_EIGENSYSTEM_H
