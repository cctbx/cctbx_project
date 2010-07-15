#ifndef SCITBX_LBFGSB_RAW_H
#define SCITBX_LBFGSB_RAW_H

#include <scitbx/lbfgs.h>
#include <scitbx/math/floating_point_epsilon.h>
#include <scitbx/array_family/ref.h>
#include <scitbx/array_family/accessors/c_grid.h>
#include <scitbx/error.h>
#include <boost/timer.hpp>

#include <cstdio>

#if !defined(SCITBX_LBFGSB_RAW_ASSERTION_FLAG)
# if (defined(__GNUC__) && __GNUC__ < 3) \
  || (defined(__ICC))
#  define SCITBX_LBFGSB_RAW_ASSERTION_FLAG 0
# else
#  define SCITBX_LBFGSB_RAW_ASSERTION_FLAG 0
# endif
#endif

namespace scitbx {
namespace lbfgsb {
//! C++ port of the raw FORTRAN interface of L-BFGS-B Version 2.1
/*! Original FORTRAN distribution:

      http://www.ece.northwestern.edu/~nocedal/lbfgsb.html

    Written by Ciyou Zhu in collaboration with R.H. Byrd, P. Lu-Chen
    and J. Nocedal.

    C++ port by Ralf W. Grosse-Kunstleve.
 */
namespace raw {

  using std::printf;

  template <typename ElementType>
  class ref2;

  //! Emulation of 1-dimensional FORTRAN arrays with offset 1.
  template <typename ElementType>
  class ref1 : public af::ref<ElementType>
  {
    public:
      ref1() {}

      ref1(af::ref<ElementType> const& r)
      :
        af::ref<ElementType>(r)
      {}

      ref1(ElementType* begin, int n)
      :
        af::ref<ElementType>(begin, n)
      {}

      ElementType&
      operator()(int i) const
      {
#if (SCITBX_LBFGSB_RAW_ASSERTION_FLAG != 0)
        SCITBX_ASSERT(i > 0);
        SCITBX_ASSERT(i <= this->size());
#endif
        return this->operator[](i-1);
      }

      ref1
      get1(int i, int n) const
      {
#if (SCITBX_LBFGSB_RAW_ASSERTION_FLAG != 0)
        SCITBX_ASSERT(i-1+n <= this->size());
#endif
        return ref1(&(this->operator()(i)), n);
      }

      ref2<ElementType>
      get2(int i, int n, int m) const
      {
#if (SCITBX_LBFGSB_RAW_ASSERTION_FLAG != 0)
        SCITBX_ASSERT(i-1+n*m <= this->size());
#endif
        return ref2<ElementType>(&(this->operator()(i)), n, m);
      }
  };

  //! Emulation of 2-dimensional FORTRAN arrays with offset 1.
  template <typename ElementType>
  class ref2 : public af::ref<ElementType, af::c_grid<2> >
  {
    public:
      typedef af::ref<ElementType, af::c_grid<2> > base_t;

      ref2() {}

      ref2(ElementType* begin, int n, int m)
      :
        af::ref<ElementType, af::c_grid<2> >(begin, af::c_grid<2>(m, n)),
        hard_size_(n*m)
      {}

      ElementType&
      operator()(int i, int j) const
      {
#if (SCITBX_LBFGSB_RAW_ASSERTION_FLAG != 0)
        SCITBX_ASSERT(i > 0);
        SCITBX_ASSERT(j > 0);
        SCITBX_ASSERT(i <= this->accessor()[1]);
        SCITBX_ASSERT(j <= this->accessor()[0]);
#endif
        std::size_t offs = this->accessor()(j-1,i-1);
#if (SCITBX_LBFGSB_RAW_ASSERTION_FLAG != 0)
        SCITBX_ASSERT(offs < hard_size_);
#endif
        return base_t::operator[](offs);
      }

      ref1<ElementType>
      get1(int i, int j, int n) const
      {
        ElementType* ref1_begin = &(this->operator()(i,j));
#if (SCITBX_LBFGSB_RAW_ASSERTION_FLAG != 0)
        SCITBX_ASSERT(ref1_begin - this->begin() + n <= this->size());
#endif
        return ref1<ElementType>(ref1_begin, n);
      }

      ref2
      get2(int i, int j, int n, int m) const
      {
        ElementType* ref2_begin = &(this->operator()(i,j));
#if (SCITBX_LBFGSB_RAW_ASSERTION_FLAG != 0)
        std::size_t begin_shift = ref2_begin - this->begin();
        SCITBX_ASSERT(begin_shift + n*m <= this->size());
        SCITBX_ASSERT(begin_shift <= hard_size_);
#endif
        return ref2(ref2_begin, n, m);
      }

      ref2
      get2soft(int i, int j, int n, int m) const
      {
        ElementType* ref2_begin = &(this->operator()(i,j));
        std::size_t begin_shift = ref2_begin - this->begin();
#if (SCITBX_LBFGSB_RAW_ASSERTION_FLAG != 0)
        SCITBX_ASSERT(begin_shift <= hard_size_);
#endif
        ref2 result(ref2_begin, n, m);
        result.hard_size_ = hard_size_ - begin_shift;
        return result;
      }

    protected:
      std::size_t hard_size_;
  };

  //! Emulation of the intrinsic function max for three arguments.
  template <typename T>
  inline
  T
  max3(T const& v0, T const& v1, T const& v2)
  {
    return std::max(std::max(v0,v1),v2);
  }

  //! Current user time for the process.
  template <typename FloatType>
  void
  timer(FloatType& ttime)
  {
    static boost::timer timer_;
    ttime = static_cast<FloatType>(timer_.elapsed());
  }

  //! Emulation of write statement with implicit loop.
  template <typename ElementType>
  void
  write_ref1(
    std::string const& label,
    ref1<ElementType> const& a,
    const char* fmt=" %11.4E")
  {
    printf("\n%s", label.c_str());
    for(int i=1;i<=a.size();i++) {
      if ((i-1)%6 == 0) {
        if (i != 1) {
          printf("\n");
          for(int j=0;j<label.size();j++) {
            printf(" ");
          }
        }
      }
      printf(fmt, a(i));
    }
    printf("\n");
  }

  //===============    L-BFGS-B (version 2.1)   ==========================

  //! copies a vector, x, to a vector, y.
  /*! uses unrolled loops for increments equal to one.
      jack dongarra, linpack, 3/11/78.
   */
  template <typename FloatType>
  void
  dcopy(
    int const& n,
    ref1<FloatType> const& dx,
    int const& incx,
    ref1<FloatType> const& dy,
    int const& incy)
  {
    if(n<=0)return;
    if(incx==1&&incy==1)goto lbl_20;
    { // scope for variables
      // code for unequal increments or equal increments
      // not equal to 1
      int ix = 1;
      int iy = 1;
      if(incx<0)ix = (-n+1)*incx + 1;
      if(incy<0)iy = (-n+1)*incy + 1;
      for(int i=1;i<=n;i++) {
        dy(iy) = dx(ix);
        ix = ix + incx;
        iy = iy + incy;
      }
      return;
    }
    // code for both increments equal to 1
    // clean-up loop
    lbl_20:
    int m = n % 7;
    if(m == 0) goto lbl_40;
    for(int i=1;i<=m;i++) {
      dy(i) = dx(i);
    }
    if(n < 7) return;
    lbl_40:
    int mp1 = m + 1;
    for(int i=mp1;i<=n;i+=7) {
      dy(i) = dx(i);
      dy(i + 1) = dx(i + 1);
      dy(i + 2) = dx(i + 2);
      dy(i + 3) = dx(i + 3);
      dy(i + 4) = dx(i + 4);
      dy(i + 5) = dx(i + 5);
      dy(i + 6) = dx(i + 6);
    }
  }

  //! scales a vector by a constant.
  /*! uses unrolled loops for increment equal to one.
      jack dongarra, linpack, 3/11/78.
      modified 3/93 to return if incx .le. 0.
   */
  template <typename FloatType>
  void
  dscal(
    int const& n,
    FloatType const& da,
    ref1<FloatType> const& dx,
    int const& incx)
  {
    if (n<=0 || incx<=0) return;
    if (incx==1) goto lbl_20;
    { // scope for variables
      // code for increment not equal to 1
      int nincx = n*incx;
      for(int i=1;i<=nincx;i+=incx) {
        dx(i) = da*dx(i);
      }
      return;
    }
    // code for increment equal to 1
    // clean-up loop
    lbl_20:
    int m = n % 5;
    if (m == 0) goto lbl_40;
    for(int i=1;i<=m;i++) {
      dx(i) = da*dx(i);
    }
    if(n < 5) return;
    lbl_40:
    for(int i=m+1;i<=n;i+=5) {
      dx(i) = da*dx(i);
      dx(i + 1) = da*dx(i + 1);
      dx(i + 2) = da*dx(i + 2);
      dx(i + 3) = da*dx(i + 3);
      dx(i + 4) = da*dx(i + 4);
    }
  }

  //! constant times a vector plus a vector.
  template <typename FloatType, typename SizeType>
  inline
  void daxpy(
    SizeType const& n,
    FloatType const& da,
    const FloatType* dx,
    FloatType* dy)
  {
    lbfgs::detail::daxpy(
      n, da, dx, SizeType(0), SizeType(1), dy, SizeType(0), SizeType(1));
  }

  /*! \brief This subroutine computes the infinity norm of the projected
       gradient.
   */
  /*! NEOS, November 1994. (Latest revision June 1996.)
      Optimization Technology Center.
      Argonne National Laboratory and Northwestern University.
      Written by
                         Ciyou Zhu
      in collaboration with R.H. Byrd, P. Lu-Chen and J. Nocedal.
   */
  template <typename FloatType>
  void
  projgr(
    int const& n,
    ref1<FloatType> const& l,
    ref1<FloatType> const& u,
    ref1<int> const& nbd,
    ref1<FloatType> const& x,
    ref1<FloatType> const& g,
    FloatType& sbgnrm)
  {
    FloatType zero = 0;
    sbgnrm = zero;
    for(int i=1;i<=n;i++) {
      FloatType gi = g(i);
      if (nbd(i) != 0) {
        if (gi < zero) {
          if (nbd(i) >= 2) gi = std::max((x(i)-u(i)),gi);
        }
        else {
          if (nbd(i) <= 2) gi = std::min((x(i)-l(i)),gi);
        }
      }
      sbgnrm = std::max(sbgnrm,fn::absolute(gi));
    }
  }

  /*! \brief Solution of t * x = b.
   */
  /*! dtrsl solves systems of the form

                    t * x = b
      or
                    trans(t) * x = b

      where t is a triangular matrix of order n. here trans(t)
      denotes the transpose of the matrix t.

      on entry

          t         double precision(ldt,n)
                    t contains the matrix of the system. the zero
                    elements of the matrix are not referenced, and
                    the corresponding elements of the array can be
                    used to store other information.

          ldt       integer
                    ldt is the leading dimension of the array t.

          n         integer
                    n is the order of the system.

          b         double precision(n).
                    b contains the right hand side of the system.

          job       integer
                    job specifies what kind of system is to be solved.
                    if job is

                         00   solve t*x=b, t lower triangular,
                         01   solve t*x=b, t upper triangular,
                         10   solve trans(t)*x=b, t lower triangular,
                         11   solve trans(t)*x=b, t upper triangular.

      on return

          b         b contains the solution, if info .eq. 0.
                    otherwise b is unaltered.

          info      integer
                    info contains zero if the system is nonsingular.
                    otherwise info contains the index of
                    the first zero diagonal element of t.

      linpack. this version dated 08/14/78 .
      g. w. stewart, university of maryland, argonne national lab.

      subroutines and functions

      blas daxpy,ddot
      fortran mod
   */
  template <typename FloatType>
  void
  dtrsl(
    ref2<FloatType> const& t,
    int const&
#if (SCITBX_LBFGSB_RAW_ASSERTION_FLAG != 0)
    ldt
#endif
    ,
    int const& n,
    ref1<FloatType> const& b,
    int const& job,
    int& info)
  {
#if (SCITBX_LBFGSB_RAW_ASSERTION_FLAG != 0)
    SCITBX_ASSERT(t.size() == ldt*n);
    SCITBX_ASSERT(b.size() == n);
#endif
    FloatType zero = 0;
    // check for zero diagonal elements.
    for(info=1;info<=n;info++) {
      if (t(info,info) == zero) return;
    }
    info = 0;
    // determine the task and go to it.
    int case_ = 1;
    if (job % 10 != 0) case_ = 2;
    if ((job % 100)/10 != 0) case_ = case_ + 2;
    if (case_ == 1) {
      // solve t*x=b for t lower triangular
      b(1) = b(1)/t(1,1);
      if (n < 2) return;
      for(int j=2;j<=n;j++) {
        FloatType temp = -b(j-1);
        daxpy(n-j+1,temp,&t(j,j-1),&b(j));
        b(j) = b(j)/t(j,j);
      }
    }
    else if (case_ == 2) {
      // solve t*x=b for t upper triangular.
      b(n) = b(n)/t(n,n);
      if (n < 2) return;
      for(int jj=2;jj<=n;jj++) {
        int j = n - jj + 1;
        FloatType temp = -b(j+1);
        daxpy(j,temp,&t(1,j+1),&b(1));
        b(j) = b(j)/t(j,j);
      }
    }
    else if (case_ == 3) {
      // solve trans(t)*x=b for t lower triangular.
      b(n) = b(n)/t(n,n);
      if (n < 2) return;
      for(int jj=2;jj<=n;jj++) {
        int j = n - jj + 1;
        b(j) = b(j) - lbfgs::detail::ddot(jj-1,&t(j+1,j),&b(j+1));
        b(j) = b(j)/t(j,j);
      }
    }
    else {
      // solve trans(t)*x=b for t upper triangular.
      b(1) = b(1)/t(1,1);
      if (n < 2) return;
      for(int j=2;j<=n;j++) {
        b(j) = b(j) - lbfgs::detail::ddot(j-1,&t(1,j),&b(1));
        b(j) = b(j)/t(j,j);
      }
    }
  }

  /*! \brief Product of the 2m x 2m middle matrix
      in the compact L-BFGS formula of B and a 2m vector v.
   */
  /*! This subroutine computes the product of the 2m x 2m middle matrix
        in the compact L-BFGS formula of B and a 2m vector v;
        it returns the product in p.

      m is an integer variable.
        On entry m is the maximum number of variable metric corrections
          used to define the limited memory matrix.
        On exit m is unchanged.

      sy is a double precision array of dimension m x m.
        On entry sy specifies the matrix S'Y.
        On exit sy is unchanged.

      wt is a double precision array of dimension m x m.
        On entry wt specifies the upper triangular matrix J' which is
          the Cholesky factor of (thetaS'S+LD^(-1)L').
        On exit wt is unchanged.

      col is an integer variable.
        On entry col specifies the number of s-vectors (or y-vectors)
          stored in the compact L-BFGS formula.
        On exit col is unchanged.

      v is a double precision array of dimension 2col.
        On entry v specifies vector v.
        On exit v is unchanged.

      p is a double precision array of dimension 2col.
        On entry p is unspecified.
        On exit p is the product Mv.

      info is an integer variable.
        On entry info is unspecified.
        On exit info = 0       for normal return,
                     = nonzero for abnormal return when the system
                                 to be solved by dtrsl is singular.

      Subprograms called:

        Linpack ... dtrsl.

                            *  *  *

      NEOS, November 1994. (Latest revision June 1996.)
      Optimization Technology Center.
      Argonne National Laboratory and Northwestern University.
      Written by
                         Ciyou Zhu
      in collaboration with R.H. Byrd, P. Lu-Chen and J. Nocedal.
   */
  template <typename FloatType>
  void
  bmv(
    int const& m,
    ref2<FloatType> const& sy,
    ref2<FloatType> const& wt,
    int const& col,
    ref1<FloatType> const& v,
    ref1<FloatType> const& p,
    int& info)
  {
#if (SCITBX_LBFGSB_RAW_ASSERTION_FLAG != 0)
    SCITBX_ASSERT(sy.size() == m*m);
    SCITBX_ASSERT(wt.size() == m*m);
    SCITBX_ASSERT(v.size() == 2*col);
    SCITBX_ASSERT(p.size() == 2*col);
#endif
    if (col == 0) return;
    // PART I: solve [  D^(1/2)      O ] [ p1 ] = [ v1 ]
    //               [ -L*D^(-1/2)   J ] [ p2 ]   [ v2 ].
    // solve Jp2=v2+LD^(-1)v1.
    p(col + 1) = v(col + 1);
    for(int i=2;i<=col;i++) {
      int i2 = col + i;
      FloatType sum = 0;
      for(int k=1;k<=i-1;k++) {
        sum = sum + sy(i,k)*v(k)/sy(k,k);
      }
      p(i2) = v(i2) + sum;
    }
    // Solve the triangular system
    dtrsl(wt.get2(1,1,m,col),m,col,p.get1(col+1,col),11,info);
    if (info != 0) return;
    // solve D^(1/2)p1=v1.
    for(int i=1;i<=col;i++) {
      p(i) = v(i)/std::sqrt(sy(i,i));
    }
    // PART II: solve [ -D^(1/2)   D^(-1/2)*L'  ] [ p1 ] = [ p1 ]
    //                [  0         J'           ] [ p2 ]   [ p2 ].
    // solve J^Tp2=p2.
    dtrsl(wt.get2(1,1,m,col),m,col,p.get1(col+1,col),01,info);
    if (info != 0) return;
    // compute p1=-D^(-1/2)(p1-D^(-1/2)L'p2)
    //           =-D^(-1/2)p1+D^(-1)L'p2.
    for(int i=1;i<=col;i++) {
      p(i) = -p(i)/std::sqrt(sy(i,i));
    }
    for(int i=1;i<=col;i++) {
      FloatType sum = 0;
      for(int k=i+1;k<=col;k++) {
        sum = sum + sy(k,i)*p(col+k)/sy(i,i);
      }
      p(i) = p(i) + sum;
    }
  }

  /*! \brief Sorting of t.
   */
  /*! This subroutine sorts out the least element of t, and puts the
        remaining elements of t in a heap.

      n is an integer variable.
        On entry n is the dimension of the arrays t and iorder.
        On exit n is unchanged.

      t is a double precision array of dimension n.
        On entry t stores the elements to be sorted,
        On exit t(n) stores the least elements of t, and t(1) to t(n-1)
          stores the remaining elements in the form of a heap.

      iorder is an integer array of dimension n.
        On entry iorder(i) is the index of t(i).
        On exit iorder(i) is still the index of t(i), but iorder may be
          permuted in accordance with t.

      iheap is an integer variable specifying the task.
        On entry iheap should be set as follows:
          iheap .eq. 0 if t(1) to t(n) is not in the form of a heap,
          iheap .ne. 0 if otherwise.
        On exit iheap is unchanged.

      References:
        Algorithm 232 of CACM (J. W. J. Williams): HEAPSORT.

                            *  *  *

      NEOS, November 1994. (Latest revision June 1996.)
      Optimization Technology Center.
      Argonne National Laboratory and Northwestern University.
      Written by
                         Ciyou Zhu
      in collaboration with R.H. Byrd, P. Lu-Chen and J. Nocedal.
   */
  template <typename FloatType>
  void
  hpsolb(
    int const& n,
    ref1<FloatType> const& t,
    ref1<int> const& iorder,
    int const& iheap)
  {
#if (SCITBX_LBFGSB_RAW_ASSERTION_FLAG != 0)
    SCITBX_ASSERT(t.size() == n);
    SCITBX_ASSERT(iorder.size() == n);
#endif
    if (iheap == 0) {
      // Rearrange the elements t(1) to t(n) to form a heap.
      for(int k=2;k<=n;k++) {
        FloatType ddum  = t(k);
        int indxin = iorder(k);
        // Add ddum to the heap.
        int i = k;
        // lbl_10:
        while (i>1) {
          int j = i/2;
          if (ddum < t(j)) {
            t(i) = t(j);
            iorder(i) = iorder(j);
            i = j;
            // goto lbl_10;
          }
          else {
            break;
          }
        }
        t(i) = ddum;
        iorder(i) = indxin;
      }
    }
    // Assign to 'out' the value of t(1), the least member of the heap,
    // and rearrange the remaining members to form a heap as
    // elements 1 to n-1 of t.
    if (n > 1) {
      int i = 1;
      FloatType out = t(1);
      int indxou = iorder(1);
      FloatType ddum = t(n);
      int indxin = iorder(n);
      // Restore the heap
      while (true) { // lbl_30:
        int j = i+i;
        if (j <= n-1) {
          if (t(j+1) < t(j)) j = j+1;
          if (t(j) < ddum ) {
            t(i) = t(j);
            iorder(i) = iorder(j);
            i = j;
            // goto lbl_30;
            continue;
          }
        }
        break;
      }
      t(i) = ddum;
      iorder(i) = indxin;
      // Put the least member in t(n).
      t(n) = out;
      iorder(n) = indxou;
    }
  }

  /*! \brief Check validity of input data.
   */
  /*! This subroutine checks the validity of the input data.

                            *  *  *

      NEOS, November 1994. (Latest revision June 1996.)
      Optimization Technology Center.
      Argonne National Laboratory and Northwestern University.
      Written by
                         Ciyou Zhu
      in collaboration with R.H. Byrd, P. Lu-Chen and J. Nocedal.
   */
  template <typename FloatType>
  void
  errclb(
    int const& n,
    int const& m,
    FloatType const& factr,
    ref1<FloatType> const& l,
    ref1<FloatType> const& u,
    ref1<int> const& nbd,
    std::string& task,
    int& info,
    int& k)
  {
    FloatType zero = 0;
    // Check the input arguments for errors.
    if (n <= 0) task = "ERROR: N .LE. 0";
    if (m <= 0) task = "ERROR: M .LE. 0";
    if (factr < zero) task = "ERROR: FACTR .LT. 0";
    // Check the validity of the arrays nbd(i), u(i), and l(i).
    for(int i=1;i<=n;i++) {
      if (nbd(i) < 0 || nbd(i) > 3) {
        // return
        task = "ERROR: INVALID NBD";
        info = -6;
        k = i;
      }
      if (nbd(i) == 2) {
        if (l(i) > u(i)) {
          // return
          task = "ERROR: NO FEASIBLE SOLUTION";
          info = -7;
          k = i;
        }
      }
    }
  }

  /*! \brief Printing function 1.
   */
  /*! This subroutine prints the input data, initial point, upper and
        lower bounds of each variable, machine precision, as well as
        the headings of the output.

                            *  *  *

      NEOS, November 1994. (Latest revision June 1996.)
      Optimization Technology Center.
      Argonne National Laboratory and Northwestern University.
      Written by
                         Ciyou Zhu
      in collaboration with R.H. Byrd, P. Lu-Chen and J. Nocedal.
   */
  template <typename FloatType>
  void
  prn1lb(
    int const& n,
    int const& m,
    ref1<FloatType> const& l,
    ref1<FloatType> const& u,
    ref1<FloatType> const& x,
    int const& iprint,
    int const& /*itfile*/,
    FloatType const& epsmch)
  {
    if (iprint >= 0) {
      printf("RUNNING THE L-BFGS-B CODE\n");
      printf("\n");
      printf("           * * *\n");
      printf("\n");
      printf("Machine precision =%10.3E\n", epsmch);
      printf(" N = %12d    M = %12d\n", n, m);
      if (iprint >= 1) {
        // iterate.dat
        printf("RUNNING THE L-BFGS-B CODE\n");
        printf("\n");
        printf("it    = iteration number\n");
        printf("nf    = number of function evaluations\n");
        printf(
          "nint  = number of segments explored during the Cauchy search\n");
        printf(
          "nact  = number of active bounds at the generalized Cauchy point\n");
        printf(
          "sub   = manner in which the subspace minimization terminated:\n");
        printf("        con = converged, bnd = a bound was reached\n");
        printf("itls  = number of iterations performed in the line search\n");
        printf("stepl = step length used\n");
        printf("tstep = norm of the displacement (total step)\n");
        printf("projg = norm of the projected gradient\n");
        printf("f     = function value\n");
        printf("\n");
        printf("           * * *\n");
        printf("\n");
        printf("Machine precision =%10.3E\n", epsmch);
        printf(" N = %12d    M = %12d\n", n, m);
        printf("\n");
        printf("   it   nf  nint  nact  sub  itls  stepl"
               "    tstep     projg        f\n");
        if (true || iprint > 100) {
          write_ref1(" L =", l);
          write_ref1("X0 =", x);
          write_ref1(" U =", u);
        }
      }
    }
  }

  /*! \brief Printing function 2.
   */
  /*! This subroutine prints out new information after a successful
        line search.

                            *  *  *

      NEOS, November 1994. (Latest revision June 1996.)
      Optimization Technology Center.
      Argonne National Laboratory and Northwestern University.
      Written by
                         Ciyou Zhu
      in collaboration with R.H. Byrd, P. Lu-Chen and J. Nocedal.
   */
  template <typename FloatType>
  void
  prn2lb(
    int const& /*n*/,
    ref1<FloatType> const& x,
    FloatType const& f,
    ref1<FloatType> const& g,
    int const& iprint,
    int const& /*itfile*/,
    int const& iter,
    int const& nfgv,
    int const& nact,
    FloatType const& sbgnrm,
    int const& nint,
    std::string& word,
    int const& iword,
    int const& iback,
    FloatType const& stp,
    FloatType const& xstep)
  {
    // 'word' records the status of subspace solutions.
    if (iword == 0) {
      // the subspace minimization converged.
      word = "con";
    }
    else if (iword == 1) {
      // the subspace minimization stopped at a bound.
      word = "bnd";
    }
    else if (iword == 5) {
      // the truncated Newton step has been used.
      word = "TNT";
    }
    else {
      word = "---";
    }
    static const char*
      fmt_2001 = "\nAt iterate%5d    f= %12.5E    |proj g|= %12.5E\n";
    if (iprint >= 99) {
      printf(" LINE SEARCH%12d times; norm of step =  %.15G\n",
             iback, xstep);
      printf(fmt_2001, iter,f,sbgnrm);
      if (iprint > 100) {
        write_ref1(" X =", x);
        write_ref1(" G =", g);
      }
    }
    else if (iprint > 0) {
      int imod = iter % iprint;
      if (imod == 0) {
        printf(fmt_2001, iter,f,sbgnrm);
      }
    }
    if (iprint >= 1) {
      // iterate.dat
      printf(" %4d %4d %5d %5d  %-3.3s %4d  %7.1E  %7.1E %10.3E %10.3E\n",
             iter,nfgv,nint,nact,word.c_str(),iback,stp,xstep,sbgnrm,f);
    }
  }

  /*! \brief Printing function 3.
   */
  /*! This subroutine prints out information when either a built-in
        convergence test is satisfied or when an error message is
        generated.

        *  *  *

      NEOS, November 1994. (Latest revision June 1996.)
      Optimization Technology Center.
      Argonne National Laboratory and Northwestern University.
      Written by
                         Ciyou Zhu
      in collaboration with R.H. Byrd, P. Lu-Chen and J. Nocedal.
   */
  template <typename FloatType>
  void
  prn3lb(
    int const& n,
    ref1<FloatType> const& x,
    FloatType const& f,
    std::string const& task,
    int const& iprint,
    int const& info,
    int const& /*itfile*/,
    int const& iter,
    int const& nfgv,
    int const& nintol,
    int const& nskip,
    int const& nact,
    FloatType const& sbgnrm,
    FloatType const& time,
    int const& nint,
    std::string const& word,
    int const& iback,
    FloatType const& stp,
    FloatType const& xstep,
    int const& k,
    FloatType const& cachyt,
    FloatType const& sbtime,
    FloatType const& lnscht)
  {
    static const char* fmt_3002 =
      " %4d %4d %5d %5d  %-3.3s %4d  %7.1E  %7.1E      -          -\n";
    static const char* fmt_3003 =
      "\n"
      "           * * *\n\n"
      "Tit   = total number of iterations\n"
      "Tnf   = total number of function evaluations\n"
      "Tnint = total number of segments explored during Cauchy searches\n"
      "Skip  = number of BFGS updates skipped\n"
      "Nact  = number of active bounds at final generalized Cauchy point\n"
      "Projg = norm of the final projected gradient\n"
      "F     = final function value\n\n"
      "           * * *\n";
    static const char* fmt_3004 =
      "\n   N   Tit  Tnf  Tnint  Skip  Nact     Projg        F\n";
    static const char* fmt_3005 =
      "%5d %4d %4d %6d  %4d %5d  %10.3E  %10.3E\n";
    static const char* fmt_3007 =
      "\n"
      " Cauchy                time%10.3E seconds.\n"
      " Subspace minimization time%10.3E seconds.\n"
      " Line search           time%10.3E seconds.\n";
    static const char* fmt_3008 =
      "\n Total User time%10.3E seconds.\n\n";
    static const char* fmt_3009 =
      "\n%-60.60s\n";
    static const char* fmt_9011 =
      "\n Matrix in 1st Cholesky factorization in formk is not Pos. Def.\n";
    static const char* fmt_9012 =
      "\n Matrix in 2st Cholesky factorization in formk is not Pos. Def.\n";
    static const char* fmt_9013 =
      "\n Matrix in the Cholesky factorization in formt is not Pos. Def.\n";
    static const char* fmt_9014 =
      "\n"
      " Derivative >= 0, backtracking line search impossible.\n"
      "   Previous x, f and g restored.\n"
      " Possible causes: 1 error in function or gradient evaluation;\n"
      "                  2 rounding errors dominate computation.\n";
    static const char* fmt_9015 =
      "\n"
      " Warning:  more than 10 function and gradient\n"
      "   evaluations in the last line search.  Termination\n"
      "   may possibly be caused by a bad search direction.\n";
    static const char* fmt_9018 =
      "\n The triangular system is singular.\n";
    static const char* fmt_9019 =
      "\n"
      " Line search cannot locate an adequate point after 20 function\n"
      "  and gradient evaluations.  Previous x, f and g restored.\n"
      " Possible causes: 1 error in function or gradient evaluation;\n"
      "                  2 rounding error dominate computation.\n";
    if (task.substr(0,5) == "ERROR") goto lbl_999;
    if (iprint >= 0) {
      printf("%s", fmt_3003);
      printf("%s", fmt_3004);
      printf(fmt_3005, n,iter,nfgv,nintol,nskip,nact,sbgnrm,f);
      if (iprint >= 100) {
        write_ref1(" X =", x);
      }
      if (iprint >= 1) printf("  F = %.15G\n", f);
    }
    lbl_999:
    if (iprint >= 0) {
      printf(fmt_3009, task.c_str());
      if (info != 0) {
        if (info == -1) printf("%s", fmt_9011);
        if (info == -2) printf("%s", fmt_9012);
        if (info == -3) printf("%s", fmt_9013);
        if (info == -4) printf("%s", fmt_9014);
        if (info == -5) printf("%s", fmt_9015);
        if (info == -6) {
          printf("  Input nbd(%12d) is invalid.\n", k);
        }
        if (info == -7) {
          printf("  l(%12d) > u(%12d).  No feasible solution.\n", k, k);
        }
        if (info == -8) printf("%s", fmt_9018);
        if (info == -9) printf("%s", fmt_9019);
      }
      if (iprint >= 1) printf(fmt_3007, cachyt,sbtime,lnscht);
      printf(fmt_3008, time);
      if (iprint >= 1) {
        if (info == -4 || info == -9) {
          printf(fmt_3002,
                 iter,nfgv,nint,nact,word.c_str(),iback,stp,xstep); // itfile
        }
        printf(fmt_3009, task.c_str()); // itfile
        if (info != 0) {
          if (info == -1) printf("%s", fmt_9011); // itfile
          if (info == -2) printf("%s", fmt_9012); // itfile
          if (info == -3) printf("%s", fmt_9013); // itfile
          if (info == -4) printf("%s", fmt_9014); // itfile
          if (info == -5) printf("%s", fmt_9015); // itfile
          if (info == -8) printf("%s", fmt_9018); // itfile
          if (info == -9) printf("%s", fmt_9019); // itfile
        }
        printf(fmt_3008, time); // itfile
      }
    }
  }

  /*! \brief Initializes iwhere and projects the initial x to
      the feasible set if necessary.
   */
  /*! This subroutine initializes iwhere and projects the initial x to
        the feasible set if necessary.

      iwhere is an integer array of dimension n.
        On entry iwhere is unspecified.
        On exit iwhere(i)=-1  if x(i) has no bounds
                          3   if l(i)=u(i)
                          0   otherwise.
        In cauchy, iwhere is given finer gradations.

                            *  *  *

      NEOS, November 1994. (Latest revision June 1996.)
      Optimization Technology Center.
      Argonne National Laboratory and Northwestern University.
      Written by
                         Ciyou Zhu
      in collaboration with R.H. Byrd, P. Lu-Chen and J. Nocedal.
   */
  template <typename FloatType>
  void
  active(
    int const& n,
    ref1<FloatType> const& l,
    ref1<FloatType> const& u,
    ref1<int> const& nbd,
    ref1<FloatType> const& x,
    ref1<int> const& iwhere,
    int const& iprint,
    bool& prjctd,
    bool& cnstnd,
    bool& boxed)
  {
    FloatType zero = 0;
    // Initialize nbdd, prjctd, cnstnd and boxed.
    int nbdd = 0;
    prjctd = false;
    cnstnd = false;
    boxed = true;
    // Project the initial x to the easible set if necessary.
    for(int i=1;i<=n;i++) {
      if (nbd(i) > 0) {
        if (nbd(i) <= 2 && x(i) <= l(i)) {
          if (x(i) < l(i)) {
            prjctd = true;
            x(i) = l(i);
          }
          nbdd = nbdd + 1;
        }
        else if (nbd(i) >= 2 && x(i) >= u(i)) {
          if (x(i) > u(i)) {
            prjctd = true;
            x(i) = u(i);
          }
          nbdd = nbdd + 1;
        }
      }
    }
    // Initialize iwhere and assign values to cnstnd and boxed.
    for(int i=1;i<=n;i++) {
      if (nbd(i) != 2) boxed = false;
      if (nbd(i) == 0) {
        // this variable is always free
        iwhere(i) = -1;
        // otherwise set x(i)=mid(x(i), u(i), l(i)).
      }
      else {
        cnstnd = true;
        if (nbd(i) == 2 && u(i) - l(i) <= zero) {
          // this variable is always fixed
          iwhere(i) = 3;
        }
        else {
          iwhere(i) = 0;
        }
      }
    }
    if (iprint >= 0) {
      if (prjctd) {
        printf(
          " The initial X is infeasible.  Restart with its projection.\n");
      }
      if (!cnstnd) {
        printf(" This problem is unconstrained.\n");
      }
    }
    if (iprint > 0) {
      printf("\nAt X0 %9d variables are exactly at the bounds\n", nbdd);
    }
  }

  //! Fragment from subroutine cauchy.
  template <typename FloatType>
  bool
  cauchy_loop(
    // arguments to cauchy
    int const& n,
    ref1<FloatType> const& x,
    ref1<FloatType> const& l,
    ref1<FloatType> const& u,
    ref1<int> const& iorder,
    ref1<int> const& iwhere,
    ref1<FloatType> const& t,
    ref1<FloatType> const& d,
    ref1<FloatType> const& xcp,
    int const& m,
    ref2<FloatType> const& wy,
    ref2<FloatType> const& ws,
    ref2<FloatType> const& sy,
    ref2<FloatType> const& wt,
    FloatType const& theta,
    int const& col,
    int const& head,
    ref1<FloatType> const& p,
    ref1<FloatType> const& c,
    ref1<FloatType> const& wbp,
    ref1<FloatType> const& v,
    int& nint,
    int const& iprint,
    int& info,
    FloatType const& epsmch,
    // local to cauchy
    FloatType const& bkmin,
    bool const& bnded,
    int const& col2,
    FloatType& dtm,
    FloatType& f1,
    FloatType& f2,
    FloatType& f2_org,
    int const& ibkmin,
    int const& nbreak,
    FloatType& tsum)
  {
    FloatType zero = 0;
    int nleft = nbreak;
    int iter = 1;
    FloatType tj = zero;
    //------------------ the beginning of the loop -------------------------
    while (true) { // lbl_777:
      // Find the next smallest breakpoint;
      // compute dt = t(nleft) - t(nleft + 1).
      FloatType tj0 = tj;
      int ibp = 0; // uninitialized
      if (iter == 1) {
        // Since we already have the smallest breakpoint we need not do
        // heapsort yet. Often only one breakpoint is used and the
        // cost of heapsort is avoided.
        tj = bkmin;
        ibp = iorder(ibkmin);
      }
      else {
        if (iter == 2) {
          // Replace the already used smallest breakpoint with the
          // breakpoint numbered nbreak > nlast, before heapsort call.
          if (ibkmin != nbreak) {
            t(ibkmin) = t(nbreak);
            iorder(ibkmin) = iorder(nbreak);
          }
          // Update heap structure of breakpoints
          // (if iter=2, initialize heap).
        }
        hpsolb(nleft,t.get1(1,nleft),iorder.get1(1,nleft),iter-2);
        tj = t(nleft);
        ibp = iorder(nleft);
      }
      FloatType dt = tj - tj0;
      if (dt != zero && iprint >= 100) {
        printf("\nPiece    %3d --f1, f2 at start point  %11.4E %11.4E\n",
               nint,f1,f2);
        printf("Distance to the next break point =  %11.4E\n", dt);
        printf("Distance to the stationary point =  %11.4E\n", dtm);
      }
      // If a minimizer is within this interval, locate the GCP and return.
      if (dtm < dt) return false; // goto lbl_888;
      // Otherwise fix one variable and
      // reset the corresponding component of d to zero.
      tsum = tsum + dt;
      nleft = nleft - 1;
      iter = iter + 1;
      FloatType dibp = d(ibp);
      d(ibp) = zero;
      FloatType zibp = 0; // uninitialized
      if (dibp > zero) {
        zibp = u(ibp) - x(ibp);
        xcp(ibp) = u(ibp);
        iwhere(ibp) = 2;
      }
      else {
        zibp = l(ibp) - x(ibp);
        xcp(ibp) = l(ibp);
        iwhere(ibp) = 1;
      }
      if (iprint >= 100) {
        printf(" Variable  %12d  is fixed.\n", ibp);
      }
      if (nleft == 0 && nbreak == n) {
        // all n variables are fixed,
        // return with xcp as GCP.
        dtm = dt;
        return true; // goto lbl_999;
      }
      // Update the derivative information.
      nint = nint + 1;
      FloatType dibp2 = dibp*dibp;
      // Update f1 and f2.
      // temporarily set f1 and f2 for col=0.
      f1 = f1 + dt*f2 + dibp2 - theta*dibp*zibp;
      f2 = f2 - theta*dibp2;
      if (col > 0) {
        // update c = c + dt*p.
        daxpy(col2,dt,p.begin(),c.begin());
        // choose wbp,
        // the row of W corresponding to the breakpoint encountered.
        int pointr = head;
        for(int j=1;j<=col;j++) {
          wbp(j) = wy(ibp,pointr);
          wbp(col + j) = theta*ws(ibp,pointr);
          pointr = pointr % m + 1;
        }
        // compute (wbp)Mc, (wbp)Mp, and (wbp)M(wbp)'.
        bmv(m,sy,wt,col,wbp.get1(1,2*col),v.get1(1,2*col),info);
        if (info != 0) return false; // goto lbl_888;
        FloatType wmc = lbfgs::detail::ddot(col2,c.begin(),v.begin());
        FloatType wmp = lbfgs::detail::ddot(col2,p.begin(),v.begin());
        FloatType wmw = lbfgs::detail::ddot(col2,wbp.begin(),v.begin());
        // update p = p - dibp*wbp.
        daxpy(col2,-dibp,wbp.begin(),p.begin());
        // complete updating f1 and f2 while col > 0.
        f1 = f1 + dibp*wmc;
        f2 = f2 + 2*dibp*wmp - dibp2*wmw;
      }
      f2 = std::max(epsmch*f2_org,f2);
      if (nleft > 0) {
        dtm = -f1/f2;
        // goto lbl_777;
        // to repeat the loop for unsearched intervals.
      }
      else if(bnded) {
        f1 = zero;
        f2 = zero;
        dtm = zero;
        break;
      }
      else {
        dtm = -f1/f2;
        break;
      }
    }
    //------------------ the end of the loop -------------------------------
    return false; // next line was lbl_888
  }

  /*! \brief Computation of the generalized Cauchy point (GCP).
   */
  /*! For given x, l, u, g (with sbgnrm > 0), and a limited memory
        BFGS matrix B defined in terms of matrices WY, WS, WT, and
        scalars head, col, and theta, this subroutine computes the
        generalized Cauchy point (GCP), defined as the first local
        minimizer of the quadratic

                   Q(x + s) = g's + 1/2 s'Bs

        along the projected gradient direction P(x-tg,l,u).
        The routine returns the GCP in xcp.

      n is an integer variable.
        On entry n is the dimension of the problem.
        On exit n is unchanged.

      x is a double precision array of dimension n.
        On entry x is the starting point for the GCP computation.
        On exit x is unchanged.

      l is a double precision array of dimension n.
        On entry l is the lower bound of x.
        On exit l is unchanged.

      u is a double precision array of dimension n.
        On entry u is the upper bound of x.
        On exit u is unchanged.

      nbd is an integer array of dimension n.
        On entry nbd represents the type of bounds imposed on the
          variables, and must be specified as follows:
          nbd(i)=0 if x(i) is unbounded,
                 1 if x(i) has only a lower bound,
                 2 if x(i) has both lower and upper bounds, and
                 3 if x(i) has only an upper bound.
        On exit nbd is unchanged.

      g is a double precision array of dimension n.
        On entry g is the gradient of f(x).  g must be a nonzero vector.
        On exit g is unchanged.

      iorder is an integer working array of dimension n.
        iorder will be used to store the breakpoints in the piecewise
        linear path and free variables encountered. On exit,
          iorder(1),...,iorder(nleft) are indices of breakpoints
                                 which have not been encountered;
          iorder(nleft+1),...,iorder(nbreak) are indices of
                                      encountered breakpoints; and
          iorder(nfree),...,iorder(n) are indices of variables which
                  have no bound constraits along the search direction.

      iwhere is an integer array of dimension n.
        On entry iwhere indicates only the permanently fixed (iwhere=3)
        or free (iwhere= -1) components of x.
        On exit iwhere records the status of the current x variables.
        iwhere(i)=-3  if x(i) is free and has bounds, but is not moved
                  0   if x(i) is free and has bounds, and is moved
                  1   if x(i) is fixed at l(i), and l(i) .ne. u(i)
                  2   if x(i) is fixed at u(i), and u(i) .ne. l(i)
                  3   if x(i) is always fixed, i.e.,  u(i)=x(i)=l(i)
                  -1  if x(i) is always free, i.e., it has no bounds.

      t is a double precision working array of dimension n.
        t will be used to store the break points.

      d is a double precision array of dimension n used to store
        the Cauchy direction P(x-tg)-x.

      xcp is a double precision array of dimension n used to return the
        GCP on exit.

      m is an integer variable.
        On entry m is the maximum number of variable metric corrections
          used to define the limited memory matrix.
        On exit m is unchanged.

      ws, wy, sy, and wt are double precision arrays.
        On entry they store information that defines the
                              limited memory BFGS matrix:
          ws(n,m) stores S, a set of s-vectors;
          wy(n,m) stores Y, a set of y-vectors;
          sy(m,m) stores S'Y;
          wt(m,m) stores the
                  Cholesky factorization of (theta*S'S+LD^(-1)L').
        On exit these arrays are unchanged.

      theta is a double precision variable.
        On entry theta is the scaling factor specifying B_0 = theta I.
        On exit theta is unchanged.

      col is an integer variable.
        On entry col is the actual number of variable metric
          corrections stored so far.
        On exit col is unchanged.

      head is an integer variable.
        On entry head is the location of the first s-vector (or y-vector)
          in S (or Y).
        On exit col is unchanged.

      p is a double precision working array of dimension 2m.
        p will be used to store the vector p = W^(T)d.

      c is a double precision working array of dimension 2m.
        c will be used to store the vector c = W^(T)(xcp-x).

      wbp is a double precision working array of dimension 2m.
        wbp will be used to store the row of W corresponding
          to a breakpoint.

      v is a double precision working array of dimension 2m.

      nint is an integer variable.
        On exit nint records the number of quadratic segments explored
          in searching for the GCP.

      sg and yg are double precision arrays of dimension m.
        On entry sg  and yg store S'g and Y'g correspondingly.
        On exit they are unchanged.

      iprint is an INTEGER variable that must be set by the user.
        It controls the frequency and type of output generated:
         iprint<0    no output is generated;
         iprint=0    print only one line at the last iteration;
         0<iprint<99 print also f and |proj g| every iprint iterations;
         iprint=99   print details of every iteration except n-vectors;
         iprint=100  print also the changes of active set and final x;
         iprint>100  print details of every iteration including x and g;
        When iprint > 0, the file iterate.dat will be created to
                         summarize the iteration.

      sbgnrm is a double precision variable.
        On entry sbgnrm is the norm of the projected gradient at x.
        On exit sbgnrm is unchanged.

      info is an integer variable.
        On entry info is 0.
        On exit info = 0       for normal return,
                     = nonzero for abnormal return when the the system
                               used in routine bmv is singular.

      Subprograms called:

        L-BFGS-B Library ... hpsolb, bmv.

        Linpack ... dscal dcopy, daxpy.

      References:

        [1] R. H. Byrd, P. Lu, J. Nocedal and C. Zhu, ``A limited
        memory algorithm for bound constrained optimization'',
        SIAM J. Scientific Computing 16 (1995), no. 5, pp. 1190--1208.

        [2] C. Zhu, R.H. Byrd, P. Lu, J. Nocedal, ``L-BFGS-B: FORTRAN
        Subroutines for Large Scale Bound Constrained Optimization''
        Tech. Report, NAM-11, EECS Department, Northwestern University,
        1994.

        (Postscript files of these papers are available via anonymous
         ftp to eecs.nwu.edu in the directory pub/lbfgs/lbfgs_bcm.)

                            *  *  *

      NEOS, November 1994. (Latest revision June 1996.)
      Optimization Technology Center.
      Argonne National Laboratory and Northwestern University.
      Written by
                         Ciyou Zhu
      in collaboration with R.H. Byrd, P. Lu-Chen and J. Nocedal.

   */
  template <typename FloatType>
  void
  cauchy(
    int const& n,
    ref1<FloatType> const& x,
    ref1<FloatType> const& l,
    ref1<FloatType> const& u,
    ref1<int> const& nbd,
    ref1<FloatType> const& g,
    ref1<int> const& iorder,
    ref1<int> const& iwhere,
    ref1<FloatType> const& t,
    ref1<FloatType> const& d,
    ref1<FloatType> const& xcp,
    int const& m,
    ref2<FloatType> const& wy,
    ref2<FloatType> const& ws,
    ref2<FloatType> const& sy,
    ref2<FloatType> const& wt,
    FloatType const& theta,
    int const& col,
    int const& head,
    ref1<FloatType> const& p,
    ref1<FloatType> const& c,
    ref1<FloatType> const& wbp,
    ref1<FloatType> const& v,
    int& nint,
    ref1<FloatType> const&
#if (SCITBX_LBFGSB_RAW_ASSERTION_FLAG != 0)
    sg
#endif
    ,
    ref1<FloatType> const&
#if (SCITBX_LBFGSB_RAW_ASSERTION_FLAG != 0)
    yg
#endif
    ,
    int const& iprint,
    FloatType const& sbgnrm,
    int& info,
    FloatType const& epsmch)
  {
#if (SCITBX_LBFGSB_RAW_ASSERTION_FLAG != 0)
    SCITBX_ASSERT(iorder.size() == n);
    SCITBX_ASSERT(iwhere.size() == n);
    SCITBX_ASSERT(t.size() == n);
    SCITBX_ASSERT(d.size() == n);
    SCITBX_ASSERT(xcp.size() == n);
    SCITBX_ASSERT(wy.size() == n*col);
    SCITBX_ASSERT(ws.size() == n*col);
    SCITBX_ASSERT(sy.size() == m*m);
    SCITBX_ASSERT(wt.size() == m*m);
    SCITBX_ASSERT(p.size() == 2*m);
    SCITBX_ASSERT(c.size() == 2*m);
    SCITBX_ASSERT(wbp.size() == 2*m);
    SCITBX_ASSERT(v.size() == 2*m);
    SCITBX_ASSERT(sg.size() == m);
    SCITBX_ASSERT(yg.size() == m);
#endif
    FloatType one = 1;
    FloatType zero = 0;
    // Check the status of the variables, reset iwhere(i) if necessary;
    // the Cauchy direction d and the breakpoints t; initialize
    // derivative f1 and the vector p = W'd (for theta = 1).
    if (sbgnrm <= zero) {
      if (iprint >= 0) {
        printf(" Subgnorm = 0.  GCP = X.\n");
      }
      dcopy(n,x,1,xcp,1);
      return;
    }
    bool bnded = true;
    int nfree = n + 1;
    int nbreak = 0;
    int ibkmin = 0;
    FloatType bkmin = zero;
    int col2 = 2*col;
    FloatType f1 = zero;
    if (iprint >= 99) {
      printf("\n---------------- CAUCHY entered-------------------\n");
    }
    // We set p to zero and build it up as we determine d.
    for(int i=1;i<=col2;i++) {
      p(i) = zero;
    }
    // In the following loop we determine for each variable its bound
    // status and its breakpoint, and update p accordingly.
    // Smallest breakpoint is identified.
    for(int i=1;i<=n;i++) {
      FloatType neggi = -g(i);
      FloatType tl = 0; // uninitialized
      FloatType tu = 0; // uninitialized
      if (iwhere(i) != 3 && iwhere(i) != -1) {
        // if x(i) is not a constant and has bounds,
        // compute the difference between x(i) and its bounds.
        if (nbd(i) <= 2) tl = x(i) - l(i);
        if (nbd(i) >= 2) tu = u(i) - x(i);

        // If a variable is close enough to a bound
        // we treat it as at bound.
        bool xlower = nbd(i) <= 2 && tl <= zero;
        bool xupper = nbd(i) >= 2 && tu <= zero;

        // reset iwhere(i).
        iwhere(i) = 0;
        if (xlower) {
          if (neggi <= zero) iwhere(i) = 1;
        }
        else if (xupper) {
          if (neggi >= zero) iwhere(i) = 2;
        }
        else {
          if (fn::absolute(neggi) <= zero) iwhere(i) = -3;
        }
      }
      int pointr = head;
      if (iwhere(i) != 0 && iwhere(i) != -1) {
        d(i) = zero;
      }
      else {
        d(i) = neggi;
        f1 = f1 - neggi*neggi;
        // calculate p := p - W'e_i* (g_i).
        for(int j=1;j<=col;j++) {
          p(j) = p(j) +  wy(i,pointr)* neggi;
          p(col + j) = p(col + j) + ws(i,pointr)*neggi;
          pointr = pointr % m + 1;
        }
        if (nbd(i) <= 2 && nbd(i) != 0 && neggi < zero) {
          // x(i) + d(i) is bounded; compute t(i).
          nbreak = nbreak + 1;
          iorder(nbreak) = i;
          t(nbreak) = tl/(-neggi);
          if (nbreak == 1 || t(nbreak) < bkmin) {
            bkmin = t(nbreak);
            ibkmin = nbreak;
          }
        }
        else if (nbd(i) >= 2 && neggi > zero) {
          // x(i) + d(i) is bounded; compute t(i).
          nbreak = nbreak + 1;
          iorder(nbreak) = i;
          t(nbreak) = tu/neggi;
          if (nbreak == 1 || t(nbreak) < bkmin) {
            bkmin = t(nbreak);
            ibkmin = nbreak;
          }
        }
        else {
          // x(i) + d(i) is not bounded.
          nfree = nfree - 1;
          iorder(nfree) = i;
          if (fn::absolute(neggi) > zero) bnded = false;
        }
      }
    }
    // The indices of the nonzero components of d are now stored
    // in iorder(1),...,iorder(nbreak) and iorder(nfree),...,iorder(n).
    // The smallest of the nbreak breakpoints is in t(ibkmin)=bkmin.
    if (theta != one) {
      // complete the initialization of p for theta not= one.
      dscal(col,theta,p.get1(col+1,col),1);
    }
    // Initialize GCP xcp = x.
    dcopy(n,x,1,xcp,1);
    if (nbreak == 0 && nfree == n + 1) {
      // is a zero vector, return with the initial xcp as GCP.
      if (iprint > 100) {
        printf("Cauchy X =  ");
        write_ref1("    ", xcp);
      }
      return;
    }
    // Initialize c = W'(xcp - x) = 0.
    for(int j=1;j<=col2;j++) {
      c(j) = zero;
    }
    // Initialize derivative f2.
    FloatType f2 = -theta*f1;
    FloatType f2_org = f2;
    if (col > 0) {
      bmv(m,sy,wt,col,p.get1(1,2*col),v.get1(1,2*col),info);
      if (info != 0) return;
      f2 = f2 - lbfgs::detail::ddot(col2,v.begin(),p.begin());
    }
    FloatType dtm = -f1/f2;
    FloatType tsum = zero;
    nint = 1;
    if (iprint >= 99) {
      printf(" There are %12d  breakpoints \n", nbreak);
    }
    // If there are no breakpoints, locate the GCP and return.
    bool skip_lbl_888 = false;
    if (nbreak != 0) {
      skip_lbl_888 = cauchy_loop(
        n, x, l, u, iorder, iwhere, t, d, xcp,
        m, wy, ws, sy, wt, theta, col, head, p, c, wbp,
        v, nint, iprint, info, epsmch,
        bkmin, bnded, col2, dtm, f1, f2, f2_org, ibkmin, nbreak, tsum);
      // goto lbl_999;
    }
    if (!skip_lbl_888) {
      // lbl_888:
      if (iprint >= 99) {
        printf(" \n");
        printf(" GCP found in this segment\n");
        printf("Piece    %3d --f1, f2 at start point  %11.4E %11.4E\n",
               nint,f1,f2);
        printf("Distance to the stationary point =  %11.4E\n", dtm);
      }
      if (dtm <= zero) dtm = zero;
      tsum = tsum + dtm;
      // Move free variables (i.e., the ones w/o breakpoints) and
      // the variables whose breakpoints haven't been reached.
      daxpy(n,tsum,d.begin(),xcp.begin());
    }
    // lbl_999:
    // Update c = c + dtm*p = W'(x^c - x)
    // which will be used in computing r = Z'(B(x^c - x) + g).
    if (col > 0) daxpy(col2,dtm,p.begin(),c.begin());
    if (iprint > 100) {
      printf("Cauchy X =  ");
      write_ref1("    ", xcp);
    }
    if (iprint >= 99) {
      printf("\n---------------- exit CAUCHY----------------------\n\n");
    }
  }

  /*! \brief Determination of the index set of free and active variables
      at the GCP.
   */
  /*! This subroutine counts the entering and leaving variables when
        iter > 0, and finds the index set of free and active variables
        at the GCP.

      cnstnd is a logical variable indicating whether bounds are present

      index is an integer array of dimension n
        for i=1,...,nfree, index(i) are the indices of free variables
        for i=nfree+1,...,n, index(i) are the indices of bound variables
        On entry after the first iteration, index gives
          the free variables at the previous iteration.
        On exit it gives the free variables based on the determination
          in cauchy using the array iwhere.

      indx2 is an integer array of dimension n
        On entry indx2 is unspecified.
        On exit with iter>0, indx2 indicates which variables
           have changed status since the previous iteration.
        For i= 1,...,nenter, indx2(i) have changed from bound to free.
        For i= ileave+1,...,n, indx2(i) have changed from free to bound.

                            *  *  *

      NEOS, November 1994. (Latest revision June 1996.)
      Optimization Technology Center.
      Argonne National Laboratory and Northwestern University.
      Written by
                         Ciyou Zhu
      in collaboration with R.H. Byrd, P. Lu-Chen and J. Nocedal.
   */
  inline
  void
  freev(
    int const& n,
    int& nfree,
    ref1<int> const& index,
    int& nenter,
    int& ileave,
    ref1<int> const& indx2,
    ref1<int> const& iwhere,
    bool& wrk,
    bool const& updatd,
    bool const& cnstnd,
    int const& iprint,
    int const& iter)
  {
#if (SCITBX_LBFGSB_RAW_ASSERTION_FLAG != 0)
    SCITBX_ASSERT(index.size() == n);
    SCITBX_ASSERT(indx2.size() == n);
    SCITBX_ASSERT(iwhere.size() == n);
#endif
    nenter = 0;
    ileave = n + 1;
    if (iter > 0 && cnstnd) {
      // count the entering and leaving variables.
      for(int i=1;i<=nfree;i++) {
        int k = index(i);
        if (iwhere(k) > 0) {
          ileave = ileave - 1;
          indx2(ileave) = k;
          if (iprint >= 100) {
            printf(" Variable %12d leaves the set of free variables\n", k);
          }
        }
      }
      for(int i=1+nfree;i<=n;i++) {
        int k = index(i);
        if (iwhere(k) <= 0) {
          nenter = nenter + 1;
          indx2(nenter) = k;
          if (iprint >= 100) {
            printf(" Variable %12d enters the set of free variables\n", k);
          }
        }
      }
      if (iprint >= 99) {
        printf("%12d variables leave; %12d variables enter\n",
               n+1-ileave, nenter);
      }
    }
    wrk = (ileave < n+1) || (nenter > 0) || updatd;
    // Find the index set of free and active variables at the GCP.
    nfree = 0;
    int iact = n + 1;
    for(int i=1;i<=n;i++) {
      if (iwhere(i) <= 0) {
        nfree = nfree + 1;
        index(nfree) = i;
      }
      else {
        iact = iact - 1;
        index(iact) = i;
      }
    }
    if (iprint >= 99) {
      printf("%12d variables are free at GCP %12d\n",
             nfree, iter + 1);
    }
  }

  /*! \brief Factorization of a symmetric positive definite matrix.
   */
  /*! dpofa factors a double precision symmetric positive definite
      matrix.

      dpofa is usually called by dpoco, but it can be called
      directly with a saving in time if  rcond  is not needed.
      (time for dpoco) = (1 + 18/n)*(time for dpofa) .

      on entry

         a       double precision(lda, n)
                 the symmetric matrix to be factored.  only the
                 diagonal and upper triangle are used.

         lda     integer
                 the leading dimension of the array  a .

         n       integer
                 the order of the matrix  a .

      on return

         a       an upper triangular matrix  r  so that  a = trans(r)*r
                 where  trans(r)  is the transpose.
                 the strict lower triangle is unaltered.
                 if  info .ne. 0 , the factorization is not complete.

         info    integer
                 = 0  for normal return.
                 = k  signals an error condition.  the leading minor
                      of order  k  is not positive definite.

      linpack.  this version dated 08/14/78 .
      cleve moler, university of new mexico, argonne national lab.

      subroutines and functions

      blas ddot
      fortran sqrt
   */
  template <typename FloatType>
  void
  dpofa(
    ref2<FloatType> const& a,
    int const&
#if (SCITBX_LBFGSB_RAW_ASSERTION_FLAG != 0)
    lda
#endif
    ,
    int const& n,
    int& info)
  {
#if (SCITBX_LBFGSB_RAW_ASSERTION_FLAG != 0)
    SCITBX_ASSERT(a.size() == lda*n);
#endif
    for(int j=1;j<=n;j++) {
      info = j;
      FloatType s = 0;
      int jm1 = j - 1;
      if (jm1 < 1) goto lbl_20;
      for(int k=1;k<=jm1;k++) {
        FloatType t = a(k,j) - lbfgs::detail::ddot(k-1,&a(1,k),&a(1,j));
        t = t/a(k,k);
        a(k,j) = t;
        s = s + t*t;
      }
      lbl_20:
      s = a(j,j) - s;
      if (s <= 0) return;
      a(j,j) = std::sqrt(s);
    }
    info = 0;
  }

  /*! \brief LEL^T factorization.
   */
  /*! This subroutine forms  the LEL^T factorization of the indefinite

        matrix    K = [-D -Y'ZZ'Y/theta     L_a'-R_z'  ]
                      [L_a -R_z           theta*S'AA'S ]
                                                     where E = [-I  0]
                                                               [ 0  I]
      The matrix K can be shown to be equal to the matrix M^[-1]N
        occurring in section 5.1 of [1], as well as to the matrix
        Mbar^[-1] Nbar in section 5.3.

      n is an integer variable.
        On entry n is the dimension of the problem.
        On exit n is unchanged.

      nsub is an integer variable
        On entry nsub is the number of subspace variables in free set.
        On exit nsub is not changed.

      ind is an integer array of dimension nsub.
        On entry ind specifies the indices of subspace variables.
        On exit ind is unchanged.

      nenter is an integer variable.
        On entry nenter is the number of variables entering the
          free set.
        On exit nenter is unchanged.

      ileave is an integer variable.
        On entry indx2(ileave),...,indx2(n) are the variables leaving
          the free set.
        On exit ileave is unchanged.

      indx2 is an integer array of dimension n.
        On entry indx2(1),...,indx2(nenter) are the variables entering
          the free set, while indx2(ileave),...,indx2(n) are the
          variables leaving the free set.
        On exit indx2 is unchanged.

      iupdat is an integer variable.
        On entry iupdat is the total number of BFGS updates made so far.
        On exit iupdat is unchanged.

      updatd is a logical variable.
        On entry 'updatd' is true if the L-BFGS matrix is updatd.
        On exit 'updatd' is unchanged.

      wn is a double precision array of dimension 2m x 2m.
        On entry wn is unspecified.
        On exit the upper triangle of wn stores the LEL^T factorization
          of the 2*col x 2*col indefinite matrix
                      [-D -Y'ZZ'Y/theta     L_a'-R_z'  ]
                      [L_a -R_z           theta*S'AA'S ]

      wn1 is a double precision array of dimension 2m x 2m.
        On entry wn1 stores the lower triangular part of
                      [Y' ZZ'Y   L_a'+R_z']
                      [L_a+R_z   S'AA'S   ]
          in the previous iteration.
        On exit wn1 stores the corresponding updated matrices.
        The purpose of wn1 is just to store these inner products
        so they can be easily updated and inserted into wn.

      m is an integer variable.
        On entry m is the maximum number of variable metric corrections
          used to define the limited memory matrix.
        On exit m is unchanged.

      ws, wy, sy, and wtyy are double precision arrays;
      theta is a double precision variable;
      col is an integer variable;
      head is an integer variable.
        On entry they store the information defining the
                                           limited memory BFGS matrix:
          ws(n,m) stores S, a set of s-vectors;
          wy(n,m) stores Y, a set of y-vectors;
          sy(m,m) stores S'Y;
          wtyy(m,m) stores the Cholesky factorization
                                    of (theta*S'S+LD^(-1)L')
          theta is the scaling factor specifying B_0 = theta I;
          col is the number of variable metric corrections stored;
          head is the location of the 1st s- (or y-) vector in S (or Y).
        On exit they are unchanged.

      info is an integer variable.
        On entry info is unspecified.
        On exit info =  0 for normal return;
                     = -1 when the 1st Cholesky factorization failed;
                     = -2 when the 2st Cholesky factorization failed.

      Subprograms called:

        Linpack ... dcopy, dpofa, dtrsl.

      References:
        [1] R. H. Byrd, P. Lu, J. Nocedal and C. Zhu, ``A limited
        memory algorithm for bound constrained optimization'',
        SIAM J. Scientific Computing 16 (1995), no. 5, pp. 1190--1208.

        [2] C. Zhu, R.H. Byrd, P. Lu, J. Nocedal, ``L-BFGS-B: a
        limited memory FORTRAN code for solving bound constrained
        optimization problems'', Tech. Report, NAM-11, EECS Department,
        Northwestern University, 1994.

        (Postscript files of these papers are available via anonymous
         ftp to eecs.nwu.edu in the directory pub/lbfgs/lbfgs_bcm.)

                            *  *  *

      NEOS, November 1994. (Latest revision June 1996.)
      Optimization Technology Center.
      Argonne National Laboratory and Northwestern University.
      Written by
                         Ciyou Zhu
      in collaboration with R.H. Byrd, P. Lu-Chen and J. Nocedal.
   */
  template <typename FloatType>
  void
  formk(
    int const& n,
    int const& nsub,
    ref1<int> const& ind,
    int const& nenter,
    int const& ileave,
    ref1<int> const& indx2,
    int const& iupdat,
    bool const& updatd,
    ref2<FloatType> const& wn,
    ref2<FloatType> const& wn1,
    int const& m,
    ref2<FloatType> const& ws,
    ref2<FloatType> const& wy,
    ref2<FloatType> const& sy,
    FloatType const& theta,
    int const& col,
    int const& head,
    int& info)
  {
#if (SCITBX_LBFGSB_RAW_ASSERTION_FLAG != 0)
    SCITBX_ASSERT(ind.size() == n);
    SCITBX_ASSERT(indx2.size() == n);
    SCITBX_ASSERT(wn.size() == 2*m*2*m);
    SCITBX_ASSERT(wn1.size() == 2*m*2*m);
    SCITBX_ASSERT(ws.size() == n*m);
    SCITBX_ASSERT(wy.size() == n*m);
    SCITBX_ASSERT(sy.size() == m*m);
#endif
    // Form the lower triangular part of
    //   WN1 = [Y' ZZ'Y   L_a'+R_z']
    //         [L_a+R_z   S'AA'S   ]
    // where L_a is the strictly lower triangular part of S'AA'Y
    //   R_z is the upper triangular part of S'ZZ'Y.
    int upcl = 0; // uninitialized
    if (updatd) {
      if (iupdat > m) {
        // shift old part of WN1.
        for(int jy=1;jy<=m-1;jy++) {
          int js = m + jy;
          dcopy(m-jy,wn1.get1(jy+1,jy+1,m-jy),1,wn1.get1(jy,jy,m-jy),1);
          dcopy(m-jy,wn1.get1(js+1,js+1,m-jy),1,wn1.get1(js,js,m-jy),1);
          dcopy(m-1,wn1.get1(m+2,jy+1,m-1),1,wn1.get1(m+1,jy,m-1),1);
        }
      }
      // put new rows in blocks (1,1), (2,1) and (2,2).
      int pbegin = 1;
      int pend = nsub;
      int dbegin = nsub + 1;
      int dend = n;
      int iy = col;
      int is = m + col;
      int ipntr = head + col - 1;
      if (ipntr > m) ipntr = ipntr - m;
      int jpntr = head;
      int jy;
      for(jy=1;jy<=col;jy++) {
        int js = m + jy;
        FloatType temp1 = 0;
        FloatType temp2 = 0;
        FloatType temp3 = 0;
        // compute element jy of row 'col' of Y'ZZ'Y
        for(int k=pbegin;k<=pend;k++) {
          int k1 = ind(k);
          temp1 = temp1 + wy(k1,ipntr)*wy(k1,jpntr);
        }
        // compute elements jy of row 'col' of L_a and S'AA'S
        for(int k=dbegin;k<=dend;k++) {
          int k1 = ind(k);
          temp2 = temp2 + ws(k1,ipntr)*ws(k1,jpntr);
          temp3 = temp3 + ws(k1,ipntr)*wy(k1,jpntr);
        }
        wn1(iy,jy) = temp1;
        wn1(is,js) = temp2;
        wn1(is,jy) = temp3;
        jpntr = jpntr % m + 1;
      }
      // put new column in block (2,1).
      jy = col;
      jpntr = head + col - 1;
      if (jpntr > m) jpntr = jpntr - m;
      ipntr = head;
      for(int i=1;i<=col;i++) {
        is = m + i;
        FloatType temp3 = 0;
        // compute element i of column 'col' of R_z
        for(int k=pbegin;k<=pend;k++) {
          int k1 = ind(k);
          temp3 = temp3 + ws(k1,ipntr)*wy(k1,jpntr);
        }
        ipntr = ipntr % m + 1;
        wn1(is,jy) = temp3;
      }
      upcl = col - 1;
    }
    else {
      upcl = col;
    }
    // modify the old parts in blocks (1,1) and (2,2) due to changes
    // in the set of free variables.
    int ipntr = head;
    for(int iy=1;iy<=upcl;iy++) {
      int is = m + iy;
      int jpntr = head;
      for(int jy=1;jy<=iy;jy++) {
        int js = m + jy;
        FloatType temp1 = 0;
        FloatType temp2 = 0;
        FloatType temp3 = 0;
        FloatType temp4 = 0;
        for(int k=1;k<=nenter;k++) {
          int k1 = indx2(k);
          temp1 = temp1 + wy(k1,ipntr)*wy(k1,jpntr);
          temp2 = temp2 + ws(k1,ipntr)*ws(k1,jpntr);
        }
        for(int k=ileave;k<=n;k++) {
          int k1 = indx2(k);
          temp3 = temp3 + wy(k1,ipntr)*wy(k1,jpntr);
          temp4 = temp4 + ws(k1,ipntr)*ws(k1,jpntr);
        }
        wn1(iy,jy) = wn1(iy,jy) + temp1 - temp3;
        wn1(is,js) = wn1(is,js) - temp2 + temp4;
        jpntr = jpntr % m + 1;
      }
      ipntr = ipntr % m + 1;
    }
    // modify the old parts in block (2,1).
    ipntr = head;
    for(int is=m+1;is<=m+upcl;is++) {
      int jpntr = head;
      for(int jy=1;jy<=upcl;jy++) {
        FloatType temp1 = 0;
        FloatType temp3 = 0;
        for(int k=1;k<=nenter;k++) {
          int k1 = indx2(k);
          temp1 = temp1 + ws(k1,ipntr)*wy(k1,jpntr);
        }
        for(int k=ileave;k<=n;k++) {
          int k1 = indx2(k);
          temp3 = temp3 + ws(k1,ipntr)*wy(k1,jpntr);
        }
        if (is <= jy + m) {
          wn1(is,jy) = wn1(is,jy) + temp1 - temp3;
        }
        else {
          wn1(is,jy) = wn1(is,jy) - temp1 + temp3;
        }
        jpntr = jpntr % m + 1;
      }
      ipntr = ipntr % m + 1;
    }
    // Form the upper triangle of WN = [D+Y' ZZ'Y/theta   -L_a'+R_z' ]
    //                                 [-L_a +R_z        S'AA'S*theta]
    int m2 = 2*m;
    for(int iy=1;iy<=col;iy++) {
      int is = col + iy;
      int is1 = m + iy;
      for(int jy=1;jy<=iy;jy++) {
        int js = col + jy;
        int js1 = m + jy;
        wn(jy,iy) = wn1(iy,jy)/theta;
        wn(js,is) = wn1(is1,js1)*theta;
      }
      for(int jy=1;jy<=iy-1;jy++) {
        wn(jy,is) = -wn1(is1,jy);
      }
      for(int jy=iy;jy<=col;jy++) {
        wn(jy,is) = wn1(is1,jy);
      }
      wn(iy,iy) = wn(iy,iy) + sy(iy,iy);
    }
    // Form the upper triangle of WN= [  LL'            L^-1(-L_a'+R_z')]
    //                                [(-L_a +R_z)L'^-1   S'AA'S*theta  ]
    // first Cholesky factor (1,1) block of wn to get LL'
    // with L' stored in the upper triangle of wn.
    dpofa(wn.get2(1,1,m2,col),m2,col,info);
    if (info != 0) {
      info = -1;
      return;
    }
    // then form L^-1(-L_a'+R_z') in the (1,2) block.
    int col2 = 2*col;
    for(int js=col+1;js<=col2;js++) {
      dtrsl(wn.get2(1,1,m2,col),m2,col,wn.get1(1,js,col),11,info);
    }
    // Form S'AA'S*theta + (L^-1(-L_a'+R_z'))'L^-1(-L_a'+R_z') in the
    // upper triangle of (2,2) block of wn.
    for(int is=col+1;is<=col2;is++) {
      for(int js=is;js<=col2;js++) {
        wn(is,js) = wn(is,js) + lbfgs::detail::ddot(col,&wn(1,is),&wn(1,js));
      }
    }
    // Cholesky factorization of (2,2) block of wn.
    dpofa(wn.get2soft(col+1,col+1,m2,col),m2,col,info);
    if (info != 0) {
      info = -2;
      return;
    }
  }

  /*! \brief Computation of -Z'B(xcp-xk)-Z'g.
   */
  /*! This subroutine computes r=-Z'B(xcp-xk)-Z'g by using
        wa(2m+1)=W'(xcp-x) from subroutine cauchy.

      Subprograms called:

        L-BFGS-B Library ... bmv.

                            *  *  *

      NEOS, November 1994. (Latest revision June 1996.)
      Optimization Technology Center.
      Argonne National Laboratory and Northwestern University.
      Written by
                         Ciyou Zhu
      in collaboration with R.H. Byrd, P. Lu-Chen and J. Nocedal.
   */
  template <typename FloatType>
  void
  cmprlb(
    int const& n,
    int const& m,
    ref1<FloatType> const& x,
    ref1<FloatType> const& g,
    ref2<FloatType> const& ws,
    ref2<FloatType> const& wy,
    ref2<FloatType> const& sy,
    ref2<FloatType> const& wt,
    ref1<FloatType> const& z,
    ref1<FloatType> const& r,
    ref1<FloatType> const& wa,
    ref1<int> const& index,
    FloatType const& theta,
    int const& col,
    int const& head,
    int const& nfree,
    bool const& cnstnd,
    int& info)
  {
#if (SCITBX_LBFGSB_RAW_ASSERTION_FLAG != 0)
    SCITBX_ASSERT(ws.size() == n*m);
    SCITBX_ASSERT(wy.size() == n*m);
    SCITBX_ASSERT(sy.size() == m*m);
    SCITBX_ASSERT(wt.size() == m*m);
    SCITBX_ASSERT(z.size() == n);
    SCITBX_ASSERT(r.size() == n);
    SCITBX_ASSERT(wa.size() == 4*m);
    SCITBX_ASSERT(index.size() == n);
#endif
    if (!cnstnd && col > 0) {
      for(int i=1;i<=n;i++) {
        r(i) = -g(i);
      }
    }
    else {
      for(int i=1;i<=nfree;i++) {
        int k = index(i);
        r(i) = -theta*(z(k) - x(k)) - g(k);
      }
      bmv(m,sy,wt,col,wa.get1(2*m+1,2*col),wa.get1(1,2*col),info);
      if (info != 0) {
        info = -8;
        return;
      }
      int pointr = head;
      for(int j=1;j<=col;j++) {
        FloatType a1 = wa(j);
        FloatType a2 = theta*wa(col + j);
        for(int i=1;i<=nfree;i++) {
          int k = index(i);
          r(i) = r(i) + wy(k,pointr)*a1 + ws(k,pointr)*a2;
        }
        pointr = pointr % m + 1;
      }
    }
  }

  /*! \brief Approximate solution of the subspace problem.
   */
  /*! Given xcp, l, u, r, an index set that specifies
        the active set at xcp, and an l-BFGS matrix B
        (in terms of WY, WS, SY, WT, head, col, and theta),
        this subroutine computes an approximate solution
        of the subspace problem

        (P)   min Q(x) = r'(x-xcp) + 1/2 (x-xcp)' B (x-xcp)

              subject to l<=x<=u
                        x_i=xcp_i for all i in A(xcp)

        along the subspace unconstrained Newton direction

           d = -(Z'BZ)^(-1) r.

        The formula for the Newton direction, given the L-BFGS matrix
        and the Sherman-Morrison formula, is

           d = (1/theta)r + (1/theta*2) Z'WK^(-1)W'Z r.

        where
                  K = [-D -Y'ZZ'Y/theta     L_a'-R_z'  ]
                      [L_a -R_z           theta*S'AA'S ]

      Note that this procedure for computing d differs
      from that described in [1]. One can show that the matrix K is
      equal to the matrix M^[-1]N in that paper.

      n is an integer variable.
        On entry n is the dimension of the problem.
        On exit n is unchanged.

      m is an integer variable.
        On entry m is the maximum number of variable metric corrections
          used to define the limited memory matrix.
        On exit m is unchanged.

      nsub is an integer variable.
        On entry nsub is the number of free variables.
        On exit nsub is unchanged.

      ind is an integer array of dimension nsub.
        On entry ind specifies the coordinate indices of free variables.
        On exit ind is unchanged.

      l is a double precision array of dimension n.
        On entry l is the lower bound of x.
        On exit l is unchanged.

      u is a double precision array of dimension n.
        On entry u is the upper bound of x.
        On exit u is unchanged.

      nbd is a integer array of dimension n.
        On entry nbd represents the type of bounds imposed on the
          variables, and must be specified as follows:
          nbd(i)=0 if x(i) is unbounded,
                 1 if x(i) has only a lower bound,
                 2 if x(i) has both lower and upper bounds, and
                 3 if x(i) has only an upper bound.
        On exit nbd is unchanged.

      x is a double precision array of dimension n.
        On entry x specifies the Cauchy point xcp.
        On exit x(i) is the minimizer of Q over the subspace of
                                                         free variables.

      d is a double precision array of dimension n.
        On entry d is the reduced gradient of Q at xcp.
        On exit d is the Newton direction of Q.

      ws and wy are double precision arrays;
      theta is a double precision variable;
      col is an integer variable;
      head is an integer variable.
        On entry they store the information defining the
                                           limited memory BFGS matrix:
          ws(n,m) stores S, a set of s-vectors;
          wy(n,m) stores Y, a set of y-vectors;
          theta is the scaling factor specifying B_0 = theta I;
          col is the number of variable metric corrections stored;
          head is the location of the 1st s- (or y-) vector in S (or Y).
        On exit they are unchanged.

      iword is an integer variable.
        On entry iword is unspecified.
        On exit iword specifies the status of the subspace solution.
          iword = 0 if the solution is in the box,
                  1 if some bound is encountered.

      wv is a double precision working array of dimension 2m.

      wn is a double precision array of dimension 2m x 2m.
        On entry the upper triangle of wn stores the LEL^T factorization
          of the indefinite matrix

               K = [-D -Y'ZZ'Y/theta     L_a'-R_z'  ]
                   [L_a -R_z           theta*S'AA'S ]
                                                     where E = [-I  0]
                                                               [ 0  I]
        On exit wn is unchanged.

      iprint is an INTEGER variable that must be set by the user.
        It controls the frequency and type of output generated:
         iprint<0    no output is generated;
         iprint=0    print only one line at the last iteration;
         0<iprint<99 print also f and |proj g| every iprint iterations;
         iprint=99   print details of every iteration except n-vectors;
         iprint=100  print also the changes of active set and final x;
         iprint>100  print details of every iteration including x and g;
        When iprint > 0, the file iterate.dat will be created to
                         summarize the iteration.

      info is an integer variable.
        On entry info is unspecified.
        On exit info = 0       for normal return,
                     = nonzero for abnormal return
                                   when the matrix K is ill-conditioned.

      Subprograms called:

        Linpack dtrsl.

      References:

        [1] R. H. Byrd, P. Lu, J. Nocedal and C. Zhu, ``A limited
        memory algorithm for bound constrained optimization'',
        SIAM J. Scientific Computing 16 (1995), no. 5, pp. 1190--1208.

                            *  *  *

      NEOS, November 1994. (Latest revision June 1996.)
      Optimization Technology Center.
      Argonne National Laboratory and Northwestern University.
      Written by
                         Ciyou Zhu
      in collaboration with R.H. Byrd, P. Lu-Chen and J. Nocedal.
   */
  template <typename FloatType>
  void
  subsm(
    int const&
#if (SCITBX_LBFGSB_RAW_ASSERTION_FLAG != 0)
    n
#endif
    ,
    int const& m,
    int const& nsub,
    ref1<int> const& ind,
    ref1<FloatType> const& l,
    ref1<FloatType> const& u,
    ref1<int> const& nbd,
    ref1<FloatType> const& x,
    ref1<FloatType> const& d,
    ref2<FloatType> const& ws,
    ref2<FloatType> const& wy,
    FloatType const& theta,
    int const& col,
    int const& head,
    int& iword,
    ref1<FloatType> const& wv,
    ref2<FloatType> const& wn,
    int const& iprint,
    int& info)
  {
#if (SCITBX_LBFGSB_RAW_ASSERTION_FLAG != 0)
    SCITBX_ASSERT(ind.size() == nsub);
    SCITBX_ASSERT(d.size() == n);
    SCITBX_ASSERT(ws.size() == n*m);
    SCITBX_ASSERT(wy.size() == n*m);
    SCITBX_ASSERT(wv.size() == 2*m);
    SCITBX_ASSERT(wn.size() == 2*m*2*m);
#endif
    if (nsub <= 0) return;
    if (iprint >= 99) {
      printf("\n----------------SUBSM entered-----------------\n\n");
    }
    // Compute wv = W'Zd.
    FloatType zero = 0;
    FloatType one = 1;
    int pointr = head;
    for(int i=1;i<=col;i++) {
      FloatType temp1 = 0;
      FloatType temp2 = 0;
      for(int j=1;j<=nsub;j++) {
        int k = ind(j);
        temp1 = temp1 + wy(k,pointr)*d(j);
        temp2 = temp2 + ws(k,pointr)*d(j);
      }
      wv(i) = temp1;
      wv(col + i) = theta*temp2;
      pointr = pointr % m + 1;
    }
    // Compute wv:=K^(-1)wv.
    int m2 = 2*m;
    int col2 = 2*col;
    dtrsl(wn.get2(1,1,m2,col2),m2,col2,wv.get1(1,col2),11,info);
    if (info != 0) return;
    for(int i=1;i<=col;i++) {
      wv(i) = -wv(i);
    }
    dtrsl(wn.get2(1,1,m2,col2),m2,col2,wv.get1(1,col2),01,info);
    if (info != 0) return;
    // Compute d = (1/theta)d + (1/theta**2)Z'W wv.
    pointr = head;
    for(int jy=1;jy<=col;jy++) {
      int js = col + jy;
      for(int i=1;i<=nsub;i++) {
        int k = ind(i);
        d(i) = d(i) + wy(k,pointr)*wv(jy)/theta
          + ws(k,pointr)*wv(js);
      }
      pointr = pointr % m + 1;
    }
    for(int i=1;i<=nsub;i++) {
      d(i) = d(i)/theta;
    }
    // Backtrack to the feasible region.
    FloatType alpha = one;
    FloatType temp1 = alpha;
    int ibd = 0; // uninitialized
    for(int i=1;i<=nsub;i++) {
      int k = ind(i);
      FloatType dk = d(i);
      if (nbd(k) != 0) {
        bool temp1_updated = false;
        if (dk < zero && nbd(k) <= 2) {
          FloatType temp2 = l(k) - x(k);
          if (temp2 >= zero) {
            temp1 = zero;
            temp1_updated = true;
          }
          else if (dk*alpha < temp2) {
            temp1 = temp2/dk;
            temp1_updated = true;
          }
        }
        else if (dk > zero && nbd(k) >= 2) {
          FloatType temp2 = u(k) - x(k);
          if (temp2 <= zero) {
            temp1 = zero;
            temp1_updated = true;
          }
          else if (dk*alpha > temp2) {
            temp1 = temp2/dk;
            temp1_updated = true;
          }
        }
        if (temp1_updated && temp1 < alpha) {
          alpha = temp1;
          ibd = i;
        }
      }
    }
    if (alpha < one) {
      FloatType dk = d(ibd);
      int k = ind(ibd);
      if (dk > zero) {
        x(k) = u(k);
        d(ibd) = zero;
      }
      else if (dk < zero) {
        x(k) = l(k);
        d(ibd) = zero;
      }
    }
    for(int i=1;i<=nsub;i++) {
      int k = ind(i);
      x(k) = x(k) + alpha*d(i);
    }
    if (iprint >= 99) {
      if (alpha < one) {
        printf("ALPHA = %7.5f backtrack to the BOX\n", alpha);
      }
      else {
        printf(" SM solution inside the box\n");
      }
      if (iprint >100) {
        printf("Subspace solution X =  ");
        write_ref1("    ", x);
      }
    }
    if (alpha < one) {
      iword = 1;
    }
    else {
      iword = 0;
    }
    if (iprint >= 99) {
      printf("\n----------------exit SUBSM --------------------\n\n");
    }
  }

  /*! \brief Computation of a safeguarded step for a search procedure.
   */
  /*! This subroutine computes a safeguarded step for a search
      procedure and updates an interval that contains a step that
      satisfies a sufficient decrease and a curvature condition.

      The parameter stx contains the step with the least function
      value. If brackt is set to .true. then a minimizer has
      been bracketed in an interval with endpoints stx and sty.
      The parameter stp contains the current step.
      The subroutine assumes that if brackt is set to .true. then

            min(stx,sty) < stp < max(stx,sty),

      and that the derivative at stx is negative in the direction
      of the step.

      The subroutine statement is

        subroutine dcstep(stx,fx,dx,sty,fy,dy,stp,fp,dp,brackt,
                          stpmin,stpmax)

      where

        stx is a double precision variable.
          On entry stx is the best step obtained so far and is an
             endpoint of the interval that contains the minimizer.
          On exit stx is the updated best step.

        fx is a double precision variable.
          On entry fx is the function at stx.
          On exit fx is the function at stx.

        dx is a double precision variable.
          On entry dx is the derivative of the function at
             stx. The derivative must be negative in the direction of
             the step, that is, dx and stp - stx must have opposite
             signs.
          On exit dx is the derivative of the function at stx.

        sty is a double precision variable.
          On entry sty is the second endpoint of the interval that
             contains the minimizer.
          On exit sty is the updated endpoint of the interval that
             contains the minimizer.

        fy is a double precision variable.
          On entry fy is the function at sty.
          On exit fy is the function at sty.

        dy is a double precision variable.
          On entry dy is the derivative of the function at sty.
          On exit dy is the derivative of the function at the exit sty.

        stp is a double precision variable.
          On entry stp is the current step. If brackt is set to .true.
             then on input stp must be between stx and sty.
          On exit stp is a new trial step.

        fp is a double precision variable.
          On entry fp is the function at stp
          On exit fp is unchanged.

        dp is a double precision variable.
          On entry dp is the the derivative of the function at stp.
          On exit dp is unchanged.

        brackt is an logical variable.
          On entry brackt specifies if a minimizer has been bracketed.
             Initially brackt must be set to .false.
          On exit brackt specifies if a minimizer has been bracketed.
             When a minimizer is bracketed brackt is set to .true.

        stpmin is a double precision variable.
          On entry stpmin is a lower bound for the step.
          On exit stpmin is unchanged.

        stpmax is a double precision variable.
          On entry stpmax is an upper bound for the step.
          On exit stpmax is unchanged.

      MINPACK-1 Project. June 1983
      Argonne National Laboratory.
      Jorge J. More' and David J. Thuente.

      MINPACK-2 Project. October 1993.
      Argonne National Laboratory and University of Minnesota.
         Brett M. Averick and Jorge J. More'.
   */
  template <typename FloatType>
  void
  dcstep(
    FloatType& stx,
    FloatType& fx,
    FloatType& dx,
    FloatType& sty,
    FloatType& fy,
    FloatType& dy,
    FloatType& stp,
    FloatType const& fp,
    FloatType const& dp,
    bool& brackt,
    FloatType const& stpmin,
    FloatType const& stpmax)
  {
    FloatType zero=0.0e0;
    FloatType p66=0.66e0;
    FloatType two=2.0e0;
    FloatType three=3.0e0;
    FloatType sgnd = dp*(dx/fn::absolute(dx));
    FloatType stpf = 0; // uninitialized
    if (fp > fx) {
      // First case: A higher function value. The minimum is bracketed.
      // If the cubic step is closer to stx than the quadratic step, the
      // cubic step is taken, otherwise the average of the cubic and
      // quadratic steps is taken.
      FloatType theta = three*(fx - fp)/(stp - stx) + dx + dp;
      FloatType s=max3(fn::absolute(theta),fn::absolute(dx),fn::absolute(dp));
      FloatType theta_s = theta/s;
      FloatType gamma = s*std::sqrt(theta_s*theta_s - (dx/s)*(dp/s));
      if (stp < stx) gamma = -gamma;
      FloatType p = (gamma - dx) + theta;
      FloatType q = ((gamma - dx) + gamma) + dp;
      FloatType r = p/q;
      FloatType stpc = stx + r*(stp - stx);
      FloatType stpq = stx + ((dx/((fx-fp)/(stp-stx) + dx))/two)*(stp - stx);
      if (fn::absolute(stpc-stx) < fn::absolute(stpq-stx)) {
        stpf = stpc;
      }
      else {
        stpf = stpc + (stpq - stpc)/two;
      }
      brackt = true;
    }
    else if (sgnd < zero) {
      // Second case: A lower function value and derivatives of opposite
      // sign. The minimum is bracketed. If the cubic step is farther from
      // stp than the secant step, the cubic step is taken, otherwise the
      // secant step is taken.
      FloatType theta = three*(fx - fp)/(stp - stx) + dx + dp;
      FloatType s=max3(fn::absolute(theta),fn::absolute(dx),fn::absolute(dp));
      FloatType theta_s = theta/s;
      FloatType gamma = s*std::sqrt(theta_s*theta_s - (dx/s)*(dp/s));
      if (stp > stx) gamma = -gamma;
      FloatType p = (gamma - dp) + theta;
      FloatType q = ((gamma - dp) + gamma) + dx;
      FloatType r = p/q;
      FloatType stpc = stp + r*(stx - stp);
      FloatType stpq = stp + (dp/(dp - dx))*(stx - stp);
      if (fn::absolute(stpc-stp) > fn::absolute(stpq-stp)) {
        stpf = stpc;
      }
      else {
        stpf = stpq;
      }
      brackt = true;
    }
    else if (fn::absolute(dp) < fn::absolute(dx)) {
      // Third case: A lower function value, derivatives of the same sign,
      // and the magnitude of the derivative decreases.
      // The cubic step is computed only if the cubic tends to infinity
      // in the direction of the step or if the minimum of the cubic
      // is beyond stp. Otherwise the cubic step is defined to be the
      // secant step.
      FloatType theta = three*(fx - fp)/(stp - stx) + dx + dp;
      FloatType s=max3(fn::absolute(theta),fn::absolute(dx),fn::absolute(dp));
      // The case gamma = 0 only arises if the cubic does not tend
      // to infinity in the direction of the step.
      FloatType theta_s = theta/s;
      FloatType
        gamma = s*std::sqrt(std::max(zero,theta_s*theta_s-(dx/s)*(dp/s)));
      if (stp > stx) gamma = -gamma;
      FloatType p = (gamma - dp) + theta;
      FloatType q = (gamma + (dx - dp)) + gamma;
      FloatType r = p/q;
      FloatType stpc;
      if (r < zero && gamma != zero) {
        stpc = stp + r*(stx - stp);
      }
      else if (stp > stx) {
        stpc = stpmax;
      }
      else {
        stpc = stpmin;
      }
      FloatType stpq = stp + (dp/(dp - dx))*(stx - stp);
      if (brackt) {
        // A minimizer has been bracketed. If the cubic step is
        // closer to stp than the secant step, the cubic step is
        // taken, otherwise the secant step is taken.
        if (fn::absolute(stpc-stp) < fn::absolute(stpq-stp)) {
          stpf = stpc;
        }
        else {
          stpf = stpq;
        }
        if (stp > stx) {
          stpf = std::min(stp+p66*(sty-stp),stpf);
        }
        else {
          stpf = std::max(stp+p66*(sty-stp),stpf);
        }
      }
      else {
        // A minimizer has not been bracketed. If the cubic step is
        // farther from stp than the secant step, the cubic step is
        // taken, otherwise the secant step is taken.
        if (fn::absolute(stpc-stp) > fn::absolute(stpq-stp)) {
          stpf = stpc;
        }
        else {
          stpf = stpq;
        }
        stpf = std::min(stpmax,stpf);
        stpf = std::max(stpmin,stpf);
      }
    }
    else {
      // Fourth case: A lower function value, derivatives of the same sign,
      // and the magnitude of the derivative does not decrease. If the
      // minimum is not bracketed, the step is either stpmin or stpmax,
      // otherwise the cubic step is taken.
      if (brackt) {
        FloatType theta = three*(fp - fy)/(sty - stp) + dy + dp;
        FloatType
          s = max3(fn::absolute(theta),fn::absolute(dy),fn::absolute(dp));
        FloatType theta_s = theta/s;
        FloatType gamma = s*std::sqrt(theta_s*theta_s - (dy/s)*(dp/s));
        if (stp > sty) gamma = -gamma;
        FloatType p = (gamma - dp) + theta;
        FloatType q = ((gamma - dp) + gamma) + dy;
        FloatType r = p/q;
        FloatType stpc = stp + r*(sty - stp);
        stpf = stpc;
      }
      else if (stp > stx) {
        stpf = stpmax;
      }
      else {
        stpf = stpmin;
      }
    }
    // Update the interval which contains a minimizer.
    if (fp > fx) {
      sty = stp;
      fy = fp;
      dy = dp;
    }
    else {
      if (sgnd < zero) {
        sty = stx;
        fy = fx;
        dy = dx;
      }
      stx = stp;
      fx = fp;
      dx = dp;
    }
    // Compute the new step.
    stp = stpf;
  }

  /*! \brief Determination of a step that satisfies a sufficient
       decrease condition and a curvature condition.
   */
  /*! This subroutine finds a step that satisfies a sufficient
      decrease condition and a curvature condition.

      Each call of the subroutine updates an interval with
      endpoints stx and sty. The interval is initially chosen
      so that it contains a minimizer of the modified function

            psi(stp) = f(stp) - f(0) - ftol*stp*f'(0).

      If psi(stp) <= 0 and f'(stp) >= 0 for some step, then the
      interval is chosen so that it contains a minimizer of f.

      The algorithm is designed to find a step that satisfies
      the sufficient decrease condition

            f(stp) <= f(0) + ftol*stp*f'(0),

      and the curvature condition

            abs(f'(stp)) <= gtol*abs(f'(0)).

      If ftol is less than gtol and if, for example, the function
      is bounded below, then there is always a step which satisfies
      both conditions.

      If no step can be found that satisfies both conditions, then
      the algorithm stops with a warning. In this case stp only
      satisfies the sufficient decrease condition.

      A typical invocation of dcsrch has the following outline:

      task = 'START'
      continue
         call dcsrch( ... )
         if (task .eq. 'FG') then
            Evaluate the function and the gradient at stp
            goto 10
            end if

      NOTE: The user must no alter work arrays between calls.

      The subroutine statement is

         subroutine dcsrch(f,g,stp,ftol,gtol,xtol,stpmin,stpmax,
                           task,isave,dsave)
      where

        f is a double precision variable.
          On initial entry f is the value of the function at 0.
             On subsequent entries f is the value of the
             function at stp.
          On exit f is the value of the function at stp.

        g is a double precision variable.
          On initial entry g is the derivative of the function at 0.
             On subsequent entries g is the derivative of the
             function at stp.
          On exit g is the derivative of the function at stp.

        stp is a double precision variable.
          On entry stp is the current estimate of a satisfactory
             step. On initial entry, a positive initial estimate
             must be provided.
          On exit stp is the current estimate of a satisfactory step
             if task = 'FG'. If task = 'CONV' then stp satisfies
             the sufficient decrease and curvature condition.

        ftol is a double precision variable.
          On entry ftol specifies a nonnegative tolerance for the
             sufficient decrease condition.
          On exit ftol is unchanged.

        gtol is a double precision variable.
          On entry gtol specifies a nonnegative tolerance for the
             curvature condition.
          On exit gtol is unchanged.

        xtol is a double precision variable.
          On entry xtol specifies a nonnegative relative tolerance
             for an acceptable step. The subroutine exits with a
             warning if the relative difference between sty and stx
             is less than xtol.
          On exit xtol is unchanged.

        stpmin is a double precision variable.
          On entry stpmin is a nonnegative lower bound for the step.
          On exit stpmin is unchanged.

        stpmax is a double precision variable.
          On entry stpmax is a nonnegative upper bound for the step.
          On exit stpmax is unchanged.

        task is a character variable of length at least 60.
          On initial entry task must be set to 'START'.
          On exit task indicates the required action:

             If task(1:2) = 'FG' then evaluate the function and
             derivative at stp and call dcsrch again.

             If task(1:4) = 'CONV' then the search is successful.

             If task(1:4) = 'WARN' then the subroutine is not able
             to satisfy the convergence conditions. The exit value of
             stp contains the best point found during the search.

             If task(1:5) = 'ERROR' then there is an error in the
             input arguments.

          On exit with convergence, a warning or an error, the
             variable task contains additional information.

        isave is an integer work array of dimension 2.

        dsave is a double precision work array of dimension 13.

      Subprograms called

        MINPACK-2 ... dcstep

      MINPACK-1 Project. June 1983.
      Argonne National Laboratory.
      Jorge J. More' and David J. Thuente.

      MINPACK-2 Project. October 1993.
      Argonne National Laboratory and University of Minnesota.
      Brett M. Averick, Richard G. Carter, and Jorge J. More'.
   */
  template <typename FloatType>
  void
  dcsrch(
    FloatType const& f,
    FloatType const& g,
    FloatType& stp,
    FloatType const& ftol,
    FloatType const& gtol,
    FloatType const& xtol,
    FloatType const& stpmin,
    FloatType const& stpmax,
    std::string& task,
    ref1<int> const& isave,
    ref1<FloatType> const& dsave)
  {
#if (SCITBX_LBFGSB_RAW_ASSERTION_FLAG != 0)
    SCITBX_ASSERT(isave.size() == 2);
    SCITBX_ASSERT(dsave.size() == 13);
#endif
    FloatType zero=0.0e0;
    FloatType p5=0.5e0;
    FloatType p66=0.66e0;
    FloatType xtrapl=1.1e0;
    FloatType xtrapu=4.0e0;
    // Initialization block.
    bool brackt = false; // uninitialized
    int stage = 0; // uninitialized
    FloatType finit = 0; // uninitialized
    FloatType ginit = 0; // uninitialized
    FloatType gtest = 0; // uninitialized
    FloatType width = 0; // uninitialized
    FloatType width1 = 0; // uninitialized
    FloatType stx = 0; // uninitialized
    FloatType fx = 0; // uninitialized
    FloatType gx = 0; // uninitialized
    FloatType sty = 0; // uninitialized
    FloatType fy = 0; // uninitialized
    FloatType gy = 0; // uninitialized
    FloatType stmin = 0; // uninitialized
    FloatType stmax = 0; // uninitialized
    while (true) { // enable use of break instead of goto lbl_1000
      if (task.substr(0,5) == "START") {
        // Check the input arguments for errors.
        if (stp < stpmin) task = "ERROR: STP .LT. STPMIN";
        if (stp > stpmax) task = "ERROR: STP .GT. STPMAX";
        if (g >= zero) task = "ERROR: INITIAL G .GE. ZERO";
        if (ftol < zero) task = "ERROR: FTOL .LT. ZERO";
        if (gtol < zero) task = "ERROR: GTOL .LT. ZERO";
        if (xtol < zero) task = "ERROR: XTOL .LT. ZERO";
        if (stpmin < zero) task = "ERROR: STPMIN .LT. ZERO";
        if (stpmax < stpmin) task = "ERROR: STPMAX .LT. STPMIN";
        // Exit if there are errors on input.
        if (task.substr(0,5) == "ERROR") return;
        // Initialize local variables.
        brackt = false;
        stage = 1;
        finit = f;
        ginit = g;
        gtest = ftol*ginit;
        width = stpmax - stpmin;
        width1 = width/p5;
        // The variables stx, fx, gx contain the values of the step,
        // function, and derivative at the best step.
        // The variables sty, fy, gy contain the value of the step,
        // function, and derivative at sty.
        // The variables stp, f, g contain the values of the step,
        // function, and derivative at stp.
        stx = zero;
        fx = finit;
        gx = ginit;
        sty = zero;
        fy = finit;
        gy = ginit;
        stmin = zero;
        stmax = stp + xtrapu*stp;
        task = "FG";
        break; // goto lbl_1000;
      }
      else {
        // Restore local variables.
        if (isave(1) == 1) {
          brackt = true;
        }
        else {
          brackt = false;
        }
        stage = isave(2);
        ginit = dsave(1);
        gtest = dsave(2);
        gx = dsave(3);
        gy = dsave(4);
        finit = dsave(5);
        fx = dsave(6);
        fy = dsave(7);
        stx = dsave(8);
        sty = dsave(9);
        stmin = dsave(10);
        stmax = dsave(11);
        width = dsave(12);
        width1 = dsave(13);
      }
      // If psi(stp) <= 0 and f'(stp) >= 0 for some step, then the
      // algorithm enters the second stage.
      FloatType ftest = finit + stp*gtest;
      if (stage == 1 && f <= ftest && g >= zero) {
        stage = 2;
      }
      // Test for warnings.
      if (brackt && (stp <= stmin || stp >= stmax)) {
        task = "WARNING: ROUNDING ERRORS PREVENT PROGRESS";
      }
      if (brackt && stmax - stmin <= xtol*stmax) {
        task = "WARNING: XTOL TEST SATISFIED";
      }
      if (stp == stpmax && f <= ftest && g <= gtest) {
        task = "WARNING: STP = STPMAX";
      }
      if (stp == stpmin && (f > ftest || g >= gtest)) {
        task = "WARNING: STP = STPMIN";
      }
      if (stp == stx) {
        task = "WARNING: STP = STX";
      }
      // Test for convergence.
      if (f <= ftest && fn::absolute(g) <= gtol*(-ginit)) {
        task = "CONVERGENCE";
      }
      // Test for termination.
      if (   task.substr(0,4) == "WARN"
          || task.substr(0,4) == "CONV") break; // goto lbl_1000;
      // A modified function is used to predict the step during the
      // first stage if a lower function value has been obtained but
      // the decrease is not sufficient.
      if (stage == 1 && f <= fx && f > ftest) {
        // Define the modified function and derivative values.
        FloatType fm = f - stp*gtest;
        FloatType fxm = fx - stx*gtest;
        FloatType fym = fy - sty*gtest;
        FloatType gm = g - gtest;
        FloatType gxm = gx - gtest;
        FloatType gym = gy - gtest;
        // Call dcstep to update stx, sty, and to compute the new step.
        dcstep(stx,fxm,gxm,sty,fym,gym,stp,fm,gm,brackt,stmin,stmax);
        // Reset the function and derivative values for f.
        fx = fxm + stx*gtest;
        fy = fym + sty*gtest;
        gx = gxm + gtest;
        gy = gym + gtest;
      }
      else {
        // Call dcstep to update stx, sty, and to compute the new step.
        dcstep(stx,fx,gx,sty,fy,gy,stp,f,g,brackt,stmin,stmax);
      }
      // Decide if a bisection step is needed.
      if (brackt) {
        if (fn::absolute(sty-stx) >= p66*width1) stp = stx + p5*(sty - stx);
        width1 = width;
        width = fn::absolute(sty-stx);
      }
      // Set the minimum and maximum steps allowed for stp.
      if (brackt) {
        stmin = std::min(stx,sty);
        stmax = std::max(stx,sty);
      }
      else {
        stmin = stp + xtrapl*(stp - stx);
        stmax = stp + xtrapu*(stp - stx);
      }
      // Force the step to be within the bounds stpmax and stpmin.
      stp = std::max(stp,stpmin);
      stp = std::min(stp,stpmax);
      // If further progress is not possible, let stp be the best
      // point obtained during the search.
      if ((brackt && (stp <= stmin || stp >= stmax))
          || (brackt && stmax-stmin <= xtol*stmax)) stp = stx;
      // Obtain another function and derivative.
      task = "FG";
      break;
    }
    // lbl_1000:
    // Save local variables.
    if (brackt) {
      isave(1) = 1;
    }
    else {
      isave(1) = 0;
    }
    isave(2) = stage;
    dsave(1) =  ginit;
    dsave(2) =  gtest;
    dsave(3) =  gx;
    dsave(4) =  gy;
    dsave(5) =  finit;
    dsave(6) =  fx;
    dsave(7) =  fy;
    dsave(8) =  stx;
    dsave(9) =  sty;
    dsave(10) = stmin;
    dsave(11) = stmax;
    dsave(12) = width;
    dsave(13) = width1;
  }

  /*! \brief Line search.
   */
  /*! This subroutine calls subroutine dcsrch from the Minpack2 library
        to perform the line search.  Subroutine dscrch is safeguarded so
        that all trial points lie within the feasible region.

      Subprograms called:

        Minpack2 Library ... dcsrch.

        Linpack ... dtrsl, ddot.

                            *  *  *

      NEOS, November 1994. (Latest revision June 1996.)
      Optimization Technology Center.
      Argonne National Laboratory and Northwestern University.
      Written by
                         Ciyou Zhu
      in collaboration with R.H. Byrd, P. Lu-Chen and J. Nocedal.
   */
  template <typename FloatType>
  void
  lnsrlb(
    int const& n,
    ref1<FloatType> const& l,
    ref1<FloatType> const& u,
    ref1<int> const& nbd,
    ref1<FloatType> const& x,
    FloatType const& f,
    FloatType& fold,
    FloatType& gd,
    FloatType& gdold,
    ref1<FloatType> const& g,
    ref1<FloatType> const& d,
    ref1<FloatType> const& r,
    ref1<FloatType> const& t,
    ref1<FloatType> const& z,
    FloatType& stp,
    FloatType& dnorm,
    FloatType& dtd,
    FloatType& xstep,
    FloatType& stpmx,
    int const& iter,
    int& ifun,
    int& iback,
    int& nfgv,
    int& info,
    std::string& task,
    bool const& boxed,
    bool const& cnstnd,
    std::string& csave,
    ref1<int> const& isave,
    ref1<FloatType> const& dsave,
    bool enable_stp_init)
  {
#if (SCITBX_LBFGSB_RAW_ASSERTION_FLAG != 0)
    SCITBX_ASSERT(d.size() == n);
    SCITBX_ASSERT(r.size() == n);
    SCITBX_ASSERT(t.size() == n);
    SCITBX_ASSERT(z.size() == n);
    SCITBX_ASSERT(isave.size() == 2);
    SCITBX_ASSERT(dsave.size() == 13);
#endif
    FloatType one=1.0e0;
    FloatType zero=0.0e0;
    FloatType big=1.0e+10;
    FloatType ftol=1.0e-3;
    FloatType gtol=0.9e0;
    FloatType xtol=0.1e0;
    if (task.substr(0,5) == "FG_LN") goto lbl_556;
    if (task.substr(0,9) == "STP_INIT:") goto lbl_after_initial_stp;
    dtd = lbfgs::detail::ddot(n,d.begin(),d.begin());
    dnorm = std::sqrt(dtd);
    // Determine the maximum step length.
    stpmx = big;
    if (cnstnd) {
      if (iter == 0) {
        stpmx = one;
      }
      else {
        for(int i=1;i<=n;i++) {
          FloatType a1 = d(i);
          if (nbd(i) != 0) {
            if (a1 < zero && nbd(i) <= 2) {
              FloatType a2 = l(i) - x(i);
              if (a2 >= zero) {
                stpmx = zero;
              }
              else if (a1*stpmx < a2) {
                stpmx = a2/a1;
              }
            }
            else if (a1 > zero && nbd(i) >= 2) {
              FloatType a2 = u(i) - x(i);
              if (a2 <= zero) {
                stpmx = zero;
              }
              else if (a1*stpmx > a2) {
                stpmx = a2/a1;
              }
            }
          }
        }
      }
    }
    if (iter == 0 && ! boxed) {
      stp = std::min(one/dnorm, stpmx);
    }
    else {
      stp = one;
    }
    if (enable_stp_init) {
      task = "STP_INIT:" + task;
      return;
    }
    lbl_after_initial_stp:
    if (enable_stp_init) {
      task = task.substr(9,task.size());
    }
    dcopy(n,x,1,t,1);
    dcopy(n,g,1,r,1);
    fold = f;
    ifun = 0;
    iback = 0;
    csave = "START";
    lbl_556:
    gd = lbfgs::detail::ddot(n,g.begin(),d.begin());
    if (ifun == 0) {
      gdold=gd;
      if (gd >= zero) {
        // the directional derivative >=0.
        // Line search is impossible.
        info = -4;
        return;
      }
    }
    dcsrch(f,gd,stp,ftol,gtol,xtol,zero,stpmx,csave,isave,dsave);
    xstep = stp*dnorm;
    if (csave.substr(0,4) != "CONV" && csave.substr(0,4) != "WARN") {
      task = "FG_LNSRCH";
      ifun = ifun + 1;
      nfgv = nfgv + 1;
      iback = ifun - 1;
      if (stp == one) {
        dcopy(n,z,1,x,1);
      }
      else {
        for(int i=1;i<=n;i++) {
          x(i) = stp*d(i) + t(i);
        }
      }
    }
    else {
      task = "NEW_X";
    }
  }

  /*! \brief Updates matrices WS and WY, and forms the middle matrix in B.
   */
  /*! This subroutine updates matrices WS and WY, and forms the
        middle matrix in B.

      Subprograms called:

        Linpack ... dcopy, ddot.

                            *  *  *

      NEOS, November 1994. (Latest revision June 1996.)
      Optimization Technology Center.
      Argonne National Laboratory and Northwestern University.
      Written by
                         Ciyou Zhu
      in collaboration with R.H. Byrd, P. Lu-Chen and J. Nocedal.
   */
  template <typename FloatType>
  void
  matupd(
    int const& n,
    int const& m,
    ref2<FloatType> const& ws,
    ref2<FloatType> const& wy,
    ref2<FloatType> const& sy,
    ref2<FloatType> const& ss,
    ref1<FloatType> const& d,
    ref1<FloatType> const& r,
    int& itail,
    int const& iupdat,
    int& col,
    int& head,
    FloatType& theta,
    FloatType const& rr,
    FloatType const& dr,
    FloatType const& stp,
    FloatType const& dtd)
  {
#if (SCITBX_LBFGSB_RAW_ASSERTION_FLAG != 0)
    SCITBX_ASSERT(d.size() == n);
    SCITBX_ASSERT(r.size() == n);
    SCITBX_ASSERT(ws.size() == n*m);
    SCITBX_ASSERT(wy.size() == n*m);
    SCITBX_ASSERT(sy.size() == m*m);
    SCITBX_ASSERT(ss.size() == m*m);
#endif
    FloatType one=1;
    // Set pointers for matrices WS and WY.
    if (iupdat <= m) {
      col = iupdat;
      itail = (head+iupdat-2) % m + 1;
    }
    else {
      itail = itail % m + 1;
      head = head % m + 1;
    }
    // Update matrices WS and WY.
    dcopy(n,d,1,ws.get1(1,itail,n),1);
    dcopy(n,r,1,wy.get1(1,itail,n),1);
    // Set theta=yy/ys.
    theta = rr/dr;
    // Form the middle matrix in B.
    // update the upper triangle of SS,
    // and the lower triangle of SY:
    if (iupdat > m) {
      // move old information
      for(int j=1;j<=col-1;j++) {
        dcopy(j,ss.get1(2,j+1,j),1,ss.get1(1,j,j),1);
        dcopy(col-j,sy.get1(j+1,j+1,col-j),1,sy.get1(j,j,col-j),1);
      }
    }
    // add new information: the last row of SY
    // and the last column of SS:
    int pointr = head;
    for(int j=1;j<=col-1;j++) {
      sy(col,j) = lbfgs::detail::ddot(n,d.begin(),&wy(1,pointr));
      ss(j,col) = lbfgs::detail::ddot(n,&ws(1,pointr),d.begin());
      pointr = pointr % m + 1;
    }
    if (stp == one) {
      ss(col,col) = dtd;
    }
    else {
      ss(col,col) = stp*stp*dtd;
    }
    sy(col,col) = dr;
  }

  /*! \brief Forms the upper half of the pos. def. and symm. T.
   */
  /*! This subroutine forms the upper half of the pos. def. and symm.
        T = theta*SS + L*D^(-1)*L', stores T in the upper triangle
        of the array wt, and performs the Cholesky factorization of T
        to produce J*J', with J' stored in the upper triangle of wt.

      Subprograms called:

        Linpack ... dpofa.

                            *  *  *

      NEOS, November 1994. (Latest revision June 1996.)
      Optimization Technology Center.
      Argonne National Laboratory and Northwestern University.
      Written by
                         Ciyou Zhu
      in collaboration with R.H. Byrd, P. Lu-Chen and J. Nocedal.
   */
  template <typename FloatType>
  void
  formt(
    int const& m,
    ref2<FloatType> const& wt,
    ref2<FloatType> const& sy,
    ref2<FloatType> const& ss,
    int const& col,
    FloatType const& theta,
    int& info)
  {
#if (SCITBX_LBFGSB_RAW_ASSERTION_FLAG != 0)
    SCITBX_ASSERT(wt.size() == m*m);
    SCITBX_ASSERT(sy.size() == m*m);
    SCITBX_ASSERT(ss.size() == m*m);
#endif
    // Form the upper half of  T = theta*SS + L*D^(-1)*L',
    // store T in the upper triangle of the array wt.
    for(int j=1;j<=col;j++) {
      wt(1,j) = theta*ss(1,j);
    }
    for(int i=2;i<=col;i++) {
      for(int j=i;j<=col;j++) {
        int k1 = std::min(i,j) - 1;
        FloatType ddum = 0;
        for(int k=1;k<=k1;k++) {
          ddum = ddum + sy(i,k)*sy(j,k)/sy(k,k);
        }
        wt(i,j) = ddum + theta*ss(i,j);
      }
    }
    // Cholesky factorize T to J*J' with
    // J' stored in the upper triangle of wt.
    dpofa(wt.get2(1,1,m,col),m,col,info);
    if (info != 0) {
      info = -3;
    }
  }

  /*! \brief Solution of bound constrained optimization problems.
   */
  /*! This subroutine solves bound constrained optimization problems by
        using the compact formula of the limited memory BFGS updates.

      n is an integer variable.
        On entry n is the number of variables.
        On exit n is unchanged.

      m is an integer variable.
        On entry m is the maximum number of variable metric
           corrections allowed in the limited memory matrix.
        On exit m is unchanged.

      x is a double precision array of dimension n.
        On entry x is an approximation to the solution.
        On exit x is the current approximation.

      l is a double precision array of dimension n.
        On entry l is the lower bound of x.
        On exit l is unchanged.

      u is a double precision array of dimension n.
        On entry u is the upper bound of x.
        On exit u is unchanged.

      nbd is an integer array of dimension n.
        On entry nbd represents the type of bounds imposed on the
          variables, and must be specified as follows:
          nbd(i)=0 if x(i) is unbounded,
                 1 if x(i) has only a lower bound,
                 2 if x(i) has both lower and upper bounds,
                 3 if x(i) has only an upper bound.
        On exit nbd is unchanged.

      f is a double precision variable.
        On first entry f is unspecified.
        On final exit f is the value of the function at x.

      g is a double precision array of dimension n.
        On first entry g is unspecified.
        On final exit g is the value of the gradient at x.

      factr is a double precision variable.
        On entry factr >= 0 is specified by the user.  The iteration
          will stop when

          (f^k - f^{k+1})/max{|f^k|,|f^{k+1}|,1} <= factr*epsmch

          where epsmch is the machine precision, which is automatically
          generated by the code.
        On exit factr is unchanged.

      pgtol is a double precision variable.
        On entry pgtol >= 0 is specified by the user.  The iteration
          will stop when

                  max{|proj g_i | i = 1, ..., n} <= pgtol

          where pg_i is the ith component of the projected gradient.
        On exit pgtol is unchanged.

      ws, wy, sy, and wt are double precision working arrays used to
        store the following information defining the limited memory
           BFGS matrix:
           ws, of dimension n x m, stores S, the matrix of s-vectors;
           wy, of dimension n x m, stores Y, the matrix of y-vectors;
           sy, of dimension m x m, stores S'Y;
           ss, of dimension m x m, stores S'S;
           yy, of dimension m x m, stores Y'Y;
           wt, of dimension m x m, stores the Cholesky factorization
                                   of (theta*S'S+LD^(-1)L'); see eq.
                                   (2.26) in [3].

      wn is a double precision working array of dimension 2m x 2m
        used to store the LEL^T factorization of the indefinite matrix
                  K = [-D -Y'ZZ'Y/theta     L_a'-R_z'  ]
                      [L_a -R_z           theta*S'AA'S ]

        where     E = [-I  0]
                      [ 0  I]

      snd is a double precision working array of dimension 2m x 2m
        used to store the lower triangular part of
                  N = [Y' ZZ'Y   L_a'+R_z']
                      [L_a +R_z  S'AA'S   ]

      z(n),r(n),d(n),t(n),wa(8*m) are double precision working arrays.
        z is used at different times to store the Cauchy point and
        the Newton point.

      sg(m),sgo(m),yg(m),ygo(m) are double precision working arrays.

      index is an integer working array of dimension n.
        In subroutine freev, index is used to store the free and fixed
           variables at the Generalized Cauchy Point (GCP).

      iwhere is an integer working array of dimension n used to record
        the status of the vector x for GCP computation.
        iwhere(i)=0 or -3 if x(i) is free and has bounds,
                  1       if x(i) is fixed at l(i), and l(i) .ne. u(i)
                  2       if x(i) is fixed at u(i), and u(i) .ne. l(i)
                  3       if x(i) is always fixed, i.e.,  u(i)=x(i)=l(i)
                 -1       if x(i) is always free, i.e., no bounds on it.

      indx2 is an integer working array of dimension n.
        Within subroutine cauchy, indx2 corresponds to the array iorder.
        In subroutine freev, a list of variables entering and leaving
        the free set is stored in indx2, and it is passed on to
        subroutine formk with this information.

      task is a working string of characters of length 60 indicating
        the current job when entering and leaving this subroutine.

      iprint is an INTEGER variable that must be set by the user.
        It controls the frequency and type of output generated:
         iprint<0    no output is generated;
         iprint=0    print only one line at the last iteration;
         0<iprint<99 print also f and |proj g| every iprint iterations;
         iprint=99   print details of every iteration except n-vectors;
         iprint=100  print also the changes of active set and final x;
         iprint>100  print details of every iteration including x and g;
        When iprint > 0, the file iterate.dat will be created to
                         summarize the iteration.

      csave is a working string of characters of length 60.

      lsave is a logical working array of dimension 4.

      isave is an integer working array of dimension 23.

      dsave is a double precision working array of dimension 29.

      Subprograms called

        L-BFGS-B Library ... cauchy, subsm, lnsrlb, formk,

         errclb, prn1lb, prn2lb, prn3lb, active, projgr,

         freev, cmprlb, matupd, formt.

        Minpack2 Library ... timer, dpmeps.

        Linpack Library ... dcopy, ddot.

      References:

        [1] R. H. Byrd, P. Lu, J. Nocedal and C. Zhu, ``A limited
        memory algorithm for bound constrained optimization'',
        SIAM J. Scientific Computing 16 (1995), no. 5, pp. 1190--1208.

        [2] C. Zhu, R.H. Byrd, P. Lu, J. Nocedal, ``L-BFGS-B: FORTRAN
        Subroutines for Large Scale Bound Constrained Optimization''
        Tech. Report, NAM-11, EECS Department, Northwestern University,
        1994.

        [3] R. Byrd, J. Nocedal and R. Schnabel "Representations of
        Quasi-Newton Matrices and their use in Limited Memory Methods'',
        Mathematical Programming 63 (1994), no. 4, pp. 129-156.

        (Postscript files of these papers are available via anonymous
         ftp to eecs.nwu.edu in the directory pub/lbfgs/lbfgs_bcm.)

                            *  *  *

      NEOS, November 1994. (Latest revision June 1996.)
      Optimization Technology Center.
      Argonne National Laboratory and Northwestern University.
      Written by
                         Ciyou Zhu
      in collaboration with R.H. Byrd, P. Lu-Chen and J. Nocedal.

   */
  template <typename FloatType>
  void
  mainlb(
    int const& n,
    int const& m,
    ref1<FloatType> const& x,
    ref1<FloatType> const& l,
    ref1<FloatType> const& u,
    ref1<int> const& nbd,
    FloatType& f,
    ref1<FloatType> const& g,
    FloatType const& factr,
    FloatType const& pgtol,
    ref2<FloatType> const& ws,
    ref2<FloatType> const& wy,
    ref2<FloatType> const& sy,
    ref2<FloatType> const& ss,
    ref2<FloatType> const&
#if (SCITBX_LBFGSB_RAW_ASSERTION_FLAG != 0)
    yy
#endif
    ,
    ref2<FloatType> const& wt,
    ref2<FloatType> const& wn,
    ref2<FloatType> const& snd,
    ref1<FloatType> const& z,
    ref1<FloatType> const& r,
    ref1<FloatType> const& d,
    ref1<FloatType> const& t,
    ref1<FloatType> const& wa,
    ref1<FloatType> const& sg,
    ref1<FloatType> const&
#if (SCITBX_LBFGSB_RAW_ASSERTION_FLAG != 0)
    sgo
#endif
    ,
    ref1<FloatType> const& yg,
    ref1<FloatType> const&
#if (SCITBX_LBFGSB_RAW_ASSERTION_FLAG != 0)
    ygo
#endif
    ,
    ref1<int> const& index,
    ref1<int> const& iwhere,
    ref1<int> const& indx2,
    std::string& task,
    int const& iprint,
    std::string& csave,
    ref1<bool> const& lsave,
    ref1<int> const& isave,
    ref1<FloatType> const& dsave,
    bool enable_stp_init)
  {
    SCITBX_ASSERT(x.size() == n);
    SCITBX_ASSERT(l.size() == n);
    SCITBX_ASSERT(u.size() == n);
    SCITBX_ASSERT(nbd.size() == n);
    SCITBX_ASSERT(g.size() == n);
#if (SCITBX_LBFGSB_RAW_ASSERTION_FLAG != 0)
    SCITBX_ASSERT(ws.size() == n*m);
    SCITBX_ASSERT(wy.size() == n*m);
    SCITBX_ASSERT(sy.size() == m*m);
    SCITBX_ASSERT(ss.size() == m*m);
    SCITBX_ASSERT(yy.size() == m*m);
    SCITBX_ASSERT(wt.size() == m*m);
    SCITBX_ASSERT(wn.size() == 2*m*2*m);
    SCITBX_ASSERT(snd.size() == 2*m*2*m);
    SCITBX_ASSERT(z.size() == n);
    SCITBX_ASSERT(r.size() == n);
    SCITBX_ASSERT(d.size() == n);
    SCITBX_ASSERT(t.size() == n);
    SCITBX_ASSERT(wa.size() == 8*m);
    SCITBX_ASSERT(sg.size() == m);
    SCITBX_ASSERT(sgo.size() == m);
    SCITBX_ASSERT(yg.size() == m);
    SCITBX_ASSERT(ygo.size() == m);
    SCITBX_ASSERT(index.size() == n);
    SCITBX_ASSERT(iwhere.size() == n);
    SCITBX_ASSERT(indx2.size() == n);
    SCITBX_ASSERT(lsave.size() == 4);
    SCITBX_ASSERT(isave.size() == 23);
    SCITBX_ASSERT(dsave.size() == 29);
#endif
    FloatType zero = 0;
    FloatType one = 1;
    // begin variables in [lid]save arrays
    bool prjctd = false; // uninitialized
    bool cnstnd = false; // uninitialized
    bool boxed = false; // uninitialized
    bool updatd = false; // uninitialized
    int nintol = 0; // uninitialized
    int itfile = 0; // uninitialized
    int iback = 0; // uninitialized
    int nskip = 0; // uninitialized
    int head = 0; // uninitialized
    int col = 0; // uninitialized
    int itail = 0; // uninitialized
    int iter = 0; // uninitialized
    int iupdat = 0; // uninitialized
    int nint = 0; // uninitialized
    int nfgv = 0; // uninitialized
    int info = 0; // uninitialized
    int ifun = 0; // uninitialized
    int iword = 0; // uninitialized
    int nfree = 0; // uninitialized
    int nact = 0; // uninitialized
    int ileave = 0; // uninitialized
    int nenter = 0; // uninitialized
    FloatType theta = 0; // uninitialized
    FloatType fold = 0; // uninitialized
    FloatType tol = 0; // uninitialized
    FloatType dnorm = 0; // uninitialized
    FloatType epsmch = 0; // uninitialized
    FloatType cpu1 = 0; // uninitialized
    FloatType cachyt = 0; // uninitialized
    FloatType sbtime = 0; // uninitialized
    FloatType lnscht = 0; // uninitialized
    FloatType time1 = 0; // uninitialized
    FloatType gd = 0; // uninitialized
    FloatType stpmx = 0; // uninitialized
    FloatType sbgnrm = 0; // uninitialized
    FloatType stp = 0; // uninitialized
    FloatType gdold = 0; // uninitialized
    FloatType dtd = 0; // uninitialized
    // end variables in [lid]save arrays
    bool wrk = false; // uninitialized
    int k = 0; // uninitialized
    FloatType cpu2 = 0; // uninitialized
    FloatType dr = 0; // uninitialized
    FloatType rr = 0; // uninitialized
    FloatType xstep = 0; // uninitialized
    std::string word; // uninitialized
    if (task.substr(0,5) == "START") {
      timer(time1);
      // Generate the current machine precision.
      epsmch = scitbx::math::floating_point_epsilon<FloatType>::get();
      // Initialize counters and scalars when task='START'.
      // for the limited memory BFGS matrices:
      col    = 0;
      head   = 1;
      theta  = one;
      iupdat = 0;
      updatd = false;
      // for operation counts:
      iter   = 0;
      nfgv   = 0;
      nint   = 0;
      nintol = 0;
      nskip  = 0;
      nfree  = n;
      // for stopping tolerance:
      tol = factr*epsmch;
      // for measuring running time:
      cachyt = 0;
      sbtime = 0;
      lnscht = 0;
      // 'word' records the status of subspace solutions.
      word = "---";
      // 'info' records the termination information.
      info = 0;
      if (iprint >= 1) {
        // open a summary file 'iterate.dat'
        //PR open (8, file = 'iterate.dat', status = 'unknown');
        itfile = 8;
      }
      // Check the input arguments for errors.
      errclb(n,m,factr,l,u,nbd,task,info,k);
      if (task.substr(0,5) == "ERROR") {
        prn3lb(n,x,f,task,iprint,info,itfile,
               iter,nfgv,nintol,nskip,nact,sbgnrm,
               zero,nint,word,iback,stp,xstep,k,
               cachyt,sbtime,lnscht);
        return;
      }
      prn1lb(n,m,l,u,x,iprint,itfile,epsmch);
      // Initialize iwhere & project x onto the feasible set.
      active(n,l,u,nbd,x,iwhere,iprint,prjctd,cnstnd,boxed);
    }
    else {
      // restore local variables.
      prjctd = lsave(1);
      cnstnd = lsave(2);
      boxed  = lsave(3);
      updatd = lsave(4);
      nintol = isave(1);
      itfile = isave(3);
      iback  = isave(4);
      nskip  = isave(5);
      head   = isave(6);
      col    = isave(7);
      itail  = isave(8);
      iter   = isave(9);
      iupdat = isave(10);
      nint   = isave(12);
      nfgv   = isave(13);
      info   = isave(14);
      ifun   = isave(15);
      iword  = isave(16);
      nfree  = isave(17);
      nact   = isave(18);
      ileave = isave(19);
      nenter = isave(20);
      theta  = dsave(1);
      fold   = dsave(2);
      tol    = dsave(3);
      dnorm  = dsave(4);
      epsmch = dsave(5);
      cpu1   = dsave(6);
      cachyt = dsave(7);
      sbtime = dsave(8);
      lnscht = dsave(9);
      time1  = dsave(10);
      gd     = dsave(11);
      stpmx  = dsave(12);
      sbgnrm = dsave(13);
      stp    = dsave(14);
      gdold  = dsave(15);
      dtd    = dsave(16);
      // After returning from the driver go to the point where execution
      // is to resume.
      if (task.substr(0,9) == "STP_INIT:") goto lbl_666;
      if (task.substr(0,5) == "FG_LN") goto lbl_666;
      if (task.substr(0,5) == "NEW_X") goto lbl_777;
      if (task.substr(0,5) == "FG_ST") goto lbl_111;
      if (task.substr(0,4) == "STOP") {
        if (task.substr(6,3) == "CPU") {
          // restore the previous iterate.
          dcopy(n,t,1,x,1);
          dcopy(n,r,1,g,1);
          f = fold;
        }
        goto lbl_999;
      }
    }
    // Compute f0 and g0.
    task = "FG_START";
    // return to the driver to calculate f and g; reenter at 111.
    goto lbl_1000;
    lbl_111:
    nfgv = 1;
    // Compute the infinity norm of the (-) projected gradient.
    projgr(n,l,u,nbd,x,g,sbgnrm);
    if (iprint >= 1) {
      printf("\nAt iterate%5d    f= %12.5E    |proj g|= %12.5E\n",
             iter,f,sbgnrm);
      printf(
        " %4d %4d     -     -   -     -     -        -    %10.3E %10.3E\n",
        iter,nfgv,sbgnrm,f);
    }
    if (sbgnrm <= pgtol) {
      // terminate the algorithm.
      task = "CONVERGENCE: NORM OF PROJECTED GRADIENT <= PGTOL";
      goto lbl_999;
    }
    //----------------- the beginning of the loop --------------------------
    lbl_222:
    if (iprint >= 99) {
      printf("\n\nITERATION %5d\n", iter + 1);
    }
    iword = -1;
    if (!cnstnd && col > 0) {
      // skip the search for GCP.
      dcopy(n,x,1,z,1);
      wrk = updatd;
      nint = 0;
      goto lbl_333;
    }
    //cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    //
    //    Compute the Generalized Cauchy Point (GCP).
    //
    //cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    timer(cpu1);
    cauchy(
           n,x,l,u,nbd,g,indx2,iwhere,t,d,z,
           m,wy.get2(1,1,n,col),ws.get2(1,1,n,col),sy,wt,theta,col,head,
           wa.get1(1, 2*m),
           wa.get1(2*m+1, 2*m),
           wa.get1(4*m+1, 2*m),
           wa.get1(6*m+1, 2*m),
           nint,sg,yg,iprint,sbgnrm,info,epsmch);
    static const char* fmt_1005 =
      "\n"
      " Singular triangular system detected;\n"
      "   refresh the lbfgs memory and restart the iteration.\n";
    if (info != 0) {
      // singular triangular system detected; refresh the lbfgs memory.
      if(iprint >= 1) printf("%s", fmt_1005);
      info   = 0;
      col    = 0;
      head   = 1;
      theta  = one;
      iupdat = 0;
      updatd = false;
      timer(cpu2);
      cachyt = cachyt + cpu2 - cpu1;
      goto lbl_222;
    }
    timer(cpu2);
    cachyt = cachyt + cpu2 - cpu1;
    nintol = nintol + nint;
    // Count the entering and leaving variables for iter > 0;
    // find the index set of free and active variables at the GCP.
    freev(n,nfree,index,nenter,ileave,indx2,
          iwhere,wrk,updatd,cnstnd,iprint,iter);
    nact = n - nfree;
    lbl_333:
    // If there are no free variables or B=theta*I, then
    // skip the subspace minimization.
    if (nfree == 0 || col == 0) goto lbl_555;
    //cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    //
    //    Subspace minimization.
    //
    //cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    timer(cpu1);
    // Form  the LEL^T factorization of the indefinite
    //   matrix    K = [-D -Y'ZZ'Y/theta     L_a'-R_z'  ]
    //                 [L_a -R_z           theta*S'AA'S ]
    //   where     E = [-I  0]
    //                 [ 0  I]
    if (wrk) formk(n,nfree,index,nenter,ileave,indx2,iupdat,
                   updatd,wn,snd,m,ws,wy,sy,theta,col,head,info);
    if (info != 0) {
      // nonpositive definiteness in Cholesky factorization;
      // refresh the lbfgs memory and restart the iteration.
      if(iprint >= 1) {
        printf("\n"
          " Nonpositive definiteness in Cholesky factorization in formk;\n"
          "   refresh the lbfgs memory and restart the iteration.\n");
      }
      info   = 0;
      col    = 0;
      head   = 1;
      theta  = one;
      iupdat = 0;
      updatd = false;
      timer(cpu2);
      sbtime = sbtime + cpu2 - cpu1;
      goto lbl_222;
    }
    // compute r=-Z'B(xcp-xk)-Z'g (using wa(2m+1)=W'(xcp-x)
    // from 'cauchy').
    cmprlb(n,m,x,g,ws,wy,sy,wt,z,r,wa.get1(1,4*m),index,
           theta,col,head,nfree,cnstnd,info);
    if (info == 0) {
      // call the direct method.
      subsm(n,m,nfree,index.get1(1,nfree),l,u,nbd,z,r,ws,wy,theta,
            col,head,iword,wa.get1(1,2*m),wn,iprint,info);
    }
    if (info != 0) {
      // singular triangular system detected;
      // refresh the lbfgs memory and restart the iteration.
      if(iprint >= 1)  printf("%s", fmt_1005);
      info   = 0;
      col    = 0;
      head   = 1;
      theta  = one;
      iupdat = 0;
      updatd = false;
      timer(cpu2);
      sbtime = sbtime + cpu2 - cpu1;
      goto lbl_222;
    }
    timer(cpu2);
    sbtime = sbtime + cpu2 - cpu1;
    lbl_555:
    //cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    //
    //    Line search and optimality tests.
    //
    //cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    // Generate the search direction d:=z-x.
    for(int i=1;i<=n;i++) {
      d(i) = z(i) - x(i);
    }
    timer(cpu1);
    lbl_666:
    lnsrlb(n,l,u,nbd,x,f,fold,gd,gdold,g,d,r,t,z,stp,dnorm,
           dtd,xstep,stpmx,iter,ifun,iback,nfgv,info,task,
           boxed,cnstnd,csave,isave.get1(22,2),dsave.get1(17,13),
           enable_stp_init);
    if (info != 0 || iback >= 20) {
      // restore the previous iterate.
      dcopy(n,t,1,x,1);
      dcopy(n,r,1,g,1);
      f = fold;
      if (col == 0) {
        // abnormal termination.
        if (info == 0) {
          info = -9;
          // restore the actual number of f and g evaluations etc.
          nfgv = nfgv - 1;
          ifun = ifun - 1;
          iback = iback - 1;
        }
        task = "ABNORMAL_TERMINATION_IN_LNSRCH";
        iter = iter + 1;
        goto lbl_999;
      }
      else {
        // refresh the lbfgs memory and restart the iteration.
        if(iprint >= 1) {
          printf("\n"
                 " Bad direction in the line search;\n"
                 "   refresh the lbfgs memory and restart the iteration.\n");
        }
        if (info == 0) nfgv = nfgv - 1;
        info   = 0;
        col    = 0;
        head   = 1;
        theta  = one;
        iupdat = 0;
        updatd = false;
        task   = "RESTART_FROM_LNSRCH";
        timer(cpu2);
        lnscht = lnscht + cpu2 - cpu1;
        goto lbl_222;
      }
    }
    else if (task.substr(0,5) == "FG_LN") {
      // return to the driver for calculating f and g; reenter at 666.
      goto lbl_1000;
    }
    else if (task.substr(0,9) == "STP_INIT:") {
      goto lbl_1000;
    }
    else {
      // calculate and print out the quantities related to the new X.
      timer(cpu2);
      lnscht = lnscht + cpu2 - cpu1;
      iter = iter + 1;
      // Compute the infinity norm of the projected (-)gradient.
      projgr(n,l,u,nbd,x,g,sbgnrm);
      // Print iteration information.
      prn2lb(n,x,f,g,iprint,itfile,iter,nfgv,nact,
             sbgnrm,nint,word,iword,iback,stp,xstep);
      goto lbl_1000;
    }
    lbl_777:
    { // scope for variables
      // Test for termination.
      if (sbgnrm <= pgtol) {
        // terminate the algorithm.
        task = "CONVERGENCE: NORM OF PROJECTED GRADIENT <= PGTOL";
        goto lbl_999;
      }
      FloatType ddum = max3(fn::absolute(fold), fn::absolute(f), one);
      if ((fold - f) <= tol*ddum) {
        // terminate the algorithm.
        task = "CONVERGENCE: REL_REDUCTION_OF_F <= FACTR*EPSMCH";
        if (iback >= 10) info = -5;
        // i.e., to issue a warning if iback>10 in the line search.
        goto lbl_999;
      }
      // Compute d=newx-oldx, r=newg-oldg, rr=y'y and dr=y's.
      for(int i=1;i<=n;i++) {
        r(i) = g(i) - r(i);
      }
      rr = lbfgs::detail::ddot(n,r.begin(),r.begin());
      if (stp == one) {
        dr = gd - gdold;
        ddum = -gdold;
      }
      else {
        dr = (gd - gdold)*stp;
        dscal(n,stp,d,1);
        ddum = -gdold*stp;
      }
      if (dr <= epsmch*ddum) {
        // skip the L-BFGS update.
        nskip = nskip + 1;
        updatd = false;
        if (iprint >= 1) {
          printf("  ys=%10.3E  -gs=%10.3E BFGS update SKIPPED\n",
                 dr, ddum);
        }
        goto lbl_222; // goto lbl_888;
      }
    }
    //cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    //
    //    Update the L-BFGS matrix.
    //
    //cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    updatd = true;
    iupdat = iupdat + 1;
    // Update matrices WS and WY and form the middle matrix in B.
    matupd(n,m,ws,wy,sy,ss,d,r,itail,
           iupdat,col,head,theta,rr,dr,stp,dtd);
    // Form the upper half of the pds T = theta*SS + L*D^(-1)*L';
    // Store T in the upper triangular of the array wt;
    // Cholesky factorize T to J*J' with
    // J' stored in the upper triangular of wt.
    formt(m,wt,sy,ss,col,theta,info);
    if (info != 0) {
      // nonpositive definiteness in Cholesky factorization;
      // refresh the lbfgs memory and restart the iteration.
      if(iprint >= 1) {
        printf("\n"
          " Nonpositive definiteness in Cholesky factorization in formt;\n"
          "   refresh the lbfgs memory and restart the iteration.\n");
      }
      info = 0;
      col = 0;
      head = 1;
      theta = one;
      iupdat = 0;
      updatd = false;
      goto lbl_222;
    }
    // Now the inverse of the middle matrix in B is
    //   [  D^(1/2)      O ] [ -D^(1/2)  D^(-1/2)*L' ]
    //   [ -L*D^(-1/2)   J ] [  0        J'          ]
    // lbl_888:
    //-------------------- the end of the loop -----------------------------
    goto lbl_222;
    lbl_999:
    { // scope for variables
      FloatType time2;
      timer(time2);
      FloatType time = time2 - time1;
      prn3lb(n,x,f,task,iprint,info,itfile,
             iter,nfgv,nintol,nskip,nact,sbgnrm,
             time,nint,word,iback,stp,xstep,k,
             cachyt,sbtime,lnscht);
    }
    lbl_1000:
    // Save local variables.
    lsave(1)  = prjctd;
    lsave(2)  = cnstnd;
    lsave(3)  = boxed;
    lsave(4)  = updatd;
    isave(1)  = nintol;
    isave(3)  = itfile;
    isave(4)  = iback;
    isave(5)  = nskip;
    isave(6)  = head;
    isave(7)  = col;
    isave(8)  = itail;
    isave(9)  = iter;
    isave(10) = iupdat;
    isave(12) = nint;
    isave(13) = nfgv;
    isave(14) = info;
    isave(15) = ifun;
    isave(16) = iword;
    isave(17) = nfree;
    isave(18) = nact;
    isave(19) = ileave;
    isave(20) = nenter;
    dsave(1)  = theta;
    dsave(2)  = fold;
    dsave(3)  = tol;
    dsave(4)  = dnorm;
    dsave(5)  = epsmch;
    dsave(6)  = cpu1;
    dsave(7)  = cachyt;
    dsave(8)  = sbtime;
    dsave(9)  = lnscht;
    dsave(10) = time1;
    dsave(11) = gd;
    dsave(12) = stpmx;
    dsave(13) = sbgnrm;
    dsave(14) = stp;
    dsave(15) = gdold;
    dsave(16) = dtd;
  }

  /*! \brief Main L-BFGS-B interface function.
   */
  /*! This subroutine partitions the working arrays wa and iwa, and
        then uses the limited memory BFGS method to solve the bound
        constrained optimization problem by calling mainlb.
        (The direct method will be used in the subspace minimization.)

      n is an integer variable.
        On entry n is the dimension of the problem.
        On exit n is unchanged.

      m is an integer variable.
        On entry m is the maximum number of variable metric corrections
          used to define the limited memory matrix.
        On exit m is unchanged.

      x is a double precision array of dimension n.
        On entry x is an approximation to the solution.
        On exit x is the current approximation.

      l is a double precision array of dimension n.
        On entry l is the lower bound on x.
        On exit l is unchanged.

      u is a double precision array of dimension n.
        On entry u is the upper bound on x.
        On exit u is unchanged.

      nbd is an integer array of dimension n.
        On entry nbd represents the type of bounds imposed on the
          variables, and must be specified as follows:
          nbd(i)=0 if x(i) is unbounded,
                 1 if x(i) has only a lower bound,
                 2 if x(i) has both lower and upper bounds, and
                 3 if x(i) has only an upper bound.
        On exit nbd is unchanged.

      f is a double precision variable.
        On first entry f is unspecified.
        On final exit f is the value of the function at x.

      g is a double precision array of dimension n.
        On first entry g is unspecified.
        On final exit g is the value of the gradient at x.

      factr is a double precision variable.
        On entry factr >= 0 is specified by the user.  The iteration
          will stop when

          (f^k - f^{k+1})/max{|f^k|,|f^{k+1}|,1} <= factr*epsmch

          where epsmch is the machine precision, which is automatically
          generated by the code. Typical values for factr: 1.d+12 for
          low accuracy; 1.d+7 for moderate accuracy; 1.d+1 for extremely
          high accuracy.
        On exit factr is unchanged.

      pgtol is a double precision variable.
        On entry pgtol >= 0 is specified by the user.  The iteration
          will stop when

                  max{|proj g_i | i = 1, ..., n} <= pgtol

          where pg_i is the ith component of the projected gradient.
        On exit pgtol is unchanged.

      wa is a double precision working array of length
        (2mmax + 4)nmax + 12mmax^2 + 12mmax.

      iwa is an integer working array of length 3nmax.

      task is a working string of characters of length 60 indicating
        the current job when entering and quitting this subroutine.

      iprint is an integer variable that must be set by the user.
        It controls the frequency and type of output generated:
         iprint<0    no output is generated;
         iprint=0    print only one line at the last iteration;
         0<iprint<99 print also f and |proj g| every iprint iterations;
         iprint=99   print details of every iteration except n-vectors;
         iprint=100  print also the changes of active set and final x;
         iprint>100  print details of every iteration including x and g;
        When iprint > 0, the file iterate.dat will be created to
                         summarize the iteration.

      csave is a working string of characters of length 60.

      lsave is a logical working array of dimension 4.
        On exit with 'task' = NEW_X, the following information is
                                                              available:
          If lsave(1) = .true.  then  the initial X has been replaced by
                                      its projection in the feasible set;
          If lsave(2) = .true.  then  the problem is constrained;
          If lsave(3) = .true.  then  each variable has upper and lower
                                      bounds;

      isave is an integer working array of dimension 44.
        On exit with 'task' = NEW_X, the following information is
                                                              available:
          isave(22) = the total number of intervals explored in the
                          search of Cauchy points;
          isave(26) = the total number of skipped BFGS updates before
                          the current iteration;
          isave(30) = the number of current iteration;
          isave(31) = the total number of BFGS updates prior the current
                          iteration;
          isave(33) = the number of intervals explored in the search of
                          Cauchy point in the current iteration;
          isave(34) = the total number of function and gradient
                          evaluations;
          isave(36) = the number of function value or gradient
                                   evaluations in the current iteration;
          if isave(37) = 0  then the subspace argmin is within the box;
          if isave(37) = 1  then the subspace argmin is beyond the box;
          isave(38) = the number of free variables in the current
                          iteration;
          isave(39) = the number of active constraints in the current
                          iteration;
          n + 1 - isave(40) = the number of variables leaving the set of
                            active constraints in the current iteration;
          isave(41) = the number of variables entering the set of active
                          constraints in the current iteration.

      dsave is a double precision working array of dimension 29.
        On exit with 'task' = NEW_X, the following information is
                                                              available:
          dsave(1) = current 'theta' in the BFGS matrix;
          dsave(2) = f(x) in the previous iteration;
          dsave(3) = factr*epsmch;
          dsave(4) = 2-norm of the line search direction vector;
          dsave(5) = the machine precision epsmch generated by the code;
          dsave(7) = the accumulated time spent on searching for
                                                          Cauchy points;
          dsave(8) = the accumulated time spent on
                                                  subspace minimization;
          dsave(9) = the accumulated time spent on line search;
          dsave(11) = the slope of the line search function at
                                   the current point of line search;
          dsave(12) = the maximum relative step length imposed in
                                                            line search;
          dsave(13) = the infinity norm of the projected gradient;
          dsave(14) = the relative step length in the line search;
          dsave(15) = the slope of the line search function at
                                  the starting point of the line search;
          dsave(16) = the square of the 2-norm of the line search
                                                       direction vector.

      Subprograms called:

        L-BFGS-B Library ... mainlb.

      References:

        [1] R. H. Byrd, P. Lu, J. Nocedal and C. Zhu, ``A limited
        memory algorithm for bound constrained optimization'',
        SIAM J. Scientific Computing 16 (1995), no. 5, pp. 1190--1208.

        [2] C. Zhu, R.H. Byrd, P. Lu, J. Nocedal, ``L-BFGS-B: a
        limited memory FORTRAN code for solving bound constrained
        optimization problems'', Tech. Report, NAM-11, EECS Department,
        Northwestern University, 1994.

        (Postscript files of these papers are available via anonymous
         ftp to eecs.nwu.edu in the directory pub/lbfgs/lbfgs_bcm.)

                            *  *  *

      NEOS, November 1994. (Latest revision June 1996.)
      Optimization Technology Center.
      Argonne National Laboratory and Northwestern University.
      Written by
                         Ciyou Zhu
      in collaboration with R.H. Byrd, P. Lu-Chen and J. Nocedal.
   */
  template <typename FloatType>
  void
  setulb(
    int const& n,
    int const& m,
    ref1<FloatType> const& x,
    ref1<FloatType> const& l,
    ref1<FloatType> const& u,
    ref1<int> const& nbd,
    FloatType& f,
    ref1<FloatType> const& g,
    FloatType const& factr,
    FloatType const& pgtol,
    ref1<FloatType> const& wa,
    ref1<int> const& iwa,
    std::string& task,
    int const& iprint,
    std::string& csave,
    ref1<bool> const& lsave,
    ref1<int> const& isave,
    ref1<FloatType> const& dsave,
    bool enable_stp_init=false)
  {
#if (SCITBX_LBFGSB_RAW_ASSERTION_FLAG != 0)
    SCITBX_ASSERT(wa.size() == 2*m*n+4*n+12*m*m+12*m);
    SCITBX_ASSERT(iwa.size() == 3*n);
    SCITBX_ASSERT(lsave.size() == 4);
    SCITBX_ASSERT(isave.size() == 44);
    SCITBX_ASSERT(dsave.size() == 29);
#endif
    if (task.substr(0,5) == "START") {
      isave(1)  = m*n;
      isave(2)  = m*m;
      isave(3)  = 4*m*m;
      isave(4)  = 1;
      isave(5)  = isave(4)  + isave(1);
      isave(6)  = isave(5)  + isave(1);
      isave(7)  = isave(6)  + isave(2);
      isave(8)  = isave(7)  + isave(2);
      isave(9)  = isave(8)  + isave(2);
      isave(10) = isave(9)  + isave(2);
      isave(11) = isave(10) + isave(3);
      isave(12) = isave(11) + isave(3);
      isave(13) = isave(12) + n;
      isave(14) = isave(13) + n;
      isave(15) = isave(14) + n;
      isave(16) = isave(15) + n;
      isave(17) = isave(16) + 8*m;
      isave(18) = isave(17) + m;
      isave(19) = isave(18) + m;
      isave(20) = isave(19) + m;
    }
    //int l1   = isave(1); // unused
    //int l2   = isave(2); // unused
    //int l3   = isave(3); // unused
    int lws  = isave(4);
    int lwy  = isave(5);
    int lsy  = isave(6);
    int lss  = isave(7);
    int lyy  = isave(8);
    int lwt  = isave(9);
    int lwn  = isave(10);
    int lsnd = isave(11);
    int lz   = isave(12);
    int lr   = isave(13);
    int ld   = isave(14);
    int lt   = isave(15);
    int lwa  = isave(16);
    int lsg  = isave(17);
    int lsgo = isave(18);
    int lyg  = isave(19);
    int lygo = isave(20);
    mainlb(n,m,x,l,u,nbd,f,g,factr,pgtol,
      wa.get2(lws, n, m),
      wa.get2(lwy, n, m),
      wa.get2(lsy, m, m),
      wa.get2(lss, m, m),
      wa.get2(lyy, m, m),
      wa.get2(lwt, m, m),
      wa.get2(lwn, 2*m, 2*m),
      wa.get2(lsnd, 2*m, 2*m),
      wa.get1(lz, n),
      wa.get1(lr, n),
      wa.get1(ld, n),
      wa.get1(lt, n),
      wa.get1(lwa, 8*m),
      wa.get1(lsg, m),
      wa.get1(lsgo, m),
      wa.get1(lyg, m),
      wa.get1(lygo, m),
      iwa.get1(1, n),
      iwa.get1(n+1, n),
      iwa.get1(2*n+1, n),
      task,iprint,csave,lsave,
      ref1<int>(&isave(22), 23),
      dsave,
      enable_stp_init);
  }

}}} // namspace scitbx::lbfgsb::raw

#endif // SCITBX_LBFGSB_RAW_H
