#ifndef SCITBX_MATRIX_GIVENS_H
#define SCITBX_MATRIX_GIVENS_H

#include <scitbx/array_family/accessors/mat_grid.h>
#include <scitbx/array_family/shared.h>
#include <boost/type_traits.hpp>
#include <cmath>

namespace scitbx { namespace matrix { namespace givens {

/// Givens rotation
/**
  R = [ c -s ]
      [ s  c ]
  with c^2 + s^2 = 1.

  R^T is applied to the left of a column vector [ x0 x1 ]^t
  whereas R is applied to the right of a row vector [ x0 x1 ], resulting in
  both case in
    x0 :=  c x0 + s x1
    x1 := -s x0 + c x1
  for a transformation in-place.

  Reference (beware sign convention for the sine):
    [1] Golub and Van Loan, section 5.1.8 and following
    [2] Ake Bjorck, Numerical Methods for Least-Squares Problems
    [3] LAPACK
*/
template <typename FloatType>
struct rotation
{
  typedef FloatType scalar_t;

  scalar_t c, s;

  /// Unitialized
  rotation()
  {}

  /// Construct the rotation with the given cosine and sine
  rotation(scalar_t cosine, scalar_t sine)
  : c(cosine), s(sine)
  {}

  /// Transform [ x0 x1 ] in-place
  void apply(scalar_t &x0, scalar_t &x1) {
    scalar_t y0 =  c*x0 + s*x1;
    scalar_t y1 = -s*x0 + c*x1;
    x0 = y0;
    x1 = y1;
  }

  /// Transform [ x0 x1 ] assuming that x1 = 0 on entry (in-place)
  void apply_assuming_null_x1(scalar_t &x0, scalar_t &x1) {
    x1 = -s*x0;
    x0 =  c*x0;
  }

  /// Transform [ x0 x1 ] assuming that x0 = 0 on entry (in-place)
  void apply_assuming_null_x0(scalar_t &x0, scalar_t &x1) {
    x0 = s*x1;
    x1 = c*x1;
  }

  /// Set this to the rotation of [ x0 x1 ] zeroing x1 and apply it to [ x0 x1 ]
  /** Reference: algorithm 5.1.3 */
  void zero_x1(scalar_t &x0, scalar_t &x1) {
    scalar_t a=x0, b=x1;
    if (b == 0) {
      c = 1;
      s = 0;
    }
    else if (a == 0) {
      c = 0;
      s = 1;
      x0 = b;
    }
    else {
      if (std::abs(b) > std::abs(a)) {
        scalar_t tau = a/b;
        scalar_t t = std::sqrt(1+tau*tau);
        s = 1/t;
        c = s*tau;
        x0 = t*b;
      }
      else {
        scalar_t tau = b/a;
        scalar_t t = std::sqrt(1+tau*tau);
        c = 1/t;
        s = c*tau;
        x0 = t*a;
      }
    }
    x1 = 0;
  }

  /// Assuming that on entry y1 = 0, setup a rotation s.t. on exit x1 = 0.
  /** The rotation is the applied in place to [ x0 x1 ] and [ y0 y1 ] */
  void chase_nonzero_from_x1_to_y1(scalar_t &x0, scalar_t &y0,
                                   scalar_t &x1, scalar_t &y1)
  {
    zero_x1(x0, x1);
    apply_assuming_null_x1(y0, y1);
  }

  /// Assuming that on entry y0 = 0, setup a rotation s.t. on exit x1 = 0.
  /** The rotation is the applied in place to [ x0 x1 ] and [ y0 y1 ] */
  void chase_nonzero_from_x1_to_y0(scalar_t &x0, scalar_t &y0,
                                   scalar_t &x1, scalar_t &y1)
  {
    zero_x1(x0, x1);
    apply_assuming_null_x0(y0, y1);
  }

  /// Setup a rotation s.t. on exit x1 = 0.
  /** The rotation is then applied in-place to [ x0 x1 ] and [ y0 y1 ] */
  void chase_nonzero_from_x1_off(scalar_t &x0, scalar_t &y0,
                                 scalar_t &x1, scalar_t &y1)
  {
    zero_x1(x0, x1);
    apply(y0, y1);
  }

  /// Assuming that on entry z0 = 0, setup a rotation s.t. on exit x1 = 0.
  /** The rotation is then applied in-place to [ x0 x1 ], [ y0 y1 ]
   and [ z0 z1 ]
   */
  void chase_nonzero_from_x1_to_z0(scalar_t &x0, scalar_t &y0, scalar_t &z0,
                                   scalar_t &x1, scalar_t &y1, scalar_t &z1)
  {
    chase_nonzero_from_x1_off(x0, y0,
                              x1, y1);
    apply_assuming_null_x0(z0, z1);
  }

  /// Perform u = g u in-place, where g is this rotation
  void apply_on_left(af::ref<scalar_t, af::mat_grid> const &u, int k0, int k1) {
    for (int j=0; j < u.n_columns(); ++j) apply(u(k0, j), u(k1, j));
  }

  /// Perform v = v g in-place, where g is this rotation
  void apply_on_right(af::ref<scalar_t, af::mat_grid> const &v, int k0, int k1) {
    for (int i=0; i < v.n_rows(); ++i) apply(v(i, k0), v(i, k1));
  }
};


/// The pair of rotations used in the implicit zero-shift QR iteration
/// devised by Demmel and Kahan for the SVD of bidiagonal matrices
/** This code is equivalent to the second part of the loop in the implicit
zero-shift algorithm of [1], followed by the first part (for a generic iteration,
i.e. not the first one which is special). That way, we don't need the variables
oldcs and oldsn.

Reference: [1] James Demmel and W. Kahan.
               Accurate singular values of bidiagonal matrices.
               SIAM J. Sci. Stat. Comput, 11:873–912, 1990.
*/
template <typename FloatType>
struct demmel_kahan_rotations
{
  typedef FloatType scalar_t;

  rotation<scalar_t> g1, g2;

  /// Assuming f0 and t being zero on entry, zero out z by applying g1 then g2
  /** Transformation in-place. Two cases:
        - g1 applied on the left and g2 on the right for a block
           [  d0  f0     ]
           [   z  d1 f1  ]
           [      t  d2  ]
        - g1 applied on the right and g2 on the left for a block
           [  d2  f1     ]
           [   t  d1 f0  ]
           [      z  d0  ]

      If the conventions used in givens::rotation were to be changed, this code
      would need to be changed too.
  */
  void chase_non_zero_from_z_to_t(scalar_t &d0, scalar_t &f0,
                                  scalar_t &z , scalar_t &d1, scalar_t &f1,
                                                scalar_t &t , scalar_t &d2)
  {
    g1.zero_x1(d0, z);
    g2.zero_x1(d1, f1);
    f0 = g1.s *d1;
    d1 = g1.c *d1;
    t  = g2.s*d2;
    d2 = g2.c*d2;
  }
};


/// Product of Given rotations \f$G_1, G_2, ... ,G_n\f$
/** rotations have contiguous indices (r, r+1), (r+1, r+2), ... (downward case)
    or (s, s-1), (s-1, s-2), ... (upward case)
*/
template <typename FloatType>
struct product
{
  typedef FloatType scalar_t;

  af::shared<rotation<scalar_t> > terms_;
  af::ref<rotation<scalar_t> > terms;
  int end;
  bool effective;

  /// Construct an empty product
  /** - \c effective controls whether \c apply_on_right and \c apply_on_left
        do anything at all
      - \c buffer_size is the size of the buffer reserved to store the rotations
        in the product
  */
  product(int buffer_size, bool effective_=true)
    : terms_(buffer_size, af::init_functor_null<rotation<scalar_t> >()),
      terms(terms_.ref()),
      end(0),
      effective(effective_)
  {}

  /// Append g to the list \f$G_1, G_2, ..., G_n\f$
  /** Buffer overflow is not guarded against. */
  void multiply_by(rotation<scalar_t> const &g) { terms[end++] = g; }

  /// \f$ U = G_n G_{n-1} ... G_1 U\f$ in-place
  void apply_downward_on_left(af::ref<scalar_t, af::mat_grid> const &u, int r) {
    if (effective) {
      for (int i=0; i<end; ++i) terms[i].apply_on_left(u, r+i, r+i+1);
    }
    end = 0;
  }

  /// \f$ V = V G_1 G_2 ... G_n\f$ in-place
  void apply_downward_on_right(af::ref<scalar_t, af::mat_grid> const &v, int r) {
    if (effective) {
      for (int i=0; i<end; ++i) terms[i].apply_on_right(v, r+i, r+i+1);
    }
    end = 0;
  }

  /// \f$ U = G_n G_{n-1} ... G_1 U\f$ in-place
  void apply_upward_on_left(af::ref<scalar_t, af::mat_grid> const &u, int s) {
    if (effective) {
      for (int i=0; i<end; ++i) terms[i].apply_on_left(u, s-i, s-i-1);
    }
    end = 0;
  }

  /// \f$ V = V G_1 G_2 ... G_n\f$ in-place
  void apply_upward_on_right(af::ref<scalar_t, af::mat_grid> const &v, int s) {
    if (effective) {
      for (int i=0; i<end; ++i) terms[i].apply_on_right(v, s-i, s-i-1);
    }
    end = 0;
  }
};

}}} // scitbx::matrix::givens

#if !defined(BOOST_NO_TEMPLATE_PARTIAL_SPECIALIZATION)

namespace boost {
  template <typename FloatType>
  struct has_trivial_destructor<scitbx::matrix::givens::rotation<FloatType> > {
    static const bool value = ::boost::has_trivial_destructor<FloatType>::value;
  };
}

#endif


#endif // GUARD
