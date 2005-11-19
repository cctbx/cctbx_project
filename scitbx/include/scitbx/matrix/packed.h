#ifndef SCITBX_MATRIX_PACKED_H
#define SCITBX_MATRIX_PACKED_H

#include <scitbx/array_family/versa.h>
#include <scitbx/array_family/shared.h>

namespace scitbx { namespace matrix {

  class packed_u_accessor
  {
    public:
      typedef unsigned index_value_type;
      typedef af::tiny_plain<unsigned, 2> index_type;

      explicit
      packed_u_accessor(unsigned n_) : n(n_) {}

      std::size_t
      size_1d() const { return n*(n+1)/2; }

      unsigned
      operator()(unsigned i, unsigned j) const
      {
        return i*(n-1)-i*(i-1)/2+j;
      }

      unsigned n;
  };

  inline
  unsigned
  symmetric_n_from_packed_size(std::size_t packed_size)
  {
    unsigned n = static_cast<unsigned>(
      (std::sqrt(1.0+8.0*static_cast<double>(packed_size))-1.0)/2.0 + 0.5);
    SCITBX_ASSERT(n*(n+1)/2 == packed_size);
    return n;
  }

  template <typename FloatType>
  af::shared<FloatType>
  upper_triangle_as_packed_u(
    af::const_ref<FloatType, af::c_grid<2> > const& a)
  {
    SCITBX_ASSERT(a.accessor().is_square());
    typename af::c_grid<2>::index_value_type n = a.accessor()[0];
    af::shared<FloatType> result(
      n*(n+1)/2, af::init_functor_null<FloatType>());
    FloatType *r = result.begin();
    std::size_t ij = 0;
    for(unsigned i=0;i<n;i++) {
      ij += i;
      for(unsigned j=i;j<n;j++) {
        *r++ = a[ij++];
      }
    }
    return result;
  }

  template <typename FloatType>
  af::versa<FloatType, af::c_grid<2> >
  packed_u_as_upper_triangle(
    af::const_ref<FloatType> const& a)
  {
    unsigned n = symmetric_n_from_packed_size(a.size());
    af::versa<FloatType, af::c_grid<2> > result(
      af::c_grid<2>(n,n), af::init_functor_null<FloatType>());
    FloatType *r = result.begin();
    std::size_t i_a = 0;
    std::size_t ij = 0;
    for(unsigned i=0;i<n;i++) {
      unsigned j = 0;
      for (;j<i;j++) r[ij++] = 0;
      for (;j<n;j++) r[ij++] = a[i_a++];
    }
    return result;
  }

  template <typename FloatType>
  af::shared<FloatType>
  lower_triangle_as_packed_l(
    af::const_ref<FloatType, af::c_grid<2> > const& a)
  {
    SCITBX_ASSERT(a.accessor().is_square());
    typename af::c_grid<2>::index_value_type n = a.accessor()[0];
    af::shared<FloatType> result(
      n*(n+1)/2, af::init_functor_null<FloatType>());
    FloatType *r = result.begin();
    std::size_t i0 = 0;
    for(unsigned i=0;i<n;i++) {
      std::size_t ij = i0;
      for(unsigned j=0;j<=i;j++) {
        *r++ = a[ij++];
      }
      i0 += n;
    }
    return result;
  }

  template <typename FloatType>
  af::versa<FloatType, af::c_grid<2> >
  packed_l_as_lower_triangle(
    af::const_ref<FloatType> const& a)
  {
    unsigned n = symmetric_n_from_packed_size(a.size());
    af::versa<FloatType, af::c_grid<2> > result(
      af::c_grid<2>(n,n), af::init_functor_null<FloatType>());
    FloatType *r = result.begin();
    std::size_t i_a = 0;
    std::size_t ij = 0;
    for(unsigned i=0;i<n;i++) {
      unsigned j = 0;
      for (;j<=i;j++) r[ij++] = a[i_a++];
      for (;j<n;j++) r[ij++] = 0;
    }
    return result;
  }

  template <typename FloatType>
  af::shared<FloatType>
  symmetric_as_packed_u(
    af::const_ref<FloatType, af::c_grid<2> > const& a,
    FloatType const& relative_eps=1.e-12)
  {
    SCITBX_ASSERT(a.accessor().is_square());
    typename af::c_grid<2>::index_value_type n = a.accessor()[0];
    af::shared<FloatType> result(
      n*(n+1)/2, af::init_functor_null<FloatType>());
    bool use_eps;
    FloatType eps;
    if (relative_eps <= 0 || n == 0) {
      use_eps = false;
    }
    else {
      eps = relative_eps * af::max_absolute(a);
      use_eps = true;
    }
    FloatType *r = result.begin();
    std::size_t ij = 0;
    for(unsigned i=0;i<n;i++) {
      ij += i;
      std::size_t jnpi = ij+n;
      *r++ = a[ij++];
      for(unsigned j=i+1;j<n;j++,jnpi+=n) {
        FloatType const& a_ij = a[ij++];
        FloatType ave = (a_ij + a[jnpi]) / 2;
        if (use_eps && fn::absolute(a_ij - ave) > eps) {
          throw std::runtime_error(
            "symmetric_as_packed_u(): matrix is not symmetric.");
        }
        *r++ = ave;
      }
    }
    return result;
  }

  template <typename FloatType>
  af::shared<FloatType>
  symmetric_as_packed_l(
    af::const_ref<FloatType, af::c_grid<2> > const& a,
    FloatType const& relative_eps=1.e-12)
  {
    SCITBX_ASSERT(a.accessor().is_square());
    typename af::c_grid<2>::index_value_type n = a.accessor()[0];
    af::shared<FloatType> result(
      n*(n+1)/2, af::init_functor_null<FloatType>());
    bool use_eps;
    FloatType eps;
    if (relative_eps <= 0 || n == 0) {
      use_eps = false;
    }
    else {
      eps = relative_eps * af::max_absolute(a);
      use_eps = true;
    }
    FloatType *r = result.begin();
    std::size_t i0 = 0;
    for(unsigned i=0;i<n;i++,i0+=n) {
      std::size_t ij = i0;
      std::size_t jnpi = i;
      for(unsigned j=0;j<i;j++,jnpi+=n) {
        FloatType const& a_ij = a[ij++];
        FloatType ave = (a_ij + a[jnpi]) / 2;
        if (use_eps && fn::absolute(a_ij - ave) > eps) {
          throw std::runtime_error(
            "symmetric_as_packed_l(): matrix is not symmetric.");
        }
        *r++ = ave;
      }
      *r++ = a[ij];
    }
    return result;
  }

  template <typename FloatType>
  af::versa<FloatType, af::c_grid<2> >
  packed_u_as_symmetric(
    af::const_ref<FloatType> const& a)
  {
    unsigned n = symmetric_n_from_packed_size(a.size());
    af::versa<FloatType, af::c_grid<2> > result(
      af::c_grid<2>(n,n), af::init_functor_null<FloatType>());
    FloatType *r = result.begin();
    std::size_t i_a = 0;
    std::size_t ij = 0;
    for(unsigned i=0;i<n;i++) {
      ij += i;
      std::size_t jnpi = ij+n;
      r[ij++] = a[i_a++];
      for(unsigned j=i+1;j<n;j++,jnpi+=n) {
        r[ij++] = r[jnpi] = a[i_a++];
      }
    }
    return result;
  }

  template <typename FloatType>
  af::versa<FloatType, af::c_grid<2> >
  packed_l_as_symmetric(
    af::const_ref<FloatType> const& a)
  {
    unsigned n = symmetric_n_from_packed_size(a.size());
    af::versa<FloatType, af::c_grid<2> > result(
      af::c_grid<2>(n,n), af::init_functor_null<FloatType>());
    FloatType *r = result.begin();
    std::size_t i_a = 0;
    std::size_t i0 = 0;
    for(unsigned i=0;i<n;i++,i0+=n) {
      std::size_t ij = i0;
      std::size_t jnpi = i;
      for(unsigned j=0;j<i;j++,jnpi+=n) {
        r[ij++] = r[jnpi] = a[i_a++];
      }
      r[ij] = a[i_a++];
    }
    return result;
  }

}} // namespace scitbx::matrix

#endif // SCITBX_MATRIX_PACKED_H
