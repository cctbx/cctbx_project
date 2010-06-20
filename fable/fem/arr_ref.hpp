#ifndef FEM_ARR_REF_HPP
#define FEM_ARR_REF_HPP

#include <fem/dimension.hpp>
#include <fem/zero.hpp>
#include <algorithm>
#include <cstring>

namespace fem {

  template <size_t Ndims>
  struct arr_ref_dims;

  template <>
  struct arr_ref_dims<1> : dim_data<1>
  {
    size_t
    index_1d(
      ssize_t i1) const
    {
      return i1 - this->origin[0];
    }
  };

  template <>
  struct arr_ref_dims<2> : dim_data<2>
  {
    size_t
    index_1d(
      ssize_t i1,
      ssize_t i2) const
    {
      return
          (i2 - this->origin[1]) * this->all[0]
        + (i1 - this->origin[0]);
    }

  };

  template <>
  struct arr_ref_dims<3> : dim_data<3>
  {
    size_t
    index_1d(
      ssize_t i1,
      ssize_t i2,
      ssize_t i3) const
    {
      return
          ((i3 - this->origin[2])  * this->all[1]
        +  (i2 - this->origin[1])) * this->all[0]
        +  (i1 - this->origin[0]);
    }
  };

  template <>
  struct arr_ref_dims<4> : dim_data<4>
  {
    size_t
    index_1d(
      ssize_t i1,
      ssize_t i2,
      ssize_t i3,
      ssize_t i4) const
    {
      return
          (((i4 - this->origin[3])  * this->all[2]
        +   (i3 - this->origin[2])) * this->all[1]
        +   (i2 - this->origin[1])) * this->all[0]
        +   (i1 - this->origin[0]);
    }
  };

  template <>
  struct arr_ref_dims<5> : dim_data<5>
  {
    size_t
    index_1d(
      ssize_t i1,
      ssize_t i2,
      ssize_t i3,
      ssize_t i4,
      ssize_t i5) const
    {
      return
          ((((i5 - this->origin[4])  * this->all[3]
        +    (i4 - this->origin[3])) * this->all[2]
        +    (i3 - this->origin[2])) * this->all[1]
        +    (i2 - this->origin[1])) * this->all[0]
        +    (i1 - this->origin[0]);
    }
  };

  template <>
  struct arr_ref_dims<6> : dim_data<6>
  {
    size_t
    index_1d(
      ssize_t i1,
      ssize_t i2,
      ssize_t i3,
      ssize_t i4,
      ssize_t i5,
      ssize_t i6) const
    {
      return
          (((((i6 - this->origin[5])  * this->all[4]
        +     (i5 - this->origin[4])) * this->all[3]
        +     (i4 - this->origin[3])) * this->all[2]
        +     (i3 - this->origin[2])) * this->all[1]
        +     (i2 - this->origin[1])) * this->all[0]
        +     (i1 - this->origin[0]);
    }
  };

  template <typename T, size_t Ndims=1>
  struct arr_cref : arr_ref_dims<Ndims>
  {
    protected:
      T const* elems_;

      arr_cref() {}

      public:

    arr_cref(
      T const& val)
    :
      elems_(&val)
    {}

    template <size_t OtherNdims>
    arr_cref(
      arr_cref<T, OtherNdims> const& other)
    {
      elems_ = other.begin();
    }

    T const*
    begin() const { return elems_; }

    template <size_t BufferNdims>
    void
    operator()(
      dim_data<BufferNdims> const& dims)
    {
      this->copy_origin_all(dims);
    }

    T const&
    operator()(
      ssize_t i1) const
    {
      return elems_[this->index_1d(i1)];
    }

    T const&
    operator()(
      ssize_t i1,
      ssize_t i2) const
    {
      return elems_[this->index_1d(i1, i2)];
    }

    T const&
    operator()(
      ssize_t i1,
      ssize_t i2,
      ssize_t i3) const
    {
      return elems_[this->index_1d(i1, i2, i3)];
    }

    T const&
    operator()(
      ssize_t i1,
      ssize_t i2,
      ssize_t i3,
      ssize_t i4) const
    {
      return elems_[this->index_1d(i1, i2, i3, i4)];
    }

    T const&
    operator()(
      ssize_t i1,
      ssize_t i2,
      ssize_t i3,
      ssize_t i4,
      ssize_t i5) const
    {
      return elems_[this->index_1d(i1, i2, i3, i4, i5)];
    }

    T const&
    operator()(
      ssize_t i1,
      ssize_t i2,
      ssize_t i3,
      ssize_t i4,
      ssize_t i5,
      ssize_t i6) const
    {
      return elems_[this->index_1d(i1, i2, i3, i4, i5, i6)];
    }
  };

  template <typename T, size_t Ndims=1>
  struct arr_ref : arr_cref<T, Ndims>
  {
    protected:
      arr_ref() {}

      public:

    arr_ref(
      T& val)
    :
      arr_cref<T, Ndims>(val)
    {}

    template <size_t OtherNdims>
    arr_ref(
      arr_ref<T, OtherNdims> const& other)
    :
      arr_cref<T, Ndims>(other)
    {}

    template <size_t BufferNdims>
    arr_ref(
      T& val,
      dim_data<BufferNdims> const& dims)
    :
      arr_cref<T, Ndims>(val)
    {
      (*this)(dims);
    }

    template <size_t BufferNdims>
    arr_ref(
      T& val,
      dim_data<BufferNdims> const& dims,
      no_fill0_type const&)
    :
      arr_cref<T, Ndims>(val)
    {
      (*this)(dims);
    }

    template <size_t BufferNdims>
    arr_ref(
      T& val,
      dim_data<BufferNdims> const& dims,
      fill0_type const&)
    :
      arr_cref<T, Ndims>(val)
    {
      (*this)(dims);
      std::memset(this->begin(), 0, this->size_1d() * sizeof(T));
    }

    T*
    begin() const { return const_cast<T*>(this->elems_); }

    operator
    T&() const { return *(this->begin()); }

    template <size_t BufferNdims>
    void
    operator()(
      dim_data<BufferNdims> const& dims)
    {
      arr_cref<T, Ndims>::operator()(dims);
    }

    T&
    operator()(
      ssize_t i1) const
    {
      return this->begin()[this->index_1d(i1)];
    }

    T&
    operator()(
      ssize_t i1,
      ssize_t i2) const
    {
      return this->begin()[this->index_1d(i1, i2)];
    }

    T&
    operator()(
      ssize_t i1,
      ssize_t i2,
      ssize_t i3) const
    {
      return this->begin()[this->index_1d(i1, i2, i3)];
    }

    T&
    operator()(
      ssize_t i1,
      ssize_t i2,
      ssize_t i3,
      ssize_t i4) const
    {
      return this->begin()[this->index_1d(i1, i2, i3, i4)];
    }

    T&
    operator()(
      ssize_t i1,
      ssize_t i2,
      ssize_t i3,
      ssize_t i4,
      ssize_t i5) const
    {
      return this->begin()[this->index_1d(i1, i2, i3, i4, i5)];
    }

    T&
    operator()(
      ssize_t i1,
      ssize_t i2,
      ssize_t i3,
      ssize_t i4,
      ssize_t i5,
      ssize_t i6) const
    {
      return this->begin()[this->index_1d(i1, i2, i3, i4, i5, i6)];
    }
  };

} // namespace fem

#endif // GUARD
