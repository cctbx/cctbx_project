#ifndef FEM_ARR_HPP
#define FEM_ARR_HPP

#include <fem/arr_ref.hpp>

namespace fem {

  template <typename T, size_t Ndims=1>
  struct arr : arr_ref<T, Ndims>
  {
    private:
      arr(
        arr const&);

      arr&
      operator=(
        arr const&);

      public:

    template <size_t BufferNdims>
    explicit
    arr(
      dim_data<BufferNdims> const& dims)
    :
      arr_ref<T, Ndims>(*(new T[dims.size_1d(Ndims)]), dims)
    {}

    template <size_t BufferNdims>
    arr(
      dim_data<BufferNdims> const& dims,
      no_fill0_type const&)
    :
      arr_ref<T, Ndims>(*(new T[dims.size_1d(Ndims)]), dims, no_fill0)
    {}

    template <size_t BufferNdims>
    arr(
      dim_data<BufferNdims> const& dims,
      fill0_type const&)
    :
      arr_ref<T, Ndims>(*(new T[dims.size_1d(Ndims)]), dims, fill0)
    {}

    ~arr()
    {
      delete[] this->elems_;
    }

    operator
    T&() { return *this->begin(); }
  };

  template <typename T, size_t Ndims=1>
  struct cmn_arr
  {
    arr_ref_dims<Ndims> const& dims_;
    T* elems_;

    private:
      cmn_arr(
        cmn_arr const&);

      cmn_arr&
      operator=(
        cmn_arr const&);

      public:

    explicit
    cmn_arr(
      arr_ref_dims<Ndims> const& dims)
    :
      dims_(dims),
      elems_(*(new T[dims.size_1d()]))
    {}

    cmn_arr(
      arr_ref_dims<Ndims> const& dims,
      fill0_type const&)
    :
      dims_(dims),
      elems_(new T[dims.size_1d()])
    {
      std::memset(elems_, 0, dims.size_1d() * sizeof(T));
    }

    ~cmn_arr()
    {
      delete[] elems_;
    }

    template <size_t OtherNdims>
    operator
    arr_ref<T, OtherNdims>() const
    {
      return arr_ref<T, OtherNdims>(*elems_);
    }

    T&
    operator()(
      ssize_t i1,
      ssize_t i2,
      ssize_t i3) const
    {
      return elems_[dims_.index_1d(i1, i2, i3)];
    }
  };

} // namespace fem

#endif // GUARD
