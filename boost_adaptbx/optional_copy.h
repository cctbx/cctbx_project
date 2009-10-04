#ifndef BOOST_ADAPTBX_OPTIONAL_COPY_H
#define BOOST_ADAPTBX_OPTIONAL_COPY_H

namespace boost_adaptbx {

  //! Optional allocation of value with new, with value-copy semantics.
  template <typename ValueType>
  class optional_copy
  {
    public:
      typedef ValueType value_type;

    protected:
      ValueType* ptr_;

    public:
      optional_copy() : ptr_(0) {}

      ~optional_copy()
      {
        if (ptr_ != 0) delete ptr_;
      }

      optional_copy(optional_copy const& other)
      :
        ptr_(other.ptr_ == 0 ? 0 : new ValueType(*(other.ptr_)))
      {}

      optional_copy&
      operator=(optional_copy const& other)
      {
        if (ptr_ != 0) delete ptr_;
        ptr_ = (other.ptr_ == 0 ? 0 : new ValueType(*(other.ptr_)));
        return *this;
      }

      explicit
      optional_copy(ValueType const& value)
      :
        ptr_(new ValueType(value))
      {}

      optional_copy&
      operator=(ValueType const& value)
      {
        if (ptr_ != 0) delete ptr_;
        ptr_ = new ValueType(value);
        return *this;
      }

      void
      release()
      {
        if (ptr_ != 0) delete ptr_;
        ptr_ = 0;
      }

      operator bool() const { return (ptr_ != 0); }

      ValueType*
      get()       { return ptr_; }

      ValueType const*
      get() const { return ptr_; }

      ValueType*
      operator->()       { return ptr_; }

      ValueType const*
      operator->() const { return ptr_; }

      ValueType&
      operator*()       { return *ptr_; }

      ValueType const&
      operator*() const { return *ptr_; }

      typename ValueType::value_type&
      operator[](typename ValueType::size_type const& i)
      {
        return (*ptr_)[i];
      }

      typename ValueType::value_type const&
      operator[](typename ValueType::size_type const& i) const
      {
        return (*ptr_)[i];
      }
  };

} // namespace boost_adaptbx

#endif // GUARD
