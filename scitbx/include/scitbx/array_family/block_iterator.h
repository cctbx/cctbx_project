#ifndef SCITBX_ARRAY_FAMILY_BLOCK_ITERATOR_H
#define SCITBX_ARRAY_FAMILY_BLOCK_ITERATOR_H

#include <scitbx/array_family/ref.h>
#include <scitbx/error.h>

namespace scitbx { namespace af {

  template <typename ElementType>
  class const_block_iterator
  {
    public:
      const_block_iterator(
        af::const_ref<ElementType> const& array,
        std::string const& error_message)
      :
        array_(array),
        error_message_(error_message),
        i_(0)
      {}

      const ElementType*
      operator()(std::size_t block_size)
      {
        if (i_ + block_size > array_.size()) {
          throw error(error_message_);
        }
        const ElementType* result = &array_[i_];
        i_ += block_size;
        return result;
      }

      ElementType const&
      operator()() { return *(*this)(1); }

      bool
      is_at_end() const { return i_ == array_.size(); }

    private:
      af::const_ref<ElementType> array_;
      std::string error_message_;
      std::size_t i_;
  };

  template <typename ElementType>
  class block_iterator
  {
    public:
      block_iterator(
        af::ref<ElementType> const& array,
        std::string const& error_message)
      :
        array_(array),
        error_message_(error_message),
        i_(0)
      {}

      ElementType*
      operator()(std::size_t block_size)
      {
        if (i_ + block_size > array_.size()) {
          throw error(error_message_);
        }
        ElementType* result = &array_[i_];
        i_ += block_size;
        return result;
      }

      ElementType&
      operator()() { return *(*this)(1); }

      bool
      is_at_end() const { return i_ == array_.size(); }

    private:
      af::ref<ElementType> array_;
      std::string error_message_;
      std::size_t i_;
  };

}} // namespace scitbx::af

#endif // SCITBX_ARRAY_FAMILY_BLOCK_ITERATOR_H
