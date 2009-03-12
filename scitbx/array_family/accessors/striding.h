#ifndef SCITBX_ARRAY_FAMILY_ACCESSORS_STRIDING_H
#define SCITBX_ARRAY_FAMILY_ACCESSORS_STRIDING_H

#include <cstddef>

namespace scitbx { namespace af {

  class striding_linear_accessor
  {
    public:
      typedef std::size_t index_type;
      struct index_value_type {};

      striding_linear_accessor() : size_(0), stride_(0) {}

      striding_linear_accessor(std::size_t n, std::size_t stride)
        : size_(n), stride_(stride)
      {}

      std::size_t size_1d() const { return size_; }

      std::size_t operator()(std::size_t i) const { return i*stride_; }

    protected:
      std::size_t size_, stride_;
  };

}} // namespace scitbx::af

#endif // GUARD
