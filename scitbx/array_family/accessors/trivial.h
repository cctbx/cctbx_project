#ifndef SCITBX_ARRAY_FAMILY_ACCESSORS_TRIVIAL_H
#define SCITBX_ARRAY_FAMILY_ACCESSORS_TRIVIAL_H

#include <cstddef>

namespace scitbx { namespace af {

  class trivial_accessor
  {
    public:
      typedef std::size_t index_type;
      struct index_value_type {};

      trivial_accessor() : size_(0) {}

      trivial_accessor(std::size_t const& n) : size_(n) {}

      std::size_t size_1d() const { return size_; }

      std::size_t operator()(std::size_t i) const { return i; }

    protected:
      std::size_t size_;
  };

}} // namespace scitbx::af

#endif // SCITBX_ARRAY_FAMILY_ACCESSORS_TRIVIAL_H
