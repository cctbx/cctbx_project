#include <scitbx/array_family/ref.h>
#include <scitbx/serialization/double_buffered.h>

namespace scitbx { namespace af { namespace boost_python {
namespace pickle_double_buffered {

  struct to_string : scitbx::serialization::double_buffered::to_string
  {
    using scitbx::serialization::double_buffered::to_string::operator<<;

    template <typename ElementType,
              typename AccessorType>
    to_string& operator<<(af::const_ref<ElementType, AccessorType> const& val)
    {
      for(std::size_t i=0;i<val.size();i++) {
        *this << val[i];
      }
      return *this;
    }
  };

  struct from_string : scitbx::serialization::double_buffered::from_string
  {
    from_string(const char* str_ptr)
    : scitbx::serialization::double_buffered::from_string(str_ptr)
    {}

    using scitbx::serialization::double_buffered::from_string::operator>>;

    template <typename ElementType>
    from_string& operator>>(ref<ElementType> val)
    {
      for(std::size_t i=0;i<val.size();i++) {
        *this >> val[i];
      }
      return *this;
    }
  };

}}}} // namespace scitbx::af::boost_python::pickle_double_buffered
