#ifndef SCITBX_SERIALIZATION_DOUBLE_BUFFERED_H
#define SCITBX_SERIALIZATION_DOUBLE_BUFFERED_H

#include <scitbx/serialization/single_buffered.h>
#include <scitbx/type_holder.h>

namespace scitbx { namespace serialization { namespace double_buffered {

  struct to_string
  {
    std::string buffer;

    to_string& operator<<(bool const& val)
    {
      if (val) buffer += "1";
      else     buffer += "0";
      return *this;
    }

    to_string& operator<<(short const& val)
    {
      char buf[64];
      buffer.append(buf, base_256::to_string(buf, val));
      return *this;
    }

    to_string& operator<<(unsigned short const& val)
    {
      char buf[64];
      buffer.append(buf, base_256::to_string(buf, val));
      return *this;
    }

    to_string& operator<<(int const& val)
    {
      char buf[64];
      buffer.append(buf, base_256::to_string(buf, val));
      return *this;
    }

    to_string& operator<<(unsigned int const& val)
    {
      char buf[64];
      buffer.append(buf, base_256::to_string(buf, val));
      return *this;
    }

    to_string& operator<<(long const& val)
    {
      char buf[64];
      buffer.append(buf, base_256::to_string(buf, val));
      return *this;
    }

    to_string& operator<<(unsigned long const& val)
    {
      char buf[64];
      buffer.append(buf, base_256::to_string(buf, val));
      return *this;
    }

    to_string& operator<<(long long const& val)
    {
      char buf[64];
      buffer.append(buf, base_256::to_string(buf, val));
      return *this;
    }

    to_string& operator<<(unsigned long long const& val)
    {
      char buf[64];
      buffer.append(buf, base_256::to_string(buf, val));
      return *this;
    }

    to_string& operator<<(float const& val)
    {
      char buf[64];
      buffer.append(buf, base_256::to_string(buf, val));
      return *this;
    }

    to_string& operator<<(double const& val)
    {
      char buf[64];
      buffer.append(buf, base_256::to_string(buf, val));
      return *this;
    }

    template <typename FloatType>
    to_string& operator<<(std::complex<FloatType> const& val)
    {
      return *this << val.real() << val.imag();
    }

    to_string& operator<<(std::string const& val)
    {
      *this << val.size();
      buffer += val;
      return *this;
    }

    to_string& operator<<(const char* val)
    {
      *this << std::string(val);
      return *this;
    }

  };

  struct from_string
  {
    const char* str_ptr;

    from_string()
    : str_ptr(0)
    {
    }

    from_string(const char* str_ptr_)
    : str_ptr(str_ptr_)
    {
      SCITBX_ASSERT(str_ptr != 0);
    }

    void assert_end() const
    {
      SCITBX_ASSERT(*str_ptr == 0);
    }

    template <typename ValueType>
    ValueType get_value(type_holder<ValueType>)
    {
      single_buffered::from_string<ValueType> proxy(str_ptr);
      str_ptr = proxy.end;
      return proxy.value;
    }

    from_string& operator>>(std::string& val)
    {
      val = get_value(type_holder<std::string>());
      return *this;
    }

    from_string& operator>>(bool& val)
    {
      val = get_value(type_holder<bool>());
      return *this;
    }

    from_string& operator>>(short& val)
    {
      val = get_value(type_holder<short>());
      return *this;
    }

    from_string& operator>>(unsigned short& val)
    {
      val = get_value(type_holder<unsigned short>());
      return *this;
    }

    from_string& operator>>(int& val)
    {
      val = get_value(type_holder<int>());
      return *this;
    }

    from_string& operator>>(unsigned int& val)
    {
      val = get_value(type_holder<unsigned int>());
      return *this;
    }

    from_string& operator>>(long& val)
    {
      val = get_value(type_holder<long>());
      return *this;
    }

    from_string& operator>>(unsigned long& val)
    {
      val = get_value(type_holder<unsigned long>());
      return *this;
    }

    from_string& operator>>(long long& val)
    {
      val = get_value(type_holder<long long>());
      return *this;
    }

    from_string& operator>>(unsigned long long& val)
    {
      val = get_value(type_holder<unsigned long long>());
      return *this;
    }

    from_string& operator>>(float& val)
    {
      val = get_value(type_holder<float>());
      return *this;
    }

    from_string& operator>>(double& val)
    {
      val = get_value(type_holder<double>());
      return *this;
    }

    template <typename FloatType>
    from_string& operator>>(std::complex<FloatType>& val)
    {
      val = get_value(type_holder<std::complex<FloatType> >());
      return *this;
    }

  };

}}} // namespace scitbx::serialization::double_buffered

#endif // SCITBX_SERIALIZATION_DOUBLE_BUFFERED_H
