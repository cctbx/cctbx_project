#include <scitbx/error.h>
#include <cstdio>
#include <cstring>
#include <string>
#include <vector>

namespace scitbx { namespace af { namespace boost_python {

  // True if format_string is one printf directive and nothing else: "%",
  // optional flags, width and precision, and one of the given conversion
  // characters.
  inline bool
  is_printf_directive(std::string const& format_string, const char* conversions)
  {
    std::size_t n = format_string.size();
    if (n < 2 || format_string[0] != '%') return false;
    std::size_t i = 1;
    while (i < n && format_string[i] != '\0'
           && std::strchr("-+ #0", format_string[i]) != 0) i++;
    while (i < n && format_string[i] >= '0' && format_string[i] <= '9') i++;
    if (i < n && format_string[i] == '.') {
      i++;
      while (i < n && format_string[i] >= '0' && format_string[i] <= '9') i++;
    }
    return i + 1 == n && format_string[i] != '\0'
        && std::strchr(conversions, format_string[i]) != 0;
  }

  // snprintf of one value with a directive is_printf_directive accepted.
  template <typename ValueType>
  std::string
  snprintf_string(const char* directive, ValueType value)
  {
    char buffer[64];
    int n = std::snprintf(buffer, sizeof(buffer), directive, value);
    SCITBX_ASSERT(n >= 0);
    if (static_cast<std::size_t>(n) < sizeof(buffer)) {
      return std::string(buffer, n);
    }
    std::vector<char> big(n + 1);
    std::snprintf(&big[0], big.size(), directive, value);
    return std::string(&big[0], n);
  }

  template <typename ElementType, typename UnsignedType>
  boost::python::object
  add_selected_unsigned_a(
    boost::python::object const& self,
    af::const_ref<UnsignedType> const& indices,
    af::const_ref<ElementType> const& values)
  {
    boost::python::extract<af::ref<ElementType> > a_proxy(self);
    af::ref<ElementType> a = a_proxy();
    SCITBX_ASSERT(indices.size() == values.size());
    for(std::size_t i=0;i<indices.size();i++) {
      SCITBX_ASSERT(indices[i] < a.size());
      a[indices[i]] += values[i];
    }
    return self;
  }

  template <typename ElementType, typename UnsignedType>
  boost::python::object
  add_selected_unsigned_s(
    boost::python::object const& self,
    af::const_ref<UnsignedType> const& indices,
    ElementType const& value)
  {
    boost::python::extract<af::ref<ElementType> > a_proxy(self);
    af::ref<ElementType> a = a_proxy();
    for(std::size_t i=0;i<indices.size();i++) {
      SCITBX_ASSERT(indices[i] < a.size());
      a[indices[i]] += value;
    }
    return self;
  }

  template <typename ElementType>
  boost::python::object
  add_selected_bool_a(
    boost::python::object const& self,
    af::const_ref<bool> const& flags,
    af::const_ref<ElementType> const& values)
  {
    boost::python::extract<af::ref<ElementType> > a_proxy(self);
    af::ref<ElementType> a = a_proxy();
    SCITBX_ASSERT(a.size() == flags.size());
    if (values.size() == flags.size()) {
      ElementType* ai = a.begin();
      const bool* fi = flags.begin();
      const ElementType* ni = values.begin();
      const ElementType* ne = values.end();
      while (ni != ne) {
        if (*fi++) *ai += *ni;
        ai++;
        ni++;
      }
    }
    else {
      std::size_t i_value = 0;
      for(std::size_t i=0;i<flags.size();i++) {
        if (flags[i]) {
          SCITBX_ASSERT(i_value < values.size());
          a[i] += values[i_value];
          i_value++;
        }
      }
      SCITBX_ASSERT(i_value == values.size());
    }
    return self;
  }


  template <typename ElementType>
  boost::python::object
  add_selected_bool_s(
    boost::python::object const& self,
    af::const_ref<bool, flex_grid<> > const& flags,
    ElementType const& value)
  {
    boost::python::extract<af::ref<ElementType, flex_grid<> > > a_proxy(self);
    af::ref<ElementType, flex_grid<> > a = a_proxy();
    SCITBX_ASSERT(a.accessor() == flags.accessor());
    for(std::size_t i=0;i<flags.size();i++) {
      if (flags[i]) a[i] += value;
    }
    return self;
  }

}}} // namespace scitbx::af::boost_python
