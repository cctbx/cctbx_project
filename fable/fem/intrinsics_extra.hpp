#ifndef FEM_INTRINSICS_EXTRA_HPP
#define FEM_INTRINSICS_EXTRA_HPP

#include <fem/str.hpp>
#include <cstdio>
#include <cstdlib>
#include <ctime>

namespace fem {

  inline
  void
  getenv(
    str_cref key,
    str_ref result)
  {
    std::string k = utils::strip_leading_and_trailing_blank_padding(key);
    char* v = std::getenv(k.c_str());
    result = (v == 0 ? " " : v);
  }

  inline
  void
  date(
    str_ref result)
  {
    static const char* months[] = {
      "Jan", "Feb", "Mar", "Apr", "May", "Jun",
      "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"};
    std::time_t now = std::time(0);
    std::tm const* tm(std::localtime(&now));
    str<10> buf;
    std::sprintf(buf.elems, "%02d-%s-%02d",
      tm->tm_mday,
      months[tm->tm_mon],
      tm->tm_year % 100);
    result = buf;
  }

  inline
  void
  time(
    str_ref result)
  {
    std::time_t now = std::time(0);
    std::tm const* tm(std::localtime(&now));
    str<10> buf;
    std::sprintf(buf.elems, "%02d:%02d:%02d",
      tm->tm_hour,
      tm->tm_min,
      tm->tm_sec);
    result = buf;
  }

} // namespace fem

#endif // GUARD
