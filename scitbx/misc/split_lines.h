#ifndef SCITBX_MISC_SPLIT_LINES_H
#define SCITBX_MISC_SPLIT_LINES_H

#include <scitbx/array_family/shared.h>
#include <boost/shared_array.hpp>
#include <string>

namespace scitbx { namespace misc {

  // Based on Python-2.4.2/Objects/stringobject.c string_splitlines()
  inline
  af::shared<std::string>
  split_lines(
    const char* data,
    std::size_t size,
    bool keep_ends=false,
    bool count_lines_first=true)
  {
    af::shared<std::string> result;
    std::size_t n_lines = 0;
    for(unsigned i_pass=(count_lines_first ? 0 : 1);i_pass<2;i_pass++) {
      int i = 0;
      int j = 0;
      while (i < size) {
        while (i < size && data[i] != '\n' && data[i] != '\r') i++;
        int eol = i;
        if (i < size) {
          if (data[i] == '\r' && i+1 < size && data[i+1] == '\n') {
            i += 2;
          }
          else {
            i++;
          }
          if (keep_ends) {
            eol = i;
          }
        }
        if (i_pass == 0) {
          n_lines++;
        }
        else {
          result.push_back(std::string(data+j, data+eol));
        }
        j = i;
      }
      if (j < size) {
        if (i_pass == 0) {
          n_lines++;
        }
        else {
          result.push_back(std::string(data+j, data+size));
        }
      }
      if (i_pass == 0) {
        result.reserve(n_lines);
      }
    }
    return result;
  }

  struct char_buffer
  {
    boost::shared_array<char> data;
    std::size_t size;

    char_buffer() {}

    char_buffer(std::size_t size_)
    :
      data(new char[size_]),
      size(size_)
    {}
  };

  inline
  af::shared<std::string>
  split_lines(
    char_buffer const& buffer,
    bool keep_ends=false,
    bool count_lines_first=true)
  {
    return split_lines(
      buffer.data.get(), buffer.size, keep_ends, count_lines_first);
  }

  inline
  af::shared<std::string>
  split_lines(
    std::string const& buffer,
    bool keep_ends=false,
    bool count_lines_first=true)
  {
    return split_lines(
      buffer.data(), buffer.size(), keep_ends, count_lines_first);
  }

}} // namespace scitbx::misc

#endif // SCITBX_MISC_SPLIT_LINES_H
