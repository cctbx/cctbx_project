#ifndef SCITBX_MISC_SPLIT_LINES_H
#define SCITBX_MISC_SPLIT_LINES_H

#include <scitbx/array_family/shared.h>
#include <string>

namespace scitbx { namespace misc {

  // Based on Python-2.4.2/Objects/stringobject.c string_splitlines()
  af::shared<std::string>
  split_lines(
    const char* data,
    std::size_t len,
    bool keep_ends=false,
    bool count_lines_first=true)
  {
    af::shared<std::string> result;
    std::size_t n_lines = 0;
    for(unsigned i_pass=(count_lines_first ? 0 : 1);i_pass<2;i_pass++) {
      int i = 0;
      int j = 0;
      while (i < len) {
        while (i < len && data[i] != '\n' && data[i] != '\r') i++;
        int eol = i;
        if (i < len) {
          if (data[i] == '\r' && i+1 < len && data[i+1] == '\n') {
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
      if (j < len) {
        if (i_pass == 0) {
          n_lines++;
        }
        else {
          result.push_back(std::string(data+j, data+len));
        }
      }
      if (i_pass == 0) {
        result.reserve(n_lines);
      }
    }
    return result;
  }

}} // namespace scitbx::misc

#endif // SCITBX_MISC_SPLIT_LINES_H
