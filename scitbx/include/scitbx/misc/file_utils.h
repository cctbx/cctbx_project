#ifndef SCITBX_MISC_FILE_UTILS_H
#define SCITBX_MISC_FILE_UTILS_H

#include <scitbx/misc/split_lines.h>
#include <boost/shared_ptr.hpp>
#include <string>
#include <stdexcept>
#include <stdio.h>
#include <sys/types.h>
#include <sys/stat.h>

#if defined(__ALPHA) && defined(__DECCXX)
extern "C" { extern int fileno __((FILE *)); }
#endif

#if defined(BOOST_MSVC)
#define SCITBX_MISC_FILE_UTILS_FILENO       _fileno
#define SCITBX_MISC_FILE_UTILS_STRUCT_STAT  struct _stat
#define SCITBX_MISC_FILE_UTILS_FSTAT        _fstat
#else
#define SCITBX_MISC_FILE_UTILS_FILENO       fileno
#define SCITBX_MISC_FILE_UTILS_STRUCT_STAT  struct stat
#define SCITBX_MISC_FILE_UTILS_FSTAT        fstat
#endif

namespace scitbx { namespace misc {

  inline
  char_buffer
  file_to_char_buffer(std::string const& file_name, bool binary=true)
  {
    FILE* fp = fopen(file_name.c_str(), (binary ? "rb" : "r"));
    if (!fp) {
      throw std::runtime_error(
        "Cannot open file for reading: \"" + file_name + "\"");
    }
    boost::shared_ptr<FILE> fp_holder(fp, fclose);
    int fn = SCITBX_MISC_FILE_UTILS_FILENO(fp);
    if (fn < 0) {
      throw std::runtime_error(
        "fileno() failed for open file: \"" + file_name + "\"");
    }
    SCITBX_MISC_FILE_UTILS_STRUCT_STAT stat_buffer;
    if (SCITBX_MISC_FILE_UTILS_FSTAT(fn, &stat_buffer) != 0) {
      throw std::runtime_error(
        "fstat() failed for open file: \"" + file_name + "\"");
    }
    char_buffer result(static_cast<std::size_t>(stat_buffer.st_size));
    std::size_t n = fread(
      static_cast<void*>(result.data.get()), 1U, result.size, fp);
    if (n != result.size || ferror(fp)) {
      throw std::runtime_error(
        "Error reading file: \"" + file_name + "\"");
    }
    return result;
  }

  inline
  af::shared<std::string>
  file_to_lines(
    std::string const& file_name,
    bool binary=true,
    bool keep_ends=false,
    bool count_lines_first=true)
  {
    return split_lines(
      file_to_char_buffer(file_name, binary),
      keep_ends,
      count_lines_first);
  }

}} // namespace scitbx::misc

#endif // SCITBX_MISC_FILE_UTILS_H
