#ifndef FEM_UTILS_PATH_HPP
#define FEM_UTILS_PATH_HPP

#include <cstdio>
#include <sys/types.h>
#include <sys/stat.h>
#include <errno.h>

#if defined(_MSC_VER)

#include <io.h> // for _chsize

#define FEM_UTILS_PATH_STRUCT_STAT struct _stat
#define FEM_UTILS_PATH_STAT _stat
#define FEM_UTILS_PATH_FTRUNCATE _chsize

#else

#include <unistd.h> // for ftruncate

#define FEM_UTILS_PATH_STRUCT_STAT struct stat
#define FEM_UTILS_PATH_STAT stat
#define FEM_UTILS_PATH_FTRUNCATE ftruncate

#endif

namespace fem { namespace utils { namespace path {

  inline
  bool
  exists(
    char const* path)
  {
    FEM_UTILS_PATH_STRUCT_STAT buf;
    int stat_result = FEM_UTILS_PATH_STAT(path, &buf);
    return (stat_result == 0 || errno != ENOENT);
  }

  inline
  bool
  truncate_file_at_current_position(
    std::FILE* fp)
  {
    long curr_pos = std::ftell(fp);
    if (curr_pos < 0) return false;
    if (FEM_UTILS_PATH_FTRUNCATE(fileno(fp), curr_pos) != 0) return false;
    std::fflush(fp);
    if (std::fseek(fp, 0L, SEEK_END) != 0) return false;
    return true;
  }

}}} // namespace fem::utils::path

#endif // GUARD
