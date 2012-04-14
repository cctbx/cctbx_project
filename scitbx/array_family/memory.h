#ifndef SCITBX_ARRAY_FAMILY_MEMORY_H
#define SCITBX_ARRAY_FAMILY_MEMORY_H

/* + On some platforms std::malloc is 16-byte aligned:
 - glibc 2.8 onwards on 64 bit system:
 http://www.gnu.org/s/libc/manual/html_node/Aligned-Memory-Blocks.html
 the version requirement is lifted from the Eigen library
 (c.f. comments in src/Core/util/Memory.h)
 - MacOS X
 - 64-bit Windows
 + POSIX provides posix_memalign if the Advanced Real Time option group is present,
 which is detected by looking at _POSIX_ADVISORY_INFO e.g.
 as per http://pubs.opengroup.org/onlinepubs/009695399/basedefs/xbd_chap02.html
 POSIX itself is detected with _XOPEN_SOURCE as per the same document.

 */
#if (                                                                 \
defined(__GLIBC__)                                              \
&& ((__GLIBC__>= 2 && __GLIBC_MINOR__ >= 8) || __GLIBC__>2)     \
&& defined(__LP64__)                                            \
)                                                                 \
||                                                                \
defined(__APPLE__)                                                \
||                                                                \
defined(_WIN64)

#define SCITBX_AF_HAS_ALIGNED_MALLOC 1
#define SCITBX_AF_USE_STD_FOR_ALIGNED_MALLOC 1

#elif defined(_XOPEN_SOURCE) && _XOPEN_SOURCE >= 600 && defined(_POSIX_ADVISORY_INFO)

#define SCITBX_AF_HAS_ALIGNED_MALLOC 1
#define SCITBX_AF_USE_POSIX_FOR_ALIGNED_MALLOC 1

#elif defined(_MSC_VER)

#define SCITBX_AF_HAS_ALIGNED_MALLOC 1
#define SCITBX_AF_USE_WIN32_FOR_ALIGNED_MALLOC 1

#endif



#ifdef SCITBX_AF_USE_STD_FOR_ALIGNED_MALLOC
#include <cstdlib>
#endif

#ifdef SCITBX_AF_USE_POSIX_FOR_ALIGNED_MALLOC
#include <stdlib.h>
#include <new>
#endif

namespace scitbx { namespace af {

  /// Like std::malloc but the returned pointer is aligned on 16-byte boundaries
  /** This is essential for efficient SSE 2 code on Intel processors */
  void *aligned_malloc(std::size_t n);

  /// Free memory allocated by aligned_malloc
  void aligned_free(void *p);


#ifdef SCITBX_AF_USE_STD_FOR_ALIGNED_MALLOC
  inline
  void *aligned_malloc(std::size_t n) {
    return std::malloc(n);
  }

  inline
  void aligned_free(void *p) {
    std::free(p);
  }
#endif

#ifdef SCITBX_AF_USE_POSIX_FOR_ALIGNED_MALLOC
  inline
  void *aligned_malloc(std::size_t n) {
    void *result;
    if(posix_memalign(&result, 16, n)) throw std::bad_alloc();
    return result;
  }

  inline
  void aligned_free(void *p) {
    std::free(p);
  }
#endif

#ifdef SCITBX_AF_USE_WIN32_FOR_ALIGNED_MALLOC
  inline
  void *aligned_malloc(std::size_t n) {
    return _aligned_malloc(n, 16);
  }

  inline
  void aligned_free(void *p) {
    _aligned_free(p);
  }
#endif

}}



#endif
