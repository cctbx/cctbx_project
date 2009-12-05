// Implementations limited to just what was needed at the time.
// Please expand as necessary.

#ifndef GLTBX_SPECIAL_WRAPPER_SUPPORT_H
#define GLTBX_SPECIAL_WRAPPER_SUPPORT_H

#include <gltbx/error.h>
#include <boost/python/ssize_t.hpp>

namespace gltbx { namespace gl {

  inline
  boost::python::ssize_t
  glReadPixels_pixels_expected_size(
    GLsizei width,
    GLsizei height,
    GLenum format,
    GLenum type)
  {
    if (format == GL_RGB && (type == GL_BYTE || type == GL_UNSIGNED_BYTE)) {
      return 3 * width * height;
    }
    throw GLTBX_ERROR("Sorry: incomplete glReadPixels() support.");
  }

}} // namespace gltbx::gl

#endif // GUARD
