#ifndef GLTBX_ERROR_UTILS_H
#define GLTBX_ERROR_UTILS_H

#include <gltbx/include_opengl.h>
#include <gltbx/error_utils.h>
#include <string>
#include <stdexcept>
#include <cstdio>

namespace gltbx {

  inline
  std::string
  opengl_error_string(GLenum code)
  {
    std::string result("OpenGL: ");
    const GLubyte* p = gluErrorString(code);
    while(*p) result += *p++;
    return result;
  }

  inline
  void
  handle_error()
  {
    GLenum code = glGetError();
    bool ok = false;
    for(unsigned i=0;i<1000;i++) {
      if (glGetError() == GL_NO_ERROR) {
        ok = true;
        break;
      }
    }
    if (!ok && glGetError() != GL_INVALID_OPERATION) {
      throw std::runtime_error("OpenGL: corrupt glGetError()");
    }
    if (code != GL_NO_ERROR) {
      throw std::runtime_error(opengl_error_string(code));
    }
  }

  namespace detail {

    inline
    std::string
    compose_error_message(
      const char* file,
      long line,
      std::string const& msg)
    {
      char buf[64];
      std::snprintf(buf, sizeof(buf), "%ld", line);
      std::string result = std::string(file) + "(" + buf + ")";
      if (msg.size()) result += std::string(": ") + msg;
      return result;
    }

  }

} // namespace gltbx

#endif // GUARD
