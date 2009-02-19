#ifndef GLTBX_UTIL_H
#define GLTBX_UTIL_H

#include <gltbx/include_opengl.h>

namespace gltbx {

  struct scope_push_matrix
  {
    scope_push_matrix() { glPushMatrix(); }

    ~scope_push_matrix() { glPopMatrix(); }

  };

} // gltbx

#define GLTBX_SCOPE_PUSH_MATRIX gltbx::scope_push_matrix sentry

#endif // GUARD
