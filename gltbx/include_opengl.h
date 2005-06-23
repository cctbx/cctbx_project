#if defined(_MSC_VER)
# include <windows.h>
#endif
#if defined(__APPLE_CC__)
# include <OpenGL/gl.h>
# include <OpenGL/glu.h>
#else
# include <GL/gl.h>
# include <GL/glu.h>
#endif

#ifdef __cplusplus
extern "C" {
#endif
#if defined(_MSC_VER)
  typedef void (__stdcall *glu_function_pointer)();
#elif defined(__APPLE_CC__)
  typedef GLvoid (*glu_function_pointer)(...);
#else
  typedef GLvoid (*glu_function_pointer)();
#endif
#ifdef __cplusplus
}
#endif
