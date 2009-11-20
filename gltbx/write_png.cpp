
// Basic PNG output.  This is mostly copied from example.c (which is
// public domain) in the main libpng distribution.
//
// I've also used the PyMOL code as a guide to adapting this, since it
// clearly uses the exact same source and has some of the comments left
// but has been tailored to PyMOL's uses.  Copyright appended at bottom.
//

/*
A* -------------------------------------------------------------------
B* This file contains source code for the PyMOL computer program
C* copyright 1998-2000 by Warren Lyford Delano of DeLano Scientific.
D* -------------------------------------------------------------------
E* It is unlawful to modify or remove this copyright notice.
F* -------------------------------------------------------------------
G* Please see the accompanying LICENSE file for further information.
H* -------------------------------------------------------------------
I* Additional authors of this source file include:
-*
-*
-*
Z* -------------------------------------------------------------------
*/

#include <gltbx/include_opengl.h>
#include <gltbx/error.h>
#include <boost/scoped_array.hpp>
#include <boost/shared_ptr.hpp>
#include <cstdio>
#include <cstdlib>

#ifdef GLTBX_HAVE_LIBPNG
extern "C" {
#include <png.h>
}
#endif

namespace gltbx { namespace viewer_utils {

bool write_png (const char* file_name, int width, int height) {
#ifdef GLTBX_HAVE_LIBPNG

  int i;
  int bit_depth = 8;
  int bytes_per_pixel = 3;

  png_structp png_ptr;
  png_infop info_ptr;
  png_colorp palette;

  boost::shared_ptr<std::FILE> fp(std::fopen(file_name, "wb"), std::fclose);
  if (fp.get() == 0) {
    throw std::runtime_error("Failed to open PNG file for writing.");
  }

  png_ptr = png_create_write_struct(PNG_LIBPNG_VER_STRING, 0, 0, 0);

  if (png_ptr == 0) {
    throw std::runtime_error("Couldn't allocate PNG output struct.");
  }

  info_ptr = png_create_info_struct(png_ptr);
  if (info_ptr == 0) {
    png_destroy_write_struct(&png_ptr,  png_infopp_NULL);
    throw std::runtime_error("Couldn't allocate PNG info struct.");
  }

  if (setjmp(png_jmpbuf(png_ptr))) {
    /* If we get here, we had a problem reading the file */
    png_destroy_write_struct(&png_ptr, &info_ptr);
    throw std::runtime_error("Couldn't read PNG output.");
  }

  png_init_io(png_ptr, fp.get());

  /* We don't really need an alpha channel - these images will be used for
   * later viewing or printing, not any application where transparency
   * might be valuable.
   */
  png_set_IHDR(png_ptr, info_ptr, width, height, bit_depth, PNG_COLOR_TYPE_RGB,
    PNG_INTERLACE_NONE, PNG_COMPRESSION_TYPE_BASE, PNG_FILTER_TYPE_BASE);
  png_write_info(png_ptr, info_ptr);

  /* As with writing TIFF files, we have to output rows in reverse
   * relative to how OpenGL stores them.
   */
  boost::scoped_array<png_bytep> row_pointers(new png_bytep[height]);
  boost::scoped_array<GLubyte> image(new GLubyte[width*height*3]);
  glPixelStorei(GL_PACK_ALIGNMENT, 1);
  glReadPixels(0, 0, width, height, GL_RGB, GL_UNSIGNED_BYTE, image.get());
  GLubyte* p = image.get();
  for (i = height - 1; i >= 0; i--) {
    row_pointers[i] = p;
    p += width * sizeof(GLubyte) * 3;
  }

  png_write_image(png_ptr, row_pointers.get());
  png_write_end(png_ptr, info_ptr);
  png_destroy_write_struct(&png_ptr, &info_ptr);

  return true;
#else
  return false;
#endif
}

} // namespace viewer_utils
} // namespace gltbx
