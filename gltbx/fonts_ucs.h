#ifndef GLTBX_FONTS_UCS_H
#define GLTBX_FONTS_UCS_H

#include <gltbx/error.h>
#include <boost/scoped_array.hpp>
#include <map>
#include <string>
#include <cstring>

namespace gltbx { namespace fonts { namespace ucs {

  struct encoding_range {
    unsigned begin;
    unsigned count;
  };

  struct bitmap_font_record {
    const char *short_name;
    const char *full_name;
    unsigned width;
    unsigned height;
    int xorig;
    int yorig;
    unsigned number_of_chars;
    unsigned raw_bitmaps_size;
    const unsigned char* raw_bitmaps;
    const encoding_range* encoding_ranges;
  };

  extern bitmap_font_record bitmap_8x13;
  extern bitmap_font_record bitmap_9x15;
  extern bitmap_font_record bitmap_10x20;

  template <typename UnsignedInt2Type>
  class bitmap
  {
    public:
      typedef UnsignedInt2Type unsigned_int2_type;

      bitmap() {}

      bitmap(const char* short_name)
      :
        have_call_lists_base(false),
        call_lists_base(0)
      {
        GLTBX_ASSERT(sizeof(UnsignedInt2Type) == 2);
        if      (std::strcmp(short_name, "8x13") == 0) {
          font_record = &bitmap_8x13;
        }
        else if (std::strcmp(short_name, "9x15") == 0) {
          font_record = &bitmap_9x15;
        }
        else if (std::strcmp(short_name, "10x20") == 0) {
          font_record = &bitmap_10x20;
        }
        else {
          throw std::runtime_error(
            std::string("Unknown bitmap font: ") + short_name);
        }
      }

      const char*
      short_name() const { return font_record->short_name; }

      const char*
      full_name() const { return font_record->full_name; }

      unsigned
      width() const { return font_record->width; }

      unsigned
      height() const { return font_record->height; }

      int
      xorig() const { return font_record->xorig; }

      int
      yorig() const { return font_record->yorig; }

      void
      setup_call_lists() const
      {
        if (have_call_lists_base) return;
        unsigned n_chars = font_record->number_of_chars;
        GLTBX_ASSERT(font_record->raw_bitmaps_size % n_chars == 0);
        unsigned bytes_per_char = font_record->raw_bitmaps_size / n_chars;
        call_lists_base = glGenLists(n_chars);
        have_call_lists_base = true;
        glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
        const unsigned char* raw_bitmap = font_record->raw_bitmaps;
        unsigned i_char;
        for(i_char=0;i_char<n_chars;i_char++) {
          glNewList(i_char+call_lists_base, GL_COMPILE);
          glBitmap(
            font_record->width,
            font_record->height,
            font_record->xorig,
            font_record->yorig,
            font_record->width,
            0,
            raw_bitmap);
          glEndList();
          raw_bitmap += bytes_per_char;
        }
        i_char = 0;
        for(const encoding_range*
            er = font_record->encoding_ranges;
            er->count != 0;
            er++) {
          GLTBX_ASSERT(i_char+er->count <= n_chars);
          UnsignedInt2Type encoding = static_cast<UnsignedInt2Type>(er->begin);
          for(unsigned i=0;i<er->count;i++) {
            encoding_to_bitmap_indices[encoding+i] = i_char++;
          }
        }
        GLTBX_ASSERT(i_char == n_chars);
      }

      GLuint
      bitmap_index(UnsignedInt2Type const& encoding) const
      {
        GLTBX_ASSERT(have_call_lists_base);
        typename std::map<UnsignedInt2Type, unsigned>::const_iterator
          pair = encoding_to_bitmap_indices.find(encoding);
        if (pair == encoding_to_bitmap_indices.end()) return 0;
        return static_cast<GLuint>(pair->second);
      }

      void
      render_bitmap_indices(unsigned size, const GLuint* bitmap_indices) const
      {
        glPushAttrib(GL_LIST_BIT);
        glListBase(call_lists_base);
        glCallLists(size, GL_UNSIGNED_INT, bitmap_indices);
        glPopAttrib();
      }

      void
      render_encodings(unsigned size, const UnsignedInt2Type* encodings) const
      {
        boost::scoped_array<GLuint> bitmap_indices(new GLuint[size]);
        GLuint* bi = bitmap_indices.get();
        for(unsigned i_char=0;i_char<size;i_char++) {
          *bi++ = bitmap_index(encodings[i_char]);
        }
        render_bitmap_indices(size, bitmap_indices.get());
      }

      void
      render_string(std::string const& string) const
      {
        boost::scoped_array<GLuint> bitmap_indices(new GLuint[string.size()]);
        GLuint* bi = bitmap_indices.get();
        for(unsigned i_char=0;i_char<string.size();i_char++) {
          char c = string[i_char];
          UnsignedInt2Type encoding = static_cast<UnsignedInt2Type>(
            *(reinterpret_cast<const unsigned char*>(&c)));
          *bi++ = bitmap_index(encoding);
        }
        render_bitmap_indices(string.size(), bitmap_indices.get());
      }

      // XXX TODO: void glDeleteLists ( GLuint list, GLsizei range )

    protected:
      bitmap_font_record* font_record;
      mutable std::map<UnsignedInt2Type, unsigned> encoding_to_bitmap_indices;
      mutable bool have_call_lists_base;
      mutable GLuint call_lists_base;
  };

}}} // namespace gltbx::fonts::ucs

#endif // GUARD
