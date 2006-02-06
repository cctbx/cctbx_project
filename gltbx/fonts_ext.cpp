#include <boost/python/module.hpp>
#include <boost/python/class.hpp>
#include <boost/python/args.hpp>

#include <gltbx/fonts_ucs.h>

#if !defined(BOOST_NO_STD_WSTRING) && defined(Py_USING_UNICODE)
# define GLTBX_USING_UNICODE
#endif

namespace gltbx { namespace fonts {

namespace {

  struct ucs_bitmap_wrappers
  {
    typedef ucs::bitmap<unsigned short> w_t;

#if defined(GLTBX_USING_UNICODE)
    static void
    render_wstring(w_t const& self, std::wstring const& string)
    {
      boost::scoped_array<GLuint> bitmap_indices(new GLuint[string.size()]);
      GLuint* bi = bitmap_indices.get();
      for(unsigned i_char=0;i_char<string.size();i_char++) {
        w_t::unsigned_int2_type encoding
          = static_cast<w_t::unsigned_int2_type>(string[i_char]);
        *bi++ = self.bitmap_index(encoding);
      }
      self.render_bitmap_indices(string.size(), bitmap_indices.get());
    }
#endif

    static void
    wrap()
    {
      using namespace boost::python;
      class_<w_t>("ucs_bitmap", no_init)
        .def(init<const char*>((arg_("short_name"))))
        .def("short_name", &w_t::short_name)
        .def("full_name", &w_t::full_name)
        .def("width", &w_t::width)
        .def("height", &w_t::height)
        .def("xorig", &w_t::xorig)
        .def("yorig", &w_t::yorig)
        .def("setup_call_lists", &w_t::setup_call_lists)
#if defined(GLTBX_USING_UNICODE)
        .def("render_string", render_wstring, (arg_("string")))
#endif
        .def("render_string", &w_t::render_string, (arg_("string")))
      ;
    }
  };

}

  void
  init_module()
  {
    ucs_bitmap_wrappers::wrap();
  }

}} // namespace gltbx::fonts

BOOST_PYTHON_MODULE(gltbx_fonts_ext)
{
  gltbx::fonts::init_module();
}
