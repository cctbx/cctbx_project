#include <cctbx/boost_python/flex_fwd.h>
#include <cctbx/uctbx.h>
#include <scitbx/array_family/flex_types.h>
#include <boost/python/class.hpp>
#include <boost/python/args.hpp>
#include <boost/python/module.hpp>
#include <fstream>
#include <cstdlib>
#include <cstdio>

namespace af = scitbx::af;

namespace {

  template <unsigned Width>
  struct format_e
  {
    static void
    throw_error()
    {
      throw scitbx::error("Floating-point value too large for format.");
    }

    format_e(const char* fmt, double val)
    {
#if !defined(BOOST_MSVC)
      s = buf;
      std::sprintf(buf, fmt, val);
      if (*(s + Width)) throw_error();
#else
      s = buf + 1;
      std::sprintf(s, fmt, val);
      char* p = s + Width;
      if (*p) {
        p++;
        if (*p) throw_error();
      }
      else {
        s--;
        *s = ' ';
      }
      if (*(p-3) != '0') throw_error();
      *(p-3) = *(p-2);
      *(p-2) = *(p-1);
      *(p-1) = '\0';
#endif
    }

    char buf[32];
    char* s;
  };

  class XplorMap
  {
    public:
      XplorMap() {}

      af::versa<double, af::flex_grid<> >
      ReadXplorMap(
        std::string const& file_name,
        std::size_t n_header_lines,
        af::flex_grid<> const& grid)
      {
        SCITBX_ASSERT(grid.nd() == 3);
        SCITBX_ASSERT(grid.all().all_gt(0));
        af::versa<double, af::flex_grid<> > map(grid, 0);
        std::ifstream cin(file_name.c_str());
        std::string line;
        for (std::size_t i=0;i<n_header_lines;i++) {
          std::getline(cin, line);
        }
        af::ref<double, af::c_grid<3> > map_ref(
          map.begin(),
          af::c_grid<3>(af::adapt(map.accessor().all())));
        for(std::size_t iz=0;iz<map_ref.accessor()[2];iz++) {
          std::getline(cin, line); // reads section number
          std::size_t i_fld = 6;
          for(std::size_t iy=0;iy<map_ref.accessor()[1];iy++) {
            for(std::size_t ix=0;ix<map_ref.accessor()[0];ix++) {
              if (i_fld == 6) {
                std::getline(cin, line);
                i_fld = 0;
              }
              map_ref(ix,iy,iz) = std::atof(line.substr(i_fld*12,12).c_str());
              i_fld++;
            }
          }
        }
        cin.close();
        return map;
      }

      void
      WriteXplorMap(
        std::string const& file_name,
        cctbx::uctbx::unit_cell const& unit_cell,
        af::const_ref<double, af::flex_grid<> > const& map,
        double average,
        double standard_deviation)
      {
        SCITBX_ASSERT(map.accessor().nd() == 3);
        SCITBX_ASSERT(map.accessor().all().all_gt(0));
        SCITBX_ASSERT(!map.accessor().is_padded());
        FILE* fh = fopen(file_name.c_str(), "ab");
        SCITBX_ASSERT(fh != 0);
        for(std::size_t i=0;i<6;i++) {
          fprintf(fh, "%s",
            format_e<12>("%12.5E", unit_cell.parameters()[i]).s);
        }
        fprintf(fh, "\n");
        fprintf(fh, "ZYX\n");
        af::const_ref<double, af::c_grid<3> > map_ref(
          map.begin(),
          af::c_grid<3>(af::adapt(map.accessor().all())));
        for(std::size_t iz=0;iz<map_ref.accessor()[2];iz++) {
          fprintf(fh, "%8lu\n", static_cast<unsigned long>(iz));
          std::size_t i_fld = 0;
          for(std::size_t iy=0;iy<map_ref.accessor()[1];iy++) {
            for(std::size_t ix=0;ix<map_ref.accessor()[0];ix++) {
              fprintf(fh, "%s", format_e<12>("%12.5E", map_ref(ix,iy,iz)).s);
              i_fld++;
              if (i_fld == 6) {
                fprintf(fh, "\n");
                i_fld = 0;
              }
            }
          }
          if (i_fld > 0) {
            fprintf(fh, "\n");
          }
        }
        fprintf(fh, "   -9999\n");
        fprintf(fh, "%s%s\n",
          format_e<12>("%12.4E", average).s,
          format_e<12>("%12.4E", standard_deviation).s);
        fclose(fh);
      }
  };

} // namespace <anonymous>

BOOST_PYTHON_MODULE(iotbx_xplor_ext)
{
  using namespace boost::python;
  typedef boost::python::arg arg_; // gcc 2.96 workaround
  class_<XplorMap>("XplorMap", init<>())
    .def("ReadXplorMap", &XplorMap::ReadXplorMap,
      (arg_("file_name"), arg_("n_header_lines"), arg_("grid")))
    .def("WriteXplorMap", &XplorMap::WriteXplorMap,
      (arg_("file_name"), arg_("unit_cell"), arg_("map"),
       arg_("average"), arg_("standard_deviation")))
  ;
}
