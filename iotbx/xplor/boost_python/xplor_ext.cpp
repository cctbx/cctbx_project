#include <cctbx/boost_python/flex_fwd.h>
#include <cctbx/maptbx/accessors/c_grid_padded_p1.h>
#include <cctbx/uctbx.h>
#include <scitbx/array_family/flex_types.h>
#include <boost/python/module.hpp>
#include <boost/python/class.hpp>
#include <boost/python/args.hpp>
#include <boost/python/def.hpp>
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

  class map_reader
  {
    public:
      map_reader() {}

      map_reader(
        std::string const& file_name,
        std::size_t n_header_lines,
        af::flex_grid<> const& grid)
      :
        data(grid, 0)
      {
        SCITBX_ASSERT(grid.nd() == 3);
        SCITBX_ASSERT(grid.all().all_gt(0));
        std::ifstream cin(file_name.c_str());
        std::string line;
        for (std::size_t i=0;i<n_header_lines;i++) {
          std::getline(cin, line);
        }
        af::ref<double, af::c_grid<3> > data_ref(
          data.begin(),
          af::c_grid<3>(af::adapt(data.accessor().all())));
        for(std::size_t iz=0;iz<data_ref.accessor()[2];iz++) {
          std::getline(cin, line); // reads section number
          std::size_t i_fld = 6;
          for(std::size_t iy=0;iy<data_ref.accessor()[1];iy++) {
            for(std::size_t ix=0;ix<data_ref.accessor()[0];ix++) {
              if (i_fld == 6) {
                std::getline(cin, line);
                i_fld = 0;
              }
              data_ref(ix,iy,iz) = std::atof(line.substr(i_fld*12,12).c_str());
              i_fld++;
            }
          }
        }
        std::getline(cin, line);
        if (line.size() == 0) {
          average = -1;
          standard_deviation = -1;
        }
        else {
          int expected_9999 = std::atoi(line.substr(0,8).c_str());
          SCITBX_ASSERT(expected_9999 == -9999);
          std::getline(cin, line);
          average = std::atof(line.substr(0,12).c_str());
          standard_deviation = std::atof(line.substr(12,12).c_str());
        }
        cin.close();
      }

      af::versa<double, af::flex_grid<> > data;
      double average;
      double standard_deviation;
  };

  FILE*
  write_head(
    std::string const& file_name,
    cctbx::uctbx::unit_cell const& unit_cell)
  {
    FILE* fh = fopen(file_name.c_str(), "ab");
    SCITBX_ASSERT(fh != 0);
    for(std::size_t i=0;i<6;i++) {
      fprintf(fh, "%s",
        format_e<12>("%12.5E", unit_cell.parameters()[i]).s);
    }
    fprintf(fh, "\n");
    fprintf(fh, "ZYX\n");
    return fh;
  }

  void
  write_tail(
    FILE* fh,
    double average,
    double standard_deviation)
  {
    fprintf(fh, "   -9999\n");
    fprintf(fh, "%s%s\n",
      format_e<12>("%12.4E", average).s,
      format_e<12>("%12.4E", standard_deviation).s);
    fclose(fh);
  }

  void
  map_writer_box(
    std::string const& file_name,
    cctbx::uctbx::unit_cell const& unit_cell,
    af::const_ref<double, af::flex_grid<> > const& data,
    double average,
    double standard_deviation)
  {
    SCITBX_ASSERT(data.accessor().nd() == 3);
    SCITBX_ASSERT(data.accessor().all().all_gt(0));
    SCITBX_ASSERT(!data.accessor().is_padded());
    FILE* fh = write_head(file_name, unit_cell);
    af::const_ref<double, af::c_grid<3> > data_ref(
      data.begin(),
      af::c_grid<3>(af::adapt(data.accessor().all())));
    for(std::size_t iz=0;iz<data_ref.accessor()[2];iz++) {
      fprintf(fh, "%8lu\n", static_cast<unsigned long>(iz));
      int i_fld = 0;
      for(std::size_t iy=0;iy<data_ref.accessor()[1];iy++) {
        for(std::size_t ix=0;ix<data_ref.accessor()[0];ix++) {
          fprintf(fh, "%s", format_e<12>("%12.5E", data_ref(ix,iy,iz)).s);
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
    write_tail(fh, average, standard_deviation);
  }

  void
  map_writer_p1_cell(
    std::string const& file_name,
    cctbx::uctbx::unit_cell const& unit_cell,
    af::int3 const& gridding_first,
    af::int3 const& gridding_last,
    af::const_ref<double, cctbx::maptbx::c_grid_padded_p1<3> > const& data,
    double average,
    double standard_deviation)
  {
    FILE* fh = write_head(file_name, unit_cell);
    unsigned i_section = 0;
    for(int iz=gridding_first[2];iz<=gridding_last[2];iz++,i_section++) {
      fprintf(fh, "%8u\n", i_section);
      int i_fld = 0;
      for(int iy=gridding_first[1];iy<=gridding_last[1];iy++) {
        for(int ix=gridding_first[0];ix<=gridding_last[0];ix++) {
          fprintf(fh, "%s", format_e<12>("%12.5E", data(ix,iy,iz)).s);
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
    write_tail(fh, average, standard_deviation);
  }

} // namespace <anonymous>

BOOST_PYTHON_MODULE(iotbx_xplor_ext)
{
  using namespace boost::python;
  typedef boost::python::arg arg_; // gcc 2.96 workaround

  class_<map_reader>("map_reader", no_init)
    .def(init<std::string const&, std::size_t, af::flex_grid<> const&>(
      (arg_("file_name"), arg_("n_header_lines"), arg_("grid"))))
    .def_readonly("data", &map_reader::data)
    .def_readonly("average", &map_reader::average)
    .def_readonly("standard_deviation", &map_reader::standard_deviation)
  ;

  def("map_writer", map_writer_box,
    (arg_("file_name"), arg_("unit_cell"),
     arg_("data"),
     arg_("average"), arg_("standard_deviation")));
  def("map_writer", map_writer_p1_cell,
    (arg_("file_name"), arg_("unit_cell"),
     arg_("gridding_first"), arg_("gridding_last"), arg_("data"),
     arg_("average"), arg_("standard_deviation")));
}
