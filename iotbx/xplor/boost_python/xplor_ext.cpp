#include <cctbx/boost_python/flex_fwd.h>
#include <cctbx/uctbx.h>
#include <scitbx/array_family/flex_types.h>
#include <boost/python/class.hpp>
#include <boost/python/args.hpp>
#include <boost/python/module.hpp>
#include <boost/smart_ptr.hpp>
#include <fstream>
#include <string>
#include <cstdlib>
#include <cstdio>

namespace af = scitbx::af;

namespace {

class ScientificFormatter {
  //Necessary because Microsoft Visual C Runtime formats printf incorrectly
  //%12.5E becomes -d.dddddE+000
  //Usage:  ScientificFormatter pretty("%12.5E");
  //Some drawbacks to this approach:
  //1.  Need to use the get method:  fprintf(FH,"%s",pretty(2.0).get());
  //2.  Takes 0.3 sec per 100000 calls on tru64 (additional time over printf)

  //XXX No attempt has been made to optimize this! Certainly it can be improved.
private:
  char oformat[10];
  char iformat[10];
  int digits, frac;
  bool stdprintf;
  const char pad;
public:
  ScientificFormatter(std::string const& s);
  boost::shared_array<char> operator() (double d);
};

ScientificFormatter::ScientificFormatter(std::string const& s): pad(' ') {
    //parse the format %12.5E becomes TOKEN digits TOKEN decimal LETTER
    //could re-implement this with the boost tokenizer
    SCITBX_ASSERT(s.find("%") == 0);

    int dot = s.find("."); std::string(s,dot)/* essentially an assertion */;
    int letter = s.find("E"); std::string(s,letter);

    digits   = std::atoi(std::string(s,1,dot-1).c_str());
    frac = std::atoi(std::string(s,dot+1,letter-dot).c_str());

    SCITBX_ASSERT(frac+7 <= digits); // otherwise format too small to
                                     // hold the frac
    sprintf (iformat, "%c%d.%dE\0",'%',frac+8,frac);
    sprintf (oformat, "%s\0",s.c_str());
    SCITBX_ASSERT(std::string(iformat).size() < 8);
    SCITBX_ASSERT(std::string(oformat).size() < 8);

    //Now determine if we use standard or MSVC printf format
    boost::shared_array<char> t(new char[frac+9]); //ie., "-d.dddddE+00\0"
    sprintf (t.get(),iformat,-1.0);
    if (t[frac+5]=='0') {stdprintf=false;} else {stdprintf=true;}
  }

boost::shared_array<char> ScientificFormatter::operator() (double d) {
    boost::shared_array<char> b(new char[frac+10]);
    if (stdprintf) {
      sprintf(b.get(),oformat,d);
      return b;
    }
    boost::shared_array<char> a(new char[frac+10]);
    sprintf(a.get(),iformat,d);
    assert (*(a.get()+frac+5)=='0'); //Prevent floating overflow or underflow
    for (int x=0; x<digits-frac-7; ++x){ //forward padding
      std::strncpy(b.get()+x, &pad, 1);}
    std::strncpy(b.get()+digits-frac-7, a.get(), frac+5);
    std::strncpy(b.get()+digits-2, a.get()+frac+6, 3);
    return b;
  }

  class XplorMap
  {
    private:
      ScientificFormatter pretty;
      ScientificFormatter pretty4;

    public:
      XplorMap();

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
          fprintf(fh, "%s", pretty(unit_cell.parameters()[i]).get());
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
              fprintf(fh, "%s", pretty(map_ref(ix,iy,iz)).get());
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
          pretty4(average).get(), pretty4(standard_deviation).get());
        fclose(fh);
      }
  };

  XplorMap::XplorMap():pretty("%12.5E"),pretty4("%12.4E"){}

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
