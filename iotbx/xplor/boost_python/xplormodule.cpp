#include <cctbx/boost_python/flex_fwd.h>
#include <cctbx/uctbx.h>
#include <scitbx/array_family/flex_types.h>
#include <boost/python/class.hpp>
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


class XplorMap {
private:
  ScientificFormatter pretty;
  ScientificFormatter pretty4;
public:
  XplorMap();
  af::versa<double,af::c_grid<3> >
  ReadXplorMap(std::string const&, int, af::int9 const& sec);
  void WriteXplorMap(cctbx::uctbx::unit_cell const&,
                     af::const_ref<double,af::c_grid<3> > const&,
                     double, double, std::string const&);
};

XplorMap::XplorMap():pretty("%12.5E"),pretty4("%12.4E"){}

af::versa<double,af::c_grid<3> >
XplorMap::ReadXplorMap(std::string const& filename,
                       int headers,
                       af::int9 const& sec)
{
  std::ifstream cin(filename.c_str());
  std::string line;
  for (std::size_t i = 0; i < headers; i++) {
    std::getline(cin,line);
  }
  // dump everything to a 3-D array
  int nX = sec[2] - sec[1] + 1;
  int nY = sec[5] - sec[4] + 1;
  int nZ = sec[8] - sec[7] + 1;
  SCITBX_ASSERT(nX > 0);
  SCITBX_ASSERT(nY > 0);
  SCITBX_ASSERT(nZ > 0);
  af::versa<double,af::c_grid<3> > m(af::c_grid<3>(nX,nY,nZ));
  af::ref<double, af::c_grid<3> > mref=m.ref();
  for (int section = 0; section < nZ; section++) {
    std::getline(cin,line); //reads section number
    int x=0;
    int y=0;
    for (int rect = 0, fast = 0; rect < (nX*nY); rect++,fast++) {
      if (fast%6==0)
        {    std::getline(cin,line); fast=0;      }
      mref(x,y,section) = std::atof(line.substr(fast*12,12).c_str());
      ++x;
      if (x==nX) {x=0; ++y;}
    }
  }

  cin.close();
  return m;
}

void XplorMap::WriteXplorMap(cctbx::uctbx::unit_cell const& uc,
                             af::const_ref<double,af::c_grid<3> > const& data,
                             double average,
                             double stddev,
                             std::string const& outputfile)
{
  FILE* fh = fopen(outputfile.c_str(),"ab");
  SCITBX_ASSERT(fh != 0);
  //Unit Cell
  for(std::size_t i=0;i<6;i++) {
    fprintf(fh, "%s", pretty(uc.parameters()[i]).get());
  }
  fprintf(fh, "\n");
  //Data
  fprintf(fh,"ZYX\n");
  int xsize = data.accessor()[0];
  int ysize = data.accessor()[1];
  for (int z = 0; z< data.accessor()[2]; z++){
    int eol=0;
    int x=0;
    int y=0;
    fprintf(fh,"%8d\n",z);
    for (int xy = 0; xy < xsize*ysize; ++xy){
        fprintf(fh,"%s",pretty(data(x,y,z)).get());
        ++x;
        if (x==xsize) {x = 0;++y;}
        ++eol;
        if (eol % 6 == 0) {fprintf(fh,"\n");}
    }
    if (eol % 6 != 0) {fprintf(fh,"\n");}
  }
  //summary
  fprintf(fh,"   -9999\n");
  fprintf(fh, "%s%s\n", pretty4(0.0).get(), pretty4(1.0).get());
  fclose(fh);
}

} // namespace <anonymous>

BOOST_PYTHON_MODULE(xplor)
{
    using namespace boost::python;
    class_<XplorMap>("XplorMap", init<>())
      .def("ReadXplorMap", &XplorMap::ReadXplorMap)
      .def("WriteXplorMap", &XplorMap::WriteXplorMap)
    ;
}
