#include <string>
#include <cstdlib>
#include <fstream>
#include <iostream>

# include <boost/python/def.hpp>
# include <boost/python/class.hpp>
# include <boost/python/extract.hpp>
# include <boost/python/list.hpp>
# include <boost/python/module.hpp>

#include <scitbx/array_family/flex_types.h>
#include <scitbx/array_family/boost_python/flex_fwd.h>
#include <cctbx/uctbx.h>
namespace af = scitbx::af;

#include <cstdio>
#include <string>
#include <cassert>
#include <boost/smart_ptr.hpp>

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
public:
  ScientificFormatter(std::string s);
  boost::shared_array<char> operator() (double d);
};

ScientificFormatter::ScientificFormatter(std::string s){
    //parse the format %12.5E becomes TOKEN digits TOKEN decimal LETTER
    //could re-implement this with the boost tokenizer
    assert (s.find("%") == 0);

    int dot = s.find("."); std::string(s,dot)/* essentially an assertion */;
    int letter = s.find("E"); std::string(s,letter);

    digits   = std::atoi(std::string(s,1,dot-1).c_str());
    frac = std::atoi(std::string(s,dot+1,letter-dot).c_str());

    assert (frac+7 <= digits); //otherwise format too small to hold the frac
    sprintf (iformat, "%c%d.%dE\0",'%',frac+8,frac);
    sprintf (oformat, "%s\0",s.c_str());

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
    std::strncpy(b.get(), a.get(), digits-2);
    std::strncpy(b.get()+digits-2, a.get()+frac+6, 3);
    return b;
  }


class XplorMap {
private:
  ScientificFormatter pretty;
  ScientificFormatter pretty4;
public:
  XplorMap();
  af::flex_double ReadXplorMap(const std::string&,int,boost::python::list);
  void WriteXplorMap(cctbx::uctbx::unit_cell, af::flex_double, 
                     double,double,std::string);
};


XplorMap::XplorMap():pretty("%12.5E"),pretty4("%12.4E"){}

af::flex_double XplorMap::ReadXplorMap(const std::string& filename,
                             int headers,
                             boost::python::list sec
                             ) {
  std::ifstream cin(filename.c_str());
  std::string line;
  for (std::size_t i = 0; i < headers; i++) {
    getline(cin,line);
  }
  // dump everything to a 3-D array
  int nX,nY,nZ;
  nX = boost::python::extract<int>(sec[2]) - boost::python::extract<int>(sec[1]) + 1;
  nY = boost::python::extract<int>(sec[5]) - boost::python::extract<int>(sec[4]) + 1;
  nZ = boost::python::extract<int>(sec[8]) - boost::python::extract<int>(sec[7]) + 1;

  af::flex_double m(af::flex_grid<>(nZ,nY,nX)); //arrays have first dimension slow
  double* e = m.begin();
  for (int section = 0; section < nZ; section++) {
    getline(cin,line); //reads section number
    for (int rect = 0, fast = 0; rect < (nX*nY); rect++,fast++) {
      if (fast%6==0) 
        {    getline(cin,line); fast=0;      }
      *e = std::atof(line.substr(fast*12,12).c_str());
      ++e;      
    }
  }
    
  cin.close();
  return m;
}

void XplorMap::WriteXplorMap(cctbx::uctbx::unit_cell uc,
                             af::flex_double data,
                             double average, 
                             double stddev,
                             std::string outputfile) {
  FILE* fh = fopen(outputfile.c_str(),"a");
  //Unit Cell
  fprintf(fh,"%s%s%s%s%s%s\n", pretty(uc.parameters()[0]).get(), 
                               pretty(uc.parameters()[1]).get(), 
                               pretty(uc.parameters()[2]).get(), 
                               pretty(uc.parameters()[3]).get(), 
                               pretty(uc.parameters()[4]).get(), 
                               pretty(uc.parameters()[5]).get() );
  //Data
  fprintf(fh,"ZYX\n");
  double* e = data.begin();
  for (int z = 0; z< data.accessor().focus()[0]; z++){
    int eol=0;
    fprintf(fh,"%8d\n",z);
    for (int xy = 0; xy < data.accessor().focus()[1]*data.accessor().focus()[2]; ++xy){
        fprintf(fh,"%s",pretty(*e).get());
        ++e;
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

#include <cctbx/boost_python/flex_fwd.h>
#include <scitbx/boost_python/utils.h>
using namespace boost::python;

BOOST_PYTHON_MODULE(xplor)
{
    class_<XplorMap>("XplorMap", init<>())
      .def("ReadXplorMap", &XplorMap::ReadXplorMap)
      .def("WriteXplorMap", &XplorMap::WriteXplorMap)
    ;
}
