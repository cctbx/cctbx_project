#include <string>
#include <cstdlib>
#include <fstream>
#include <iostream>

# include <boost/python/def.hpp>
# include <boost/python/extract.hpp>
# include <boost/python/list.hpp>
# include <boost/python/module.hpp>

#include <scitbx/array_family/flex_types.h>
#include <scitbx/array_family/boost_python/flex_fwd.h>
namespace af = scitbx::af;

af::flex_double ReadXplorMap(const std::string& filename,
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

void WriteXplorMap(std::string outputfile) {
/*  FILE* fh = fopen(outputfile.c_str(),"w");

  MillerIndex realsize = realGrid();

  ScientificFormatter pretty("%12.5E");

  //Unit Cell
  const cctbx::uctbx::uc_params& p = unitcell().getParameters();

  fprintf(fh,"%s%s%s%s%s%s\n", pretty(p.Len(0)).get(), 
                               pretty(p.Len(1)).get(), 
                               pretty(p.Len(2)).get(), 
                               pretty(p.Ang(0)).get(), 
                               pretty(p.Ang(1)).get(), 
                               pretty(p.Ang(2)).get() );

  fprintf(fh,"ZYX\n");

  int xp, yp, zp;

  //sections: wrap the output around so that O displays a full unit cell
  for (int z = 0; z<= realsize[2]; z++){
    int eol=0;
    fprintf(fh,"%8d\n",z);
    for (int y = 0; y<= realsize[1]; y++){
      for (int x = 0; x<= realsize[0]; x++){
        xp = (x != realsize[0])? x :0;
        yp = (y != realsize[1])? y :0;
        zp = (z != realsize[2])? z :0;
        //fprintf(fh,"%12.5E", (realData(xp,yp,zp)-average_)/stddev_);
        //fprintf(fh,"%12.5E", (realData(xp,yp,zp)-0.0)/1.0);
        try {
        fprintf(fh,"%s",pretty((realData(xp,yp,zp)-average_)/stddev_).get());
        } catch (...) {
        cout << "Caught assertion while printing "<<xp<<" "<<yp<<" "<<zp<<endl;
        printf("%12.5E", (realData(xp,yp,zp)-average_)/stddev_);
        throw;}
        eol++;
        if (eol % 6 == 0) {fprintf(fh,"\n");}
      }
    }
    if (eol % 6 != 0) {fprintf(fh,"\n");}
  }

  //summary
  fprintf(fh,"   -9999\n");
  //  fprintf(fh,"%12.4E %12.4E \n",average(),stddev());
  ScientificFormatter pretty4("%12.4E");  
  fprintf(fh, "%s%s\n", pretty4(0.0).get(), pretty4(1.0).get());
  fclose(fh);
  */
}


#include <cctbx/boost_python/flex_fwd.h>
#include <scitbx/boost_python/utils.h>
using namespace boost::python;

BOOST_PYTHON_MODULE(xplor)
{
   def("ReadXplorMap", ReadXplorMap);
}
