#ifndef CBFLIB_BASIC_AD_H
#define CBFLIB_BASIC_AD_H
#include <string>
#include <vector>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <exception>
#include <scitbx/array_family/flex_types.h>

namespace iotbx {
  namespace detectors {
class Error : public std::exception {
private:
  std::string s;
public:
  inline Error(std::string s):s(s){/*SCITBX_EXAMINE(s);*/}
  inline virtual const char* what() const throw() {return s.c_str();}
  inline virtual ~Error() throw() {}
};
class OpenFileError : public std::exception {
private:
  std::string s;
public:
  inline OpenFileError(std::string s, FILE * stream):s(s){
    /*SCITBX_EXAMINE(s);*/
    std::fclose(stream);}
  inline virtual const char* what() const throw() {return s.c_str();}
  inline virtual ~OpenFileError() throw() {}
};

  }//namespace detectors
}//namespace iotbx

#endif
