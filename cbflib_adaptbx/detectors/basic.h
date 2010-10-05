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
  inline Error(std::string s):s(s){}
  inline virtual const char* what() const throw() {return s.c_str();}
  inline virtual ~Error() throw() {}
};

  }//namespace detectors
}//namespace iotbx

#endif
