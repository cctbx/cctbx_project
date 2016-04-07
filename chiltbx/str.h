/*

  (c) Copyright 2005, Jacob N. Smith, Erik Mckee, Texas Agricultural
  Experiment Station & Texas Engineering Experiment Station.

*/
#ifndef CHILTBX_STRING_HPP
#define CHILTBX_STRING_HPP

#include<string>

namespace chiltbx {

#ifndef CHILTBX_DEFAULT_CHARACTER_TYPE
# define CHILTBX_DEFAULT_CHARACTER_TYPE char
#endif //CHILTBX_DEFAULT_CHARACTER_TYPE

typedef std::basic_string<CHILTBX_DEFAULT_CHARACTER_TYPE> string;

}// end chiltbx namespace

#endif//CHILTBX_STRING_HPP
