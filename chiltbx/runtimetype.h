/*

  (c) Copyright 2005, Jacob N. Smith, Erik Mckee, Texas Agricultural
  Experiment Station & Texas Engineering Experiment Station.

*/
#ifndef CHILTBX_RUNTIMETYPE_H
#define CHILTBX_RUNTIMETYPE_H

#include<typeinfo>
#include<chiltbx/str.h>

namespace chiltbx {

typedef chiltbx::string ns_string;

// as simple as it seems
// just allows the user to override the type
// for use with serialization dictionaries

template < typename T > struct runtimetype {
  static ns_string name () {
    return ns_string(typeid(T).name());
  }
  static ns_string name ( T const& value ) {
    return ns_string(typeid(value).name());
  }
};

}// end chiltbx namespace

#endif//CHILTBX_RUNTIMETYPE_H
