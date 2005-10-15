/*

  (c) Copyright 2005, Jacob N. Smith, Erik Mckee, Texas Agricultural
  Experiment Station & Texas Engineering Experiment Station.

*/
#include<cstring>
#include<iostream>
#include<chiltbx/grammar/standard_parsers.h>

int main ( int argc, char *argv[] ) {

        if ( argc<2 )
                return 0;

        char *beg = argv[1];
        char *end = beg + std::strlen(argv[1]);

        using namespace chiltbx::grammar;

        natural<16,char,long long> n;
        result<long long> R = n.parse<void>(beg,end);

        if ( R.truth )
                std::cout << R.value << std::endl;
        else
                std::cout << "Failed" << std::endl;

        return 0;
}
