/*
  (c) Copyright 2005, Jacob N. Smith, Erik Mckee, Texas Agricultural
  Experiment Station & Texas Engineering Experiment Station.
*/

#include<iostream>
#include<cstring>
#include<chiltbx/grammar/standard_parsers.h>
#include<chiltbx/chil.h>

struct foo {
        int operator () ( int ) const {
                return 5001;
        }
};

typedef int(foo::*bar)( int ) const;

int main ( int argc, char *argv[] ) {

        bar B = &foo::operator();

        char *statement = "include \t\t\t   foo\n";

        if ( argc>1 )
                statement = argv[1];

        char *beg = statement;
        char *end = beg + std::strlen(statement);

        using namespace chiltbx::grammar;

        chiltbx::chil::parser<> P(false);
        P.parse(beg,end);

}
