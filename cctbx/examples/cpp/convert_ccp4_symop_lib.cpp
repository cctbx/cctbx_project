/*
Example by Kevin Cowtan, 2001. Public Domain.
This little program reads the CCP4 symmetry library file, and interprets
the symops. It then returns the Hall code for the spacegroup. Used to
create a new library file to allow compatible spacegroup naming between
CCP4 and cctbx.

usage: convert_ccp4_symop_lib < symop.lib
*/
#include <iostream>
#include <cctbx/sgtbx/symbols.h>
#include <cctbx/sgtbx/groups.h>

int main()
{
  using namespace std;
  string line;
  int nspgrp,nsym;
  while ( !cin.eof() ) {
    cin >> nspgrp >> nsym;          // read spacegroup number and nsym
    if ( cin.eof() ) break;
    getline( cin, line );           // and the rest of the line
    cout << nspgrp << " " << nsym << line << "\n";         // print it all
    sgtbx::SpaceGroup SgOps;             // now interpret the symops
    for (int i = 0; i < nsym; i++) {
      getline( cin, line );                     // get the i'th symop
      // cout << line << "\n";
      SgOps.expandSMx( sgtbx::RTMx(line) );     // and interpret
    }
    cout << SgOps.BuildHallSymbol() << "\n";    // now produce the sg symbol
    cout << SgOps.BuildLookupSymbol() << "\n";
  }
  return 0;
}
