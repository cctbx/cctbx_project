/*
Example by Kevin Cowtan, 2001. Public Domain.
This little program reads the CCP4 symmetry library file, and interprets
the symops. It then returns the Hall code for the spacegroup. Used to
create a new library file to allow compatible spacegroup naming between
CCP4 and cctbx.

usage: convert_ccp4_symop_lib < symop.lib
*/
#include <cctbx/sgtbx/space_group_type.h>
#include <iostream>

int main()
{
  using cctbx::sgtbx::space_group_type;
  using cctbx::sgtbx::space_group;
  using cctbx::sgtbx::rt_mx;
  using namespace std;
  string line;
  int nspgrp,nsym;
  while ( !cin.eof() ) {
    cin >> nspgrp >> nsym;          // read spacegroup number and nsym
    if ( cin.eof() ) break;
    getline( cin, line );           // and the rest of the line
    cout << nspgrp << " " << nsym << line << "\n";         // print it all
    space_group sg; // now interpret the symops
    for (int i = 0; i < nsym; i++) {
      getline( cin, line );                     // get the i'th symop
      // cout << line << "\n";
      sg.expand_smx( rt_mx(line) );     // and interpret
    }
    space_group_type sg_type(sg);
    cout << sg_type.hall_symbol() << "\n";    // now produce the sg symbol
    cout << sg_type.lookup_symbol() << "\n";
  }
  return 0;
}
