// temp comment out #include <cmtzlib.h>
#include <iostream>
// temp comment out #include "iotbx/cppmtz.h"
//  print "Unitcell:", tuple(p.UnitCell(0))

//  print "Columns:", tuple(p.columns())

void tester(){
/* temp comment out
  iotbx::mtz::Mtz p("test.mtz");
  p.UnitCell(0);
  int N = p.size();
  //std::cout<<p.SpaceGroup()<<std::endl;
  iotbx::mtz::Column H = p.getColumn("H");
  iotbx::mtz::Column K = p.getColumn("K");
  iotbx::mtz::Column L = p.getColumn("L");
  iotbx::mtz::Column I = p.getColumn("I");
  iotbx::mtz::Column SIGI = p.getColumn("SIGI");
  iotbx::mtz::Column XDET = p.getColumn("XDET");
  iotbx::mtz::Column YDET = p.getColumn("YDET");
  iotbx::mtz::Column BATCH = p.getColumn("BATCH");

  temp comment out */
  /*for (int i=0; i< p.size(); i++) {
    if (H.lookup(i)==K.lookup(i) && K.lookup(i)==L.lookup(i)) {
      std::cout<<H.lookup(i)<<" "<<K.lookup(i)<<" "<<L.lookup(i);
      std::cout<<" "<<I.lookup(i)<<" "<<SIGI.lookup(i)<<std::endl;
    }
  }
  */

}

int main(){
  int i;
  for ( i=0; i<10000; i++){tester();std::cout<<i<<std::endl;}
}
