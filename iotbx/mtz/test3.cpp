#include <iostream>
#include <iotbx/cppmtz.h>

void tester(){
  iotbx::mtz::Mtz p("phase_noanom.mtz");
  
  iotbx::mtz::Crystal cryst = p.columnToCrystal("FP");
  cctbx::uctbx::unit_cell ux = cryst.UnitCell();
  std::cout<<ux.parameters()[0]<<std::endl;
  std::cout<<ux.parameters()[1]<<std::endl;
  std::cout<<ux.parameters()[2]<<std::endl;
  std::cout<<ux.parameters()[3]<<std::endl;
  std::cout<<ux.parameters()[4]<<std::endl;
  std::cout<<ux.parameters()[5]<<std::endl;
  std::cout<<p.SpaceGroup()<<std::endl;
}

int main(){
  int i;
  for ( i=0; i<1; i++){tester();std::cout<<i<<std::endl;}
}
