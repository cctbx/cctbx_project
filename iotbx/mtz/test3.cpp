#include <iostream>
#include <iotbx/cppmtz.h>

void tester(){
  iotbx::mtz::Mtz p("phase_noanom.mtz");
  
  iotbx::mtz::Crystal cryst = p.columnToCrystal("FP");
  af::tiny<double,6> ux = cryst.UnitCell();
  std::cout<<ux[0]<<std::endl;
  std::cout<<ux[1]<<std::endl;
  std::cout<<ux[2]<<std::endl;
  std::cout<<ux[3]<<std::endl;
  std::cout<<ux[4]<<std::endl;
  std::cout<<ux[5]<<std::endl;
  std::cout<<p.SpaceGroup()<<std::endl;
}

int main(){
  int i;
  for ( i=0; i<1; i++){tester();std::cout<<i<<std::endl;}
}
