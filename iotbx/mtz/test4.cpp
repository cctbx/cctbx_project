#include <cmtzlib.h>
#include <iostream>
#include "iotbx/cppmtz.h"

void tester(){
  iotbx::mtz::Foo p;
  std::cout<<p.value()[0]<<std::endl;
}
int main(){
  for ( int i=0; i<100000; i++){tester();std::cout<<i<<std::endl;}
}
