#include <cctbx/basic/dynamic.h>
#include <iostream>
using std::cout;
using std::endl;
using cctbx::array_family::detail::handle;

struct twodouble{
  double a[2];
  double& operator[](int i) {return a[i];}
};

//may need to go up to 12 for alloc failure
int alloc_fail = 4;

void f() {  // test of raw handle
  try {
    int m;
    for (int i = 0,m=10; i < alloc_fail; i++,m*=10) { 
      std::cout<<"class handle allocate "<<m<<std::endl;

      handle<char> a(m);
      std::string s = "hello";
      for (int z=0; z<5; ++z) {a[z]=s[z];}
      for (int z=0; z<5; ++z) {std::cout<<a[z];} std::cout<<a.size()<<std::endl;
 
      handle<char> b(m);
      s = "bye..";
      for (int z=0; z<5; ++z) {b[z]=s[z];}
      for (int z=0; z<5; ++z) {std::cout<<b[z];} std::cout<<b.size()<<std::endl;
 
      handle<char> c(a);
      for (int z=0; z<5; ++z) {std::cout<<c[z];} std::cout<<c.size()<<std::endl;

      handle<char> d = a;
      for (int z=0; z<5; ++z) {std::cout<<d[z];} std::cout<<d.size()<<std::endl;
    
      char* pc = new char[5];
      for (int z=0; z<5; ++z) {pc[z]=s[z];}
      handle<char>::Courier courier(pc,5);
      handle<char> e(courier);
      for (int z=0; z<5; ++z) {std::cout<<e[z];} std::cout<<e.size()<<std::endl;

    }
  } catch (const std::exception& e) {
    std::cout<<e.what()<<std::endl;
  }
}

using cctbx::array_family::dynamic;
void g() {  // test of dynamic
  try {
    int m;
    for (int i = 0,m=10; i < alloc_fail; i++,m*=10) {
      std::cout<<"class dynamic allocate "<<m<<std::endl;

      dynamic<double> a(m);
      for (int z=0; z<m; ++z) {a[z]=z;}
      for (int z=0; z<5; ++z) {std::cout<<a[z];} 
      std::cout<<"a size "<<a.size()<<" "<<a.handle().use_count()<<std::endl;
      std::cout<<"size is"<<sizeof(a.element_size())<<std::endl;

      dynamic<double> b(m);
      for (int z=0; z<m; ++z) {b[z]=z+5;}
      for (int z=0; z<5; ++z) {std::cout<<b[z];}
      std::cout<<"b size "<<b.size()<<" "<<b.handle().use_count()<<std::endl;

      dynamic<double> c(a);
      for (int z=0; z<5; ++z) {std::cout<<c[z];}
      std::cout<<"c size "<<c.size()<<" "<<c.handle().use_count()<<std::endl;

      dynamic<double> d = a;
      for (int z=0; z<5; ++z) {std::cout<<d[z];}
      std::cout<<"d size "<<d.size()<<" "<<d.handle().use_count()<<std::endl;
      
      double* pd = new double[5];
      for (int z=0; z<5; ++z) {pd[z]=z;}
      handle<double>::Courier courier(pd,5);
      dynamic<double> e(courier);
      for (int z=0; z<5; ++z) {std::cout<<e[z];}
      std::cout<<" "<<e.size()<<" "<<e.handle().use_count()<<std::endl;

      a.resize(m*10);
      for (int z=0; z<5; ++z) {std::cout<<a[z];}
      std::cout<<" "<<a.size()<<" "<<a.handle().use_count()<<std::endl;
      for (int z=0; z<5; ++z) {std::cout<<d[z];}
      std::cout<<" "<<d.size()<<" "<<d.handle().use_count()<<std::endl;

      dynamic<float> f(b);
      for (int z=0; z<5; ++z) {std::cout<<f[z];}
      std::cout<<"f size "<<f.size()<<" "<<f.handle().use_count()<<std::endl;

      dynamic<float> g(b.size());
      g = b;
      for (int z=0; z<5; ++z) {std::cout<<g[z];}
      std::cout<<"g size "<<g.size()<<" "<<f.handle().use_count()<<std::endl;

      dynamic<twodouble> h(a.handle());
      for (int z=0; z<5; ++z) {std::cout<<h[z][0];}
      std::cout<<" "<<h.size()<<" "<<h.handle().use_count()<<std::endl;
     
    }
  } catch (const std::exception& e) {
    std::cout<<e.what()<<std::endl;
  }
}


int main() {
  for (int i = 0; i < 100; ++i) {
    cout<<"-------Iteration "<<i<<endl;
    f();g();}
}
