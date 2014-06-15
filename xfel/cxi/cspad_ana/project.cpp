#include <iostream>
#include <math.h>
typedef scitbx::af::versa<double, scitbx::af::flex_grid<> >  flex_double;

static boost::python::object
compute_projection(const flex_double& pixels,
                   const int limit,
                   const float sin,
                   const float cos)
{
//    using namespace std;
    int x=0; int y=0;
    float raw_data[limit][limit];
    shared_double oneD(limit,0);
    memset(raw_data, 0, limit*limit*sizeof(float));
    for (int i=0; i<1024; i++)
       for (int j=0; j<1024; j++)
         {
           x=int(round(j*cos+i*sin));
           y=int(round(-j*sin+i*cos))+int(round(sin*1023));
           raw_data[y][x]=pixels(i,j);
         }
    int count=0;
    for (int i=0; i<limit; i++)
      {
       count=0;
       for (int j=0; j<limit; j++)
         {
          if(raw_data[i][j]!=0)
           {
           oneD[i]=oneD[i]+raw_data[i][j];
           count++;
           }
         }
      if(count>0) oneD[i]=oneD[i]/count;
      }
         return  boost::python::object(oneD);

}
