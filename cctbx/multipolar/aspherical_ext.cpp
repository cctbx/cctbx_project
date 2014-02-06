#include "aspherical.h"

#include <boost/python.hpp>
#include <sstream>

namespace cctbx {
namespace multipolar {

namespace{
std::string int2string(int i)
{
  std::string result;
  std::stringstream ss;
  ss<< i;
  ss>> result;
  return result;
}
}

   af::shared<std::string> assign_atom_types(
            af::shared<int> atomic_numbers,
            af::shared<scitbx::vec3<double> > coordinates)
           {
              af::shared<std::string> types;
              size_t i,n;
              n = atomic_numbers.size();

              for(i=0;i<n;i++)
                 types.push_back(int2string(atomic_numbers[i]));

              return types;
           }

   void multipolar_test(){}

   void init_module(){

       //using namespace boost::python;
       boost::python::def("multipolar_test",multipolar_test);
       boost::python::def("assign_atom_types",assign_atom_types);

   }


}
} //namespace cctbx

BOOST_PYTHON_MODULE(cctbx_multipolar_ext)
{
  cctbx::multipolar::init_module();
}
