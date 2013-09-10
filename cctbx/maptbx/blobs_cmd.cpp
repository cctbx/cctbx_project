#include "cctbx/maptbx/skeletons.h"
#include "iotbx/xplor/map_reader.h"

#include <iostream>

int main(int argc, char *argv[])
{
  try {
    if( argc<2 )
    {
      std::cout << "Provide map file on command line" << std::endl;
      return 0;
    }
    std::string file_name(argv[1]);
    iotbx::xplor::map_reader rdr(file_name);
    cctbx::sgtbx::space_group p1;
    cctbx::maptbx::asymmetric_map amap(p1.type(), rdr.data.const_ref());
    return 0;
  }
  catch(std::exception &err)
  {
    std::cerr << "ERROR! " << err.what() << std::endl;
  }
  return -1;
}
