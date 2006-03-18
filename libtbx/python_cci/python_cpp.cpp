#include "Python.h"
#include <iostream>
#include <stdexcept>

int
main(int argc, char **argv)
{
  int result = 0;
  try {
    // ensure I/O system is initialized before any extensions are imported
    std::cout << std::flush;
    result = Py_Main(argc, argv);
  }
  catch (std::exception const& exc) {
    std::cerr << "python: C++ exception: " << exc.what()
              << std::endl << std::flush;
    return 1;
  }
  catch (...) {
    std::cerr << "python: unknown C++ exception"
              << std::endl << std::flush;
    return 2;
  }
  return result;
}
