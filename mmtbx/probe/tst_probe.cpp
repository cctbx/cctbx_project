// Copyright(c) 2021, Richardson Lab at Duke
// Licensed under the Apache 2 license

#include "Scoring.h"
#include "DotSpheres.h"
#include "SpatialQuery.h"
#include <iostream>

int main(int argc, const char* argv[])
{
  std::string ret;

  // Test DotSpheres module.
  ret = molprobity::probe::DotSpheres_test();
  if (!ret.empty()) {
    std::cerr << "Error: " << ret << std::endl;
    return 1;
  }

  // Test SpatialQuery module.
  ret = molprobity::probe::SpatialQuery_test();
  if (!ret.empty()) {
    std::cerr << "Error: " << ret << std::endl;
    return 1;
  }

  // Test Scoring module.
  ret = molprobity::probe::Scoring_test();
  if (!ret.empty()) {
    std::cerr << "Error: " << ret << std::endl;
    return 1;
  }

  std::cout << "Success!" << std::endl;
  return 0;
}
