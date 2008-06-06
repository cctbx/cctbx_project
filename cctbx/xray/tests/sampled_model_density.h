#ifndef CCTBX_XRAY_TESTS_SAMPLED_MODEL_DENSITY_H
#define CCTBX_XRAY_TESTS_SAMPLED_MODEL_DENSITY_H

#include <cctbx/xray/sampled_model_density.h>
#include <scitbx/array_family/simple_io.h>
#include <iostream>

using namespace cctbx;

struct sampled_model_density_test_case
{
  xray::sampled_model_density<double> *sampled_model_density;

  sampled_model_density_test_case();

  ~sampled_model_density_test_case() { delete sampled_model_density; }
};

#endif // GUARD
