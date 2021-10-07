// Copyright(c) 2021, Richardson Lab at Duke
// Licensed under the Apache 2 license
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//     http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissionsand
// limitations under the License.

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
