// Copyright(c) 2023, Richardson Lab at Duke
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

#include "PositionReturn.h"

typedef scitbx::vec3<double> vec3;

namespace molprobity {
  namespace reduce {


std::string PositionReturn_test()
{
  std::string ret;

  // Fill in all entry types in the structure and verify that they
  // can be read back.
  PositionReturn pr;

  iotbx::pdb::hierarchy::atom atom;
  atom.set_name(" CB ");
  pr.atoms.push_back(atom);
  if (atom.data->name != pr.atoms.back().data->name) {
    return "molprobity::reduce::PositionReturn_test() failed: atom test failed";
  }

  molprobity::probe::ExtraAtomInfo ei;
  ei.setIsAcceptor(true);
  scitbx::af::shared<molprobity::probe::ExtraAtomInfo> eis;
  eis.push_back(ei);
  pr.extraInfos.push_back(eis);
  if (ei != pr.extraInfos.front().front()) {
    return "molprobity::reduce::PositionReturn_test() failed: extraAtomInfo test failed";
  };

  bool dm = true;
  scitbx::af::shared<bool> dms;
  dms.push_back(dm);
  pr.deleteMes.push_back(dms);
  if (dm != pr.deleteMes.front().front()) {
    return "molprobity::reduce::PositionReturn_test() failed: deleteMe test failed";
  };

  double pe = 1.0;
  pr.preferenceEnergies.push_back(pe);
  if (pe != pr.preferenceEnergies.front()) {
    return "molprobity::reduce::PositionReturn_test() failed: preferenceEnergy test failed";
  };

  // All tests passed.
  return "";
}


  } // end namespace reduce
} // end namespace molprobity
