##################################################################################
# This is a test program to validate that the Python wrapping of Reduce worked.
#

#                Copyright 2021  Richardson Lab at Duke University
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

from __future__ import print_function, nested_scopes, generators, division
from __future__ import absolute_import
from mmtbx.reduce import Movers
from mmtbx.reduce import InteractionGraph
from mmtbx.reduce import Optimizers

def RunReduceTests():

  #========================================================================
  # Call the test functions for all of our files.

  print('Testing Movers objects')
  ret = Movers.Test()
  if len(ret) != 0:
    print("Failure:",ret)
  assert (len(ret) == 0)

  print('Testing InteractionGraph objects')
  ret = InteractionGraph.Test()
  if len(ret) != 0:
    print("Failure:",ret)
  assert (len(ret) == 0)

  print('Testing Optimizers')
  ret = Optimizers.Test()
  if len(ret) != 0:
    print("Failure:",ret)
  assert (len(ret) == 0)

  return ret

if __name__ == '__main__':

  ret = RunReduceTests()
  if len(ret) == 0:
    print('Success!')
  else:
    print(ret)

  assert (len(ret) == 0)
