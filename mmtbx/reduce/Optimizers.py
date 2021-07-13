##################################################################################
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

import Movers
import InteractionGraph

##################################################################################
# This is a set of functions that implement placement and optimization of
# Reduce's "Movers".

##################################################################################
# Test function to verify that all functions behave properly.

def Test():
  """Test function for all functions provided above.
  :returns Empty string on success, string describing the problem on failure.
  :returns Empty string on success, string describing the problem on failure.
  """

  # @todo
  return ""

##################################################################################
# If we're run on the command line, test our classes and functions.
if __name__ == '__main__':

  ret = Test()
  if len(ret) == 0:
    print('Success!')
  else:
    print(ret)

  assert (len(ret) == 0)
