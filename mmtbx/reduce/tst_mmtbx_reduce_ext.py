##################################################################################
# This is a test program to validate that the Python wrapping of Reduce worked.
#
#                Copyright 2023  Richardson Lab at Duke University
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

# These other imports are needed within the Scons build environment.
from scitbx.array_family import flex    # import dependency
from iotbx import pdb                   # import dependency
import boost_adaptbx.boost.python as bp
bp.import_ext("mmtbx_probe_ext")
bp.import_ext("mmtbx_reduce_ext")
import mmtbx_reduce_ext as reduceExt

#========================================================================
# Call the test functions for the libraries we test.

ret = reduceExt.PositionReturn_test()
assert len(ret) == 0, "PositionReturn_test() failed: " + ret

ret = reduceExt.InteractionGraph_test()
assert len(ret) == 0, "reduceExt.InteractionGraph_test() failed: " + ret

ret = reduceExt.Optimizers_test()
assert len(ret) == 0, "reduceExt.Optimizers_test() failed: " + ret


#========================================================================
# Now ensure that we can use the C++-wrapped classes as intended to make sure
# that the wrapping code or parameters have not changed.

#========================================================================
# Make sure we can get at the objects and their methods
pr = reduceExt.PositionReturn
# We cannot use iotbx or reduce types in this test code, so we leave them blank.
# This will build a broken description, but it does test that we can convert lists
# into the appropriate types.
pr2 = reduceExt.PositionReturn(
        # Single array of atoms
        [ ],
        # List of lists of positions.
        [ [  ] ],
        # List of lists of ExtraAtomInfos
        [ [  ] ],
        # List of lists of DeleteMes
        [ [ False ] ],
        # Single list of preference energies
        [ 0.0 ])
assert len(pr2.atoms) == 0, "Unexpected atoms in list for argument-constructed PositionReturn"
assert len(pr2.positions[0]) == 0, "Unexpected positions in list for argument-constructed PositionReturn"
assert len(pr2.extraInfos[0]) == 0, "Unexpected extra-atom info in list for argument-constructed PositionReturn"
assert pr2.deleteMes[0][0] == False, "Unexpected deleteme value in list for argument-constructed PositionReturn"
