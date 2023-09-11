##################################################################################
# This is a test program to validate that the Python wrapping of Probe worked.
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
import boost_adaptbx.boost.python as bp
bp.import_ext("mmtbx_probe_ext")
import mmtbx_probe_ext as probeext

#========================================================================
# Call the test functions for the libraries we test.

ret = probeext.DotSpheres_test()
assert len(ret) == 0, "DotSpheres_test() failed: " + ret

ret = probeext.SpatialQuery_test()
assert len(ret) == 0, "SpatialQuery_test() failed: " + ret

ret = probeext.Scoring_test()
assert len(ret) == 0, "Scoring_test() failed: " + ret

#========================================================================
# Now ensure that we can use the C++-wrapped classes as intended to make sure
# that the wrapping code or parameters have not changed.

#========================================================================
# Make sure we can get at the DotSphere objects and their methods
cache = probeext.DotSphereCache(10)
sphere1 = cache.get_sphere(1)
rad = sphere1.radius()
density = sphere1.density()

# This one requires knowing about the af::shared type, which requires
# including another module, so we don't test it here.
#dots = sphere1.dots()
