from __future__ import absolute_import, division, print_function
from crys3d.regression import tests_PhenixHKLviewer as tsthkl
import sys, os

# Using a webbrowser exercise_OSbrowser() runs HKLviewer to enact the settings in philstr2 and
# eventually asserts that the visible reflections in the browser match the indices in
# reflections2match2. Due to occasional instability of websockets on virtual machines the test
# is run in a loop until it passes but no longer than maxruns times.

# The test loads 1upp_lowres.mtz, then runs xtriage on the data, then
# tests for the visible reflections of 1upp_lowres.mtz when the sphere of reflections
# is sliced perpendicular to the twin axis detected by xtriage at layer 13 and reflections have been
# divided into 4 bins according to I_lowres values and only reflections of the highest bin with
# values above 20000 are displayed. Twinning suggests the pattern of the slice should be
# close to 2-fold symmetry

def run():
  if sys.version_info < (3, 0):
    return
  if 'linux' in sys.platform and os.environ.get("DISPLAY", None) is None:
    return # no DISPLAY on Linux on azure virtual machines
  count = 0
  while True:
    print("running %d" %count)
    # websockets employed by HKLviewer is slightly unstable on virtual machines used in CI on Azure.
    # This might yield a bogus failure of the test. If so, repeat the test at most maxruns times
    # or until it passes whichever comes first.
    if not tsthkl.runagain(tsthkl.exercise_OSbrowser,
                                    tsthkl.philstr2,
                                    tsthkl.reflections2match2,
                                    "OSbrowserXtriageTwinning"):
      break
    count +=1
    assert(count < tsthkl.maxruns)


if __name__ == '__main__':
  run()
