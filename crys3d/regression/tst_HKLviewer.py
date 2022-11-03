from __future__ import absolute_import, division, print_function

import libtbx.load_env, os.path, time
from libtbx import easy_run
from crys3d.hklviewer import hklview_frame


philstr = """
clip_plane {
  normal_vector = "K-axis (0,1,0)"
  is_assoc_real_space_vector = True
  clip_width = 2
  hkldist = -9
}
viewer {
  data_array {
    label = "FP,SIGFP"
    datatype = 'Amplitude'
  }
  show_vector = "['K-axis (0,1,0)', True]"
  fixorientation = *vector None
}
hkls {
  expand_to_p1 = True
  expand_anomalous = True
}
"""

file_name = libtbx.env.find_in_repositories(
  relative_path="iotbx/regression/data/phaser_1.mtz",
  test=os.path.isfile)

# These are the indices of visible reflections of phaser_1.mtz when the sphere of reflections 
# have been sliced with a clip plane at  k= -9
reflections2match = set(  [(-3, -9, -1), (-3, -9, -2), (-3, -9, 0), (1, -9, -1), (4, -9, -2), 
  (4, -9, -1), (1, -9, -2), (-1, -9, -4), (1, -9, -3), (-1, -9, -3), (-2, -9, -3), (1, -9, -4), 
  (-1, -9, -1), (-1, -9, -2), (-2, -9, -1), (-2, -9, -2), (0, -9, 4), (1, -9, 4), (2, -9, -4),
  (3, -9, 1), (2, -9, -3), (0, -9, 2), (3, -9, 0), (-4, -9, 2), (2, -9, -1), (2, -9, -2),
  (0, -9, 3), (3, -9, 2), (-4, -9, 0), (0, -9, 1), (-4, -9, -1), (-4, -9, 1), (0, -9, -1),
  (0, -9, -2), (-2, -9, 4), (-1, -9, 4), (3, -9, -3), (2, -9, 0), (0, -9, -4), (2, -9, 1), 
  (0, -9, -3), (2, -9, 2), (-1, -9, 0), (3, -9, -1), (3, -9, -2), (-2, -9, 0), (2, -9, 3), 
  (-2, -9, 1), (-1, -9, 1), (1, -9, 3), (-2, -9, 2), (-1, -9, 2), (-3, -9, 3), (4, -9, 0), 
  (1, -9, 2), (-2, -9, 3), (-1, -9, 3), (-3, -9, 2), (4, -9, 1), (1, -9, 1), (-3, -9, 1), (1, -9, 0)]
 )


def exercise1():
  assert os.path.isfile(file_name)
  myHKLview = hklview_frame.HKLViewFrame(verbose=0)
  myHKLview.LoadReflectionsFile(file_name)
  # slice expanded sphere of reflections with clip plane and translate to
  # the 46th plane of reflections along the K axis

  myHKLview.update_from_philstr(philstr)
  # wait a little until NGLs autoView() has completed
  time.sleep(30)
  # then copy indices of visible reflections
  refls = myHKLview.viewer.visible_hkls[:]
  # Destroying HKLViewFrame releases javascipt objects from browser
  myHKLview.__exit__()

  # check that only the following 108 reflections in reflections2match were visible
  assert set(refls) == reflections2match

    
def exercise2():
  import re
  with open("philinput.txt","w") as f:
    f.write(philstr)
  assert os.path.isfile(file_name)
  cmdargs = ["cctbx.HKLviewer",
             file_name,
             "philinput.txt",
             "verbose=frustum", # dump displayed hkls to stdout when clipplaning
             "output_filename=myoutput.log", # file with stdout, stderr from hklview_frame
             "closingtime=30", # close HKLviewer after 30 seconds
            ]

  assert ( easy_run.call(command=" ".join(cmdargs))  == 0 )
  with open("myoutput.log", "r") as f:
    mstr = f.read()
  # peruse output file for the list of displayed reflections
  match = re.findall("visible \s+ hkls\: \s* (\[ .+ \])", mstr, re.VERBOSE)
  refls = []
  if match:
    refls = eval(match[0])
  # check that only the following 108 reflections in reflections2match were visible
  assert set(refls) == reflections2match
  # tidy up
  os.remove("philinput.txt")
  os.remove("myoutput.log")


def run():
  #exercise1()
  exercise2()
  print("OK")

if __name__ == '__main__':
  run()
