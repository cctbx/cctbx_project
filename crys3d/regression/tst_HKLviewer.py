from __future__ import absolute_import, division, print_function

import libtbx.load_env, os.path, time
from libtbx import easy_run
from crys3d.hklviewer import hklview_frame


philstr = """
clip_plane {
  normal_vector = "K-axis"
  is_assoc_real_space_vector = True
  clip_width = 1.184
  hkldist = 46
  normal_vector_length_scale = -1
}
viewer {
  data_array {
    label = 'F,SIGF'
    datatype = 'Amplitude'
  }
  show_vector = "['K-axis', True]"
  fixorientation = *vector None
}
hkls {
  expand_to_p1 = True
  expand_anomalous = True
}
"""

file_name = libtbx.env.find_in_repositories(
  relative_path="phenix_regression/reflection_files/2caz.mtz",
  test=os.path.isfile)

reflections2match = set( [(-24, 46, 1), (-17, 46, 0), (-28, 46, 1), (-24, 46, 0), (-17, 46, 1),
                    (-28, 46, 0), (-24, 46, -2), (-29, 46, 0), (-22, 46, 2), (-18, 46, -1),
                    (-26, 46, 2), (-24, 46, 2), (-29, 46, 1), (-22, 46, 1), (-26, 46, 1),
                    (-22, 46, 0), (-19, 46, -1), (-23, 46, 0), (-26, 46, 0), (-27, 46, 0),
                    (-23, 46, 1), (-21, 46, 2), (-27, 46, 1), (-25, 46, 2), (-20, 46, 1),
                    (-23, 46, 2), (-25, 46, -1), (-25, 46, -2), (-20, 46, 0), (-21, 46, -1),
                    (-21, 46, -2), (-21, 46, 0), (-25, 46, 0), (-21, 46, 1), (-25, 46, 1),
                    (-20, 46, 2), (-20, 46, -1), (-27, 46, -1), (-20, 46, -2), (-23, 46, -1),
                    (-23, 46, -2), (-19, 46, 0), (-19, 46, 1), (-26, 46, -1), (-26, 46, -2),
                    (-22, 46, -1), (-29, 46, -1), (-22, 46, -2), (-18, 46, 1), (-17, 46, -1),
                    (-18, 46, 0), (-28, 46, -1), (-24, 46, -1)]
              )


def exercise1():
  myHKLview = hklview_frame.HKLViewFrame(verbose=0)
  myHKLview.LoadReflectionsFile(file_name)
  # slice expanded sphere of reflections with clip plane and translate to
  # the 46th plane of reflections along the K axis

  myHKLview.update_from_philstr(philstr)
  # wait a little until NGLs autoView() has completed
  time.sleep(2)
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
  os.remove("philinput.txt")
  os.remove("myoutput.log")


def run():
  #exercise1()
  exercise2()
  print("OK")

if __name__ == '__main__':
  run()
