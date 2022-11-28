from __future__ import absolute_import, division, print_function

import libtbx.load_env, os.path, re, os, time
from libtbx import easy_run
from crys3d.hklviewer import cmdlineframes, jsview_3d


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

datafname = libtbx.env.find_in_repositories(
  relative_path="iotbx/regression/data/phaser_1.mtz",
  test=os.path.isfile)

# These are the indices of visible reflections of phaser_1.mtz when the sphere of reflections
# have been sliced with a clip plane at k= -9
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


def check_log_file(fname):
  with open(fname, "r") as f:
    mstr = f.read()
  # check output file that reflections are reported to have been drawn
  assert re.findall(r"RenderStageObjects\(\) has drawn reflections in the browser", mstr) != []
  # peruse output file for the list of displayed reflections
  match = re.findall(r"visible \s+ hkls\: \s* (\[ .+ \])", mstr, re.VERBOSE)
  refls = []
  if match:
    refls = eval(match[0])
  # check that only the following 108 reflections in reflections2match were visible
  return set(refls) == reflections2match


def exercise1():
  assert os.path.isfile(datafname)
  outputfname = "HKLviewer1_test.log"

  with open("environ.txt","w") as mfile:
    # print environment variables to log file
    for k,v in os.environ.items():
      mfile.write( k + "=" + v + "\n")

  with open("HKLviewer_philinput.txt","w") as f:
    f.write(philstr)

  # check we can actually open a browser
  browser = "chrome"
  #browser = "firefox"
  #browser = "default"
  _, webctrl = jsview_3d.get_browser_ctrl(browser)
  assert webctrl.open("https://get.webgl.org/")
  time.sleep(10)

  cmdargs = [datafname,
            "phil_file=HKLviewer_philinput.txt",
            "verbose=4_frustum_threadingmsg", # dump displayed hkls to stdout when clipplaning as well as verbose=2
            "image_file=HKLviewer1_testimage.png",
            "UseOSBrowser=%s" %browser,
            "output_filename=" + outputfname, # file with stdout, stderr from hklview_frame
            "closing_time=40",
          ]

  assert cmdlineframes.run(cmdargs)
  assert check_log_file(outputfname)


def exercise2():
  # First delete any settings from previous HKLviewer runs that might be present on this platform
  print("Removing any previous Qsettings...")
  assert ( easy_run.fully_buffered(command="cctbx.HKLviewer remove_settings").return_code == 0 )

  print("Starting the real HKLviewer test...")

  with open("HKLviewer_philinput.txt","w") as f:
    f.write(philstr)
  assert os.path.isfile(datafname)

  outputfname = "HKLviewer2_test.log"
  if os.path.isfile(outputfname):
    os.remove(outputfname)

  cmdargs = ["cctbx.HKLviewer",
             datafname,
             "phil_file=HKLviewer_philinput.txt",
             "verbose=4_frustum_threadingmsg", # dump displayed hkls to stdout when clipplaning as well as verbose=2
             "image_file=HKLviewer2_testimage.png",
             "output_filename=" + outputfname, # file with stdout, stderr from hklview_frame
             "closing_time=90", # close HKLviewer after 25 seconds
            ]

  result = easy_run.fully_buffered(" ".join(cmdargs))
  # write terminal output to our log file
  with open(outputfname, "a") as f:
    f.write("\nstdout in terminal: \n" + "-" * 80 + "\n")
    for line in result.stdout_lines:
      f.write(line + "\n")
    f.write("\nstderr in terminal: \n" + "-" * 80 + "\n")
    for line in result.stderr_lines:
      f.write(line + "\n")
  #assert result.return_code == 0
  assert check_log_file(outputfname)




if __name__ == '__main__':
  exercise1()
  exercise2()
  print("OK")
