from __future__ import absolute_import, division, print_function

import libtbx.load_env, os.path, re, os, time, subprocess
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
hkls.expand_to_p1 = True
hkls.expand_anomalous = True
use_wireframe = True

"""


philstr2 = """
clip_plane {
  hkldist = -8
  normal_vector = "H-axis (1,0,0)"
  is_assoc_real_space_vector = True
  clip_width = 2.030853224
}
binning {
  nbins = 2
}
viewer {
  data_array {
    label = "FP,SIGFP"
    datatype = "Amplitude"
  }
  show_vector = "['K-axis (0,1,0)', True]"
  fixorientation = *vector None
}
use_wireframe = True

"""

datafname = libtbx.env.find_in_repositories(
  relative_path="iotbx/regression/data/phaser_1.mtz",
  test=os.path.isfile)

# These are the indices of visible reflections of phaser_1.mtz when the sphere of reflections
# have been sliced with a clip plane at k= -9, expanded to P1 and Friedel mates
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


# These are the indices of visible reflections of phaser_1.mtz when the sphere of reflections
# have been sliced with a clip plane at h= -8
reflections2match2 = set(  [(-8, 4, 3), (-8, 4, 2), (-8, 3, 4), (-8, 4, 1), (-8, 3, 5), (-8, 5, 4),
   (-8, 1, 1), (-8, 4, 5), (-8, 1, 2), (-8, 4, 4), (-8, 1, 3), (-8, 1, 4), (-8, 2, 1), (-8, 1, 5),
   (-8, 5, 1), (-8, 0, 7), (-8, 1, 6), (-8, 5, 2), (-8, 2, 3), (-8, 0, 6), (-8, 1, 7), (-8, 5, 3),
   (-8, 2, 2), (-8, 2, 5), (-8, 0, 4), (-8, 3, 2), (-8, 2, 4), (-8, 0, 3), (-8, 3, 3), (-8, 0, 2),
   (-8, 2, 6), (-8, 3, 1), (-8, 3, 6)]
 )





def check_log_file(fname, refls2match):
  with open(fname, "r") as f:
    mstr = f.read()
  # check output file that reflections are reported to have been drawn
  assert re.findall(r"RenderStageObjects\(\) has drawn reflections in the browser", mstr) != []
  # peruse output file for the list of displayed reflections
  match = re.findall(r"visible \s+ hkls\: \s* (\[ .+ \])", mstr, re.VERBOSE)
  refls = []
  if match:
    refls = eval(match[-1]) # use the last match of reflections in the log file
  # check that only the following 108 reflections in refls2match were visible
  setrefls = set(refls)
  if setrefls != refls2match:
    print("Indices of visible reflections:\n%s" %str(setrefls))
    print("Do not match the expected indices:\n%s" %str(refls2match))
  assert setrefls == refls2match
  print("Indices of visible reflections match the expected ones.")


def Append2LogFile(fname, res):
  # write terminal output to our log file
  with open(fname, "a") as f:
    f.write("\nstdout in terminal: \n" + "-" * 80 + "\n")
    for line in res.stdout_lines:
      f.write(line + "\n")
    f.write("\nstderr in terminal: \n" + "-" * 80 + "\n")
    for line in res.stderr_lines:
      f.write(line + "\n")


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
  #browser = "chrome"
  browser = "firefox"
  #browser = "default"
  browserpath, webctrl = jsview_3d.get_browser_ctrl(browser)
  #assert webctrl.open("https://get.webgl.org/")
  #subprocess.run('"' + browserpath + '"  https://get.webgl.org/ &', shell=True,
  #               capture_output=False, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
  #time.sleep(10)

  cmdargs = [datafname,
            "phil_file=HKLviewer_philinput.txt",
            "verbose=4_frustum_threadingmsg_orientmsg_browser", # dump displayed hkls to stdout when clipplaning as well as verbose=2
            "image_file=HKLviewer1_testimage.png",
            "UseOSBrowser=%s" %browser,
            "output_filename=" + outputfname, # file with stdout, stderr from hklview_frame
            "closing_time=140",
            #"debug=True"
          ]

  assert cmdlineframes.run(cmdargs)

  check_log_file(outputfname, reflections2match)


def exercise2():
  assert os.path.isfile(datafname)
  # First delete any settings from previous HKLviewer runs that might be present on this platform
  print("Removing any previous Qsettings...")
  remove_settings_result = easy_run.fully_buffered(command="cctbx.HKLviewer remove_settings")

  print("Starting the real HKLviewer test...")

  with open("HKLviewer_philinput.txt","w") as f:
    f.write(philstr)

  outputfname = "HKLviewer2_test.log"
  if os.path.isfile(outputfname):
    os.remove(outputfname)

  cmdargs = ["cctbx.HKLviewer",
             datafname,
             "phil_file=HKLviewer_philinput.txt",
             "verbose=4_frustum_threadingmsg_orientmsg_browser", # dump displayed hkls to stdout when clipplaning as well as verbose=2
             "image_file=HKLviewer2_testimage.png",
             "output_filename=" + outputfname, # file with stdout, stderr from hklview_frame
             "closing_time=140", # close HKLviewer after 25 seconds
             #"debug=True"
            ]

  HKLviewer_result = easy_run.fully_buffered(" ".join(cmdargs))
  # append terminal output to log file
  Append2LogFile(outputfname, remove_settings_result)
  Append2LogFile(outputfname, HKLviewer_result)
  print("retval: " + str(HKLviewer_result.return_code))
  #assert HKLviewer_result.return_code == 0
  #assert remove_settings_result.return_code == 0
  check_log_file(outputfname, reflections2match)




if __name__ == '__main__':
  exercise1()
  exercise2()
  print("OK")
