# -*- coding: utf-8 -*-
from __future__ import absolute_import, division, print_function

import libtbx.load_env, os.path, re, sys, os, time, subprocess
from crys3d.hklviewer import cmdlineframes, jsview_3d
import traceback

os.environ['PYTHONIOENCODING'] = 'UTF-8'


# The tests below uses datasets in this file here
datafname = libtbx.env.find_in_repositories(
  relative_path="iotbx/regression/data/phaser_1.mtz",
  test=os.path.isfile)

closetime = 150 # about half the maximum time each test will run
# HKLviewer uses websockets which is slightly unstable on virtual machines used on Azure.
# This might yield a bogus failure of the test. If so, repeat the test at most maxruns times
# or until it passes whichever comes first.
maxruns = 4
browser = "firefox"

philstr1 = """
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
max_reflections_in_frustum = 70

"""
# These are the indices of visible reflections of phaser_1.mtz when the sphere of reflections
# have been sliced with a clip plane at k= -9, expanded to P1 and Friedel mates
reflections2match1 = set(  [(-3, -9, -1), (-3, -9, -2), (-3, -9, 0), (1, -9, -1), (4, -9, -2),
  (4, -9, -1), (1, -9, -2), (-1, -9, -4), (1, -9, -3), (-1, -9, -3), (-2, -9, -3), (1, -9, -4),
  (-1, -9, -1), (-1, -9, -2), (-2, -9, -1), (-2, -9, -2), (0, -9, 4), (1, -9, 4), (2, -9, -4),
  (3, -9, 1), (2, -9, -3), (0, -9, 2), (3, -9, 0), (-4, -9, 2), (2, -9, -1), (2, -9, -2),
  (0, -9, 3), (3, -9, 2), (-4, -9, 0), (0, -9, 1), (-4, -9, -1), (-4, -9, 1), (0, -9, -1),
  (0, -9, -2), (-2, -9, 4), (-1, -9, 4), (3, -9, -3), (2, -9, 0), (0, -9, -4), (2, -9, 1),
  (0, -9, -3), (2, -9, 2), (-1, -9, 0), (3, -9, -1), (3, -9, -2), (-2, -9, 0), (2, -9, 3),
  (-2, -9, 1), (-1, -9, 1), (1, -9, 3), (-2, -9, 2), (-1, -9, 2), (-3, -9, 3), (4, -9, 0),
  (1, -9, 2), (-2, -9, 3), (-1, -9, 3), (-3, -9, 2), (4, -9, 1), (1, -9, 1), (-3, -9, 1), (1, -9, 0)]
 )


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
max_reflections_in_frustum = 40

"""
# These are the indices of visible reflections of phaser_1.mtz when the sphere of reflections
# have been sliced with a clip plane at h= -8
reflections2match2 = set(  [(-8, 4, 3), (-8, 4, 2), (-8, 3, 4), (-8, 4, 1), (-8, 3, 5), (-8, 5, 4),
   (-8, 1, 1), (-8, 4, 5), (-8, 1, 2), (-8, 4, 4), (-8, 1, 3), (-8, 1, 4), (-8, 2, 1), (-8, 1, 5),
   (-8, 5, 1), (-8, 0, 7), (-8, 1, 6), (-8, 5, 2), (-8, 2, 3), (-8, 0, 6), (-8, 1, 7), (-8, 5, 3),
   (-8, 2, 2), (-8, 2, 5), (-8, 0, 4), (-8, 3, 2), (-8, 2, 4), (-8, 0, 3), (-8, 3, 3), (-8, 0, 2),
   (-8, 2, 6), (-8, 3, 1), (-8, 3, 6)]
 )


# Create an array of F/SigF values, make 6 bins of reflections of equal size sorted with values and
# select only reflections from the two bins with the lowest F/SigF values. Then save those to a new file
philstr3 = """
miller_array_operation = "('newarray._data= array1.data()/array1.sigmas()\\nnewarray._sigmas = None', 'FoverSigF2', ['FOBS,SIGFOBS', 'Amplitude'], ['', ''])"
clip_plane {
  hkldist = 9
  normal_vector = "L-axis (0,0,1)"
  is_assoc_real_space_vector = True
  clip_width = 2.030853224
}
binning {
  binlabel = 'FoverSigF2'
  bin_opacity = 1 0
  bin_opacity = 1 1
  bin_opacity = 0 2
  bin_opacity = 0 3
  bin_opacity = 0 4
  bin_opacity = 0 5
  nbins = 6
}
viewer {
  data_array {
    label = "FoverSigF2"
    datatype = "Amplitude"
  }
  fixorientation = *vector None
}
hkls {
  expand_to_p1 = True
  expand_anomalous = True
}
visible_dataset_label = "LowValuesFSigF"
savefilename = "%s"
datasets_to_save = 8
max_reflections_in_frustum = 40
"""
# These are the indices of visible reflections of phaser_1.mtz of the F/SigF dataset created on the fly
# where the sphere of reflections have been sliced with a clip plane at l=9 and only reflections
# with F/SigF<=1 are displayed
reflections2match3 = set([(-5, -2, 9), (2, 1, 9), (-3, -1, 9), (-3, -2, 9), (-2, 4, 9), (-3, 0, 9),
   (-5, 0, 9), (-5, 2, 9), (-4, 1, 9), (-1, 0, 9), (-3, 2, 9), (-1, -3, 9), (-2, -3, 9), (0, -4, 9),
   (-2, 1, 9), (1, -2, 9), (-2, 3, 9), (-1, 3, 9), (2, -1, 9), (2, -2, 9), (-4, -1, 9), (-3, 1, 9),
   (0, 4, 9), (-2, -4, 9), (-2, -1, 9), (2, 2, 9), (1, 2, 9), (-2, 0, 9)]
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


def Append2LogFile(fname, souterr):
  # write terminal output to our log file
  str1 = souterr.decode().replace("\r\n", "\n") # omit \r\n line endings on Windows
  str2 = str(str1).encode(sys.stdout.encoding, errors='ignore').decode(sys.stdout.encoding)
  with open(fname, "a", encoding="utf-8") as f:
    f.write("\nstdout, stderr in terminal: \n" + "-" * 80 + "\n")
    f.write(str2 + "\n")


def exercise_OSbrowser(philstr, refl2match, prefix=""):
  assert os.path.isfile(datafname)
  outputfname = prefix + "HKLviewer.log"

  with open(prefix + "environ.txt","w") as mfile:
    # print environment variables to log file
    for k,v in os.environ.items():
      mfile.write( k + "=" + v + "\n")

  with open(prefix + "HKLviewer_philinput.txt","w") as f:
    f.write(philstr)

  # check we can actually open a browser
  browserpath, webctrl = jsview_3d.get_browser_ctrl(browser)
  #assert webctrl.open("https://get.webgl.org/")
  #subprocess.run('"' + browserpath + '"  https://get.webgl.org/ &', shell=True,
  #               capture_output=False, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
  #time.sleep(10)

  cmdargs = [datafname,
            "phil_file=%sHKLviewer_philinput.txt" %prefix,
            "verbose=2_frustum_threadingmsg_orientmsg_browser", # dump displayed hkls to stdout when clipplaning as well as verbose=2
            "image_file=%sHKLviewer.png" %prefix,
            "UseOSBrowser=%s" %browser,
            "output_filename=" + outputfname, # file with stdout, stderr from hklview_frame
            "closing_time=%d" %closetime,
          ]

  cmdlineframes.run(cmdargs)
  print("=" * 80)
  check_log_file(outputfname, refl2match)


def exerciseQtGUI(philstr, refl2match, prefix=""):
  # These flags enables QWebEngine, Qt5.15 to work on VMs used on Azure
  os.environ["QTWEBENGINE_CHROMIUM_FLAGS"] = " --disable-web-security" \
            + " --enable-webgl-software-rendering --disable-gpu-compositing" \
            + " --disable_chromium_framebuffer_multisample --use-gl=swiftshader" \
            + " --swiftshader --swiftshader-webgl --ignore-gpu-blocklist"
  with open(prefix + "environ.txt","w") as mfile:
    # print environment variables to log file
    for k,v in os.environ.items():
      mfile.write( k + "=" + v + "\n")

  assert os.path.isfile(datafname)
  # First delete any settings from previous HKLviewer runs that might be present on this platform
  print("Removing any previous Qsettings...")
  obj = subprocess.Popen("cctbx.HKLviewer remove_settings",
                         shell=True,
                         env = os.environ,
                         stdin = subprocess.PIPE,
                         stdout = subprocess.PIPE,
                         stderr = subprocess.STDOUT)
  remove_settings_result,err = obj.communicate()

  print("Starting the real HKLviewer test...")
  with open(prefix + "HKLviewer_philinput.txt","w") as f:
    f.write(philstr)

  outputfname = prefix + "HKLviewer.log"
  if os.path.isfile(outputfname):
    os.remove(outputfname)

  cmdargs = ["cctbx.HKLviewer",
             datafname,
             "phil_file=%sHKLviewer_philinput.txt" %prefix,
             "verbose=2_frustum_threadingmsg_orientmsg_browser", # dump displayed hkls to stdout when clipplaning as well as verbose=2
             "image_file=%sHKLviewer.png" %prefix,
             "output_filename=" + outputfname, # file with stdout, stderr from hklview_frame
             "closing_time=%d" %closetime, # close HKLviewer after 25 seconds
            ]

  obj = subprocess.Popen(" ".join(cmdargs),
                         shell=True,
                         env = os.environ,
                         stdin = subprocess.PIPE,
                         stdout = subprocess.PIPE,
                         stderr = subprocess.STDOUT)
  HKLviewer_result,err = obj.communicate()
  # append terminal output to log file
  Append2LogFile(outputfname, remove_settings_result)
  Append2LogFile(outputfname, HKLviewer_result)
  print("retval: " + str(obj.returncode))
  print("=" * 80)
  check_log_file(outputfname, refl2match)



def runagain(func, philstr, refl2match, name):
  try:
    func(philstr, refl2match, name)
    print("OK\n")
    return False
  except Exception as e:
    print( str(e) + traceback.format_exc(limit=10))
    return True
