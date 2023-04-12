from __future__ import absolute_import, division, print_function

import libtbx.load_env, os.path, re, os, time, subprocess
from crys3d.hklviewer import cmdlineframes, jsview_3d
import traceback



# The tests below uses datasets in this file here
datafname = libtbx.env.find_in_repositories(
  relative_path="phenix_regression/reflection_files/1upp_lowres.mtz",
  test=os.path.isfile)

closetime = 150 # about half the maximum time each test will run
# HKLviewer uses websockets which is slightly unstable on virtual machines used on Azure.
# This might yield a bogus failure of the test. If so, repeat the test at most maxruns times
# or until it passes whichever comes first.
maxruns = 4
browser = "firefox"


philstr1 = """
external_cmd = None *runXtricorder runXtriage
clip_plane {
  hkldist = 33
  normal_vector = "TNCS_xtricorder"
  clip_width = 0.1997394793
}
binning {
  scene_bin_thresholds = -0.09 0.66 0.88 1.12 1.5 2.09
  binlabel = 'TEPS'
  bin_opacity = 0 0
  bin_opacity = 0 1
  bin_opacity = 0 2
  bin_opacity = 0 3
  bin_opacity = 1 4
  bin_opacity = 0 5
  bin_opacity = 0 6
  nbins = 5
}
viewer {
  data_array {
    label = "TEPS"
    datatype = "Floating-point"
  }
  show_vector = "['TNCS_xtricorder', True]"
  fixorientation = *vector None
}
hkls {
  expand_to_p1 = True
  expand_anomalous = True
}
max_reflections_in_frustum = 70
"""
# These are the indices of visible TEPS reflections of processed with xtricorder 1upp_lowres.mtz
# when the sphere of reflections is sliced perpendicular to the TNCS vector at layer 33 in the
# TNCS modulation and reflections have been divided into 5 bins according to TNCS modulation
# values but with explicit bin threshold values and only reflections of the highest bin are displayed
reflections2match1 = set(  [(23, 1, 11), (17, -5, 17), (16, 0, 18), (19, 7, 15), (22, 4, 12), (22, 2, 12),
  (20, -4, 14), (19, 1, 15), (18, 6, 16), (21, 5, 13), (17, 1, 17), (20, 0, 14), (16, 4, 18), (16, -2, 18),
  (19, 3, 15), (20, -6, 14), (19, -3, 15), (18, 2, 16), (23, -1, 11), (17, -1, 17), (17, 5, 17),
  (16, -4, 18), (22, -2, 12), (21, -5, 13), (23, -3, 11), (20, -2, 14), (18, -6, 16), (19, -7, 15),
  (18, -2, 16), (17, 3, 17), (17, -3, 17), (20, 2, 14), (22, -4, 12), (15, 1, 19), (21, -1, 13),
  (19, 5, 15), (20, 4, 14), (18, 4, 16), (15, -1, 19), (21, 3, 13), (20, 6, 14), (23, 3, 11), (21, -3, 13),
  (16, 2, 18), (19, -1, 15), (18, -4, 16), (21, 1, 13), (19, -5, 15), (18, 0, 16)]
 )



philstr2 = """
external_cmd = None runXtricorder *runXtriage
clip_plane {
  hkldist = 20
  normal_vector = "2-fold_twin_ 0"
  clip_width = 0.5
}
binning {
  scene_bin_thresholds = 0.7 1100 2500 20000 100000
  binlabel = 'I_lowres,SIGI_lowres'
  bin_opacity = 0 0
  bin_opacity = 0 1
  bin_opacity = 0 2
  bin_opacity = 1 3
  bin_opacity = 1 4
  bin_opacity = 1 5
  nbins = 4
}
viewer {
  data_array {
    label = "I_lowres,SIGI_lowres"
    datatype = "Intensity"
  }
  show_vector = "['2-fold_twin_ 0', True]"
  user_vector {
    label = "2-fold_twin_ 0"
  }
  fixorientation = *vector None
}
hkls {
  expand_to_p1 = True
  expand_anomalous = True
}
max_reflections_in_frustum = 30
"""
# These are the indices of visible reflections of 1upp_lowres.mtz when the sphere of reflections
# is sliced perpendicular to the twin axis detected by xtriage at layer 20 and reflections have been
# divided into 4 bins according to I_lowres values with explicit bin thresholds and only reflections
# of the highest bin with values above 20000 are displayed. Twinning suggests the pattern of the
# slice is close to 2 fold symmetry.
reflections2match2 = set( [(-23, 21, 18), (-10, 34, 6), (-19, 25, 29), (-8, 36, -2), (-23, 21, -18),
  (-34, 10, 6), (-26, 18, -17), (-36, 8, 2), (-25, 19, 29), (-21, 23, 18), (-25, 19, -29), (-36, 8, -2),
  (-21, 23, -18), (-26, 18, 17), (-34, 10, -6), (-8, 36, 2), (-34, 10, -8), (-19, 25, -29),
  (-10, 34, -6), (-34, 10, 8)]
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
  with open(fname, "a") as f:
    f.write("\nstdout, stderr in terminal: \n" + "-" * 80 + "\n")
    f.write(souterr + "\n")


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
            "verbose=4_frustum_threadingmsg_orientmsg_browser", # dump displayed hkls to stdout when clipplaning as well as verbose=2
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
  out,err = obj.communicate()
  remove_settings_result = out.decode().replace("\r\n", "\n") # omit \r\n line endings on Windows

  print("Starting the real HKLviewer test...")
  with open(prefix + "HKLviewer_philinput.txt","w") as f:
    f.write(philstr)

  outputfname = prefix + "HKLviewer.log"
  if os.path.isfile(outputfname):
    os.remove(outputfname)

  cmdargs = ["cctbx.HKLviewer",
             datafname,
             "phil_file=%sHKLviewer_philinput.txt" %prefix,
             "verbose=4_frustum_threadingmsg_orientmsg_browser", # dump displayed hkls to stdout when clipplaning as well as verbose=2
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
  out,err = obj.communicate()
  HKLviewer_result = out.decode().replace("\r\n", "\n") # omit \r\n line endings on Windows
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
