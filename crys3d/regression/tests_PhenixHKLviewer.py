from __future__ import absolute_import, division, print_function

import libtbx.load_env, os.path, re, os, time, subprocess
from crys3d.hklviewer import cmdlineframes, jsview_3d
import traceback



# The tests below uses datasets in this file here
datafname = libtbx.env.find_in_repositories(
  relative_path=r"C:\Users\oeffner\Work\HKLviewerTests\1upp_lowres.mtz",
  test=os.path.isfile)

closetime = 150 # about half the maximum time each test will run
maxruns = 4 # maximum number to repeat unstable test until it passes
browser = "firefox"


philstr1 = """
external_cmd = None *runXtricorder runXtriage
clip_plane {
  hkldist = 32
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

"""
# These are the indices of visible TEPS reflections of processed with xtricorder 1upp_lowres.mtz
# when the sphere of reflections
# is sliced perpendicular to the TNCS vector at layer 32 in the TNCS modulation and reflections have been
# divided into 6 bins according to TNCS modulation values and only reflections of the two highest bins
# are displayed
reflections2match1 = set(  [(21, 5, 11), (20, 10, 12), (18, 6, 14), (15, 7, 17), (15, 9, 17), (13, 7, 19),
     (14, 8, 18), (16, 10, 16), (18, 10, 14), (17, 7, 15), (21, 7, 11), (17, 9, 15), (16, 8, 16), (23, 7, 9),
     (19, 9, 13), (22, 6, 10), (22, 8, 10), (20, 6, 12), (18, 8, 14), (23, 5, 9), (21, 9, 11), (20, 8, 12),
     (19, 7, 13), (16, 6, 16)]
 )



philstr2 = """
clip_plane {
  hkldist = 13
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

"""
# These are the indices of visible reflections of 1upp_lowres.mtz when the sphere of reflections
# is sliced perpendicular to the twin axis detected by xtriage at layer 13 and reflections have been
# divided into 4 bins according to I_lowres values and only reflections of the highest bin with
# values above 20000 are displayed. Twinning suggests the pattern of the slice should be close to 2 fold symmetry
reflections2match1 = set( [(-33, -7, -7), (-19, 7, 5), (-19, 7, -1), (-7, 19, 1), (-32, -6, 20),
   (7, 33, 7), (-7, 19, 5), (-7, 19, -1), (-32, -6, -20), (-7, 19, -5), (7, 33, -7), (-20, 6, 0),
   (6, 32, -20), (-19, 7, -5), (-33, -7, 7), (-6, 20, 0), (-19, 7, 1), (6, 32, 20)]
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
  #if "linux" in sys.platform:
  #  os.environ["DISPLAY"] = "O:O"
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
  #if "linux" in sys.platform:
  #  os.environ["DISPLAY"] = "O:O"
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
