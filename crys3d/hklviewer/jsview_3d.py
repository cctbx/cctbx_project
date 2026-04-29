# -*- coding: utf-8 -*-
from __future__ import absolute_import, division, print_function
from libtbx.math_utils import roundoff
import traceback
from crys3d.hklviewer import display2 as display
from cctbx.array_family import flex
from cctbx import miller, sgtbx
from scitbx import graphics_utils
from scitbx import matrix
import scitbx.math
from libtbx import group_args
from libtbx.test_utils import approx_equal
from libtbx.utils import Sorry, to_str
import threading, math, sys, cmath
from crys3d.hklviewer.webbrowser_messenger_py3 import WBmessenger

import os.path, time, copy, re, io, subprocess
import libtbx
import webbrowser, tempfile
from six.moves import range


class HKLviewerError(Exception):
  # Unrecoverable errors we detect such as failure to connect via websocket or no WebGL
  def __init__(self, value):
    self.value = value
  def __str__(self):
    return(repr(self.value))


def has_phil_path(philobj, *paths): # variable number of arguments
  for path in paths:
    if len([ e.path for e in philobj.all_definitions() if path in e.path.split(".") ]):
      return True
  return False


def MakeHKLscene( proc_array, foms_array, pidx, fidx, renderscale, hkls, mprint=sys.stdout.write):
  """
  Conpute the hklscene for proc_array. If it's a complex array and foms_array!=None
  then also compute an hklscene with colours of each hkl attenuated by the corresponding FOM value.
  """
  from iotbx.gui_tools.reflections import ArrayInfo
  scenemaxdata =[]
  scenemindata =[]
  scenemaxsigmas = []
  sceneminsigmas = []
  scenearrayinfos = []
  hklscenes = []
  if (hkls.expand_anomalous or hkls.expand_to_p1) \
      and not proc_array.is_unique_set_under_symmetry():
    mprint("The " + proc_array.info().label_string() + \
         " array is not symmetry unique. Expansion may lead to multiple reflections rendered on top of each other.", verbose=1)
  hkls.expand_anomalous = False # don't expand reflections in display2.py as this now is done in HKLJavaScripts.js
  hkls.expand_to_p1 = False    #  TODO: remove this functionality altogether from display2.py

  hklscene = display.scene(miller_array=proc_array, merge=None, renderscale=renderscale,
    settings=hkls, foms_array=foms_array, fullprocessarray=True, mprint=mprint)
  if not hklscene.SceneCreated:
    mprint("The " + proc_array.info().label_string() + " array was not processed")
  #import code, traceback; code.interact(local=locals(), banner="".join( traceback.format_stack(limit=10) ) )
  # cast any NAN values to 1 of the colours and radii to 0.2 before writing javascript
  else:
    hklscenes.append( hklscene)
    hklscene.colors = graphics_utils.NoNansvec3( hklscene.colors, 1.0, 1.0, 1.0)
    hklscene.radii = graphics_utils.NoNansArray( hklscene.radii, 0.2)
    fomslabel = None
    if foms_array:
      fomslabel = foms_array.info().label_string()
    arrayinfo = ArrayInfo(hklscene.work_array)
    scenemindata.append(arrayinfo.minmaxdata[0])
    scenemaxdata.append(arrayinfo.minmaxdata[1])
    sceneminsigmas.append(arrayinfo.minmaxsigs[0])
    scenemaxsigmas.append(arrayinfo.minmaxsigs[1])
    lbl = arrayinfo.labelstr
    hassigmas=True
    if math.isnan(arrayinfo.maxsigmas):
      hassigmas=False
    if fomslabel:
      lbl = arrayinfo.labelstr + " + " + fomslabel
    (dummy1, infolst, dummy2, dummy3), dummy4, dummy5 = arrayinfo.get_selected_info_columns_from_phil()
    scenearrayinfos.append([infolst, pidx, fidx, lbl, infolst[1], hassigmas])
  return (hklscenes, scenemaxdata, scenemindata, scenemaxsigmas, sceneminsigmas, scenearrayinfos)


def get_browser_ctrl(using=None):
  if using is None or using=="default":
    return "default", webbrowser.get()

  if using=="firefox":
    if sys.platform == "win32":
      browser = "C:/Program Files/Mozilla Firefox/firefox.exe"
      if not os.path.isfile(browser):
        browser = "C:/Program Files (x86)/Mozilla Firefox/firefox.exe"
    if sys.platform.startswith("darwin"):
      browser = "/Applications/Firefox.app/Contents/MacOS/firefox"
    if sys.platform == "linux":
      browser = "/usr/bin/firefox"

  if using=="chrome":
    if sys.platform == "win32":
      browser = "C:/Program Files (x86)/Google/Chrome/Application/chrome.exe"
      if not os.path.isfile(browser):
        browser = "C:/Program Files/Google/Chrome/Application/chrome.exe"
    if sys.platform.startswith("darwin"):
      browser = '"/Applications/Google Chrome.app/Contents/MacOS/Google Chrome"'
      if not os.path.isfile(browser):
        browser = '"/Applications/Google Chrome.app"'
    if sys.platform == "linux":
      #pass
      browser = "/usr/bin/firefox"

  webbrowser.register(using, None, webbrowser.BackgroundBrowser(browser))
  webctrl = webbrowser.get(using)
  assert os.path.isfile(browser)
  return browser, webctrl


lock_timeout=120 # for the sempahores. Rendering could take a while for very large file. Until that
# has been completed, geometries of the NGL stage such as clipnear, clipfar, cameraZ and bounding box
# are undefined. websocket connection could take a while on Azure pipelines


class HKLview_3d:
  def __init__ (self, *args, **kwds) :
    self.diff_phil = None
    self.params = None # first assigned in HKLViewFrame().ResetPhil()
    self.miller_array = None
    self.symops = []
    self.sg = None
    self.tooltipstrings = []
    self.tooltipstringsdict = {}
    self.d_min = None
    self.scene = None
    self.lastscene_id = None
    self.merge = False
    self.url = ""
    self.bin_labels_type_idxs = []
    self.colours = []
    self.positions = []
    self.radii2 = []
    self.spbufttips = []
    self.rot_recip_zvec = None
    self.rot_zvec = None
    self.meanradius = -1
    self.past = time.time()
    self.orientmessage = None
    self.clipNear = None
    self.clipFar = None
    self.cameraPosZ = None
    self.zoom = None
    self.cameratranslation = ( 0,0,0 )
    self.planescalarvalue =0
    self.planenormalhklvec =None
    self.renderscale = 100
    #self.angle_x_svec = 0.0
    #self.angle_y_svec = 0.0
    self.angle_z_svec = 0.0
    #self.angle_z_yzvec = 0.0
    #self.angle_y_yzvec = 0.0
    #self.angle_y_xyvec = 0.0
    self.angle_x_xyvec = 0.0
    self.vecrotmx = None
    self.currentrotvec = None
    self.unit_h_axis = None
    self.unit_k_axis = None
    self.unit_l_axis = None
    self.normal_hk = None
    self.normal_kl = None
    self.normal_lh = None
    self.normal_vecnr = -1
    self.isnewfile = False
    self.has_new_miller_array = False
    self.sleeptime = 0.01 # 0.025 # used in sleep() for ProcessBrowserMessage and elsewhere in WBmessenger
    self.binvals = []
    self.binvalsboundaries = []
    self.oldnbinvalsboundaries = None
    self.proc_arrays = []
    self.HKLscene = []
    self.HKLscenedict = {}
    self.HKLscenesdict = {}
    self.HKLscenesMaxdata = []
    self.HKLscenesMindata = []
    self.HKLscenesMaxsigmas = []
    self.HKLscenesMinsigmas = []
    self.bindata = None
    self.reciproc_scale = 1.0
    self.realspace_scale = 1.0
    self.visual_symHKLs = []
    self.visual_symmxs= []
    self.visible_hkls = [] # Populated when applying clip planes. To be examined in regression tests
    self.outsideplane_hkls = []
    self.sceneisdirty = True
    self.max_reflections_in_frustum = 0
    self.imagename = None
    self.imgdatastr = ""
    self.hkl_scenes_info = []
    self.hkl_scenes_infos = []
    self.match_valarrays = []
    self.array_info_format_tpl = []
    self.binstrs = []
    self.rotation_operators = []
    self.all_vectors = []
    self.realSpaceMag = 1
    self.recipSpaceMag = 1
    self.has_unmerged_data = False
    self.cosine = 1
    self.L = 1.0
    self.nuniqueval = 0
    self.bin_infotpls = []
    self.executing_preset_btn = False
    self.mapcoef_fom_dict = {}
    # colourmap=brg, colourpower=1, nth_power_scale_radii=nan, radiiscale=1
    self.datatypedefault = ["brg", 1.0, float('nan'), 1.0]
    self.datatypedict = { }
    self.sceneid_from_arrayid = []
    self.parent = None
    if 'parent' in kwds:
      self.parent = kwds['parent']
    self.debug = None
    if 'debug' in kwds:
      self.debug = eval( kwds['debug'])
    self.mprint = sys.stdout.write
    if 'mprint' in kwds:
      self.mprint = kwds['mprint']
    self.nbinvalsboundaries = 0
    self.websockport = 7894
    if 'websockport' in kwds:
      self.websockport = kwds['websockport']
    tempdir = tempfile.gettempdir()
    # ensure unique file name by including port number in filename
    self.hklfname = os.path.join(tempdir, "hkl_%d.htm" %self.websockport )
    if os.path.isfile(self.hklfname):
      os.remove(self.hklfname)
    if 'htmlfname' in kwds and kwds['htmlfname']:
      self.hklfname = kwds['htmlfname']
    self.hklfname = os.path.abspath( self.hklfname )
    self.isHKLviewer= "false"
    self.verbose = 1
    self.verbosebrowser = False
    if 'verbose' in kwds:
      try:
        self.verbose = eval(kwds['verbose'])
      except Exception as e:
        self.verbose = kwds['verbose']
        if "browser" in self.verbose:
          self.verbosebrowser = True
    self.send_info_to_gui = None
    if 'send_info_to_gui' in kwds:
      self.send_info_to_gui = kwds['send_info_to_gui']
      self.isHKLviewer= "true"
    if 'fileinfo' in kwds:
      return
    self.mprint('Rendering done via websocket in \"%s\"'  %self.hklfname)
    self.hklhtml = r"""
<!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 4.0 Transitional//EN">
<html><head><meta charset="utf-8" /></head>
<body>
<script>var isHKLviewer = %s; </script>
<script>var websocket_portnumber = %s; </script>
<script src="%s" type="text/javascript"></script>
<script src="%s" type="text/javascript"></script>
<script src="%s" type="text/javascript"></script>
<script src="%s" type="text/javascript"></script>
<div id="viewport" style="width:100%%; height:100%%;"></div>
</body></html>

    """
    WeblglChecklibpath = os.path.join(os.path.dirname(os.path.abspath(__file__)), "webgl_check.js")
    Html2Canvaslibpath = os.path.join(os.path.dirname(os.path.abspath(__file__)), "html2canvas.min.js")
    NGLlibpath = os.path.join(os.path.dirname(os.path.abspath(__file__)), "ngl.js")
    HKLjscriptpath = os.path.join(os.path.dirname(os.path.abspath(__file__)), "HKLJavaScripts.js")
    WeblglCheckliburl = "file:///" + WeblglChecklibpath.replace("\\","/")
    Html2Canvasliburl = "file:///" + Html2Canvaslibpath.replace("\\","/")
    NGLliburl = "file:///" + NGLlibpath.replace("\\","/")
    HKLjscripturl = "file:///" + HKLjscriptpath.replace("\\","/")
    self.htmlstr = self.hklhtml %(self.isHKLviewer,
                                  self.websockport,
                                  WeblglCheckliburl,
                                  Html2Canvasliburl,
                                  NGLliburl,
                                  HKLjscripturl)
    self.colourgradientvalues = []
    self.browserpath, self.webctrl = get_browser_ctrl()
    self.UseOSBrowser = ""
    if 'useGuiSocket' not in kwds:
      self.UseOSBrowser = "default"
    ldic=locals()
    if 'UseOSBrowser' in kwds:
      exec("UseOSBrowser = kwds['UseOSBrowser']", globals(), ldic)
      self.UseOSBrowser = ldic["UseOSBrowser"]
      self.UseOSBrowser = self.UseOSBrowser.replace("\\","/")
      #if self.UseOSBrowser != "default" and not os.path.isfile(self.UseOSBrowser):
      #  raise Sorry("Error: %s does not exist" %self.UseOSBrowser)
      self.browserpath, self.webctrl = get_browser_ctrl(self.UseOSBrowser)

    self.viewmtrx = None
    self.lastviewmtrx = None
    self.currentRotmx = matrix.identity(3)
    self.include_tooltip_lst = []
    self.mouse_moved = False
    self.webgl_OK = True
    self.HKLsceneKey = None
    self.handshakewait = 5
    if 'handshakewait' in kwds:
      self.handshakewait = eval(kwds['handshakewait'])
    self.lastmsg = "" # "Ready"
    self.use_semaphore = True
    self.clipplane_msg_sem = threading.Semaphore()
    self.mousespeed_msg_sem = threading.BoundedSemaphore()
    self.hkls_drawn_sem = threading.Semaphore()
    self.autoview_sem = threading.BoundedSemaphore()
    self.browser_connect_sem = threading.BoundedSemaphore()
    self.WBmessenger = WBmessenger(self)
    self.AddToBrowserMsgQueue = self.WBmessenger.AddToBrowserMsgQueue
    # Don't overwhelm ProcessBrowserMessage() with the flurry of tooltips that can get emitted by javascript
    # when tooltips are shown for a large dataset as this makes the HKLviewer unresponsive.
    # replace_msg_lst is a list of strings that matches the
    # particular messages we want to replace rather than append to the message queue.
    self.WBmessenger.replace_msg_lst = ["tooltip_id:"]
    self.WBmessenger.StartWebsocket()
    self.javascriptcleaned = False


  def __exit__(self, exc_type, exc_value, traceback):
    # not called unless instantiated with a "with HKLview_3d ... " statement
    self.JavaScriptCleanUp()
    self.SendInfoToGUI( { "datatype_dict": self.datatypedict } ) # so the GUI can persist these across sessions
    nwait = 0
    if self.params.viewer.scene_id is None:
      self.WBmessenger.StopWebsocket()
    while not self.WBmessenger.isterminating and nwait < 5:
      time.sleep(self.sleeptime)
      nwait += self.sleeptime
    if os.path.isfile(self.hklfname):
      os.remove(self.hklfname)
    self.mprint("Destroying HKLview_3d", 1)


  def SendInfoToGUI(self, mydict):
    if self.send_info_to_gui:
      self.send_info_to_gui( mydict )


  def GuardedAddToBrowserMsgQueue(self, semaphorename, msgtype, msg="", binary=False,
                                  funcname="", posteriorcheck=True):
    if self.use_semaphore:
      semaphore = self.__dict__.get(semaphorename, None)
      self.mprint("%s waiting for %s.acquire" %(funcname, semaphorename), verbose="threadingmsg")
      if not semaphore.acquire(blocking=True, timeout=lock_timeout):
        self.mprint("Error! Timed out waiting for %s semaphore within %s seconds" %(semaphorename,lock_timeout), verbose=1)
      self.mprint("%s got %s" %(funcname, semaphorename), verbose="threadingmsg")

    self.AddToBrowserMsgQueue( msgtype, msg, binary)
    time.sleep(0.1)
    if posteriorcheck and self.use_semaphore:
      self.mprint("%s waiting for %s.acquire again" %(funcname, semaphorename), verbose="threadingmsg")
      if not semaphore.acquire(blocking=True, timeout=lock_timeout):
        self.mprint("Error! Timed out waiting for %s semaphore within %s seconds" %(semaphorename,lock_timeout), verbose=1)
      self.mprint("%s got %s again" %(funcname, semaphorename), verbose="threadingmsg")
      semaphore.release()
      self.mprint("%s released %s" %(funcname, semaphorename), verbose="threadingmsg")


  def process_PHIL_parameters(self, diff_phil, curphilparam) :
    """
    Event handler for zmq messages from the GUI or simply for commandline interaction when
    scripting HKLviewer with python
    """
    self.params = curphilparam
    self.diff_phil = diff_phil

    if has_phil_path(diff_phil,
            "openfilename",
            "use_provided_miller_arrays",
            "spacegroup_choice",
            "using_space_subgroup",
            "camera_type",
            "miller_array_operation",
       ) or has_phil_path(diff_phil, "viewer") \
          and has_phil_path(diff_phil,
                  "show_missing",
                  "show_only_missing",
                  "show_systematic_absences",
                  "slice_axis",
                  "slice_index",
                  "sigma_color_radius",
                  "scene_id",
                  "data_array",
                  "color_scheme",
                  "color_powscale",
                  "scale",
                  "nth_power_scale_radii"
              ) :
                  self.sceneisdirty = True
                  if has_phil_path(diff_phil,
                          "spacegroup_choice",
                          "show_missing",
                          "show_only_missing",
                          "show_systematic_absences",
                          "slice_axis",
                          "slice_index",
                          "sigma_color_radius",
                          "scene_id",
                          "data_array",
                          "use_provided_miller_arrays",
                          "color_scheme",
                          "color_powscale",
                          "scale",
                          "nth_power_scale_radii"
                     ):
                        self.ConstructReciprocalSpace(scene_id=self.params.viewer.scene_id )
    msg = ""
    if self.params.viewer.scene_id is not None and \
       has_phil_path(diff_phil,
            #"scene_id",
            "show_missing",
            "show_only_missing",
            "show_systematic_absences",
            "binner_idx",
            "binlabel",
            "nbins",
        ) and not has_phil_path(diff_phil, "scene_bin_thresholds"):
      self.binvalsboundaries = []
      self.binvals, self.nuniqueval = self.calc_bin_thresholds(curphilparam.binning.binner_idx,
                                                               curphilparam.binning.nbins)
      self.sceneisdirty = True

    self.params.max_reflections_in_frustum

    if has_phil_path(diff_phil, "max_reflections_in_frustum"):
      self.max_reflections_in_frustum = self.params.max_reflections_in_frustum

    if has_phil_path(diff_phil, "sigma_color_radius"):
      self.sceneisdirty = True

    if has_phil_path(diff_phil, "scene_bin_thresholds"):
      self.sceneisdirty = True

    if has_phil_path(diff_phil, "color_scheme"):
      self.add_colour_scheme_to_dict()
      self.sceneisdirty = True

    if has_phil_path(diff_phil, "color_powscale"):
      self.add_colour_powscale_to_dict()
      self.sceneisdirty = True

    if has_phil_path(diff_phil, "nth_power_scale_radii"):
      self.add_nth_power_scale_radii_to_dict()
      self.sceneisdirty = True

    if has_phil_path(diff_phil, "scale"):
      self.add_radii_scale_to_dict()
      self.sceneisdirty = True

    if has_phil_path(diff_phil, "camera_type"):
      self.set_camera_type()

    if has_phil_path(diff_phil, "show_hkl"):
      self.show_hkl()

    if has_phil_path(diff_phil, "tooltip_data"):
      tablerow, binclude = eval(self.params.tooltip_data)
      self.include_tooltip_lst[tablerow] = binclude

    if has_phil_path(diff_phil, "background_colour"):
      self.set_background_colour()

    if has_phil_path(diff_phil, "show_tooltips"):
      self.set_show_tooltips()

    if has_phil_path(diff_phil, "tooltip_alpha"):
      self.set_tooltip_opacity()

    if has_phil_path(diff_phil, "angle_around_vector"):
      i,deg = self.rotate_around_numbered_vector()
      self.params.viewer.angle_around_vector = str([i, deg])

    if has_phil_path(diff_phil, "angle_around_XHKL_vector"):
      self.rotate_stage_around_cartesian_vector([1,0,0], self.params.viewer.angle_around_XHKL_vector)
      self.params.viewer.angle_around_XHKL_vector = None

    if has_phil_path(diff_phil, "angle_around_YHKL_vector"):
      self.rotate_stage_around_cartesian_vector([0,1,0], self.params.viewer.angle_around_YHKL_vector)
      self.params.viewer.angle_around_YHKL_vector = None

    if has_phil_path(diff_phil, "angle_around_ZHKL_vector"):
      self.rotate_stage_around_cartesian_vector([0,0,1], self.params.viewer.angle_around_ZHKL_vector)
      self.params.viewer.angle_around_ZHKL_vector = None

    if has_phil_path(diff_phil,
                      "spacegroup_choice",
                      "use_provided_miller_arrays",
                      "binning",
                      "fontsize",
                      "vector_width",
                      "data_array",
                      "miller_array_operation",
                      "mouse_sensitivity",
                      "real_space_unit_cell_scale_fraction",
                      "reciprocal_unit_cell_scale_fraction",
                      "draw_real_space_unit_cell",
                      "draw_reciprocal_unit_cell",
                      "clip_plane",
                      "show_vector",
                      "show_all_vectors",
                      "hkls",
                      "use_wireframe",
                      "viewer") and self.params.viewer.scene_id is not None:
       # any change to parameters in the master phil in display2.py
      self.scene = self.HKLscene_dict_val(self.params.viewer.scene_id).scene
      self.DrawNGLJavaScript()
      self.mprint( "Rendered %d reflections" % self.scene.points.size(), verbose=1)
      #time.sleep(25) # for debugging
      self.realSpaceMag = (self.realspace_scale - 1.0)*self.params.real_space_unit_cell_scale_fraction + 1.0
      self.recipSpaceMag = (self.reciproc_scale - 1.0)*self.params.reciprocal_unit_cell_scale_fraction + 1.0

      if has_phil_path(diff_phil, "vector_width"):
        self.SetVectorWidth(self.params.NGL.vector_width)

      if has_phil_path(diff_phil, "hkldist"):
        self.visual_symHKLs = []

      if has_phil_path(diff_phil, "show_vector",
                                  "real_space_unit_cell_scale_fraction",
                                  "reciprocal_unit_cell_scale_fraction"):
        self.show_vectors(self.params.viewer.show_vector, diff_phil)

      if has_phil_path(diff_phil, "show_all_vectors"):
        self.show_all_vectors()

      if has_phil_path(diff_phil, "normal_vector") and self.params.clip_plane.normal_vector != "":
        found = False
        for (opnr, label, order, cartvec, hklop, hkl, abc, length) in self.all_vectors:
          if self.params.clip_plane.normal_vector in label:
            found=True
        if not found:
          raise Sorry("No vector present with substring: %s" %self.params.clip_plane.normal_vector)

      self.set_volatile_params()

    if has_phil_path(diff_phil, "fontsize"):
      self.SetFontSize(self.params.NGL.fontsize)

    if has_phil_path(diff_phil, "animate_rotation_around_vector"):
      i,speed = self.animate_rotate_around_vector()
      self.params.viewer.animate_rotation_around_vector = str([i, speed])

    if self.params.viewer.scene_id is None:
      self.DrawNGLJavaScript(blankscene=True)
    return curphilparam


  def set_volatile_params(self, use_semaphore=True):
    """
    Change the view of the reflections according to whatever the values are of the volatile parameters.
    Volatile parameters are those that do not require heavy computing (like position of WebGL primitives)
    but can change the appearance of primitives instantly like opacity or clipplane position. Expansion
    in browser of coordinates to P1 are also considered volatile as this operation is very fast.
    """
    self.use_semaphore = use_semaphore
    if self.params.viewer.scene_id is not None:
      if self.isnewfile:
        self.SetDefaultOrientation()
        if not self.autoview_sem.acquire(blocking=True, timeout=lock_timeout):
          self.mprint("Error! Timed out waiting for autoview_sem semaphore within %s seconds" %lock_timeout, verbose=1)
        self.mprint("set_volatile_params got autoview_sem", verbose="threadingmsg")
        self.autoview_sem.release()
        self.mprint("set_volatile_params released clipplane_msg_sem", verbose="threadingmsg")

        while len(self.WBmessenger.clientmsgqueue):
          self.mprint("set_volatile_params sleep", verbose=1)
          time.sleep(0.2)

      if self.params.viewer.fixorientation == "vector":
        self.orient_vector_to_screen(self.currentrotvec)
      self.SetMouseSpeed(self.params.NGL.mouse_sensitivity)
      hkldist = -1
      clipwidth = None
      self.fix_orientation()
      uc = self.miller_array.unit_cell()
      if self.params.clip_plane.clip_width: # then we are clipping
        self.UseCameraZoom()
        if self.params.clip_plane.auto_clip_width: # set the default spacing between layers of reflections
          self.params.clip_plane.clip_width = 0.5*self.L # equal to half the hkl vector length
        clipwidth = self.params.clip_plane.clip_width
        hkldist = -self.params.clip_plane.hkldist * self.L *self.cosine
        self.mprint("clip plane distance from origin: %s" %hkldist)
      else:
        self.UseZoomDrag()
      infomsg = ""
      self.normal_vecnr = -1
      for (opnr, label, order, cartvec, hklop, hkl, abc, length) in self.all_vectors:
        if self.params.clip_plane.normal_vector in label and self.params.clip_plane.normal_vector != "":
          self.normal_vecnr = opnr

      if self.normal_vecnr != -1: # then we are orienting clip plane with a vector
        # cartvec can be hklvec vector in cartesian coordinates
        # or abcvec vector in cartesian coordinates
        cartvec = self.all_vectors[ self.normal_vecnr ][3]
        self.L = self.all_vectors[ self.normal_vecnr ][7]
        # Use half the length of the tncs vector to allow stepping through alternating weak and strong layers
        # of reflections in the GUI when orienting clip plane perpendicular to the tncs vector
        if "tNCS" in self.all_vectors[ self.normal_vecnr ][1]:
          self.L *= 0.5

        if self.params.clip_plane.auto_clip_width: # set the default spacing between layers of reflections
          self.params.clip_plane.clip_width = 0.5*self.L # equal to half the hkl vector length
        clipwidth = self.params.clip_plane.clip_width
        # hklvec is reciprocal vector in reciprocal coordinates.
        # First try and see if they are stored in self.all_vectors[..][5].
        # If not then convert the cartesian representation cartvec of hklvec
        # into the reciprocal coordinates
        try:
          hklvec = eval(self.all_vectors[ self.normal_vecnr ][5])
        except Exception as e:
          hklvec = list(self.reciprocal_from_real_space_vector(cartvec ))
        # Get corresponding real space vector to the hkl vector (as cartesian coordinates)
        real_space_vec = hklvec * matrix.sqr(uc.orthogonalization_matrix())
        # In the general case real_space_vec is not parallel to hklvec
        # Orient the clip plane perpendicular to real_space_vec while at the
        # same time slide clip plane along the cartvec (reciprocal vector) direction
        # in units of cartvec projected onto real_space_vec
        self.cosine, _, _ = self.project_vector1_vector2(cartvec, real_space_vec)
        hkldist = -self.params.clip_plane.hkldist * self.L *self.cosine
        self.mprint("clip plane distance from origin: %s" %hkldist)
        if self.params.clip_plane.is_assoc_real_space_vector:
          orientvector = real_space_vec
          self.mprint("clip plane perpendicular to realspace vector associated with hkl vector: %s" %str(hklvec), verbose=1)
        else:
          orientvector = cartvec
          abcvec = self.all_vectors[ self.normal_vecnr ][6]
          self.mprint("clip plane perpendicular to realspace vector: %s" %str(abcvec), verbose=1)
          infomsg = "Vector distance from origin: %d" %(self.params.clip_plane.hkldist)
          if "tNCS" in self.all_vectors[ self.normal_vecnr ][1]:
            """ Clip plane width for tncs should be around 1/4 of the tncs modulation length
            as to ensure we only get the strongest/weakest reflections between clipnear, clipfar.
            The tncs modulation length is the inverse length of the tncs vector as defined in
            HKLViewFrame.list_vectors() where the length is stored as half the length of the tncs vector
            for the sake of stepping through alternating weak and strong layers with the +/- buttons.
            So set clip plane width to 0.5*0.5/tncs-vector-length.
            """
            if self.params.clip_plane.auto_clip_width: # use the default spacing between tncs layers of reflections
              self.params.clip_plane.clip_width = 0.5*self.L
            clipwidth = self.params.clip_plane.clip_width
            # Want the radius of the sphere of reflections so get a reflection at highest resolution
            # Decrease resolution by 0.00001 to avoid machine precision errors yielding 0 reflections
            dminhkl = self.miller_array.resolution_filter(d_min=0,
                                                          d_max=(self.miller_array.d_min()+0.00001)).indices()[0]
            dmincartvec = list( dminhkl * matrix.sqr(uc.fractionalization_matrix()).transpose() )
            sphereradius = math.sqrt(dmincartvec[0]*dmincartvec[0] + dmincartvec[1]*dmincartvec[1]
                                      + dmincartvec[2]*dmincartvec[2] )
            n_tncs_layers = sphereradius*self.renderscale/self.L
            infomsg = "tNCS layer: %d out of +-%2.2f" %(self.params.clip_plane.hkldist, n_tncs_layers)

          if "-fold" in self.all_vectors[ self.normal_vecnr ][1]:
            clipwidth = self.params.clip_plane.clip_width
            infomsg = "Reflection layer %d related through %s" \
              %(self.params.clip_plane.hkldist, self.all_vectors[ self.normal_vecnr ][1])

        self.orient_vector_to_screen(orientvector)
        scalefactor = 1.0
        if self.params.clip_plane.normal_vector_length_scale > 0 and self.all_vectors[self.normal_vecnr][1] != "tNCS":
          scalefactor = self.L/self.params.clip_plane.normal_vector_length_scale
          self.L = self.params.clip_plane.normal_vector_length_scale
        # Make a string of the equation of the plane of reflections
        hklvecsqr = hklvec[0]*hklvec[0] + hklvec[1]*hklvec[1] + hklvec[2]*hklvec[2]
        if self.params.clip_plane.is_assoc_real_space_vector:
          self.planescalarvalue = self.params.clip_plane.hkldist * hklvecsqr*scalefactor
          self.planenormalhklvec = hklvec
          infomsg = "Reflections satisfying: %s*h + %s*k + %s*l = %s" \
            %(roundoff(hklvec[0],4), roundoff(hklvec[1],4), roundoff(hklvec[2],4), roundoff(self.planescalarvalue))
        self.cosine, _, _ = self.project_vector1_vector2(cartvec, real_space_vec)
      # show equation or info in the browser
      self.AddToBrowserMsgQueue("PrintInformation", infomsg)
      self.ExpandInBrowser()
      retstr = ""
      if self.miller_array and self.params.binning.bin_opacity:
        bin_opacitieslst = self.params.binning.bin_opacity
        for alpha,bin in bin_opacitieslst:
          retstr += self.set_opacity(bin, alpha)
        self.SendInfoToGUI( { "bin_opacity": self.params.binning.bin_opacity } )
      self.mprint( retstr, verbose=1)
      self.DrawUnitCell()
      self.DrawReciprocalUnitCell()
      self.set_tooltip_opacity()
      self.set_show_tooltips()
      self.visualise_sym_HKLs()
      self.isnewfile = False
      self.make_clip_plane(hkldist, clipwidth)


  def set_scene(self):
    self.binvals = []
    if self.params.viewer.scene_id is None:
      return False
    self.set_miller_array(self.params.viewer.scene_id)
    if (self.miller_array is None):
      raise Sorry("No data loaded!")
    self.mprint( "Miller array %s runs from hkls: %s to %s" \
     %(self.miller_array.info().label_string(), self.miller_array.index_span().min(),
        self.miller_array.index_span().max() ) )
    self.mprint("Spacegroup: %s" %self.miller_array.space_group().info().symbol_and_number())
    return True


  def set_miller_array(self, scene_id=None, merge=None, details=""):
    if scene_id is not None:
      self.params.viewer.scene_id = scene_id
    if self.params.hkls and self.params.viewer.scene_id is not None and self.params.viewer.scene_id >= 0:
      self.miller_array = self.HKLscene_dict_val().scene.miller_array
      self.scene = self.HKLscene_dict_val().scene
    self.merge = merge
    if (self.miller_array is None):
      return
    self.inspect_arrays()
    self.GetUnitcellScales()
    self.d_min = self.miller_array.d_min()
    array_info = self.miller_array.info()
    # capture the currently selected spacegroup if not the default
    self.sg = self.proc_arrays[self.scene_id_to_array_id(self.params.viewer.scene_id)].space_group()
    self.symops = list(self.sg.all_ops())
    uc = "a=%g b=%g c=%g angles=%g,%g,%g" % self.miller_array.unit_cell().parameters()
    self.mprint( "Data: %s %s, %d reflections in space group: %s, unit Cell: %s" \
      % (array_info.label_string(), details, self.miller_array.indices().size(), \
          self.miller_array.space_group_info(), uc), verbose=0 )


  def Complex2AmplitudesPhases(self, data):
    ampls = flex.abs(data)
    phases = flex.arg(data) * 180.0/math.pi
    # purge nan values from array to avoid crash in fmod_positive()
    #b = flex.bool([bool(math.isnan(e)) for e in phases])
    # replace the nan values with an arbitrary float value
    #phases = phases.set_selected(b, 42.4242)
    phases = graphics_utils.NoNansArray( phases, 42.4242)
    # Cast negative degrees to equivalent positive degrees
    phases = flex.fmod_positive(phases, 360.0)
    return ampls, phases


  def get_rothkl_from_IDs(self, id, sym_id):
    anomalous = False
    if id < 0:
      anomalous = True
      id = abs(id) # we set anomalous ids as the negative of the original ids
    hkl = self.scene.indices[id]
    hklvec = flex.vec3_double( [(hkl[0], hkl[1], hkl[2])])
    rotmx=None
    if sym_id >= 0 and sym_id < len(self.symops):
      # symid tells which symmetry operator was used in HKLJavaScripts.js onMessage() Expand()
      rotmx = self.symops[sym_id].r()
    Rhkl = hklvec[0]
    if rotmx:
      # if a symmetry mate was clicked then deduce its hkl coordinate by
      # applying the rotation to the original hkl coordinate
      Rhkl = hklvec[0] * rotmx
    rothkl = (int(Rhkl[0]), int(Rhkl[1]), int(Rhkl[2]))
    if anomalous:
      rothkl = (int(-Rhkl[0]), int(-Rhkl[1]), int(-Rhkl[2]))
    return rothkl, hkl


  def make_visual_symHKLs(self, id, sym_id):
    symid = sym_id
    rothkl, dummy = self.get_rothkl_from_IDs(id, symid) # and use it
    if self.visual_symmxs:
      # if a list of symmetry matrices have been deduced from a selected rotation operator
      # then also compute the other symmetry mates of the current hkl
      self.visual_symHKLs = []
      for symmx,hklstr in self.visual_symmxs:
        vissymrothkl = rothkl * symmx.transpose()
        self.visual_symHKLs.append( (vissymrothkl, hklstr) )


  def GetTooltipOnTheFly(self, id, sym_id, anomalous=False):
    rothkl, hkl = self.get_rothkl_from_IDs(id, sym_id)
    spbufttip = '\'HKL: [%d,%d,%d]' %(rothkl[0], rothkl[1], rothkl[2])
    if rothkl != hkl: # then show the original hkl before P1 or anomalous expansion
      spbufttip += ', (asu): [%d,%d,%d]' %(hkl[0], hkl[1], hkl[2])
    # resolution and Angstrom character for javascript
    spbufttip += '\\ndres: %s \'+ String.fromCharCode(197) +\'' \
      %str(roundoff(self.miller_array.unit_cell().d(hkl), 2) )
    for tablerow,proc_array in enumerate(self.proc_arrays):
      if not self.include_tooltip_lst[tablerow]:
        continue
      sigvals = []
      datvals = []
      if proc_array.sigmas() is not None:
        sigvals = list( proc_array.select(proc_array.indices() == hkl).sigmas() )
      datval = None
      if hkl in proc_array.indices():
        datvals = list( proc_array.select(proc_array.indices() == hkl).data() )
      else:
        if id >= proc_array.size():
          continue
      for i,datval in enumerate(datvals):
        if proc_array.is_hendrickson_lattman_array() and math.isnan(datval[0] + datval[1] + datval[2] + datval[3]):
          continue
        if not isinstance(datval, tuple) and (math.isnan( abs(datval) ) or datval == display.inanval):
          continue
        if proc_array.is_complex_array():
          ampl = abs(datval)
          phase = cmath.phase(datval) * 180.0/math.pi
          # purge nan values from array to avoid crash in fmod_positive()
          # and replace the nan values with an arbitrary float value
          if math.isnan(phase):
            phase = 42.4242
          # Cast negative degrees to equivalent positive degrees
          phase = phase % 360.0
        spbufttip +="\\n" + proc_array.info().label_string() + ': '
        if proc_array.is_complex_array():
          spbufttip += str(roundoff(ampl, 2)) + ", " + str(roundoff(phase, 2)) + \
            "\'+ String.fromCharCode(176) +\'" # degree character for javascript
        elif sigvals:
          sigma = sigvals[i]
          spbufttip += str(roundoff(datval, 2)) + ", " + str(roundoff(sigma, 2))
        else:
          spbufttip += str(roundoff(datval, 2))
    spbufttip += '\\n\\n%d,%d' %(id, sym_id) # compared by the javascript
    spbufttip += '\''
    return spbufttip


  def get_col_fomcol(self, idx):
    if len(self.HKLscene_dict_val().arrayinfo) == 0:
      return -1, -1
    return self.HKLscene_dict_val(idx).arrayinfo[6], self.HKLscene_dict_val(idx).arrayinfo[7]


  def ConstructReciprocalSpace(self, scene_id=None):
    if len(self.proc_arrays) == 0 or scene_id is None:
      return False
    sceneid = scene_id
    self.HKLsceneKey = self.Sceneid_to_SceneKey(sceneid)
    if self.HKLsceneKey in self.HKLscenedict and not self.has_new_miller_array:
      self.HKLscene = self.HKLscenedict.get(self.HKLsceneKey, False)
      if self.HKLscene:
        self.mprint("Using cached HKL scene", verbose=1)
        return True
    if self.has_new_miller_array:
      self.inspect_arrays()
    self.mprint("Constructing HKL scenes...", verbose=1)
    idx,fdx = self.scene_id_to_array_and_foms_id(scene_id)
    fomarray = None
    if fdx >= 0: # index is -1 if idx is not paired with a FOM array
      fomarray = self.proc_arrays[fdx].deep_copy()
    (hklscenes, scenemaxdata,
      scenemindata, scenemaxsigmas,
        sceneminsigmas, scenearrayinfos
    ) = MakeHKLscene( self.proc_arrays[idx].deep_copy(), fomarray, idx, fdx,
                     self.renderscale, copy.deepcopy(self.params.hkls), self.mprint )
    for i,inf in enumerate(scenearrayinfos):
      self.mprint("%d, %s" %(idx+i+1, inf[3]), verbose=1)
      self.HKLsceneKey = self.Sceneid_to_SceneKey(sceneid)
      self.HKLscenedict[self.HKLsceneKey] =  group_args(scene = hklscenes[i],
                                                        maxdata = scenemaxdata[i],
                                                        mindata = scenemindata[i],
                                                        maxsigma = scenemaxsigmas[i],
                                                        minsigma = sceneminsigmas[i],
                                                        arrayinfo = inf)
      sceneid += 1
    self.HKLscene = self.HKLscenedict[self.HKLsceneKey].scene
    self.HKLscenesMaxdata = self.HKLscenedict[self.HKLsceneKey].maxdata
    self.HKLscenesMindata = self.HKLscenedict[self.HKLsceneKey].mindata
    self.HKLscenesMaxsigmas = self.HKLscenedict[self.HKLsceneKey].maxsigma
    self.HKLscenesMinsigmas = self.HKLscenedict[self.HKLsceneKey].minsigma
    self.hkl_scenes_info = self.HKLscenedict[self.HKLsceneKey].arrayinfo
    self.sceneisdirty = True
    self.has_new_miller_array = False
    self.mprint("Done constructing HKL scenes", verbose=1)
    return True


  def Sceneid_to_SceneKey(self, sceneid):
    return (self.params.spacegroup_choice,
                      self.params.using_space_subgroup,
                      self.params.hkls.expand_anomalous,
                      self.params.hkls.expand_to_p1,
                      self.params.hkls.slice_axis,
                      self.params.hkls.slice_index,
                      self.params.hkls.show_missing,
                      self.params.hkls.show_only_missing,
                      self.params.hkls.show_systematic_absences,
                      self.params.hkls.sigma_color_radius,
                      self.params.hkls.color_scheme,
                      self.params.hkls.color_powscale,
                      sceneid,
                      self.params.hkls.scale,
                      self.params.hkls.nth_power_scale_radii
                      )


  def HKLscene_dict_val(self, sceneid=None):
    if sceneid is None:
      sceneid = self.params.viewer.scene_id
    HKLsceneKey = self.Sceneid_to_SceneKey(sceneid)
    if not self.HKLscenedict.get(HKLsceneKey, False):
      self.ConstructReciprocalSpace(scene_id=sceneid)
    return self.HKLscenedict[HKLsceneKey]


  def inspect_arrays(self):
    self.mprint("Matching complex arrays to suitable FOM arrays", verbose=1)
    self.mapcoef_fom_dict = {}
    self.sceneid_from_arrayid = []
    self.has_unmerged_data = False

    for k,proc_array in enumerate(self.proc_arrays):
      fom_arrays_idx = []
      array_scene_ids = [(k,-1)] # using -1 to indicate not paired with a fom array
      for i,foms_array in enumerate(self.proc_arrays):
        if not proc_array.is_complex_array() or not foms_array.is_real_array():
          continue
        if proc_array.size() != foms_array.size():
          continue
        if  min(foms_array.data()) < 0.0 or flex.max(foms_array.data()) > 1.0:
          continue
        fom_arrays_idx.append( (foms_array, i) )
        array_scene_ids.append((k,i))
      self.sceneid_from_arrayid.extend( array_scene_ids)
      self.mapcoef_fom_dict[proc_array.info().label_string()] = fom_arrays_idx
      if proc_array.is_unmerged_intensity_array():
        self.has_unmerged_data = True


  def get_scene_id_from_label_or_type(self, datalabel, datatype=None):
    """ Try finding a matching sceneid to the datalabel provided. As a fallback
    try finding a sceneid for the first matching datatype regardless of its label
    """
    assert datalabel is not None
    for i,e in enumerate(self.hkl_scenes_infos):
      if e[3] == datalabel:
        return i
    if datatype is not None:
      for i,e in enumerate(self.hkl_scenes_infos):
        if e[4] == datatype:
          return i
    self.mprint("Currently no dataset with label %s or type %s" %(datalabel, datatype), verbose=1)
    return None


  def get_label_type_from_scene_id(self, sceneid):
    if sceneid is None:
      return None,None
    # Find data label and type for a particular sceneid
    assert sceneid < len(self.hkl_scenes_infos)
    datalabel = self.hkl_scenes_infos[sceneid][3]
    datatype = self.hkl_scenes_infos[sceneid][4]
    return datalabel, datatype


  def get_binner_idx_from_label(self, binlabel):
    for i,e in enumerate(self.bin_labels_type_idxs):
      if binlabel == e[0]:
        return i
    return -1


  def get_binlabel_from_binner_idx(self, idx):
    return self.bin_labels_type_idxs[idx][0]


  def scene_id_to_array_id(self, scene_id):
    for i,array_scene_id in enumerate(self.sceneid_from_arrayid):
      if scene_id == i:
        return array_scene_id[0]
    raise Sorry("scene_id, %d, is out of range" %scene_id)


  def scene_id_to_array_and_foms_id(self, scene_id):
    for i,array_scene_id in enumerate(self.sceneid_from_arrayid):
      if scene_id == i:
        return array_scene_id
    raise Sorry("scene_id, %s, is out of range" %scene_id)


  def calc_bin_thresholds(self, binner_idx, nbins):
    # make default bin thresholds if scene_bin_thresholds is not set
    binscenelabel = self.bin_labels_type_idxs[binner_idx][0]
    self.mprint("Using %s for binning" %binscenelabel)
    if binscenelabel=="Resolution":
      warray = self.HKLscene_dict_val().scene.work_array
      dres = self.HKLscene_dict_val().scene.dres
      uc = warray.unit_cell()
      indices = self.HKLscene_dict_val().scene.indices
      dmax,dmin = warray.d_max_min(d_max_is_highest_defined_if_infinite=True) # to avoid any F000 reflection
      if dmax == dmin: # say if only one reflection
        binvals = [dres[0]-0.1, dmin +0.1]
        nuniquevalues = 2
      else: # use generic binning function from cctbx
        binning = miller.binning( uc, nbins, indices, dmax, dmin )
        binvals = [ binning.bin_d_range(n)[0] for n in binning.range_all() ]
        binvals = [ e for e in binvals if e != -1.0] # delete dummy limit
        binvals = list( 1.0/flex.double(binvals) )
        nuniquevalues = len(set(list(dres)))
    elif "Singletons" in binscenelabel:
      binvals = [ -0.1, 0.1 ] # if this dataset hasn't got any singletons
      if len(set(list(self.scene.singletonsiness))) == 3: # symmetry unique anomalous data with some singletons
        binvals = [ -1.5, -0.5, 0.5, 1.5 ]
      nuniquevalues = len(binvals)
    else:
      bindata, dummy = self.get_matched_binarray(binner_idx)
      selection = flex.sort_permutation( bindata )
      bindata_sorted = bindata.select(selection)
      # First check for case where all unique values could be covered by
      # number of requested bins (e.g. multiplicity values)
      uniquevalues = list((set(list(bindata))))
      nuniquevalues = len(uniquevalues)
      if nuniquevalues <= nbins:
        uniquevalues.sort()
        binvals = [uniquevalues[0]-1]
        for ival in range(nuniquevalues):
          binvals.append(uniquevalues[ival])
        nuniquevalues = len(binvals)

      else:
        # Get binvals by dividing bindata_sorted with nbins
        # This yields approximately the same number of reflections in each bin
        binvals = [ bindata_sorted[0] ]
        nbins_used = 0
        float_data_used = 0.0
        num_per_bin = float(len(bindata_sorted))/nbins
        while nbins_used < nbins:
          index = round(float_data_used + num_per_bin) - 1
          threshold = bindata_sorted[index]
          # Handle case where there are a lot of repeated values
          float_data_used += num_per_bin
          if threshold > binvals[-1]:
            binvals.append(threshold)
            nbins_used += 1
          else:
            # Split remaining data over remaining requested bins
            num_per_bin = (float(len(bindata_sorted))-float_data_used)/(nbins-nbins_used)
        if bindata_sorted[-1] > binvals[-1]:
          binvals.append(bindata_sorted[-1])
        nuniquevalues = len(binvals)

    binvals.sort()
    self.mprint("Bin thresholds are:\n" + str(binvals), verbose=1)
    return binvals, nuniquevalues


  def UpdateBinValues(self, binner_idx, binvals = None, nuniquevalues = -1):
    if binvals:
      binvals.sort()
      self.binvals = binvals
    else: # ensure default resolution interval includes all data by avoiding rounding errors
      self.binvals = [ 1.0/(self.miller_array.d_max_min()[0]*1.001),
                       1.0/(self.miller_array.d_max_min()[1]*0.999) ]
    if nuniquevalues == -1:
      if binner_idx==0: # i.e. the resolution array of the hkls
        nuniquevalues = len(set(list( self.HKLscene_dict_val().scene.dres )))
      elif binner_idx==1:  # i.e. singletons
        binvals = [ -0.1, 0.1 ] # if this dataset hasn't got any singletons
        if len(set(list(self.scene.singletonsiness))) == 3: # symmetry unique anomalous data with some singletons
          binvals = [ -1.5, -0.5, 0.5, 1.5 ]
        nuniquevalues = len(binvals)
      else: # one of the normal datasets
        bindata, dummy = self.get_matched_binarray(binner_idx)
        nuniquevalues = len(set(list(bindata)))
    self.nuniqueval = nuniquevalues


  def get_matched_binarray(self, binner_idx):
    binscenelabel, datatype, sceneid = self.bin_labels_type_idxs[binner_idx]
    label = self.HKLscene_dict_val(sceneid).scene.work_array.info().label_string()
    if datatype == "hassigmas" and binscenelabel == "Sigmas of " + label:
      bindata = self.HKLscene_dict_val(sceneid).scene.sigmas.deep_copy()
      binvalsboundaries = [ self.HKLscene_dict_val(sceneid).minsigma - 0.1 , self.HKLscene_dict_val(sceneid).maxsigma + 0.1 ]
    elif datatype in "Map coeffs" and "Phases of " + label in binscenelabel:
      bindata = self.HKLscene_dict_val(sceneid).scene.phases.deep_copy()
      # preselect centric reflections, i.e. those with phi = 0 or 180
      binvalsboundaries = [-0.01, 0.01, 179.99, 180.01, 359.99, 360]
    elif datatype in "Map coeffs" and "Amplitudes of " + label in binscenelabel:
      bindata = self.HKLscene_dict_val(sceneid).scene.ampl.deep_copy()
      binvalsboundaries = [ self.HKLscene_dict_val(sceneid).mindata - 0.1 , self.HKLscene_dict_val(sceneid).maxdata + 0.1 ]
    else:
      bindata = self.HKLscene_dict_val(sceneid).scene.data.deep_copy()
      binvalsboundaries = [ self.HKLscene_dict_val(sceneid).mindata - 0.1 , self.HKLscene_dict_val(sceneid).maxdata + 0.1 ]
    return bindata, binvalsboundaries


  def MatchBinArrayToSceneArray(self):
    # match bindata with data or sigmas
    if self.bin_labels_type_idxs[self.params.binning.binner_idx][0] == "Resolution":
      return 1.0/self.scene.dres
    binarraydata, dummy = self.get_matched_binarray(self.params.binning.binner_idx)
    scenearraydata = self.HKLscene_dict_val().scene.data
    binlabel, _, ibinarray = self.bin_labels_type_idxs[self.params.binning.binner_idx]
    if len(set(self.HKLscene_dict_val(ibinarray).scene.indices)) < self.HKLscene_dict_val(ibinarray).scene.indices.size():
      raise Sorry("Error: The HKL indices in %s are not unique. Use a merged dataset instead!" %binlabel)
    matchindices = miller.match_multi_indices(self.HKLscene_dict_val(ibinarray).scene.indices,
                               self.HKLscene_dict_val().scene.indices )
    matched_binarray = binarraydata.select( matchindices.pairs().column(0) )
    # patch the bin array so its sequence matches the scene array
    patched_binarraydata = []
    c = 0
    for b in matchindices.pair_selection(1):
      if b:
        patched_binarraydata.append(matched_binarray[c])
        c +=1
      else:
        patched_binarraydata.append(float("nan"))
    return flex.double(patched_binarraydata)


  def OperateOn1MillerArray(self, millarr, operation):
    # lets user specify a python expression operating on millarr
    newarray = millarr.deep_copy()
    dres = newarray.unit_cell().d( newarray.indices() )
    try:
      ldic= {'dres': dres, 'array1': newarray, 'newarray': newarray }
      exec(operation, globals(), ldic)
      newarray = ldic.get("newarray", None)
      return newarray
    except Exception as e:
      raise Sorry(str(e))


  def OperateOn2MillerArrays(self, millarr1, millarr2, operation):
    # lets user specify a python expression operating on millarr1 and millarr2
    matchindices = miller.match_indices(millarr1.indices(), millarr2.indices() )
    matcharr1 = millarr1.select( matchindices.pairs().column(0) ).deep_copy()
    matcharr2 = millarr2.select( matchindices.pairs().column(1) ).deep_copy()
    dres = matcharr1.unit_cell().d( matcharr1.indices() )
    matcharr2._observation_type = millarr1._observation_type
    newarray = matcharr2.deep_copy()
    try:
      ldic= { 'dres': dres, 'array1': matcharr1, 'array2': matcharr2, 'newarray': newarray }
      exec(operation, globals(), ldic)
      newarray = ldic.get("newarray", None)
      return newarray
    except Exception as e:
      raise Sorry(str(e))


  def get_colour_map_radii_power(self):
    datatype = self.get_current_datatype_or_default_dict()
    colourscheme, colourpower, powerscale, radiiscale = \
        self.datatypedict.get( datatype, self.datatypedefault[:] )
    return colourscheme, colourpower, powerscale, radiiscale


  def add_colour_scheme_to_dict(self):
    datatype = self.get_current_datatype_or_default_dict()
    self.datatypedict[datatype][0] = self.params.hkls.color_scheme


  def add_colour_powscale_to_dict(self):
    datatype = self.get_current_datatype_or_default_dict()
    self.datatypedict[datatype][1] = self.params.hkls.color_powscale


  def add_nth_power_scale_radii_to_dict(self):
    datatype = self.get_current_datatype_or_default_dict()
    self.datatypedict[datatype][2] = self.params.hkls.nth_power_scale_radii


  def add_radii_scale_to_dict(self):
    datatype = self.get_current_datatype_or_default_dict()
    self.datatypedict[datatype][3] = self.params.hkls.scale


  def get_current_datatype_or_default_dict(self):
    datatype = self.get_current_datatype()
    if datatype is None:
      return
    if self.params.hkls.sigma_color_radius:
      datatype = datatype + "_sigmas"
    if datatype not in self.datatypedict.keys():
        # ensure individual copies of datatypedefault and not references to the same
      self.datatypedict[ datatype ] = self.datatypedefault[:]
    return datatype


  def DrawNGLJavaScript(self, blankscene=False):
    if not self.scene or not self.sceneisdirty:
      return
    if self.scene.points.size() == 0:
      blankscene = True
    if self.miller_array is None :
      self.mprint( "Select a dataset to display reflections" )
      blankscene = True
    else:
      self.mprint("Rendering reflections.", end="")

    if self.scene is not None and self.miller_array is not None: # expansion always done in the browser now
     assert (self.scene.settings.expand_to_p1==False and self.scene.settings.expand_anomalous==False)

    h_axis = flex.vec3_double([self.scene.axes[0]])
    k_axis = flex.vec3_double([self.scene.axes[1]])
    l_axis = flex.vec3_double([self.scene.axes[2]])
    self.unit_h_axis = 1.0/h_axis.norm() * h_axis
    self.unit_k_axis = 1.0/k_axis.norm() * k_axis
    self.unit_l_axis = 1.0/l_axis.norm() * l_axis
    self.unit_normal_hk = self.unit_h_axis.cross( self.unit_k_axis )
    self.unit_normal_kl = self.unit_k_axis.cross( self.unit_l_axis )
    self.unit_normal_lh = self.unit_l_axis.cross( self.unit_h_axis )
    self.normal_hk = h_axis.cross( k_axis )
    self.normal_kl = k_axis.cross( l_axis )
    self.normal_lh = l_axis.cross( h_axis )
    maxnorm = max(h_axis.norm(), max(k_axis.norm(), l_axis.norm()))
    l1 = self.renderscale * maxnorm * 1.1
    l2= self.renderscale * maxnorm * 1.15
    Hstararrowstart = roundoff( [-self.unit_h_axis[0][0]*l1, -self.unit_h_axis[0][1]*l1, -self.unit_h_axis[0][2]*l1] )
    Hstararrowend = roundoff( [self.unit_h_axis[0][0]*l1, self.unit_h_axis[0][1]*l1, self.unit_h_axis[0][2]*l1] )
    Hstararrowtxt  = roundoff( [self.unit_h_axis[0][0]*l2, self.unit_h_axis[0][1]*l2, self.unit_h_axis[0][2]*l2] )
    Kstararrowstart = roundoff( [-self.unit_k_axis[0][0]*l1, -self.unit_k_axis[0][1]*l1, -self.unit_k_axis[0][2]*l1] )
    Kstararrowend = roundoff( [self.unit_k_axis[0][0]*l1, self.unit_k_axis[0][1]*l1, self.unit_k_axis[0][2]*l1] )
    Kstararrowtxt  = roundoff( [self.unit_k_axis[0][0]*l2, self.unit_k_axis[0][1]*l2, self.unit_k_axis[0][2]*l2] )
    Lstararrowstart = roundoff( [-self.unit_l_axis[0][0]*l1, -self.unit_l_axis[0][1]*l1, -self.unit_l_axis[0][2]*l1] )
    Lstararrowend = roundoff( [self.unit_l_axis[0][0]*l1, self.unit_l_axis[0][1]*l1, self.unit_l_axis[0][2]*l1] )
    Lstararrowtxt  = roundoff( [self.unit_l_axis[0][0]*l2, self.unit_l_axis[0][1]*l2, self.unit_l_axis[0][2]*l2] )

    if not blankscene:
      self.params.hkls.color_scheme, self.params.hkls.color_powscale, self.params.hkls.nth_power_scale_radii, \
        self.params.hkls.scale = self.get_colour_map_radii_power()

      # Make colour gradient array used for drawing a bar of colours next to associated values on the rendered html
      mincolourscalar = self.HKLscene_dict_val().mindata
      maxcolourscalar = self.HKLscene_dict_val().maxdata
      if self.params.hkls.sigma_color_radius:
        mincolourscalar = self.HKLscene_dict_val().minsigma
        maxcolourscalar = self.HKLscene_dict_val().maxsigma
      span = maxcolourscalar - mincolourscalar
      ln = 120
      incr = span/ln
      colourgradarrays = []
      val = mincolourscalar
      colourscalararray = flex.double()
      colourscalararray.append( val )
      for j,sc in enumerate(range(ln)):
        val += incr
        colourscalararray.append( val )
      if self.HKLscene_dict_val().scene.miller_array.is_complex_array():
        # When displaying phases from map coefficients together with fom values
        # compute colour map chart as a function of fom and phase values (x,y axis)
        incr = 360.0/ln
        val = 0.0
        colourscalararray = flex.double()
        colourscalararray.append( val )
        for j in enumerate(range(ln)):
          val += incr
          colourscalararray.append( val )

        fomarrays = []
        COL = display.MplColorHelper(self.params.hkls.color_scheme, 0, 360)
        rgbcolarray = flex.vec3_double( [ COL.get_rgb(d)[0:3] for d in colourscalararray ] )

        if self.HKLscene_dict_val().scene.isUsingFOMs():
          fomln = 50
          fom = 1.0
          fomdecr = 1.0/(fomln-1.0)
          for j in range(fomln):
            fomarrays.append( flex.double(len(colourscalararray), fom) )
            fom -= fomdecr
          for j in range(fomln):
            arr = graphics_utils.map_to_rgb_colourmap(
                data_for_colors= colourscalararray,
                colormap= rgbcolarray,
                selection=flex.bool(colourscalararray.size(), True),
                attenuation = fomarrays[j]
              )
            colourgradarrays.append( arr*256 )
        else:
          fomln =1
          fomarrays = [1.0]
          arr = graphics_utils.map_to_rgb_colourmap(
              data_for_colors= colourscalararray,
              colormap = rgbcolarray,
              selection=flex.bool(colourscalararray.size(), True),
            )
          colourgradarrays.append(  arr*256 )
      else:
        fomln = 1
        fomarrays = [1.0]
        COL = display.MplColorHelper(self.params.hkls.color_scheme, mincolourscalar, maxcolourscalar)
        rgbcolarray = flex.vec3_double( [ COL.get_rgb(d)[0:3] for d in colourscalararray ])

        arr = graphics_utils.map_to_rgb_colourmap(
            data_for_colors= colourscalararray,
            colormap = rgbcolarray,
            selection=flex.bool(colourscalararray.size(), True),
            powscale = self.params.hkls.color_powscale
          )

        colourgradarrays.append(arr*256)
      colors = self.HKLscene_dict_val().scene.colors
      radii = self.HKLscene_dict_val().scene.radii
      self.meanradius = flex.mean(radii)

    bin_labels_type_idx = self.bin_labels_type_idxs[self.params.binning.binner_idx]
    self.mprint(".", end="")
    if blankscene:
      points = flex.vec3_double( [ ] )
      colors = flex.vec3_double( [ ] )
      radii = flex.double( [ ] )
      bin_labels_type_idx = self.bin_labels_type_idxs[0]
    else:
      points = self.scene.points

    nrefls = points.size()
    hkls = self.scene.indices
    dres = self.scene.dres
    if bin_labels_type_idx[0] =="Resolution":
      colstr = "dres"
    elif "Singletons" in bin_labels_type_idx[0]:
      colstr = "Singleton"
    else:
      if not blankscene:
        colstr = self.HKLscene_dict_val(bin_labels_type_idx[2]).scene.work_array.info().label_string()
    data = self.scene.data
    if not blankscene:
      colourlabel = self.HKLscene_dict_val().scene.colourlabel
      fomlabel = self.HKLscene_dict_val().scene.fomlabel
    #import code, traceback; code.interact(local=locals(), banner="".join( traceback.format_stack(limit=10) ) )
    assert (colors.size() == radii.size() == nrefls)
    self.colours = []
    self.positions = []
    self.radii2 = []
    self.spbufttips = []

    if not blankscene:
      if bin_labels_type_idx[0] =="Resolution":
        self.binvalsboundaries = self.binvals
        self.bindata = 1.0/self.scene.dres
      elif "Singletons" in bin_labels_type_idx[0]:
        self.binvalsboundaries = self.binvals
        self.bindata = self.scene.singletonsiness
      else:
        # get upper and lower bounds for the dataset used for binning
        self.bindata = self.MatchBinArrayToSceneArray()
        if ( len(self.binvalsboundaries)==0 or len(self.params.binning.scene_bin_thresholds) > 0):
          dummy, self.binvalsboundaries = self.get_matched_binarray(self.params.binning.binner_idx)
        # binvals derived from scene_bin_thresholds must be sorted
        # if minimum or maximum of binvals are smaller or bigger than lower or
        # upper bounds then use those values instead
        self.binvals.sort()
        vals = self.binvals[:]
        # ignoring nan values add binvalsboundaries if these are smaller or bigger than values in binvals
        nonanbinvals = [e for e in self.binvals if not math.isnan(e)]
        if nonanbinvals[0] > self.binvalsboundaries[0]:
          vals[0] = self.binvalsboundaries[0]
        if nonanbinvals[-1] < self.binvalsboundaries[1]:
          vals[-1] = self.binvalsboundaries[-1]
        # if nan values are present then sort with nan being the last value
        vals = list(set( vals)) # no duplicates
        self.binvalsboundaries = sorted(vals, key= lambda e: sys.maxsize if math.isnan(e) else e)

    self.nbinvalsboundaries = len(self.binvalsboundaries)
    # avoid resetting opacities of bins unless we change the number of bins
    if self.oldnbinvalsboundaries != self.nbinvalsboundaries and not self.executing_preset_btn:
      self.params.binning.bin_opacity = [ [1.0, e] for e in range(self.nbinvalsboundaries ) ]
    self.oldnbinvalsboundaries = self.nbinvalsboundaries
    # Un-binnable data are scene data values where there are no matching reflections in the bin data
    # Put these in a separate bin and be diligent with the book keeping!
    for ibin in range(self.nbinvalsboundaries+1): # adding the extra bin for un-binnable data
      self.colours.append([]) # colours and positions are 3 x size of data()
      self.positions.append([])
      self.radii2.append([])
      self.spbufttips.append([])

    def data2bin(d, binvalsboundaries, nbinvalsboundaries):
      for ibin, binval in enumerate(binvalsboundaries):
        if math.isnan(d): # NaN values are un-binnable. Tag them for an additional last bin
          return nbinvalsboundaries-1
        if (ibin+1) == nbinvalsboundaries:
          return ibin
        if d > binval and d <= binvalsboundaries[ibin+1]:
          return ibin
      raise Sorry("data2bin: Should never get here")

    def getprecision(v1,v2):
      diff = abs(v1-v2); precision = 1; e = 1
      while diff*e < 1.0 and diff > 0:
        e *= 10
        precision += 1
      return precision

    if nrefls > 0 and self.bindata.size() != points.size():
      raise Sorry("Not the same number of reflections in bin-data and displayed data")

    start_time = time.time()
    for i, hklstars in enumerate(points):
      # bin currently displayed data according to the values of another miller array
      ibin = data2bin( self.bindata[i], self.binvalsboundaries, self.nbinvalsboundaries )
      self.positions[ibin].extend( hklstars )
      self.colours[ibin].extend( colors[i] )
      self.radii2[ibin].append( radii[i] )
      self.spbufttips[ibin].append( i )

    elapsed_time = time.time() - start_time
    self.mprint("elapsed time: %s" %elapsed_time, verbose=1)

    if not blankscene:
      cntbin = 0
      self.binstrs = []
      self.bin_infotpls = []
      if self.nuniqueval < self.params.binning.nbins:
        self.mprint("%d bins was requested but %s data has only %d unique value(s)!" %(self.params.binning.nbins, colstr, self.nuniqueval), 0)
      for ibin in range(self.nbinvalsboundaries+1):
        mstr =""
        nreflsinbin = len(self.radii2[ibin])
        bin2 = float("nan"); bin1= float("nan") # indicates un-binned data
        if ibin == self.nbinvalsboundaries:
          mstr= "bin[%d] has %d reflections with no %s values (assigned to %2.3f)" %(cntbin, nreflsinbin, \
                  colstr, bin1)
        precision = 3
        if ibin < (self.nbinvalsboundaries-1):
          bin1 = self.binvalsboundaries[ibin]
          bin2 = self.binvalsboundaries[ibin+1]
          bin3 = bin2
          if ibin < (self.nbinvalsboundaries-2):
            bin3= self.binvalsboundaries[ibin+2]
          if colstr=="dres":
            bin1= 1.0/self.binvalsboundaries[ibin]
            bin2= 1.0/self.binvalsboundaries[ibin+1]
            if ibin < (self.nbinvalsboundaries-2):
              bin3= 1.0/self.binvalsboundaries[ibin+2]
          #calculate precision by comparing a bin value with bin value below and above it
          prec1 = getprecision(bin1, bin2)
          prec2 = prec1
          if bin2 != bin3:
            prec2 = getprecision(bin3, bin2)
          precision = max(prec1, prec2)
          # format bin values string with necessary decimal places (precision)
          binformatstr = "]%2." + str(precision) + "f; %2." + str(precision) + "f]"
          mstr= "bin[%d] has %d reflections with %s in " %(cntbin, nreflsinbin, colstr)
          mstr += binformatstr %(bin1, bin2)
        if len(self.bin_infotpls) > 0 \
         and math.isnan(self.bin_infotpls[-1][1]) \
         and math.isnan(self.bin_infotpls[-1][2])  \
         and nreflsinbin == 0:
          continue

        if nreflsinbin > 0 and not math.isnan(bin1) and not math.isnan(bin2):
          self.bin_infotpls.append( roundoff((nreflsinbin, bin1, bin2 ), precision) )
        self.binstrs.append(mstr)
        self.mprint(mstr, verbose=1)
        cntbin += 1

      if self.params.binning.bin_opacity != None:
        opqlist = self.params.binning.bin_opacity
        if len(opqlist) < self.params.binning.nbins-1:
          # an extra bin may be added when editing scene_bin_thresholds. If so, don't reset opacities to 1
          self.params.binning.bin_opacity = [ [1.0, e] for e in range(cntbin) ]

      self.params.binning.nbins = len(self.bin_infotpls)
      self.SendInfoToGUI( { "bin_opacity": self.params.binning.bin_opacity,
                            "bin_infotpls": self.bin_infotpls,
                            "binner_idx": self.params.binning.binner_idx,
                            "tooltip_opacity": self.params.NGL.tooltip_alpha
                           } )

      self.calc_rotation_axes()
      nvaluelabels = int(ln/self.params.viewer.ncolourlabels )
      colourgradstrs = []
      if self.params.hkls.sigma_color_radius:
        lst = list(colourscalararray)
        lst.reverse() # flip labels on chart when showing sigmas
        colourscalararray = flex.double(lst )

      # if displaying phases from map coefficients together with fom values then
      for g,colourgradarray in enumerate(colourgradarrays):
        self.colourgradientvalues = []
        for j,e in enumerate(colourgradarray):
          self.colourgradientvalues.append( [colourscalararray[j], e] )
        self.colourgradientvalues = roundoff( self.colourgradientvalues )
        fom = fomarrays[g]
        colourgradstr = []
        for j,val in enumerate(self.colourgradientvalues):
          vstr = "null"
          alpha = 1.0
          rgb = (int(val[1][0]), int(val[1][1]), int(val[1][2]) )
          if j % nvaluelabels == 0 or j==(len(self.colourgradientvalues)-1) :
            vstr = roundoff(val[0], 2)
          colourgradstr.append([vstr, rgb[0], rgb[1], rgb[2] ])
        colourgradstrs.append(colourgradstr)

    self.mprint(".", end="")
    if not self.WBmessenger.browserisopen:
      self.ReloadNGL()
      # if sempahore is not free then websocket failed to connect to a browser. Critical error!
      self.mprint("DrawNGLJavaScript waiting for browser_connect_sem semaphore", verbose="threadingmsg")
      if not self.browser_connect_sem.acquire(timeout= lock_timeout):
        raise HKLviewerError("Timed out connecting to a web browser after %s seconds. Ensure web browser security settings permit websocket protocol." %lock_timeout)
      self.mprint("DrawNGLJavaScript released browser_connect_sem", verbose="threadingmsg")
      self.browser_connect_sem.release()

    if self.verbosebrowser:
      self.SetBrowserDebug("true")

    if self.params.use_wireframe:
      self.UseWireFrame("true")
    else:
      self.UseWireFrame("false")

    if not blankscene: # and self.webgl_OK:
      self.RemoveStageObjects()
      self.SetFontSize(self.params.NGL.fontsize)
      self.SetVectorWidth(self.params.NGL.vector_width)
      for ibin in range(self.nbinvalsboundaries+1):
        nreflsinbin = len(self.radii2[ibin])
        self.DefineHKL_Axes(str(Hstararrowstart), str(Hstararrowend),
          str(Kstararrowstart), str(Kstararrowend),
          str(Lstararrowstart), str(Lstararrowend),
          Hstararrowtxt, Kstararrowtxt, Lstararrowtxt )
        self.SendCoordinates2Browser(self.positions[ibin], self.colours[ibin],
                                     self.radii2[ibin], self.spbufttips[ibin] )
      self.mprint(".", end="")
      self.RenderStageObjects()
      self.mprint(".", end="")
      self.MakeColourChart(colourlabel, fomlabel, colourgradstrs)
      self.GetClipPlaneDistances()
      self.SetMouseSpeed( self.params.NGL.mouse_sensitivity )
    self.sceneisdirty = False
    self.lastscene_id = self.params.viewer.scene_id
    self.SendInfoToGUI( { "CurrentDatatype": self.get_current_datatype(),
         "current_labels": self.get_current_labels() } )
    self.mprint("\nSubmitted reflections and other objects to browser for rendering.", verbose=1)


  def get_current_labels(self):
    return self.get_label_type_from_scene_id(self.params.viewer.scene_id)[0]


  def get_visible_current_miller_array(self):
    arrayidxs = []
    if self.miller_array and self.params.binning.bin_opacity:
      bin_opacitieslst = self.params.binning.bin_opacity
      for alpha,bin in bin_opacitieslst:
        ibin = int(bin)
        if ibin > self.nbinvalsboundaries:
          continue
        if alpha==1.0:
          arrayidxs.extend(self.spbufttips[ibin])
    visarray = self.miller_array.select_indices(flex.miller_index(
      # has to be a better way of doing this
                          [ self.miller_array.indices()[i] for i in arrayidxs ] ))
    return visarray.deep_copy()


  def release_all_semaphores(self):
    # avoid potential deadlock by releasing any pending sempahores
    self.clipplane_msg_sem.release()
    self.autoview_sem.release()
    self.mousespeed_msg_sem.release()
    self.hkls_drawn_sem.release()
    self.browser_connect_sem.release()
    self.mprint( "All sempahores released", verbose="threadingmsg")


  def ProcessBrowserMessage(self, message):
    # Runs in WebsocketClientMessageThread handling messages from the browser displaying our reflections
    # started in webbrowser_messenger_py3
    try:
      if sys.version_info[0] > 2:
        ustr = str
      else:
        ustr = unicode
      if isinstance(message, bytes) and isinstance(self.lastmsg, ustr) and "Imageblob" in self.lastmsg:
        self.mprint( "Saving image to %s" %self.imagename, verbose=0)
        with open( self.imagename, "wb") as imgfile:
          imgfile.write( message)
      philchanged = False
      if isinstance(message, ustr) and message != "":
        if 'Critical WebGL problem' in message:
          self.mprint(message + "\n\nCommencing initiation of protocols for invoking program termination procedures...\n", verbose=0)
          self.webgl_OK = False
          self.SendInfoToGUI( { "closing_time": True } )
        elif 'Browser.JavaScript' in message:
          self.mprint( message, verbose=1)
        elif "JavaScriptError" in message:
          self.mprint( message, verbose=0)
          self.release_all_semaphores()
        elif "Orientation" in message:
          self.ProcessOrientationMessage(message)
        elif 'WebGL' in message:
          self.mprint( message, verbose=1)
        elif "websocket" in message:
          self.mprint( message, verbose=1)
        elif "Refreshing" in message or "disconnecting" in message:
          self.mprint( message, verbose=1)
          time.sleep(self.sleeptime)
        elif "AutoViewSet" in message:
          self.set_volatile_params()
          self.mprint( message, verbose=3)
        elif "SetAutoView" in message:
          self.mprint( message, verbose=3)
        elif "AutoViewFinished_AfterRendering" in message:
          self.mprint("ProcessBrowserMessage, %s released autoview_sem" %message, verbose="threadingmsg")
          self.autoview_sem.release()
        elif "JavaScriptCleanUpDone:" in message:
          self.mprint( message, verbose=1)
          time.sleep(0.5) # time for browser to clean up
          if not self.isnewfile:
            self.WBmessenger.StopWebsocket()
        elif "Expanded rotation operator" in message:
          self.mprint( message, verbose="expansionmsg")
        elif "Expand" in message:
          self.mprint( message, verbose=2)
        elif "Connection lost" in message:
          self.mprint( message, verbose=1)
        elif "Warning!: Web browser closed unexpectedly" in message:
          self.mprint( message, verbose=1)
        elif "ToggleAnimation" in message:
          vecnr,speed = eval(self.params.viewer.animate_rotation_around_vector)
          speed = -speed # negative speed tells HKLjavascripts to pause animating
          self.params.viewer.animate_rotation_around_vector = "[%s, %s]" %(vecnr,speed)
          philchanged = True
          self.parent.SendCurrentPhilValues() # update GUI to correspond to current phil parameters
        elif "Imageblob" in message:
          self.mprint( "Image blob to be received", verbose=1)
        elif "ImageWritten" in message:
          self.mprint( "Image blob sent to CCTBX", verbose=1)
          self.mprint("ProcessBrowserMessage, ImageWritten released self.hkls_drawn_sem", verbose="threadingmsg")
          self.hkls_drawn_sem.release()
        elif "ClipPlaneDistancesSet" in message:
          if self.use_semaphore:
            self.mprint("ProcessBrowserMessage, ClipPlaneDistancesSet released clipplane_msg_sem", verbose="threadingmsg")
            self.clipplane_msg_sem.release() # as was set by make_clip_plane
        elif "ExpandedInBrowser_AfterRendering" in message:
          if self.use_semaphore:
            self.mprint("ProcessBrowserMessage, ExpandedInBrowser released clipplane_msg_sem", verbose="threadingmsg")
            self.clipplane_msg_sem.release() # as was set by make_clip_plane
        elif "ReturnClipPlaneDistances:" in message:
          datastr = message[ message.find("\n") + 1: ]
          lst = datastr.split(",")
          flst = [float(e) for e in lst[0:4]]
          self.clipNear = flst[0]
          self.clipFar = flst[1]
          self.cameraPosZ = flst[2]
          self.zoom = flst[3]
          calledby = lst[4]
          self.mprint("ReturnClipPlaneDistances(%s): cameraPosZ: %s, zoom: %s" %(calledby, self.cameraPosZ, self.zoom), verbose="orientmsg")
          if self.use_semaphore:
            if calledby == "GetClipPlaneDistances": # only unlock if requested by GetClipPlaneDistances()
              self.mprint("ProcessBrowserMessage, ReturnClipPlaneDistances released hkls_drawn_sem", verbose="threadingmsg")
              self.hkls_drawn_sem.release()
              self.mprint("ProcessBrowserMessage, ReturnClipPlaneDistances released clipplane_msg_sem", verbose="threadingmsg")
              self.clipplane_msg_sem.release()
              self.mprint("ProcessBrowserMessage, ReturnClipPlaneDistances released autoview_sem", verbose="threadingmsg")
              self.autoview_sem.release()
            if calledby == "RotateStage": # only unlock if requested by GetClipPlaneDistances()
              self.mprint("ProcessBrowserMessage, RotateStage released clipplane_msg_sem", verbose="threadingmsg")
              self.clipplane_msg_sem.release()
        elif "ReturnMouseSpeed" in message:
          datastr = message[ message.find("\n") + 1: ]
          lst = datastr.split(",")
          flst = [float(e) for e in lst]
          if flst[0] is not None and not cmath.isnan(flst[0]):
            self.params.NGL.mouse_sensitivity = flst[0]
          if self.use_semaphore:
            self.mprint("ProcessBrowserMessage, ReturnMouseSpeed released mousespeed_msg_sem", verbose="threadingmsg")
            self.mousespeed_msg_sem.release()
        elif "tooltip_id:" in message:
          ttipids = message.split("tooltip_id:")[1]
          hklid = eval(message.split("tooltip_id:")[1])[0]
          sym_id = eval(message.split("tooltip_id:")[1])[1]
          ttip = self.GetTooltipOnTheFly(hklid, sym_id)
          self.AddToBrowserMsgQueue("ShowThisTooltip", ttip)
        elif "match_hkl_id:" in message:
          hklid = eval(message.split("match_hkl_id:")[1])[0]
          sym_id = eval(message.split("match_hkl_id:")[1])[1]
          if self.sg.info().symbol_and_number() == self.miller_array.space_group().info().symbol_and_number():
            self.make_visual_symHKLs(hklid, sym_id)
            self.visualise_sym_HKLs()
            hkl = self.scene.indices[abs(hklid)]
            hklmatches = miller.match_indices(self.parent.origarrays["HKLs"], [hkl])
            orig_hkl_ids = list(hklmatches.pairs().column(0))
            self.SendInfoToGUI( { "clicked_HKL": hkl, "orig_hkl_ids": orig_hkl_ids })
        elif "onClick colour chart" in message:
          self.onClickColourChart()
        elif "SelectedBrowserDataColumnComboBox" in message:
          sceneid = int(message.split(":")[1])
          self.parent.SetScene(sceneid)
        elif "InFrustum:" in message:
          # if GetReflectionsInFrustum() finds no reflections
          # then message="InFrustum::" which crashes eval(). Avoid this
          if "InFrustum::" not in message:
            hklids = eval(message.split(":")[1])
            rotids = eval(message.split(":")[2])
            self.visible_hkls = []
            self.outsideplane_hkls = []
            for i,hklid in enumerate(hklids):
              hkl, _ = self.get_rothkl_from_IDs(hklid, rotids[i])
              self.visible_hkls.append(hkl)
              if self.normal_vecnr != -1 and self.params.clip_plane.is_assoc_real_space_vector and \
                self.planescalarvalue != (self.planenormalhklvec[0]*hkl[0] + self.planenormalhklvec[1]*hkl[1] + self.planenormalhklvec[2]*hkl[2]):
                self.outsideplane_hkls.append(hkl)
            self.visible_hkls = list(set(self.visible_hkls))
            self.outsideplane_hkls = list(set(self.outsideplane_hkls))
            self.mprint( "visible hkls: " + str(self.visible_hkls), verbose="frustum")
            if len(self.outsideplane_hkls):
              self.mprint("hkls not satisfying plane equation: " + str(self.outsideplane_hkls))
              self.mprint("Consider reducing the clip plane width on the \"Slicing\" tab")
          self.mprint( message, verbose=3)
        elif "notify_cctbx_AfterRendering" in message:
          if self.use_semaphore:
            self.mprint("ProcessBrowserMessage, notify_cctbx_AfterRendering released self.hkls_drawn_sem", verbose="threadingmsg")
            self.hkls_drawn_sem.release()
          self.GetReflectionsInFrustum()
        elif "MoveClipPlanesUp" in message:
          self.params.clip_plane.hkldist += 1
          self.visual_symHKLs = []
          self.set_volatile_params(use_semaphore=False)
          philchanged = True
        elif "MoveClipPlanesDown" in message:
          self.params.clip_plane.hkldist -= 1
          self.visual_symHKLs = []
          self.set_volatile_params(use_semaphore=False)
          philchanged = True
        elif "RenderStageObjects" in message: # reflections have been drawn
          if self.use_semaphore:
            self.mprint("RenderStageObjects() has drawn reflections in the browser", verbose=1)
            self.hkls_drawn_sem.release()
        elif "Ready " in message:
          self.mprint( message, verbose=5)
        if philchanged:
          self.parent.SendCurrentPhilValues() # update GUI to correspond to current phil parameters
    except Exception as e:
      self.mprint( to_str(e) + "\n" + traceback.format_exc(limit=10), verbose=0)
    self.lastmsg = message


  def GetCameraPosRotTrans(self, viewmtrx):
    lst = viewmtrx.split(",")
    flst = [float(e) for e in lst]
    ScaleRotMx = matrix.sqr( (flst[0], flst[4], flst[8],
                          flst[1], flst[5], flst[9],
                          flst[2], flst[6], flst[10]
                          )
    )
    cameratranslation = (flst[12], flst[13], flst[14])
    self.mprint("translation: %s" %str(roundoff(cameratranslation)), verbose="orientmsg")
    alllst = roundoff(flst)
    self.mprint("""Orientation matrix:
  %s,  %s,  %s,  %s
  %s,  %s,  %s,  %s
  %s,  %s,  %s,  %s
  %s,  %s,  %s,  %s
Distance: %s
    """ %tuple(alllst), verbose="orientmsg")
    rotdet = ScaleRotMx.determinant()
    if rotdet <= 0.0:
      self.mprint("Negative orientation matrix determinant!!", verbose=1)
      self.SetAutoView() # return old values as a fall back even if they're out of date
      return self.cameraPosZ, self.currentRotmx, self.cameratranslation
    else:
      cameradist = math.pow(rotdet, 1.0/3.0)
    self.mprint("Scale distance: %s" %roundoff(cameradist), verbose="orientmsg")
    currentRotmx = matrix.identity(3)
    if cameradist > 0.0:
      currentRotmx = ScaleRotMx/cameradist
      cameraPosZ = cameradist
    return cameraPosZ, currentRotmx, cameratranslation


  def ProcessOrientationMessage(self, message):
    if self.params.viewer.scene_id is None or self.miller_array is None:
      return
    if message.find("NaN")>=0 or message.find("undefined")>=0 or message.find("Browser.JavaScript")>=0:
      return
    msgname = message[ 0 : message.find("\n")-1]
    self.viewmtrx = message[ message.find("\n") + 1: ]
    if "OrientationBeforeReload:" in message:
      if not self.isnewfile:
        self.lastviewmtrx = self.viewmtrx
      self.isnewfile = False
    self.cameraPosZ, self.currentRotmx, self.cameratranslation = self.GetCameraPosRotTrans( self.viewmtrx)
    rotlst = roundoff(self.currentRotmx.elems, 4)
    self.mprint(msgname + """, Rotation matrix:
  %s,  %s,  %s
  %s,  %s,  %s
  %s,  %s,  %s
    """ %rotlst, verbose="orientmsg")
    uc = self.miller_array.unit_cell()
    OrtMx = matrix.sqr( uc.fractionalization_matrix() )
    InvMx = OrtMx.inverse()
    # Our local coordinate system has x-axis pointing right and z axis pointing out of the screen
    # unlike threeJS so rotate the coordinates emitted from there before presenting them
    Xvec = matrix.rec([1,0,0], n=(1,3))
    Yvec = matrix.rec([0,1,0], n=(1,3))
    Zvec = matrix.rec([0,0,1], n=(1,3))

    RotAroundYMx = matrix.sqr([-1.0,0.0,0.0, 0.0,1.0,0.0, 0.0,0.0,-1.0])
    Xhkl = list(InvMx.transpose()*self.currentRotmx.inverse()* RotAroundYMx.transpose()* Xvec.transpose())
    Yhkl = list(InvMx.transpose()*self.currentRotmx.inverse()* RotAroundYMx.transpose()* Yvec.transpose())
    Zhkl = list(InvMx.transpose()*self.currentRotmx.inverse()* RotAroundYMx.transpose()* Zvec.transpose())

    if self.debug:
      self.SendInfoToGUI( { "StatusBar": "RotMx: %s, X: %s, Y: %s, Z: %s" \
        %(str(roundoff(self.currentRotmx,4)), str(roundoff(Xhkl, 2)),
                                              str(roundoff(Yhkl, 2)),
                                              str(roundoff(Zhkl, 2))),
                           }
                         )
    else:
      self.SendInfoToGUI( { "StatusBar": "X: %s , Y: %s , Z: %s" %(str(roundoff(Xhkl, 2)),
                                                     str(roundoff(Yhkl, 2)),
                                                     str(roundoff(Zhkl, 2))),
                           } )
    if "MouseMoved_Orientation:" in message:
      self.mouse_moved = True
    if self.currentRotmx.is_r3_rotation_matrix():
      # Round off matrix elements to avoid machine imprecision errors that might cast
      # any matrix element into a number strictly larger than 1 which would
      # crash r3_rotation_matrix_as_x_y_z_angles()
      self.currentRotmx = matrix.sqr(roundoff(self.currentRotmx.elems, 8) )
      angles = self.currentRotmx.r3_rotation_matrix_as_x_y_z_angles(deg=True)
      self.mprint(msgname + ", angles: %s" %str(roundoff(angles)), verbose=3)
      z_vec = flex.vec3_double( [(0,0,1)])
      self.rot_zvec = z_vec * self.currentRotmx
      self.mprint(msgname + ", Rotated cartesian Z direction : %s" %str(roundoff(self.rot_zvec[0])), verbose=3)
      rfracmx = matrix.sqr( self.miller_array.unit_cell().reciprocal().fractionalization_matrix() )
      self.rot_recip_zvec = self.rot_zvec * rfracmx
      self.rot_recip_zvec = (1.0/self.rot_recip_zvec.norm()) * self.rot_recip_zvec
      self.mprint(msgname + ", Rotated reciprocal L direction : %s" %str(roundoff(self.rot_recip_zvec[0])), verbose=3)


  def OpenBrowser(self):
    if self.params.viewer.scene_id is not None and not self.WBmessenger.websockclient \
       and not self.WBmessenger.browserisopen or self.isnewfile:
      # don't block in case we're called again and first time failed conecting to a browser
      self.browser_connect_sem.acquire(blocking = False)
      with open(self.hklfname, "w") as f:
        f.write( self.htmlstr )
      self.url = "file:///" + os.path.abspath( self.hklfname )
      self.url = self.url.replace("\\", "/")
      self.mprint( "Writing %s and connecting to its websocket client..." %self.hklfname, verbose=1)
      # ensure websockets server starts before the webbrowser loads page with javascript websocket client
      self.mprint("OpenBrowser waiting for listening_sem.acquire", verbose="threadingmsg")
      if not self.WBmessenger.listening_sem.acquire(blocking=True, timeout=lock_timeout):
        self.mprint("Error! Timed out waiting for listening_sem semaphore within %s seconds" %lock_timeout, verbose=1)
      self.mprint("OpenBrowser got listening_sem", verbose="threadingmsg")
      self.WBmessenger.listening_sem.release()
      self.mprint("OpenBrowser released listening_sem", verbose="threadingmsg")
      time.sleep(0.5)

      if self.UseOSBrowser=="default":
        if not self.webctrl.open(self.url):
          self.mprint("Could not open the default web browser")
          return False
      if self.UseOSBrowser != "default" and self.UseOSBrowser != "":
        subprocess.run('"' + self.browserpath + '" ' + self.url + ' &',
                       shell=True,
        # the following flags ensures external browser process doesn't hang during regression tests
                       capture_output=False,  # regression test wants to capture stdout/stderr
                       stdout=subprocess.DEVNULL,
                       stderr=subprocess.DEVNULL)
      self.SendInfoToGUI({ "html_url": self.url } )
      return True
    return False


  def on_browser_connection(self):
    try:
      self.browser_connect_sem.release()
      self.mprint("on_browser_connection released browser_connect_sem", verbose="threadingmsg")
      self.WBmessenger.browserisopen = True
      self.mprint("Successfully connected to browser", verbose=1)
    except ValueError as e:
      self.mprint( "Trying to reload webpage in browser", verbose=0)
      self.ReloadNGL()
    except Exception as e:
      self.mprint( to_str(e) + "\n" + traceback.format_exc(limit=10), verbose=0)


  def GetReflectionsInFrustum(self):
    msg = str(self.max_reflections_in_frustum)
    self.AddToBrowserMsgQueue("GetReflectionsInFrustum", msg)


  def RedrawNGL(self):
    self.AddToBrowserMsgQueue("Redraw")


  def ReloadNGL(self): # expensive as javascript may be several Mbytes large
    self.mprint("Rendering JavaScript...", verbose=1)
    if not self.OpenBrowser():
      self.AddToBrowserMsgQueue("Reload")


  def set_show_tooltips(self):
    msg = "%s" %self.params.NGL.show_tooltips
    self.AddToBrowserMsgQueue("DisplayTooltips", msg)


  def set_tooltip_opacity(self):
    msg = "%f" %self.params.NGL.tooltip_alpha
    self.AddToBrowserMsgQueue("TooltipOpacity", msg)


  def set_background_colour(self, r,g,b):
    msg = "rgb(%d, %d, %d)" %(r,g,b)
    self.AddToBrowserMsgQueue("BackgroundColour", msg)


  def set_opacity(self, bin, alpha):
    if bin > self.nbinvalsboundaries-1:
      return "There are only %d bins present\n" %self.nbinvalsboundaries
    msg = "%d, %f" %(bin, alpha)
    self.AddToBrowserMsgQueue("alpha", msg)
    return "Opacity %s set on bin[%d]\n" %(alpha, bin)


  def JavaScriptCleanUp(self, ):
    self.AddToBrowserMsgQueue("JavaScriptCleanUp")


  def ExpandInBrowser(self):
    """
    Expansion of reflections stored in an assymetric unit wedge defined by the spacegroup is
    done by applying the rotation matrices defined by the spacegroup on the reflections.
    Applying these matrices on all reflections is done much faster in WebGL in the browser.
    Before sending the rotation matrices to the browser first convert them into cartesian
    coordinates.
    """
    if self.sceneisdirty:
      self.mprint( "Not expanding in browser", verbose=1)
      return
    uc = self.miller_array.unit_cell()
    OrtMx = matrix.sqr( uc.orthogonalization_matrix())
    InvMx = OrtMx.inverse()
    msgtype = "Expand"
    msg = ""
    unique_rot_ops = []
    if self.params.hkls.expand_to_p1:
      msgtype += "P1"
      unique_rot_ops = self.symops[ 0 : self.sg.order_p() ] # avoid duplicate rotation matrices
      retmsg = "Expanding to P1 in browser"
      if not self.miller_array.is_unique_set_under_symmetry():
        retmsg += "\nNot all reflections are in the same asymmetric unit in reciprocal space.\n"
        retmsg += "Some reflections might be displayed on top of one another.\n"
      self.mprint( retmsg, verbose=1)
    else:
      unique_rot_ops = [ self.symops[0] ] # No P1 expansion. So only submit the identity matrix
    if self.params.hkls.expand_anomalous and not self.miller_array.anomalous_flag():
      msgtype += "Friedel"
      self.mprint( "Expanding Friedel mates in browser", verbose=1)
    for i, symop in enumerate(unique_rot_ops):
      RotMx = matrix.sqr( symop.r().as_double())
      ortrot = (OrtMx * RotMx * InvMx).as_mat3()
      if RotMx.is_r3_identity_matrix():
        # avoid machine precision rounding errors converting 1.0 to 0.99999999..
        ortrot = (1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0)
      str_rot = str(ortrot)
      str_rot = str_rot.replace("(", "")
      str_rot = str_rot.replace(")", "")
      msg += str_rot + "\n" # add rotation matrix to end of message string
    self.GuardedAddToBrowserMsgQueue(semaphorename="clipplane_msg_sem", msgtype=msgtype, msg=msg,
                                     funcname="ExpandInBrowser",  posteriorcheck=False)
    if self.use_semaphore:
      self.clipplane_msg_sem.release()
      self.mprint("ExpandInBrowser released clipplane_msg_sem", verbose="threadingmsg")

  def draw_sphere(self, s1, s2, s3, isreciprocal=True,
                  r=0, g=0, b=0, name="", radius = 1.0, mesh=False):
    """
    Place sphere at [s1, s2, s3]  with colour r,g,b. If name=="", the creation
    is deferred until draw_sphere is eventually called with name != "". These
    spheres are then joined in the same NGL representation.
    """
    uc = self.miller_array.unit_cell()
    vec = (s1*self.renderscale, s2*self.renderscale, s3*self.renderscale)
    #svec = list(vec)
    if isreciprocal:
      # uc.reciprocal_space_vector() only takes integer miller indices so compute the cartesian coordinates
      # for floating valued miller indices with the transpose of the fractionalization matrix
      vec = list( vec * matrix.sqr(uc.fractionalization_matrix()).transpose() )
      svec = [ vec[0], vec[1], vec[2] ]
    else: # real space fractional values
      vec = list( vec * matrix.sqr(uc.orthogonalization_matrix()) )
      vscale =  1.0/self.renderscale
      # TODO: find suitable scale factor for displaying real space vector together with reciprocal vectors
      svec = [ vscale*vec[0], vscale*vec[1], vscale*vec[2] ]
    self.draw_cartesian_sphere(svec[0], svec[1], svec[2], r, g, b, name, radius, mesh)


  def draw_cartesian_sphere(self, s1, s2, s3, r=0, g=0, b=0, name="", radius = 1.0, mesh=False):
    self.mprint("cartesian sphere is at: %s" %(str(roundoff([s1, s2, s3]))), verbose=2)
    self.AddToBrowserMsgQueue("DrawSphere", "%s;; %s;; %s;; %s;; %s;; %s;; %s;; %s;; %s" \
         %(s1, s2, s3, r, g, b, radius, name, int(mesh)) )
    if name=="":
      self.mprint("deferred rendering sphere at (%s, %s, %s)" %(s1, s2, s3), verbose=2)


  def draw_vector(self, s1, s2, s3, t1, t2, t3, isreciprocal=True, label="",
                  r=0, g=0, b=0, name="", radius=1.0, labelpos=0.8, autozoom = True):
    """
    Place vector from [s1, s2, s3] to [t1, t2, t3] with colour r,g,b and label
    If name=="" creation is deferred until draw_vector is eventually called with name != ""
    These vectors are then joined in the same NGL representation
    """
    uc = self.miller_array.unit_cell()
    vec1 = (s1*self.renderscale, s2*self.renderscale, s3*self.renderscale)
    vec2 = (t1*self.renderscale, t2*self.renderscale, t3*self.renderscale)
    #svec = list(vec)
    if isreciprocal:
      # uc.reciprocal_space_vector() only takes integer miller indices so compute the cartesian coordinates
      # for floating valued miller indices with the transpose of the fractionalization matrix
      vec1 = list( vec1 * matrix.sqr(uc.fractionalization_matrix()).transpose() )
      vec2 = list( vec2 * matrix.sqr(uc.fractionalization_matrix()).transpose() )
      svec1 = [ vec1[0], vec1[1], vec1[2] ]
      svec2 = [ vec2[0], vec2[1], vec2[2] ]
    else: # real space fractional values
      vec1 = list( vec1 * matrix.sqr(uc.orthogonalization_matrix()) )
      vec2 = list( vec2 * matrix.sqr(uc.orthogonalization_matrix()) )
      vscale =  1.0/self.renderscale
      # TODO: find suitable scale factor for displaying real space vector together with reciprocal vectors
      svec1 = [ vscale*vec1[0], vscale*vec1[1], vscale*vec1[2] ]
      svec2 = [ vscale*vec2[0], vscale*vec2[1], vscale*vec2[2] ]
    self.draw_cartesian_vector(svec1[0], svec1[1], svec1[2], svec2[0], svec2[1], svec2[2],
                            label, r, g, b, name, radius, labelpos, autozoom)


  def draw_cartesian_vector(self, s1, s2, s3, t1, t2, t3, label="",
                            r=0, g=0, b=0, name="", radius = 1.0, labelpos=0.8, autozoom = True ):
    self.mprint("cartesian vector is: %s to %s" %(str(roundoff([s1, s2, s3])), str(roundoff([t1, t2, t3]))), verbose="vector")
    rad = radius #self.params.NGL.vector_width
    self.AddToBrowserMsgQueue("DrawVector", "%s;;%s;;%s;;%s;;%s;;%s;;%s;;%s;;%s;;%s;;%s;;%s;;%s;;%s" \
         %(s1, s2, s3, t1, t2, t3, r, g, b, label, name, rad, labelpos, autozoom) )
    if name=="":
      self.mprint("deferred rendering vector from (%s, %s, %s) to (%s, %s, %s)" %(s1, s2, s3, t1, t2, t3), verbose=2)


  def get_cartesian_vector_angles(self, s1, s2, s3, t1, t2, t3):
    self.mprint("in get cartesian angles", verbose="angles")
    svec = [t1-s1, t2-s2, t3-s3]
    svecnorm = math.sqrt( svec[0]*svec[0] + svec[1]*svec[1] + svec[2]*svec[2] )
    xyvec = svec[:] # deep copying
    xyvec[2] = 0.0 # projection vector of svec in the xy plane
    xyvecnorm = math.sqrt( xyvec[0]*xyvec[0] + xyvec[1]*xyvec[1] )
    output = io.StringIO()
    if xyvecnorm > 0.0 and not approx_equal(xyvecnorm/svecnorm, 0.0, out=output):
      angle_x_xyvec = math.acos( xyvec[0]/xyvecnorm )*180.0/math.pi
      angle_y_xyvec = math.acos( xyvec[1]/xyvecnorm )*180.0/math.pi
    else:
      angle_x_xyvec = 90.0
      angle_y_xyvec = 90.0
    self.mprint(output.getvalue(), verbose="angles")
    output.close()
    yzvec = svec[:]
    yzvec[0] = 0.0 # projection vector of svec in the yz plane
    yzvecnorm = math.sqrt( yzvec[1]*yzvec[1] + yzvec[2]*yzvec[2] )

    output = io.StringIO()
    if yzvecnorm > 0.0 and not approx_equal(yzvecnorm/svecnorm, 0.0, out=output):
      angle_y_yzvec = math.acos( yzvec[1]/yzvecnorm )*180.0/math.pi
      angle_z_yzvec = math.acos( yzvec[2]/yzvecnorm )*180.0/math.pi
    else:
      angle_y_yzvec = 90.0
      angle_z_yzvec = 90.0
    self.mprint(output.getvalue(), verbose="angles")
    output.close()
    angle_x_svec = math.acos( svec[0]/svecnorm )*180.0/math.pi
    angle_y_svec = math.acos( svec[1]/svecnorm )*180.0/math.pi
    angle_z_svec = math.acos( svec[2]/svecnorm )*180.0/math.pi
    if angle_y_svec > 90.0:
      angle_x_xyvec = -angle_x_xyvec
    self.mprint("angles in xy plane to x,y axis are: %s, %s" %(angle_x_xyvec, angle_y_xyvec), verbose="angles")
    self.mprint("angles in yz plane to y,z axis are: %s, %s" %(angle_y_yzvec, angle_z_yzvec), verbose="angles")
    self.mprint("angles to x,y,z axis are: %s, %s, %s" %(angle_x_svec, angle_y_svec, angle_z_svec ), verbose=2)
    return angle_x_xyvec, angle_z_svec


  def PointVectorPerpendicularToScreen(self, angle_x_xyvec, angle_z_svec):
    rotmx = self.Euler2RotMatrix(( angle_x_xyvec, angle_z_svec + 180.0, 0.0 ))
    if rotmx.determinant() < 0.99999:
      self.mprint("Rotation matrix determinant is less than 1")
      return rotmx
    self.currentRotmx = rotmx
    self.RotateMxStage(rotmx)
    return rotmx


  def PointVectorParallelToScreen(self, angle_x_xyvec, angle_z_svec):
    rotmx = self.Euler2RotMatrix(( angle_x_xyvec, angle_z_svec + 90.0, 90.0 ))
    if rotmx.determinant() < 0.99999:
      self.mprint("Rotation matrix determinant is less than 1")
      return rotmx
    self.currentRotmx = rotmx
    self.RotateMxStage(rotmx)
    return rotmx


  def GetVectorAndAngleFromRotationMx(self, rot, ma=None, rectify_improper_rotation=True):
    # rectify_improper_rotation should be False when used for Xtriage rotations as
    # they may mistakenly be seen as being improper rotations
    RotMx = matrix.sqr(rot.as_double())
    if ma==None:
      ma = self.miller_array
    uc = ma.unit_cell()
    spg = ma.space_group()
    OrtMx = matrix.sqr( uc.orthogonalization_matrix())
    InvMx = OrtMx.inverse()
    ortrotmx = (OrtMx * RotMx * InvMx)
    isProperRotation = True
    ortrot = ortrotmx.as_mat3()
    label=""
    order = 0
    if not ortrotmx.is_r3_rotation_matrix() and rectify_improper_rotation:
      isProperRotation = False
      self.mprint("""Warning! The operation '%s' is not a proper rotation
in the space group %s\nwith unit cell %s""" \
        %(rot.as_hkl(), spg.info().symbol_and_number(), str(uc) ), verbose=2)
      self.mprint("Inverse of implied rotation matrix,\n%s\nis not equal to its transpose,\n%s" \
        %(str(roundoff(ortrotmx.inverse(),4)), str(roundoff(ortrotmx.transpose(),4))), verbose=2)
      improper_vec_angle = scitbx.math.r3_rotation_axis_and_angle_from_matrix(ortrot)
      self.mprint("\nTrying to find nearest orthonormal matrix approximtion", verbose=2)
      Rmx = matrix.find_nearest_orthonormal_matrix(ortrotmx)
      self.mprint("New suggested rotation matrix is\n%s" %str(roundoff(Rmx,4)), verbose=2)
      if not Rmx.is_r3_rotation_matrix() or (Rmx - ortrotmx).norm_sq() > 0.5:
      # norm_sq of a rotation matrix should be 3. Deviating from that indicates improper rotation
        self.mprint("Failed finding an approximate rotation matrix for \"%s\" in %s" \
         %(rot.as_hkl(), spg.info().symbol_and_number()), verbose=1)
        return (0,0,0), 0.0, label, order
      ortrotmx = Rmx
    ortrot = ortrotmx.as_mat3()
    r11,r12,r13,r21,r22,r23,r31,r32,r33 = ortrot
    theta =  math.acos(roundoff((r11+r22+r33-1.0)*0.5, 10))
    rotaxis = flex.vec3_double([(0,0,0)])
    self.mprint(str(ortrot), verbose=2)
    vec_angle = scitbx.math.r3_rotation_axis_and_angle_from_matrix(ortrot)
    rotaxis = flex.vec3_double([ vec_angle.axis ])
    if not isProperRotation:
      # Divine revelation: The new proper rotation from above axis is halfway
      # of being correctly aligned so subtract it from twice the improper axis
      # to get the desired rotation axis vector
      improp_rotaxis = flex.vec3_double([ improper_vec_angle.axis ])
      rotaxis = 2*rotaxis - improp_rotaxis
      # for debugging deduce the corresponding rotation matrix from this new axis
      usedrotmx = scitbx.math.r3_rotation_axis_and_angle_as_matrix( rotaxis[0], theta )
      self.mprint("Final proper rotation matrix:\n%s" %str(roundoff(matrix.sqr(usedrotmx),4)), verbose=1)
    ## adjust the length of the rotation axes to be compatible with the sphere of reflections
    if abs(theta) > 0.0001 and rotaxis.norm() > 0.01: # avoid nullvector
      order = int(roundoff(2*math.pi/theta, 0)) # how many times to rotate before its the identity operator
      forder = roundoff(2*math.pi/theta, 2)
      label = "%s-fold" %str(order)
    return list((rotaxis)[0]), theta, label, order


  def calc_rotation_axes(self, ma=None):
    if ma is not None:
      self.sg = ma.space_group()
      self.symops = list(self.sg.all_ops())
    if self.sg:
      unique_rot_ops = self.symops[ 0 : self.sg.order_p() ] # avoid duplicate rotation matrices
      self.rotation_operators = []
      for i,op in enumerate(unique_rot_ops): # skip the last op for javascript drawing purposes
        (cartvec, a, label, order) = self.GetVectorAndAngleFromRotationMx( rot=op.r(), ma=ma )
        if label != "":
          vs = 1+len(self.rotation_operators)/20
          cartvec[0] *= vs
          cartvec[1] *= vs
          cartvec[2] *= vs
          self.mprint( str(i) + ": " + str(roundoff(cartvec)) + ", " + label, verbose=1)
          veclength = math.sqrt( cartvec[0]*cartvec[0] + cartvec[1]*cartvec[1] + cartvec[2]*cartvec[2] )
          self.rotation_operators.append( (label + "#%d"%i, order , cartvec, op.r().as_hkl(), "", "", veclength) )


  def show_all_vectors(self):
    for (opnr, label, order, cartvec, hklop, hkl, abc, length) in self.all_vectors:
      self.show_labelled_vector(self.params.viewer.show_all_vectors==1, label, order, cartvec, hklop, autozoom=False)


  def show_vector(self, val, isvisible, autozoom=True):
    # val can be either the number (zero offset) of the vector in the list of vectors
    # or the label name of the vector in the list of vectors
    mag=1
    if isinstance(val, int):
      if val >= len(self.all_vectors):
        return str([])
      (opnr, label, order, cartvec, hklop, hkl, abc, length) = self.all_vectors[val]
      scale = 1
      if len(abc) > 0:
        mag = self.realSpaceMag
      if len(hkl) > 0:
        mag = self.recipSpaceMag
      self.show_labelled_vector(isvisible, label, order, cartvec, hklop,
                                autozoom=autozoom, mag=mag)
      if not isvisible:
        self.params.viewer.show_all_vectors = 0
      return str([val, isvisible])
    if isinstance(val, str):
      for i,(opnr, label, order, cartvec, hklop, hkl, abc, length) in enumerate(self.all_vectors):
        #if val in label: # so that user_vector.label="twin" declared for a preset button will match "2-fold_mytwin"
        if val == label:
          if len(abc) > 0:
            mag = self.realSpaceMag
          if len(hkl) > 0:
            mag = self.recipSpaceMag
          self.show_labelled_vector(isvisible, label, order, cartvec, hklop,
                                    autozoom=autozoom, mag=mag)
          if not isvisible:
            self.params.viewer.show_all_vectors = 0
          return str([i, isvisible])
    raise Sorry("No vector present with label or index: %s" %val)


  def show_vectors(self, philvectors, diff_phil):
    # autozoom may cause deadlock with unreleased semaphore clipplane_msg_sem
    # if more than 1 vectors are to be drawn at once. Avoid that.
    m = re.findall("(True)", str(philvectors)) # are there more than 1 vectors to be drawn?
    doautozoom = True
    if len(m) > 1 or self.params.clip_plane.clip_width: # don't autozoom if we are clipping
      doautozoom = False

    self.visual_symmxs = []
    self.visual_symHKLs = []

    for i,ivec in enumerate(philvectors):
      try:
        [val, isvisible] = eval(ivec)
        # in case val is the label for one of the vectors let show_vector() find the
        #  corresponding number and reassign ivec to "[number, bool]"
        if has_phil_path(diff_phil, "animate_rotation_around_vector"):
          # don't zoom if also initiating animation from this set of phil parameters
          ivec = self.show_vector(val, isvisible, autozoom=False)
        else:
          val2 = isvisible and doautozoom
          #if val2:
            #if not self.autoview_sem.acquire(blocking=True, timeout=lock_timeout):
            #  self.mprint("Error! Timed out waiting for autoview_sem semaphore within %s seconds" %lock_timeout, verbose=1)
            #self.mprint("show_vectors got autoview_sem", verbose="threadingmsg")
          ivec = self.show_vector(val, isvisible, autozoom=val2)
        philvectors[i] = ivec
      except Exception as e:
        pass


  def show_labelled_vector(self, isvisible, label, order, cartvec, hklop, autozoom=True, mag=1):
    # avoid onMessage-DrawVector in HKLJavaScripts.js misinterpreting the commas in strings like "-x,z+y,-y"
    name = label + "_" + hklop.replace(",", "_")
    if isvisible:
      self.currentrotvec = cartvec # cartesian vector to display and used for aligning
      if order > 0 and hklop != "":
# if this is a rotation operator deduce the group of successive rotation matrices it belongs to
        rt = sgtbx.rt_mx(symbol= hklop, r_den=12, t_den=144)
        RotMx = matrix.sqr(rt.r().as_double() )
        self.visual_symmxs.append( (RotMx, rt.r().as_hkl()) )
        nfoldrotmx = RotMx
        nfoldrot = rt.r()
        self.visual_symmxs = []
        self.visual_symHKLs = []
        for n in range(order): # append successive rotations to self.visual_symmxs
          nfoldrotmx = RotMx * nfoldrotmx
          nfoldrot = nfoldrot.multiply( rt.r() )
          self.visual_symmxs.append( (nfoldrotmx, nfoldrot.as_hkl()) )
        # adjust the length of the rotation axes to be compatible with the sphere of reflections
        uc = self.miller_array.unit_cell()
        OrtMx = matrix.sqr( uc.orthogonalization_matrix())
        s = math.sqrt(OrtMx.transpose().norm_sq())*self.realspace_scale*mag
        self.currentrotvec = [s*cartvec[0], s*cartvec[1], s*cartvec[2]]
      self.currentrotvec = [mag*self.currentrotvec[0], mag*self.currentrotvec[1], mag*self.currentrotvec[2]]
      self.draw_cartesian_vector(0, 0, 0, self.currentrotvec[0], self.currentrotvec[1],
                                  self.currentrotvec[2], r=0.1, g=0.1,b=0.1,
                                  label=label, name=name, radius=0.8, labelpos=1.0, autozoom=autozoom)
    else:
      self.RemovePrimitives(name)
    self.RemovePrimitives("sym_HKLs") # delete other symmetry hkls from a previous rotation operator if any


  def visualise_sym_HKLs(self):
    self.RemovePrimitives("sym_HKLs")
    if len(self.visual_symHKLs):
      for i,(hkl,hklstr) in enumerate(self.visual_symHKLs):
        thkl = tuple(hkl)
        hklstr = "H,K,L: %d,%d,%d" %thkl
        if i < len(self.visual_symHKLs)-1:
          self.draw_vector(0,0,0, hkl[0],hkl[1],hkl[2], isreciprocal=True, label=hklstr, r=0.5, g=0.3, b=0.3,
                           radius=0.8, labelpos=1.0, autozoom = False)
        else: # supplying a name for the vector last graphics primitive draws them all
          self.draw_vector(0,0,0, hkl[0],hkl[1],hkl[2], isreciprocal=True, label=hklstr, name="sym_HKLs",
                           r=0.5, g=0.3, b=0.3, radius=0.8, labelpos=1.0, autozoom = False)


  def show_hkl(self, bigwireframe=True):
    """
    Draw a wireframe sphere around a reflection selected with a double click in
    the millerarraytable in the GUI
    """
    rad = self.HKLscene_dict_val().scene.max_radius*1.5
    if not bigwireframe:
      rad = self.HKLscene_dict_val().scene.min_radius*0.9
    self.RemovePrimitives("highlight_HKL")
    if self.params.viewer.show_hkl != "deselect":
      hkl = eval(self.params.viewer.show_hkl)
      if self.sg.info().symbol_and_number() == self.miller_array.space_group().info().symbol_and_number():
        self.draw_sphere(hkl[0],hkl[1],hkl[2], isreciprocal=True, name="highlight_HKL",
                          r=1, g=0.0, b=0.0, radius= rad, mesh=True)
      else:
        self.mprint("Cannot currently associate reflection in original space group with reflection in different space group.")
    self.params.viewer.show_hkl = "" # to allow clicking on the same entry in the millerarraytable


  def get_vectors_labels_from_ids(self, idvectorlst):
    labelveclst = []
    for idvecval in idvectorlst:
      try:
        id, some_val = eval(idvecval)
        if isinstance(id, int):
          for opnr, label, order, cartvec, hklop, hkl, abc, length in self.all_vectors:
            if opnr==id:
              labelveclst.append([label, some_val])
              break
        if isinstance(id, str):
          labelveclst.append([id, some_val])
      except Exception as e:
        pass
    return labelveclst


  def get_vecid_from_label(self, val ):
    vecnr = -1
    if isinstance(val, int):
      vecnr = val
    if isinstance(val, str):
      for i,(opnr, label, order, cartvec, hklop, hkl, abc, length) in enumerate(self.all_vectors):
        if val in label:
          vecnr = i
    if not (vecnr>=0 and vecnr < len(self.all_vectors)):
      raise Sorry("No vector present in file with label or index: %s" %val)
    return vecnr


  def rotate_around_numbered_vector(self):
    val, deg = eval(self.params.viewer.angle_around_vector)
    vecnr = self.get_vecid_from_label(val)
    self.rotate_components_around_cartesian_vector(self.all_vectors[vecnr][3], deg)
    return vecnr,deg


  def rotate_components_around_cartesian_vector(self, cartvec, deg):
    phi = cmath.pi*deg/180
    normR = math.sqrt(cartvec[0]*cartvec[0] + cartvec[1]*cartvec[1] + cartvec[2]*cartvec[2] )
    ux = cartvec[0]/normR
    uy = cartvec[1]/normR
    uz = cartvec[2]/normR
    self.RotateAxisComponents([ux,uy,uz], phi, True)


  def rotate_stage_around_cartesian_vector(self, cartvec, deg):
    phi = cmath.pi*deg/180
    normR = math.sqrt(cartvec[0]*cartvec[0] + cartvec[1]*cartvec[1] + cartvec[2]*cartvec[2] )
    ux = cartvec[0]/normR
    uy = cartvec[1]/normR
    uz = cartvec[2]/normR
    self.RotateAxisMx([ux,uy,uz], phi, True)


  def animate_rotate_around_vector(self):
    val, speed = eval(self.params.viewer.animate_rotation_around_vector)
    vecnr = -1
    vecnr = self.get_vecid_from_label(val)
    cartvec = self.all_vectors[vecnr][3]
    normR = math.sqrt(cartvec[0]*cartvec[0] + cartvec[1]*cartvec[1] + cartvec[2]*cartvec[2] )
    ux = cartvec[0]/normR
    uy = cartvec[1]/normR
    uz = cartvec[2]/normR
    self.AnimateRotateAxisComponents([ux,uy,uz], speed, True)
    return vecnr,speed


  def DrawUnitCell(self):
    if self.params.draw_real_space_unit_cell is False:
      self.RemovePrimitives("unitcell")
      self.mprint( "Removing real space unit cell", verbose=2)
      return
    rad = 0.8 # scale # * 0.05 #  1000/ uc.volume()
    self.draw_vector(0,0,0, self.realSpaceMag,0,0, False, label="a", r=0.5, g=0.8, b=0.8, radius=rad)
    self.draw_vector(0,0,0, 0,self.realSpaceMag,0, False, label="b", r=0.8, g=0.5, b=0.8, radius=rad)
    self.draw_vector(0,0,0, 0,0,self.realSpaceMag, False, label="c", r=0.8, g=0.8, b=0.5, radius=rad)
    self.draw_vector(self.realSpaceMag,0,0, self.realSpaceMag,self.realSpaceMag,0, False, r=0.8, g=0.5, b=0.8, radius=rad)
    self.draw_vector(0,self.realSpaceMag,0, self.realSpaceMag,self.realSpaceMag,0, False, r=0.5, g=0.8, b=0.8, radius=rad)
    self.draw_vector(0,0,self.realSpaceMag, self.realSpaceMag,0,self.realSpaceMag, False, r=0.5, g=0.8, b=0.8, radius=rad)
    self.draw_vector(0,0,self.realSpaceMag, 0,self.realSpaceMag,self.realSpaceMag, False, r=0.8, g=0.5, b=0.8, radius=rad)
    self.draw_vector(0,self.realSpaceMag,self.realSpaceMag, self.realSpaceMag,self.realSpaceMag,self.realSpaceMag, False, r=0.5, g=0.8, b=0.8, radius=rad)
    self.draw_vector(self.realSpaceMag,0,self.realSpaceMag, self.realSpaceMag,self.realSpaceMag,self.realSpaceMag, False, r=0.8, g=0.5, b=0.8, radius=rad)
    self.draw_vector(self.realSpaceMag,0,0, self.realSpaceMag,0,self.realSpaceMag, False, r=0.8, g=0.8, b=0.5, radius=rad)
    self.draw_vector(0,self.realSpaceMag,0, 0,self.realSpaceMag,self.realSpaceMag, False, r=0.8, g=0.8, b=0.5, radius=rad)
    self.draw_vector(self.realSpaceMag,self.realSpaceMag,0, self.realSpaceMag,self.realSpaceMag,self.realSpaceMag, False, r=0.8, g=0.8, b=0.5, radius=rad,
                     name="unitcell", autozoom=False)
    self.mprint( "Adding real space unit cell", verbose=1)


  def DrawReciprocalUnitCell(self):
    if self.params.draw_reciprocal_unit_cell is False:
      self.RemovePrimitives("reciprocal_unitcell")
      self.mprint( "Removing reciprocal unit cell", verbose=2)
      return
    rad = 0.8 # 0.05 * scale
    self.draw_vector(0,0,0, self.recipSpaceMag,0,0, label="a*", r=0.5, g=0.3, b=0.3, radius=rad)
    self.draw_vector(0,0,0, 0,self.recipSpaceMag,0, label="b*", r=0.3, g=0.5, b=0.3, radius=rad)
    self.draw_vector(0,0,0, 0,0,self.recipSpaceMag, label="c*", r=0.3, g=0.3, b=0.5, radius=rad)
    self.draw_vector(self.recipSpaceMag,0,0, self.recipSpaceMag,self.recipSpaceMag,0, r=0.3, g=0.5, b=0.3, radius=rad)
    self.draw_vector(0,self.recipSpaceMag,0, self.recipSpaceMag,self.recipSpaceMag,0, r=0.5, g=0.3, b=0.3, radius=rad)
    self.draw_vector(0,0,self.recipSpaceMag, self.recipSpaceMag,0,self.recipSpaceMag, r=0.5, g=0.3, b=0.3, radius=rad)
    self.draw_vector(0,0,self.recipSpaceMag, 0,self.recipSpaceMag,self.recipSpaceMag, r=0.3, g=0.5, b=0.3, radius=rad)
    self.draw_vector(0,self.recipSpaceMag,self.recipSpaceMag, self.recipSpaceMag,self.recipSpaceMag,self.recipSpaceMag, r=0.5, g=0.3, b=0.3, radius=rad)
    self.draw_vector(self.recipSpaceMag,0,self.recipSpaceMag, self.recipSpaceMag,self.recipSpaceMag,self.recipSpaceMag, r=0.3, g=0.5, b=0.3, radius=rad)
    self.draw_vector(self.recipSpaceMag,0,0, self.recipSpaceMag,0,self.recipSpaceMag, r=0.3, g=0.3, b=0.5, radius=rad)
    self.draw_vector(0,self.recipSpaceMag,0, 0,self.recipSpaceMag,self.recipSpaceMag, r=0.3, g=0.3, b=0.5, radius=rad)
    self.draw_vector(self.recipSpaceMag,self.recipSpaceMag,0, self.recipSpaceMag,self.recipSpaceMag,self.recipSpaceMag, r=0.3, g=0.3, b=0.5, radius=rad,
                     name="reciprocal_unitcell", autozoom=False)
    self.mprint( "Adding reciprocal unit cell", verbose=1)


  def GetUnitcellScales(self):
    spanmin, spanmax = ( self.miller_array.index_span().min(), self.miller_array.index_span().max())
    uc = self.miller_array.unit_cell()
    vec = (1.0, 1.0, 1.0)
    # uc.reciprocal_space_vector() only takes integer miller indices so compute
    # the cartesian coordinates for real valued miller indices with the transpose of the fractionalization matrix
    vec1 = vec * matrix.sqr(uc.fractionalization_matrix()).transpose()
    reciproc_bodydiagonal_length = vec1.length()
    reciprocspanmaxvec = spanmax * matrix.sqr(uc.fractionalization_matrix()).transpose()
    reciproc_spanmax_length = reciprocspanmaxvec.length()
    reciprocspanminvec = spanmax * matrix.sqr(uc.fractionalization_matrix()).transpose()
    reciproc_spanmin_length = reciprocspanminvec.length()
    reciprocspan_length = max(reciproc_spanmax_length, reciproc_spanmin_length)
    self.reciproc_scale = reciprocspan_length / reciproc_bodydiagonal_length
    # for real space vector
    vec2 = vec * matrix.sqr(uc.orthogonalization_matrix())
    bodydiagonal_length =  vec2.length()
    self.realspace_scale = self.renderscale * reciprocspan_length / bodydiagonal_length


  def real_space_associated_with_reciprocal_vector(self, hklvec):
    # Get corresponding real space vector (in real space units) to the hkl vector
    uc = self.miller_array.unit_cell()
    R= hklvec[0] * self.normal_kl + hklvec[1] * self.normal_lh - hklvec[2] * self.normal_hk
    return (R[0]*matrix.sqr(uc.orthogonalization_matrix()).inverse()) * self.renderscale
    #return hklvec * matrix.sqr(uc.orthogonalization_matrix())


  def reciprocal_from_real_space_vector(self, realspacevec):
    uc = self.miller_array.unit_cell()
    return matrix.sqr(uc.orthogonalization_matrix()).inverse() * realspacevec


  def reciprocal_associated_with_real_space_vector(self, hklvec):
    uc = self.miller_array.unit_cell()
    cartvec = hklvec * matrix.sqr(uc.fractionalization_matrix()).transpose()
    myhkl = matrix.sqr(uc.orthogonalization_matrix()).inverse() * cartvec * self.renderscale
    #hkl = myhkl[0] * self.normal_bc + myhkl[1] * self.normal_ca - myhkl[2] * self.normal_ab
    return myhkl



  def project_vector1_vector2(self, vec1, vec2):
    # cartesian projection of vec1 onto vec2
    L1 = math.sqrt( vec1[0]*vec1[0] + vec1[1]*vec1[1] + vec1[2]*vec1[2] )
    L2 = math.sqrt( vec2[0]*vec2[0] + vec2[1]*vec2[1] + vec2[2]*vec2[2] )
    dotproduct = vec1[0]*vec2[0] + vec1[1]*vec2[1] + vec1[2]*vec2[2]
    cosine = dotproduct/(L1*L2)
    projvec1 = (vec1[0]*cosine, vec1[1]*cosine, vec1[2]*cosine)
    projvec2 = (vec2[0]*cosine, vec2[1]*cosine, vec2[2]*cosine)
    return cosine, projvec1, projvec2


  def orient_vector_to_screen(self, cartvec):
    if cartvec is None:
      return
    angle_x_xyvec, angle_z_svec = self.get_cartesian_vector_angles(0, 0, 0,
                                                                    cartvec[0],
                                                                    cartvec[1],
                                                                    cartvec[2])
    cartveclength = math.sqrt(cartvec[0]*cartvec[0] + cartvec[1]*cartvec[1] + cartvec[2]*cartvec[2] )
    self.mprint( "cartveclength= %s" %roundoff(cartveclength), verbose=1)

    if self.params.viewer.is_parallel:
      self.PointVectorParallelToScreen(angle_x_xyvec, angle_z_svec)
    else:
      self.PointVectorPerpendicularToScreen(angle_x_xyvec, angle_z_svec)


  def fix_orientation(self):
    if self.params.viewer.fixorientation != "None":
      self.DisableMouseRotation()
    else:
      self.EnableMouseRotation()


  def make_clip_plane(self, hkldist=0.0, clipwidth=None):
    # create clip plane oriented parallel or perpendicular to abc vector
    if clipwidth is None:
      self.SetClipPlaneDistances(0, 0) # tells NGL not to do clipping
      return
    self.mprint("Applying clip plane to reflections", verbose=1)
    self.RemovePrimitives("clip_vector")
    if self.cameraPosZ is None or self.cameraPosZ == 1.0:
      self.GetClipPlaneDistances()
    if self.use_semaphore:
      if not self.clipplane_msg_sem.acquire(blocking=True, timeout=lock_timeout):
        self.mprint("Error! Timed out waiting for clipplane_msg_sem semaphore within %s seconds" %lock_timeout, verbose=1)
      self.mprint("make_clip_plane got clipplane_msg_sem", verbose="threadingmsg")
      if not self.autoview_sem.acquire(blocking=True, timeout=lock_timeout):
        self.mprint("Error! Timed out waiting for autoview_sem semaphore within %s seconds" %lock_timeout, verbose=1)
      self.mprint("make_clip_plane got autoview_sem", verbose="threadingmsg")
    halfdist = self.cameraPosZ + hkldist # self.viewer.boundingZ*0.5
    if clipwidth == 0.0:
      clipwidth = self.meanradius
    clipNear = halfdist - clipwidth # 50/self.viewer.boundingZ
    clipFar = halfdist + clipwidth  #50/self.viewer.boundingZ
    if self.use_semaphore:
      self.mprint("make_clip_plane released clipplane_msg_sem, autoview_sem", verbose="threadingmsg")
      self.clipplane_msg_sem.release()
      self.autoview_sem.release()
    self.SetClipPlaneDistances(clipNear, clipFar, -self.cameraPosZ, self.zoom)
    self.mprint("clipnear: %s, clipfar: %s, cameraZ: %s, zoom: %s" %(clipNear, clipFar, -self.cameraPosZ, self.zoom), verbose=1)


  def set_camera_type(self):
    self.AddToBrowserMsgQueue("SetCameraType", self.params.NGL.camera_type)


  def set_background_colour(self):
    self.AddToBrowserMsgQueue("BackgroundColour", self.params.NGL.background_colour)


  def get_labels_of_data_for_binning(self, arrayinfos):
    self.hkl_scenes_infos = []
    sceneid = 0
    for pidx,arrayinfo in enumerate(arrayinfos):
      hassigmas=True
      if math.isnan(arrayinfo.maxsigmas):
        hassigmas=False
      (dummy1, infolst, dummy2, dummy3), dummy4, dummy5 = arrayinfo.get_selected_info_columns_from_phil()

      fomsarrays_idx = [(None, None)]
      mextnd = self.mapcoef_fom_dict.get(infolst[0])
      if infolst[1] in ['Map coeffs'] and mextnd is not None:
        fomsarrays_idx.extend( mextnd )
      for (fomsarray, fidx) in fomsarrays_idx:
        lbl = arrayinfo.labelstr
        fomslabel = None
        if fomsarray:
          fomslabel = fomsarray.info().label_string()
          lbl = arrayinfo.labelstr + " + " + fomslabel
        self.hkl_scenes_infos.append([infolst, pidx, fidx, lbl, infolst[1], hassigmas, sceneid])
        sceneid += 1

    scenearraylabeltypes = [ (e[3], e[4], e[1], e[5], e[6]) for e in self.hkl_scenes_infos ]
    self.SendInfoToGUI({ "scene_array_label_types": scenearraylabeltypes, "NewHKLscenes" : True })

    self.bin_labels_type_idxs = []
    self.bin_labels_type_idxs.append(("Resolution",  "", -1 ))
    self.bin_labels_type_idxs.append(("Singletons (current data)", "", -1 ))
    for label,labeltype,idx,hassigmas,sceneid in scenearraylabeltypes:
      if labeltype not in  ["Map coeffs", "Map coeffs_fom", "HL coeffs"]:
        self.bin_labels_type_idxs.append((label, labeltype, sceneid))
      if hassigmas:
        self.bin_labels_type_idxs.append(("Sigmas of " + label, "hassigmas", sceneid))
      if labeltype == "Map coeffs":
        self.bin_labels_type_idxs.append(("Phases of " + label, labeltype, sceneid))
        self.bin_labels_type_idxs.append(("Amplitudes of " + label, labeltype, sceneid))

    self.SendInfoToGUI({ "bin_labels_type_idxs": self.bin_labels_type_idxs})

    self.mprint("Data can be binned according to:", verbose=1)
    for i,e in enumerate(self.bin_labels_type_idxs):
      self.mprint("%d, %s" %(i, e[0]), verbose=1)


  def get_binner_idx_from_label(self, binlabel):
    for i,e in enumerate(self.bin_labels_type_idxs):
      if binlabel == e[0]:
        return i
    return -1


  def SetFontSize(self, fontsize):
    msg = str(fontsize)
    self.AddToBrowserMsgQueue("SetFontSize", msg)


  def SetVectorWidth(self, vecwidth):
    msg = str(vecwidth)
    self.AddToBrowserMsgQueue("SetVectorWidth", msg)


  def UseWireFrame(self, iswireframe):
    msg = str(iswireframe)
    self.AddToBrowserMsgQueue("UseWireFrame", msg)


  def SetBrowserDebug(self, isdebug):
    msg = str(isdebug)
    self.AddToBrowserMsgQueue("SetBrowserDebug", msg)


  def SetMouseSpeed(self, trackspeed):
    msg = str(trackspeed)
    self.AddToBrowserMsgQueue("SetMouseSpeed", msg)
    #self.GetMouseSpeed() # TODO: fix wait time
    self.mprint("In SetMouseSpeed:\n" + "".join( traceback.format_stack(limit=4) ), verbose="stacktrace"  )


  def GetMouseSpeed(self):
    if self.use_semaphore:
      self.mprint("GetMouseSpeed waiting for mousespeed_msg_sem.acquire", verbose="threadingmsg")
      if not self.mousespeed_msg_sem.acquire(blocking=True, timeout=lock_timeout):
        self.mprint("Error! Timed out waiting for mousespeed_msg_sem semaphore within %s seconds" %lock_timeout, verbose=1)
      self.mprint("GetMouseSpeed got mousespeed_msg_sem", verbose="threadingmsg")
    self.params.NGL.mouse_sensitivity = None
    self.AddToBrowserMsgQueue("GetMouseSpeed", "")


  def SetClipPlaneDistances(self, near, far, cameraPosZ=None, zoom=None):
    if cameraPosZ is None:
      cameraPosZ = self.cameraPosZ
    if zoom is None:
      zoom= self.zoom
    msg = str(near) + ", " + str(far) + ", " + str(cameraPosZ) + ", " + str(zoom)
    self.GuardedAddToBrowserMsgQueue("clipplane_msg_sem", "SetClipPlaneDistances", msg,
                                     funcname="SetClipPlaneDistances", posteriorcheck=False)


  def GetClipPlaneDistances(self):
    if self.use_semaphore:
      self.mprint("GetClipPlaneDistances waiting for hkls_drawn_sem.acquire", verbose="threadingmsg")
      if not self.hkls_drawn_sem.acquire(timeout=lock_timeout):
        self.mprint("Error! Timed out waiting for hkls_drawn_sem semaphore within %s seconds" %lock_timeout, verbose=1)
      self.mprint("GetClipPlaneDistances got hkls_drawn_sem", verbose="threadingmsg")
      self.mprint("GetClipPlaneDistances waiting for autoview_sem.acquire", verbose="threadingmsg")
      if not self.autoview_sem.acquire(blocking=True, timeout=lock_timeout):
        self.mprint("Error! Timed out waiting for autoview_sem semaphore within %s seconds" %lock_timeout, verbose=1)
      self.mprint("GetClipPlaneDistances got autoview_sem", verbose="threadingmsg")

    self.clipNear = None
    self.clipFar = None
    self.cameraPosZ = None
    self.zoom = None # very first call will yield bogus zoom value
    self.GuardedAddToBrowserMsgQueue("clipplane_msg_sem", "GetClipPlaneDistances",
                                     funcname="GetClipPlaneDistances")



  def RemovePrimitives(self, reprname=""):
    self.AddToBrowserMsgQueue("RemovePrimitives", reprname )


  def SetAutoView(self):
    rotmx = self.Euler2RotMatrix( ( 0.0, 0.0, 0.0 ) )
    self.currentRotmx = rotmx
    self.RotateMxStage(rotmx)
    self.GuardedAddToBrowserMsgQueue("autoview_sem", "SetAutoView", funcname="SetAutoView")


  def SetDefaultOrientation(self):
    if self.params.clip_plane.clip_width and not self.isnewfile:
      # if self.params.clip_plane.clip_width
      # then we are clipping and using camerazoom instead of camera.position.z
      # Autoview used by SetDefaultOrientation will mess that up. So bail out.
      return
    self.GuardedAddToBrowserMsgQueue("autoview_sem", "SetDefaultOrientation",
                                     funcname="SetDefaultOrientation")


  def TestNewFunction(self):
    self.AddToBrowserMsgQueue("Testing")


  def MakeImage(self, filename):
    self.imagename = filename
    if self.use_semaphore:
      self.mprint("MakeImage waiting for hkls_drawn_sem semaphore", verbose="threadingmsg")
      if not self.hkls_drawn_sem.acquire(blocking=True, timeout=lock_timeout):
        self.mprint("MakeImage failed acquiring hkls_drawn_sem semaphore within %s seconds" %lock_timeout, verbose=1)
      self.mprint("MakeImage got hkls_drawn_sem semaphore", verbose="threadingmsg")
    self.AddToBrowserMsgQueue("MakeImage2", "HKLviewer.png,"+ str(sys.version_info[0]) )


  def DisableMouseRotation(self): # disable rotating with the mouse
    self.AddToBrowserMsgQueue("DisableMouseRotation")


  def EnableMouseRotation(self): # enable rotating with the mouse
    self.AddToBrowserMsgQueue("EnableMouseRotation")


  def UseCameraZoom(self): # disable zoom with the mouse use webgl camera zoom instead
    self.AddToBrowserMsgQueue("DisableZoomDrag")


  def UseZoomDrag(self): # enable zoom with the mouse
    self.AddToBrowserMsgQueue("EnableZoomDrag")


  def SimulateClick(self):
    self.AddToBrowserMsgQueue("SimulateClick")


  def ReOrientStage(self):
    if self.viewmtrx:
      self.AddToBrowserMsgQueue("SetAutoView", self.viewmtrx)


  def Euler2RotMatrix(self, eulerangles):
    eulerangles1 = eulerangles
    radangles = [e*math.pi/180.0 for e in eulerangles1]
    # scitbx is using ZYZ convention for euler angles
    # https://en.wikipedia.org/wiki/Euler_angles#Rotation_matrix
    RotMx = scitbx.math.euler_angles_as_matrix(radangles)
    return RotMx


  def RotateMxStage(self, rotmx, quietbrowser=True):
    if self.cameraPosZ is None:
      # in case HKLJavaScripts.onMessage() crashed and failed returning cameraPosZ
      self.GetClipPlaneDistances()
    if self.cameraPosZ is not None:
      if self.use_semaphore:
        self.mprint("RotateMxStage waiting for clipplane_msg_sem.acquire", verbose="threadingmsg")
        if not self.clipplane_msg_sem.acquire(blocking=True, timeout=lock_timeout):
          self.mprint("Error! Timed out waiting for clipplane_msg_sem semaphore within %s seconds" %lock_timeout, verbose=1)
        self.mprint("RotateMxStage got clipplane_msg_sem", verbose="threadingmsg")
      scaleRot = rotmx * self.cameraPosZ
      ortrot = scaleRot.as_mat3()
      str_rot = str(ortrot)
      str_rot = str_rot.replace("(", "")
      str_rot = str_rot.replace(")", "")
      str_rot = str_rot + ", " + str(self.zoom)
      msg = str_rot + ", quiet\n"
      if not quietbrowser:
        msg = msg + ", verbose\n"

      self.AddToBrowserMsgQueue("RotateStage", msg)


  def RotateAxisMx(self, vec, theta, quietbrowser=True):
    if self.cameraPosZ is None:
      return
    str_rot = str(list(vec)) + ", " + str(theta)
    str_rot = str_rot.replace("[", "")
    str_rot = str_rot.replace("]", "")
    msg = str_rot + ", quiet\n"
    if not quietbrowser:
      msg = str_rot + ", verbose\n"
    self.AddToBrowserMsgQueue("RotateAxisStage", msg)


  def RotateMxComponents(self, rotmx, quietbrowser=True):
    if self.cameraPosZ is None:
      return
    #scaleRot = rotmx * self.cameraPosZ
    ortrot = rotmx.as_mat3()
    str_rot = str(ortrot)
    str_rot = str_rot.replace("(", "")
    str_rot = str_rot.replace(")", "")
    msg = str_rot + ", quiet\n"
    if not quietbrowser:
      msg = str_rot + ", verbose\n"
    self.AddToBrowserMsgQueue("RotateComponents", msg)


  def RotateAxisComponents(self, vec, theta, quietbrowser=True):
    if self.cameraPosZ is None:
      return
    str_rot = str(list(vec)) + ", " + str(theta)
    str_rot = str_rot.replace("[", "")
    str_rot = str_rot.replace("]", "")
    msg = str_rot + ", quiet\n"
    if not quietbrowser:
      msg = str_rot + ", verbose\n"
    self.AddToBrowserMsgQueue("RotateAxisComponents", msg)


  def AnimateRotateAxisComponents(self, vec, speed, quietbrowser=True):
    if self.cameraPosZ is None:
      return
    str_rot = str(list(vec)) + ", " + str(speed)
    str_rot = str_rot.replace("[", "")
    str_rot = str_rot.replace("]", "")
    msg = str_rot + ", quiet\n"
    if not quietbrowser:
      msg = str_rot + ", verbose\n"
    self.AddToBrowserMsgQueue("AnimateRotateAxisComponents", msg)


  def RemoveStageObjects(self):
    self.AddToBrowserMsgQueue("RemoveStageObjects")


  def DefineHKL_Axes(self, Hstararrowstart, Hstararrowend, Kstararrowstart,
                     Kstararrowend, Lstararrowstart, Lstararrowend,
                     Hlabelpos, Klabelpos, Llabelpos):
    strdata = ""
    strdata += "%s\n\n" %str(Hstararrowstart)
    strdata += "%s\n\n" %str(Hstararrowend)
    strdata += "%s\n\n" %str(Kstararrowstart)
    strdata += "%s\n\n" %str(Kstararrowend)
    strdata += "%s\n\n" %str(Lstararrowstart)
    strdata += "%s\n\n" %str(Lstararrowend)
    strdata += "%s\n\n" %str(Hlabelpos)
    strdata += "%s\n\n" %str(Klabelpos)
    strdata += "%s\n\n" %str(Llabelpos)
    self.AddToBrowserMsgQueue("DefineHKL_Axes", strdata)


  def SendCoordinates2Browser(self, positions, colours, radii, ttipids ):
    # send data in binary format rather than string to browser for the sake of speed and
    # having to avoid rounding off numbers which is also slow
    self.AddToBrowserMsgQueue("AddHKLCoordinates", positions, binary=True)
    self.mprint(".", end="")
    self.AddToBrowserMsgQueue("AddHKLColours", colours, binary=True)
    self.mprint(".", end="")
    self.AddToBrowserMsgQueue("AddHKLRadii", radii, binary=True)
    self.mprint(".", end="")
    self.AddToBrowserMsgQueue("AddHKLTTipIds", ttipids, binary=True)


  def RenderStageObjects(self):
    self.GuardedAddToBrowserMsgQueue("hkls_drawn_sem", "RenderStageObjects", funcname="RenderStageObjects")


  def MakeColourChart(self, label, fomlabel, colourgradarray):
    msg = "%s\n\n%s\n\n%s" %(label, fomlabel, str(colourgradarray) )
    self.AddToBrowserMsgQueue("MakeColourChart", msg )


  def get_current_datatype(self):
    # Amplitudes, Map coeffs, weights, floating points, etc
    if self.params.viewer.scene_id is None:
      return None
    dtype = self.array_info_format_tpl[ self.scene_id_to_array_id(self.params.viewer.scene_id )][1][1]
    # if dtype is boring generic then use the name of the data column for dtype
    if dtype in ["Floating-point", "Integer"]:
      dtype = self.array_info_format_tpl[ self.scene_id_to_array_id(self.params.viewer.scene_id )][1][0]
    return dtype


  def onClickColourChart(self):
    # if running the GUI show the colour chart selection dialog
    self.SendInfoToGUI( { "ColourChart": self.params.hkls.color_scheme,
                          "ColourPowerScale": self.params.hkls.color_powscale,
                          "CurrentDatatype": self.get_current_datatype(),
                          "ShowColourMapDialog": 1
                         } )


  def MakeBrowserDataColumnComboBox(self):
    datcolstr = ""
    for i,lbl in enumerate(self.hkl_scenes_infos):
      datcolstr = datcolstr + ",".join(lbl[3]) + "\n" + str(i)
      if i < len(self.hkl_scenes_infos)-1:
        datcolstr = datcolstr + "\n\n"
    datcolstr = datcolstr + "\n\n" + str(self.params.viewer.scene_id)
    self.AddToBrowserMsgQueue("MakeBrowserDataColumnComboBox", datcolstr)




ngl_philstr = """
  mouse_sensitivity = 0.06
    .type = float
    .help = "Controls the speed of movement when adjusting view with the mouse"
  tooltip_alpha = 0.80
    .type = float
    .help = "Opacity of tooltips showing data values of reflections when clicking or hovering the mouse on reflections"
  vector_width = 5
    .type = int(value_min=1, value_max=30)
    .help = "Thickness of vectors and axes"
  fontsize = 9
    .type = int
    .help = "Font size for browser window displaying reflections"
  background_colour = 'rgb(127, 127, 127)'
    .type = str
    .help = "String of RGB colour values for the background of the browser"
  show_tooltips = none *click hover
    .type = choice
    .help = "Specifies whether tooltips for reflections should show by hovering or by clicking on a reflection" \
            "If the displayed data has a very large number of reflections it is best to select "click""
  camera_type = *orthographic perspective
    .type = choice
"""

NGLmaster_phil = libtbx.phil.parse( ngl_philstr )
NGLparams = NGLmaster_phil.fetch().extract()

def reset_NGLsettings():
  """
  Reset NGL settings to their default values as specified in the phil definition string
  """
  #global NGLmaster_phil
  #global ngl_philstr
  #global NGLparams
  NGLparams = NGLmaster_phil.fetch(source = libtbx.phil.parse( ngl_philstr) ).extract()


def NGLsettings():
  """
  Get a global phil parameters object containing some NGL settings
  """
  #global NGLparams
  return NGLparams









"""
# python2 code

from websocket_server import WebsocketServer
import threading, math
from time import sleep

nc = {}
def new_client(client, server):
  nc = client
  print "got a new client:", nc

def on_message(client, server, message):
    print message

#websocket.enableTrace(True)
server = WebsocketServer(7894, host='127.0.0.1')
server.set_fn_new_client(new_client)
server.set_fn_message_received(on_message)

wst = threading.Thread(target=server.run_forever)
wst.daemon = True
wst.start()

def LoopSendMessages():
  x = 0.0
  i=0
  while server.clients:
    nc = server.clients[0]
    x += 0.2
    alpha =  (math.cos(x) +1.0 )/2.0
    msg = u"alpha, 2, %f" %alpha
    server.send_message(server.clients[0], msg )
    r = (math.cos(x) +1.0 )/2.0
    g = (math.cos(x+1) +1.0 )/2.0
    b = (math.cos(x+2) +1.0 )/2.0
    msg = u"colour, 1, %d, %f, %f, %f" %(i,r,g,b)
    server.send_message(server.clients[0], msg )
    sleep(0.2)



"""


"""

# python3 code

# WS server example
import asyncio
import websockets

async def hello(websocket, path):
  while True:
    name = await websocket.recv()
    print(f"< {name}")
    greeting = f"Hello {name}!"
    await websocket.send(greeting)
    if name=="STOP":
      return
    await asyncio.sleep(0.2)

start_server = websockets.serve(hello, "localhost", 8765)
asyncio.get_event_loop().run_until_complete(start_server)
asyncio.get_event_loop().run_forever()



# WS client example
import asyncio
import websockets

async def hello():
  uri = "ws://localhost:8765"
  async with websockets.connect(uri) as websocket:
    while True:
      name = input("What's your name?\n" )
      await websocket.send(name)
      print(f"> {name}")
      greeting = await websocket.recv()
      print(f"< {greeting}")

asyncio.get_event_loop().run_until_complete(hello())

"""
