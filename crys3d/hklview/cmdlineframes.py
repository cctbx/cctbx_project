
from __future__ import absolute_import, division, print_function

# TODO:
#  - prompt user for missing symmetry
#  - cached scenes



r"""
Load mtz file for viewing reflections in a webbrowser using NGL.
Webbrowser input is a javascript file that together with ngl.js displays reflections
positioned in reciprocal space using a webbrowser. Usage:


from crys3d.hklview import cmdlineframes
myHKLview = cmdlineframes.HKLViewFrame(jscriptfname = r"C:\Users\oeffner\Buser\NGL_HKLviewer\myjstr.js", htmlfname = "C:\Users\oeffner\Buser\NGL_HKLviewer\myhkl.html")
myHKLview.LoadReflectionsFile(r"C:\Users\oeffner\Buser\Work\ANI_TNCS\4PA9\4pa9.tncs.mtz")
myHKLview.GetArrayInfo()
myHKLview.SetScene(0)

myHKLview.SetRadiiScale(3, nth_power_scale=1)
myHKLview.SetColourScene(3)
myHKLview.SetRadiusScene(7)

myHKLview = cmdlineframes.HKLViewFrame(jscriptfname = "myjstr.js", htmlfname = "myhkl.html", UseOSBrowser= False)

myHKLview.SetScene(3)
myHKLview.SetRadiiScale(3, nth_root_scale=1)
myHKLview.SetSceneBinThresholds([0, 3.4, 3.5, 3.6, 3.8, 4.0, 4.5, 5.0, 6.0, 8, 9,11, 15, 20, 30, 100 ], "Resolution")
for i in range(13):
  myHKLview.SetOpacity(i, 0.0)

myHKLview.SetSceneBinThresholds([0, 0.1, 1, 10], 3)
myHKLview.ExpandToP1(True)
myHKLview.SetOpacity(2, 0.1)

myHKLview.ShowSlice(True, "l", 17)
myHKLview.GetNGLstring()
myHKLview.ShowSlice(False)
myHKLview.ExpandAnomalous(True)
myHKLview.ShowMissing(True)
myHKLview.GetSpaceGroupChoices()
myHKLview.SetSpaceGroupChoice(3)

myHKLview.LoadReflectionsFile(r"C:\Users\oeffner\Buser\Phenix\phenix-installer-dev-3484-win7vc90\modules\phenix_examples\lysozyme-MRSAD\lyso2001_scala1.mtz")
myHKLview.SetCameraType("persp")

myHKLview.LoadReflectionsFile(r"C:\Users\oeffner\Buser\Experiments\CRLF3\DLS20151206CRLF3\5840-F11-X1-Hg-SAD-ONsoak\5840-F11-X1_pk_5_5_1_\xia2\dials-run\DataFiles\mx11235v49_x5840F11X1pk551_free.mtz")

myHKLview.SetScene(0)
myHKLview.SetRadiiScale(1, nth_power_scale=0.2)
myHKLview.SetScene(1)
myHKLview.SetSceneBinThresholds(1, [-20, 30, 300, 3000])


myHKLview.ShowSlice(True, "h", 20)

# P1 crystal
myHKLview.LoadReflectionsFile(r"C:\Users\oeffner\Buser\Tests\P1xtals\5aws-sf.mtz")
myHKLview.SetScene(0)



myHKLview.LoadReflectionsFile(r"C:\Users\oeffner\Buser\Phenix\phenix-installer-dev-3484-win7vc90\modules\phenix_examples\beta-blip\beta_blip_P3221.mtz")
myHKLview.SetScene(0)
myHKLview.SetSceneBinThresholds([3, 3.2, 3.4, 3.8, 4.5, 6, 9, 30])
for a in range(6):
  myHKLview.SetOpacity(a, 0.0)



# Create a small mtz file with 10 reflections and 3 miller arrays. Then load into NGL_HKL viewer

from cctbx.xray import observation_types
from cctbx.array_family import flex
from cctbx import miller
from cctbx import crystal

xs = crystal.symmetry(unit_cell=(50,50,40, 90,90,120), space_group_symbol="P3 1")
mi = flex.miller_index([ (1,-2,3), (0,0,-4), (1, 2, 3), (0, 1, 2),
                        (1, 0, 2), (-1, 1, -2), (2, -2, -2),
                        (-2, 1, 0) , (1, 0, -2), (0, 0, 2) ]
)

ma = miller.array( miller.set(xs, mi) )

ma1 = miller.array( miller.set(xs, mi), flex.double( [11.205, 6.353, 26.167, 14.94, 2.42, 24.921, 16.185, 11.798, 21.183, 4.98] ),
                   sigmas=flex.double( [13.695, 6.353, 24.921, 6.225, 11.193, 26.167, 8.715, 4.538, 27.413, 21.165] )
                   ).set_observation_type( observation_types.intensity() )
ma1.set_info(miller.array_info(source="artificial file", labels=["MyI", "SigMyI"]))

mi2 = flex.miller_index([ (1,-2,3), (0,0,-4), (1, 2, 3), (0, 1, 2),
                        (1, 0, 2), (-1, 1, -2), (2, -2, -2),
                        (-2, 1, 0) , (0, 0, 2),  ]
                        )

ma2 = miller.array(miller.set(xs, mi2),
                    flex.complex_double( [
                                -1.0 + 0.0j,
                                 -1.5 + 2.598075j,
                                 0.0 + 0.8j,
                                 1.0 + 1.7320508j,
                                 4.0 + 0.0j,
                                 0.5 - 0.866025j,
                                 0.0 - 1.0j,
                                 -2.5 - 4.330127j,
                                 -4.24264 + 4.24264j
                                 ] ) )

ma2.set_info(miller.array_info(source="artificial file", labels=["MyMap", "PhiMyMap"]))

mafom = miller.array(miller.set(xs, mi2),
                    flex.double( [0.0, 0.1, 0.25, 0.35, 0.4, 0.5, 0.6, 0.75, 1.0 ] ) )
mafom.set_info(miller.array_info(source="artificial file", labels=["FOM"]))

mi3 = flex.miller_index([ (1, 0, 2), (-1, 1, -2), (-2, 1, 0) , (1, 0, -2), (1,-2,3), (0,0,-4), (0, 1, 2), (0, 0, 2) ] )
ma3 = miller.array(miller.set(xs, mi3), data=flex.double( [22.429, 28.635, 3.328, 3.738, 24.9, 14.521, 3.738, 19.92] ) )
ma3.set_info(miller.array_info(source="artificial file", labels=["Foo"]))

mi4 = flex.miller_index([ (1,-2,3), (0,0,-4), (1, 2, 3), (0, 1, 2), (1, 0, 2), (-1, 1, -2),   (0, 0, 2) ] )
ma4 = miller.array(miller.set(xs, mi4), data=flex.std_string( ["foo", "bar", "wibble", "waffle", "mumble", "muffle","babble"] ) )
ma4.set_info(miller.array_info(source="artificial file", labels=["Bar"]))

mi5 = flex.miller_index([(1, -2, 3), (0, 0, -3), (1, 2, 3), (0, 1, 2),
                         (1, 0, 2), (2, -2, -2), (-2, 1, 0), (0, 0, 2)]
)
ma5 = miller.array(miller.set(xs, mi5), data=flex.double( [12.429, 38.635, -3.328, 3.738, 4.9, -5.521, 10.738, 19.92] ) )
ma5.set_info(miller.array_info(source="artificial file", labels=["BarFoo"]))



mtz1 = ma1.as_mtz_dataset(column_root_label="I")
mtz1.add_miller_array(ma2, column_root_label="MyMap")
mtz1.add_miller_array(ma3, column_root_label="Oink")
#mtz1.add_miller_array(ma4, column_root_label="blip")
mtz1.add_miller_array(ma5, column_root_label="bleep")
mtz1.add_miller_array(mafom, column_root_label="FOM")
mtz1.set_wavelength(1.2)
mtz1.set_name("MyTestData")
mtz1.mtz_object().write("mymtz.mtz")


from crys3d.hklview import cmdlineframes
myHKLview = cmdlineframes.HKLViewFrame(jscriptfname = "myjstr.js", high_quality=True, verbose=1)
myHKLview.LoadReflectionsFile("mymtz.mtz")
myHKLview.SetScene(0)

myHKLview.SetRadiiScale(1, nth_power_scale=0.2)
myHKLview.SetSceneBinThresholds([50, 15, 9])
myHKLview.SetOpacity(0, 0.0)
myHKLview.SetOpacity(1, 0.0)
myHKLview.SetOpacity(2, 0.0)


myHKLview.ExpandToP1(True)
myHKLview.ExpandAnomalous(True)
myHKLview.ShowMissing(True)
myHKLview.ShowSystematicAbsences(True)

myHKLview.SetColourScene(2)
myHKLview.SetRadiusScene(3)
myHKLview.SetScene(1)

myHKLview.SetScene(0)
myHKLview.SetRadiiToSigmas(True)

from crys3d.hklview import cmdlineframes
myHKLview = cmdlineframes.HKLViewFrame(jscriptfname = "myjstr.js", htmlfname = "myhkl.html")
myHKLview.LoadReflectionsFile(r"C:\Users\oeffner\Buser\Tests\MRproblem\MRproblem_15.1.mtz")

myHKLview.SetScene(2)
myHKLview.SetColoursToPhases(True)


from crys3d.hklview import cmdlineframes
myHKLview = cmdlineframes.HKLViewFrame(jscriptfname = "myjstr.js")
myHKLview.LoadReflectionsFile("3RP2_A.1.mtz")
myHKLview.SetScene(0)
myHKLview.ClipPlaneParallelToHKLplane(0, 3, 0)
myHKLview.ClipPlaneParallelToHKLplane(0, 3, 3, clipNear=-1, clipFar=1, fixorientation=False)
myHKLview.SetRadiiScale(1, nth_power_scale=0.0)
myHKLview.ExpandToP1(True)
myHKLview.ExpandAnomalous(True)
myHKLview.ShowSlice(True, "l", 25)
myHKLview.ShowMissing(True)


myHKLview.LoadReflectionsFile(r"C:\Users\oeffner\Buser\Work\TNCS\4N3E\4n3e_final.mtz")


myHKLview.LoadReflectionsFile("r4v9hsf.cif")
myHKLview.SetScene(0)

myHKLview.SetSceneBinThresholds([0,2.8, 2.86, 2.88, 2.898, 3.0, 3.03, 3.06, 3.1, 3.15, 3.2, 3.25, 3.3, 3.35, 3.4, 3.5, 3.6, 3.8, 4.0, 4.5, 5.0, 6.0, 8, 9,11, 15, 20, 30, 100 ])
for i in range(27):
  myHKLview.SetOpacity(i, 0.0)

myHKLview.SetOpacity(26, 0.0)
myHKLview.SetOpacity(26, 1.0)
myHKLview.SetOpacity(25, 1.0)
myHKLview.SetOpacity(24, 1.0)
myHKLview.SetOpacity(10, 1.0)
myHKLview.SetOpacity(14, 1.0)



from crys3d.hklview import cmdlineframes
myHKLview = cmdlineframes.HKLViewFrame(jscriptfname = "myjstr.js", htmlfname = "myhkl.html", high_quality=True, verbose=1)
myHKLview.LoadReflectionsFile(r"C:\Users\oeffner\Buser\Phenix\phenix-installer-dev-3484-win7vc90\modules\phenix_examples\beta-blip\beta_blip_P3221.mtz")
myHKLview.SetScene(0)
myHKLview.SetRadiiScale(1,0)
myHKLview.ExpandToP1(True)
myHKLview.ExpandAnomalous(True)
myHKLview.ClipPlaneParallelToHKLplane(0, 3, 3, hkldist=0)
myHKLview.RemoveNormalVectorToClipPlane()



from crys3d.hklview import cmdlineframes
myHKLview = cmdlineframes.HKLViewFrame(jscriptfname = "myjstr.js", verbose=1)
myHKLview.LoadReflectionsFile("4XVY.tncs.mtz")
myHKLview.SetScene(0)
myHKLview.SetSceneNbins(5)




"""



from cctbx.miller import display2 as display
from crys3d.hklview import jsview_3d as view_3d
from crys3d.hklview.jsview_3d import ArrayInfo
from libtbx.str_utils import format_value
from cctbx.array_family import flex
from libtbx.utils import Sorry, to_str
from libtbx import group_args
import libtbx
from cctbx import miller
import traceback
import sys, zmq, threading,  time

from six.moves import input


argn = 1
argc = len(sys.argv)

# prompt user for value if it's not on the commandline
def Inputarg(varname):
  global argn
  global argc
  if argc > 1 and argn < argc:
    myvar = sys.argv[argn]
    argn = argn + 1
    print(varname + " " + myvar)
  else:
    myvar = input(varname)
  return myvar

NOREFLDATA = "No reflection data has been selected"


class settings_window () :
  def set_index_span (self, index_span) :
    self._index_span = index_span


  def update_reflection_info (self, hkl, d_min, value) :
    print(hkl, value)
    if (hkl is None) :
      self.hkl_info.SetValue("")
      self.d_min_info.SetValue("")
      self.value_info.SetValue("")
    else :
      self.hkl_info.SetValue("%d, %d, %d" % hkl)
      d_min_str = format_value("%.3g", d_min)
      self.d_min_info.SetValue(d_min_str)
      value_str = format_value("%.3g", value, replace_none_with="---")
      self.value_info.SetValue(value_str)


  def clear_reflection_info (self) :
    self.update_reflection_info(None, None, None)


class HKLViewFrame() :
  def __init__ (self, *args, **kwds) :
    self.valid_arrays = []
    self.spacegroup_choices = []
    self.procarrays = []
    self.merge_answer = [None]
    self.dmin = -1
    self.settings = display.settings()
    self.verbose = 0
    if 'verbose' in kwds:
      self.verbose = kwds['verbose']
    kwds['settings'] = self.settings
    kwds['mprint'] = self.mprint
    self.infostr = ""
    self.guiSocketPort=None
    self.new_miller_array_operations_lst = []
    if 'useGuiSocket' in kwds:
      self.guiSocketPort = kwds['useGuiSocket']
      self.context = zmq.Context()
      self.guisocket = self.context.socket(zmq.PAIR)
      self.guisocket.connect("tcp://127.0.0.1:%s" %self.guiSocketPort )
      self.STOP = False
      print("starting socketthread")
      self.msgqueuethrd = threading.Thread(target = self.zmq_listen )
      self.msgqueuethrd.daemon = True
      self.msgqueuethrd.start()
      kwds['guisocket'] = self.guisocket
    kwds['websockport'] = self.find_free_port()
    self.viewer = view_3d.hklview_3d( **kwds )
    self.ResetPhilandViewer()
    self.NewFileLoaded = False


  def mprint(self, msg, verbose=0):
    if verbose <= self.verbose:
      if self.guiSocketPort:
        self.SendInfoToGUI( { "info": msg } )
      else:
        print(msg)


  def find_free_port(self):
    import socket
    s = socket.socket()
    s.bind(('', 0))      # Bind to a free port provided by the host.
    port = s.getsockname()[1]
    s.close()
    return port


  def zmq_listen(self):
    while not self.STOP:
      try:
        philstr = self.guisocket.recv()
        philstr = str(philstr)
        self.mprint("Received phil string:\n" + philstr, verbose=1)
        new_phil = libtbx.phil.parse(philstr)
        self.update_settings(new_phil)
        #print( "in zmq listen")
        time.sleep(0.1)
      except Exception as e:
        self.mprint( str(e) + traceback.format_exc(limit=10), verbose=1)

    del self.guisocket
    self.mprint( "Exiting zmq_listen() thread", 1)


  def ResetPhilandViewer(self, extraphil=None):
    self.master_phil = libtbx.phil.parse( masterphilstr )
    self.currentphil = self.master_phil
    if extraphil:
      self.currentphil = self.currentphil.fetch(source = extraphil)
    self.params = self.currentphil.fetch().extract()
    self.viewer.viewerparams = self.params.NGL_HKLviewer.viewer
    self.viewer.symops = []
    self.viewer.sg = None
    self.viewer.proc_arrays = []
    self.viewer.HKLscenesdict = {}
    self.viewer.sceneisdirty = True
    if self.viewer.miller_array:
      self.viewer.DrawNGLJavaScript( blankscene=True)
    self.viewer.miller_array = None


  def GetNewCurrentPhilFromString(self, philstr, oldcurrentphil):
    user_phil = libtbx.phil.parse(philstr)
    newcurrentphil = oldcurrentphil.fetch(source = user_phil)
    diffphil = oldcurrentphil.fetch_diff(source = user_phil)
    return newcurrentphil, diffphil


  def GetNewCurrentPhilFromPython(self, pyphilobj, oldcurrentphil):
    newcurrentphil = oldcurrentphil.fetch(source = pyphilobj)
    diffphil = oldcurrentphil.fetch_diff(source = pyphilobj)

    oldcolbintrshld = oldcurrentphil.extract().NGL_HKLviewer.scene_bin_thresholds
    newcolbintrshld = oldcolbintrshld
    if hasattr(pyphilobj.extract().NGL_HKLviewer, "scene_bin_thresholds"):
      newcolbintrshld = pyphilobj.extract().NGL_HKLviewer.scene_bin_thresholds
    # fetch_diff doesn't seem able to correclty spot changes
    # in the multiple scope phil object "NGL_HKLviewer.scene_bin_thresholds"
    # Must do it manually
    params = newcurrentphil.extract()
    if oldcolbintrshld != newcolbintrshld: # or old_binopacities != new_binopacities:
      params.NGL_HKLviewer.scene_bin_thresholds = newcolbintrshld
      newcurrentphil = self.master_phil.format(python_object = params)
      diffphil = self.master_phil.fetch_diff(source = newcurrentphil)
    #import code, traceback; code.interact(local=locals(), banner="".join( traceback.format_stack(limit=10) ) )
    return newcurrentphil, diffphil


  def SetCurrentPhilAsPython(self, pyphil):
    newphil = master_phil.format(python_object= pyphil)
    currphil = master_phil.fetch(source = newphil)


  def update_settings(self, new_phil=None):
    try:
      if not new_phil:
        new_phil = self.master_phil.format(python_object = self.params)
      #import code, traceback; code.interact(local=locals(), banner="".join( traceback.format_stack(limit=10) ) )
      self.currentphil, diff_phil = self.GetNewCurrentPhilFromPython(new_phil, self.currentphil)
      diff = None
      if len(diff_phil.all_definitions()) < 1:
        self.mprint( "Nothing's changed")
        return False

      diff = diff_phil.extract().NGL_HKLviewer
      self.mprint("diff phil:\n" + diff_phil.as_str(), verbose=1 )

      self.params = self.currentphil.extract()
      phl = self.params.NGL_HKLviewer

      if view_3d.has_phil_path(diff_phil, "filename"):
        self.ResetPhilandViewer(diff_phil)
        if not self.load_reflections_file(phl.filename):
          return False

      if view_3d.has_phil_path(diff_phil, "scene_id") \
       or view_3d.has_phil_path(diff_phil, "merge_data") \
       or view_3d.has_phil_path(diff_phil, "show_missing") \
       or view_3d.has_phil_path(diff_phil, "show_only_missing") \
       or view_3d.has_phil_path(diff_phil, "show_systematic_absences") \
       or view_3d.has_phil_path(diff_phil, "nbins") \
       or view_3d.has_phil_path(diff_phil, "bin_scene_label") \
       or view_3d.has_phil_path(diff_phil, "scene_bin_thresholds"):
        #import code; code.interact(local=locals(), banner="".join( traceback.format_stack(limit=10) ) )
        if self.set_scene(phl.viewer.scene_id):
          self.set_scene_bin_thresholds(binvals=phl.scene_bin_thresholds,
                                         bin_scene_label=phl.bin_scene_label,
                                         nbins=phl.nbins )

      if view_3d.has_phil_path(diff_phil, "spacegroup_choice"):
        self.set_spacegroup_choice(phl.spacegroup_choice)

      if view_3d.has_phil_path(diff_phil, "miller_array_operations"):
        self.make_new_miller_array()

      if view_3d.has_phil_path(diff_phil, "using_space_subgroup") and phl.using_space_subgroup==False:
        self.set_default_spacegroup()

      if view_3d.has_phil_path(diff_phil, "shape_primitive"):
        self.set_shape_primitive(phl.shape_primitive)

      if view_3d.has_phil_path(diff_phil, "tooltips_in_script"):
        self.viewer.script_has_tooltips = phl.tooltips_in_script

      if view_3d.has_phil_path(diff_phil, "viewer"):
        self.viewer.settings = phl.viewer
        self.settings = phl.viewer

      msg, self.params.NGL_HKLviewer = self.viewer.update_settings(diff_phil, phl)
      # parameters might have been changed. So update self.currentphil accordingly
      self.currentphil = self.master_phil.format(python_object = self.params)
      self.mprint( msg, verbose=1)
      #import code, traceback; code.interact(local=locals(), banner="".join( traceback.format_stack(limit=10) ) )
      self.NewFileLoaded = False
      if (self.viewer.miller_array is None) :
        self.mprint( NOREFLDATA, True)
        return False
      return True
    except Exception as e:
      self.mprint(to_str(e) + "\n" + traceback.format_exc(), 0)
      return False


  def update_clicked (self, index) :#hkl, d_min=None, value=None) :
    if (index is None) :
      self.settings_panel.clear_reflection_info()
    else :
      hkl, d_min, value = self.viewer.scene.get_reflection_info(index)
      self.settings_panel.update_reflection_info(hkl, d_min, value)


  def detect_Rfree(self, array):
    from iotbx.reflection_file_utils import looks_like_r_free_flags_info
    info = array.info()
    if (array.is_integer_array()) and (looks_like_r_free_flags_info(info)) :
      from iotbx.reflection_file_utils import get_r_free_flags_scores
      score_array = get_r_free_flags_scores([array], None)
      test_flag_value = score_array.test_flag_values[0]
      array = array.customized_copy(data=(array.data() == test_flag_value))
      array.set_info(info)
    return array


  def process_miller_array(self, array, merge_answer=[None]) :
    if (array is None) : return
    if (array.is_hendrickson_lattman_array()) :
      raise Sorry("Hendrickson-Lattman coefficients are currently not supported.")
    info = array.info()
    if isinstance(info, str) :
      labels = "TEST DATA"
    else :
      labels = info.label_string()
    if (array.unit_cell() is None) or (array.space_group() is None) :
      raise Sorry("No space group info is present in data")
    details = []
    self.infostr = ""
    if (not array.is_unique_set_under_symmetry() and self.params.NGL_HKLviewer.merge_data is None):
      shouldmergestr = "The data in the selected array are not symmetry-" + \
        "unique, which usually means they are unmerged (but could also be due "+ \
        "to different indexing conventions).  Do you want to merge equivalent "+ \
        "observations (preserving anomalous data if present), or view the "+ \
        "array unmodified?  (Note that if you do not merge the array, the "+ \
        "options to expand to P1 or generate Friedel pairs will be be disabled"+ \
        ", and the 2D view will only show indices present in the file, rather "+ \
        "than a full pseudo-precession view.). yes/no?\n"
      if self.useSocket:
        self.infostr = shouldmergestr
        while 1:
          philstr = self.guisocket.recv()
          philstr = str(philstr)
          self.mprint("Received phil:\n" + philstr, verbose=2)
          new_phil = libtbx.phil.parse(philstr)
          #working_phil = self.master_phil.fetch(source = new_phil)
          params = new_phil.extract().NGL_HKLviewer
          if hasattr(params, "merge_data"): # awaiting user to tick a checkbox on the gui
            self.params.NGL_HKLviewer.mergedata = params.merge_data
            break
          time.sleep(0.1)
      else:
        self.params.NGL_HKLviewer.merge_data = Inputarg(shouldmergestr).lower()[0] == "y"
        if self.params.NGL_HKLviewer.merge_data:
          details.append("merged")
        else :
          details.append("unmerged data")
          self.settings.expand_to_p1 = False
          self.settings.expand_anomalous = False
    array = self.detect_Rfree(array)
    sg = "%s" % array.space_group_info()
    uc = "a=%g b=%g c=%g angles=%g,%g,%g" % array.unit_cell().parameters()
    details_str = ""
    if (len(details) > 0) :
      details_str = "(%s)" % ", ".join(details)
    array_info = group_args(
      labels=labels,
      details_str=details_str,
      merge=self.params.NGL_HKLviewer.merge_data,
      sg=sg,
      uc=uc)
    return array, array_info


  def process_all_miller_arrays(self, col):
    #print "in process_all_miller_arrays"
    self.mprint("Processing reflection data...")
    self.procarrays = []
    if self.params.NGL_HKLviewer.merge_data == False:
      self.settings.expand_to_p1 = False
      self.settings.expand_anomalous = False
    for c,arr in enumerate(self.valid_arrays):
      procarray, procarray_info = self.process_miller_array(arr,
                                            merge_answer=self.merge_answer)
      self.procarrays.append(procarray)
      if c==col:
        array_info = procarray_info
        self.viewer.miller_array = procarray
    if col is None:
      array_info = procarray_info
    self.merge_answer = [None]
    #import code, traceback; code.interact(local=locals(), banner="".join( traceback.format_stack(limit=10) ) )
    #self.viewer.set_miller_array(col, merge=array_info.merge,
    #   details=array_info.details_str, proc_arrays=self.procarrays)
    #self.viewer.identify_suitable_fomsarrays()
    return array_info


  def set_miller_array(self, col=None) :
    #print "in set_miller_array"
    #if col >= len(self.valid_arrays):
    if col is not None and col >= len(self.viewer.hkl_scenes_info ):
      return
    array_info = self.process_all_miller_arrays(col)
    self.viewer.set_miller_array(col, merge=array_info.merge,
       details=array_info.details_str)
    self.viewer.proc_arrays = self.procarrays
    self.viewer.identify_suitable_fomsarrays()
    self.update_space_group_choices()


  def update_space_group_choices(self) :
    if self.viewer.miller_array is None or \
      self.params.NGL_HKLviewer.using_space_subgroup:
      return
    current_miller_array_idx = self.viewer.hkl_scenes_info[self.params.NGL_HKLviewer.viewer.scene_id][1]
    #matching_valid_array = self.valid_arrays[ current_miller_array_idx ]
    matching_valid_array = self.procarrays[ current_miller_array_idx ]
    from cctbx.sgtbx.subgroups import subgroups
    from cctbx import sgtbx
    sg_info = matching_valid_array.space_group_info()
    subgrs = subgroups(sg_info).groups_parent_setting()
    self.spacegroup_choices = []
    for i,subgroup in enumerate(subgrs) :
      subgroup_info = sgtbx.space_group_info(group=subgroup)
      self.spacegroup_choices.append(subgroup_info)
    if (sg_info in self.spacegroup_choices) :
      self.current_spacegroup = self.spacegroup_choices.index(sg_info)
    else :
      self.spacegroup_choices.insert(0, sg_info)
      self.current_spacegroup = sg_info
    mydict = { "spacegroups": [e.symbol_and_number() for e in self.spacegroup_choices] }
    self.SendInfoToGUI(mydict)


  def set_spacegroup_choice(self, n) :
    if (self.viewer.miller_array is None) :
      raise Sorry("No data loaded!")
    self.current_spacegroup = self.spacegroup_choices[n]
    from cctbx import crystal
    symm = crystal.symmetry(
      space_group_info= self.current_spacegroup,
      unit_cell=self.viewer.miller_array.unit_cell())
    othervalidarrays = []
    #for validarray in self.valid_arrays:
    for validarray in self.procarrays:
      # TODO: check if array is unmerged i.e. not symmetry unique
      #print "Space group casting ", validarray.info().label_string()
      arr = validarray.expand_to_p1().customized_copy(crystal_symmetry=symm)
      arr = arr.merge_equivalents().array().set_info(validarray.info())
      arr = self.detect_Rfree(arr)
      othervalidarrays.append( arr )
    self.mprint( "MERGING 2", verbose=2)
    self.viewer.proc_arrays = othervalidarrays
    self.params.NGL_HKLviewer.using_space_subgroup = True
    #self.viewer.identify_suitable_fomsarrays()
    self.viewer.set_miller_array()
    #self.viewer.DrawNGLJavaScript()
    for i,e in enumerate(self.spacegroup_choices):
      self.mprint("%d, %s" %(i,e.symbol_and_number()) , verbose=0)


  def SetSpaceGroupChoice(self, n):
    self.params.NGL_HKLviewer.spacegroup_choice = n
    self.update_settings()


  def SetDefaultSpaceGroup(self):
    self.params.NGL_HKLviewer.using_space_subgroup = False
    self.update_settings()


  def set_default_spacegroup(self):
    self.viewer.proc_arrays = self.procarrays
    self.viewer.set_miller_array()
    self.viewer.identify_suitable_fomsarrays()


  def MakeNewMillerArrayFrom(self, operation, label, arrid1, arrid2=None):
    # get list of existing new miller arrays and operations if present
    miller_array_operations_lst = []
    if self.params.NGL_HKLviewer.miller_array_operations:
      miller_array_operations_lst = eval(self.params.NGL_HKLviewer.miller_array_operations)

    miller_array_operations_lst.append( ( operation, label, arrid1, arrid2 ) )
    self.params.NGL_HKLviewer.miller_array_operations = str( miller_array_operations_lst )
    self.update_settings()


  def make_new_miller_array(self):
    miller_array_operations_lst = eval(self.params.NGL_HKLviewer.miller_array_operations)

    unique_miller_array_operations_lst = []
    for (operation, label, arrid1, arrid2) in miller_array_operations_lst:
      isunique = True
      for arr in self.procarrays:
        if label in arr.info().label_string():
          self.mprint(label + " is already labelling one of the original miller arrays")
          isunique = False
          break
      if isunique:
        unique_miller_array_operations_lst.append( (operation, label, arrid1, arrid2) )

    self.params.NGL_HKLviewer.miller_array_operations = str(unique_miller_array_operations_lst)
    from copy import deepcopy
    millarr1 = deepcopy(self.procarrays[arrid1])
    newarray = None
    if arrid2:
      millarr2 = deepcopy(self.procarrays[arrid2])
      newarray = self.viewer.OperateOn2MillerArrays(millarr1, millarr2, operation)
    else:
      newarray = self.viewer.OperateOn1MillerArray(millarr1, operation)
    if newarray:
      newarray.set_info(millarr1._info )
      newarray._info.labels = [ label ]
      procarray, procarray_info = self.process_miller_array(newarray)
      self.procarrays.append(procarray)
      self.viewer.proc_arrays = self.procarrays
      self.viewer.has_new_miller_array = True
      self.viewer.array_infostrs.append( ArrayInfo(procarray, self.mprint).infostr )
      self.viewer.array_infotpls.append( ArrayInfo(procarray, self.mprint).infotpl )
      mydict = { "array_infotpls": self.viewer.array_infotpls, "NewHKLscenes" : True, "NewMillerArray" : True}
      self.SendInfoToGUI(mydict)




  def load_reflections_file(self, file_name, set_array=True, data_only=False):
    file_name = to_str(file_name)
    if (file_name != ""):
      self.mprint("Reading file...")
      from iotbx.reflection_file_reader import any_reflection_file
      self.viewer.isnewfile = True
      #self.params.NGL_HKLviewer.mergedata = None
      self.params.NGL_HKLviewer.viewer.scene_id = None
      self.viewer.colour_scene_id = None
      self.viewer.radii_scene_id = None
      self.viewer.match_valarrays = []
      self.viewer.proc_arrays = []
      self.spacegroup_choices = []
      display.reset_settings()
      self.settings = display.settings()
      self.viewer.settings = self.params.NGL_HKLviewer.viewer
      self.viewer.mapcoef_fom_dict = {}
      try :
        hkl_file = any_reflection_file(file_name)
        arrays = hkl_file.as_miller_arrays(merge_equivalents=False,
          )#observation_type_callback=misc_dialogs.get_shelx_file_data_type)
        #arrays = f.file_server.miller_arrays
      except Exception as e :
        self.NewFileLoaded=False
        self.mprint(to_str(e))
        arrays = []

      valid_arrays = []
      self.viewer.array_infostrs = []
      self.viewer.array_infotpls = []
      for array in arrays :
        #if array.is_hendrickson_lattman_array() :
        #  continue
        #elif (data_only) :
        if (not array.is_real_array()) and (not array.is_complex_array()) \
         and (not array.is_integer_array()) and (not array.is_bool_array()) :
          self.mprint('Ignoring miller array \"%s\" of %s' \
            %(array.info().label_string(), type(array.data()[0]) ) )
          continue
        self.viewer.array_infostrs.append( ArrayInfo(array, self.mprint).infostr )
        self.viewer.array_infotpls.append( ArrayInfo(array, self.mprint).infotpl )
        valid_arrays.append(array)
      self.valid_arrays = valid_arrays
      self.mprint("Miller arrays in this file:")
      for e in self.viewer.array_infostrs:
        self.mprint("%s" %e)
      self.mprint("\n")
      self.NewFileLoaded = True
      if (len(valid_arrays) == 0):
        msg = "No arrays of the supported types in this file."
        self.mprint(msg)
        self.NewFileLoaded=False
        return False
      elif (len(valid_arrays) >= 1):
        if (set_array):
          self.set_miller_array()
        mydict = { "info": self.infostr,
                   "array_infotpls": self.viewer.array_infotpls,
                   "bin_infotpls": self.viewer.bin_infotpls,
                   "html_url": self.viewer.url,
                   "merge_data": self.params.NGL_HKLviewer.merge_data,
                   "spacegroups": [e.symbol_and_number() for e in self.spacegroup_choices],
                   "NewFileLoaded": self.NewFileLoaded
                  }
        self.SendInfoToGUI(mydict)
        return True


  def LoadReflectionsFile(self, filename):
    #import code, traceback; code.interact(local=locals(), banner="".join( traceback.format_stack(limit=10) ) )
    self.params.NGL_HKLviewer.filename = filename
    self.update_settings()


  def SetCameraType(self, camtype):
    self.params.NGL_HKLviewer.viewer.NGL.camera_type = camtype
    self.update_settings()


  def ExpandToP1(self, val, inbrowser=True):
    self.params.NGL_HKLviewer.viewer.expand_to_p1 = val
    self.params.NGL_HKLviewer.viewer.inbrowser = inbrowser
    self.update_settings()


  def ExpandAnomalous(self, val, inbrowser=True):
    self.params.NGL_HKLviewer.viewer.expand_anomalous = val
    self.params.NGL_HKLviewer.viewer.inbrowser = inbrowser
    self.update_settings()


  def ShowOnlyMissing(self, val):
    self.params.NGL_HKLviewer.viewer.show_only_missing = val
    self.update_settings()


  def ShowMissing(self, val):
    self.params.NGL_HKLviewer.viewer.show_missing = val
    self.update_settings()


  def ShowDataOverSigma(self, val):
    self.params.NGL_HKLviewer.viewer.show_data_over_sigma = val
    self.update_settings()


  def ShowSystematicAbsences(self, val):
    self.params.NGL_HKLviewer.viewer.show_systematic_absences = val
    self.update_settings()


  def ShowSlice(self, val, axis="h", index=0):
    axisstr = axis.lower()
    self.params.NGL_HKLviewer.viewer.slice_mode = val
    self.params.NGL_HKLviewer.viewer.slice_axis = axisstr
    self.params.NGL_HKLviewer.viewer.slice_index = index
    self.update_settings()


  def set_scene_bin_thresholds(self, binvals=None, bin_scene_label="Resolution", nbins=6):
    if binvals:
      binvals = list( 1.0/flex.double(binvals) )
    else:
      binvals = self.viewer.calc_bin_thresholds(bin_scene_label, nbins)
    self.viewer.UpdateBinValues( binvals )


  def SetSceneBinLabel(self, bin_scene_label="Resolution"):
    self.params.NGL_HKLviewer.bin_scene_label = bin_scene_label
    self.update_settings()


  def SetSceneNbins(self, nbins):
    self.params.NGL_HKLviewer.nbins = nbins
    self.params.NGL_HKLviewer.viewer.NGL.bin_opacities = str([ (1.0, e) for e in range(nbins) ])
    self.update_settings()


  def SetSceneBinThresholds(self, binvals=[]):
    self.params.NGL_HKLviewer.scene_bin_thresholds = binvals
    self.update_settings()


  def SetToolTipOpacity(self, val):
    self.params.NGL_HKLviewer.viewer.NGL.tooltip_alpha = val
    self.update_settings()


  def SetOpacities(self, bin_opacities):
    self.params.NGL_HKLviewer.viewer.NGL.bin_opacities = bin_opacities
    self.update_settings()


  def set_scene(self, scene_id) :
    self.viewer.binvals = []
    self.viewer.isinjected = False
    if scene_id is None:
      return False
    self.viewer.colour_scene_id = scene_id
    self.viewer.radii_scene_id = scene_id
    self.viewer.set_miller_array(scene_id)
    if (self.viewer.miller_array is None):
      raise Sorry("No data loaded!")
    self.mprint( "Miller array %s runs from hkls: %s to %s" \
     %(self.viewer.miller_array.info().label_string(), self.viewer.miller_array.index_span().min(),
        self.viewer.miller_array.index_span().max() ) )
    self.mprint("Spacegroup: %s" %self.viewer.miller_array.space_group().info().symbol_and_number())
    self.update_space_group_choices()
    return True


  def SetScene(self, scene_id):
    self.params.NGL_HKLviewer.viewer.scene_id = scene_id
    self.update_settings()


  def SetMergeData(self, val):
    self.params.NGL_HKLviewer.merge_data = val
    self.update_settings()

  def SetColourScene(self, colourcol):
    self.params.NGL_HKLviewer.viewer.colour_scene_id = colourcol
    self.update_settings()


  def SetRadiusScene(self, radiuscol):
    self.params.NGL_HKLviewer.viewer.radii_scene_id = radiuscol
    self.update_settings()


  def SetRadiiScale(self, scale, nth_power_scale = -1.0):
    """
    Scale radii. Decrease the contrast between large and small radii with nth_root_scale < 1.0
    If nth_power_scale=0.0 then all radii will have the same size regardless of data values.
    If nth_power_scale < 0.0 an automatic power will be computed ensuring the smallest radius
    is 0.1 times the maximum radius
    """
    self.params.NGL_HKLviewer.viewer.scale = scale
    self.params.NGL_HKLviewer.viewer.nth_power_scale_radii = nth_power_scale
    self.update_settings()


  def SetRadiiToSigmas(self, val):
    self.params.NGL_HKLviewer.viewer.sigma_radius = val
    self.update_settings()


  def SetColoursToSigmas(self, val):
    self.params.NGL_HKLviewer.viewer.sigma_color = val
    self.update_settings()


  def SetSqrtScaleColours(self, val):
    self.params.NGL_HKLviewer.viewer.sqrt_scale_colors = val
    self.update_settings()


  def SetColoursToPhases(self, val):
    self.params.NGL_HKLviewer.viewer.phase_color = val
    self.update_settings()


  def SetShapePrimitive(self, val):
    self.params.NGL_HKLviewer.shape_primitive = val
    self.update_settings()


  def set_shape_primitive(self, val):
    if val == "points":
      self.viewer.primitivetype = "PointBuffer"
    else:
      self.viewer.primitivetype = "sphereBuffer"
    #self.viewer.DrawNGLJavaScript()


  def ShowTooltips(self, val):
    self.params.NGL_HKLviewer.tooltips_in_script  = val
    self.update_settings()


  def NormalVectorToClipPlane(self, h, k, l, hkldist=0.0,
                           clipNear=None, clipFar=None, fixorientation=True):
    self.viewer.RemoveAllReciprocalVectors()
    self.viewer.AddVector(h, k, l)
    if fixorientation:
      self.viewer.DisableMouseRotation()
    else:
      self.viewer.EnableMouseRotation()
    self.viewer.PointVectorOut()
    if clipNear is None or clipFar is None:
      halfdist = (self.viewer.OrigClipFar - self.viewer.OrigClipNear) / 2.0
      self.viewer.GetBoundingBox()
      clipNear = halfdist - self.viewer.scene.min_dist*50/self.viewer.boundingZ
      clipFar = halfdist + self.viewer.scene.min_dist*50/self.viewer.boundingZ
    self.viewer.SetClipPlaneDistances(clipNear, clipFar, self.viewer.cameraPosZ)
    self.viewer.TranslateHKLpoints(h, k, l, hkldist)


  def ClipPlaneNormalToHKLvector(self, h, k, l, hkldist=0.0,
             clipwidth=None, fixorientation=True):
    self.params.NGL_HKLviewer.normal_clip_plane.h = h
    self.params.NGL_HKLviewer.normal_clip_plane.k = k
    self.params.NGL_HKLviewer.normal_clip_plane.l = l
    self.params.NGL_HKLviewer.normal_clip_plane.hkldist = hkldist
    self.params.NGL_HKLviewer.normal_clip_plane.clipwidth = clipwidth
    self.params.NGL_HKLviewer.viewer.NGL.fixorientation = fixorientation
    self.update_settings()


  def SetTrackBallRotateSpeed(self, trackspeed):
    self.params.NGL_HKLviewer.viewer.NGL.mouse_sensitivity = trackspeed
    self.update_settings()


  def GetTrackBallRotateSpeed(self):
    self.viewer.GetTrackBallRotateSpeed()
    while self.params.viewer.NGL.mouse_sensitivity is None:
      time.sleep(0.1)
    return self.params.viewer.NGL.mouse_sensitivity


  def GetSpaceGroupChoices(self):
    """
    return array of strings with available subgroups of the space group
    """
    if (self.viewer.miller_array is None) :
      self.mprint( NOREFLDATA)
    if self.spacegroup_choices:
      return [e.symbol_and_number() for e in self.spacegroup_choices]
    return []


  def GetHtmlURL(self):
    return self.viewer.url


  def GetNGLstring(self):
    return self.viewer.NGLscriptstr


  def GetArrayInfotpls(self):
    """
    return array of tuples with information on each miller array
    """
    return self.array_infotpls


  def GetHklScenesInfo(self):
    """
    return array of strings with information on each processed miller array
    which may have been expanded with anomalous reflections or truncated to non-anomalous reflections
    as to match the currently selected miller array
    """
    return self.viewer.hkl_scenes_info


  def GetBinInfo(self):
    """
    return array of number of hkls and bin boundaries of the bins the current miller array data has been sorted into.
    Useful when deciding which bin of reflections to make transparent
    """
    return self.viewer.binstrs


  def SendInfoToGUI(self, infodict):
    if self.guiSocketPort:
      self.guisocket.send( str(infodict).encode("utf-8") )


masterphilstr = """
NGL_HKLviewer {
  filename = None
    .type = path
  merge_data = False
    .type = bool
  miller_array_operations = ''
    .type = str
  spacegroup_choice = None
    .type = int
  using_space_subgroup = False
    .type = bool
  normal_clip_plane {
    h = 2.0
      .type = float
    k = 0
      .type = float
    l = 0
      .type = float
    hkldist = 0.0
      .type = float
    clipwidth = None
      .type = float
  }
  scene_bin_thresholds = None
    .type = float
    .multiple = True
  bin_scene_label = 'Resolution'
    .type = str
  nbins = 6
    .type = int
  shape_primitive = *'spheres' 'points'
    .type = choice
  tooltips_in_script = False
    .type = bool
  viewer {
    scene_id = None
      .type = int
    %s
    NGL {
      %s
    }
  }
}

""" %(display.philstr, view_3d.ngl_philstr)



if __name__ == '__main__':
  myHKLview = HKLViewFrame()
