
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
myHKLview.SetColumn(0)

myHKLview.SetRadiiScale(3, nth_power_scale=1)
myHKLview.SetColourColumn(3)
myHKLview.SetRadiusColumn(7)

myHKLview = cmdlineframes.HKLViewFrame(jscriptfname = "myjstr.js", htmlfname = "myhkl.html", UseOSBrowser= False)

myHKLview.SetColumn(3)
myHKLview.SetRadiiScale(3, nth_root_scale=1)
myHKLview.SetColumnBinThresholds([0, 3.4, 3.5, 3.6, 3.8, 4.0, 4.5, 5.0, 6.0, 8, 9,11, 15, 20, 30, 100 ], "Resolution")
for i in range(13):
  myHKLview.SetOpacity(i, 0.0)

myHKLview.SetColumnBinThresholds([0, 0.1, 1, 10], 3)
myHKLview.ExpandToP1(True)
myHKLview.SetOpacity(2, 0.1)

myHKLview.ShowSlice(True, "l", 17)
myHKLview.GetNGLstring()
myHKLview.ShowSlice(False)
myHKLview.ExpandAnomalous(True)
myHKLview.ShowMissing(True)
myHKLview.GetSpaceGroupChoices()
myHKLview.SetSpaceGroupChoice(3)

myHKLview.LoadReflectionsFile(r"C:\Users\oeffner\Buser\Phenix\dev-2814-working\modules\phenix_examples\lysozyme-MRSAD\lyso2001_scala1.mtz")
myHKLview.SetCameraType("persp")

myHKLview.LoadReflectionsFile(r"C:\Users\oeffner\Buser\Experiments\CRLF3\DLS20151206CRLF3\5840-F11-X1-Hg-SAD-ONsoak\5840-F11-X1_pk_5_5_1_\xia2\dials-run\DataFiles\mx11235v49_x5840F11X1pk551_free.mtz")

myHKLview.SetColumn(0)
myHKLview.SetRadiiScale(1, nth_power_scale=0.2)
myHKLview.SetColumn(1)
myHKLview.SetColumnBinThresholds(1, [-20, 30, 300, 3000])


myHKLview.ShowSlice(True, "h", 20)


myHKLview.LoadReflectionsFile(r"C:\Users\oeffner\Buser\Phenix\dev-2814-working\modules\phenix_examples\beta-blip\beta_blip_P3221.mtz")
myHKLview.SetColumn(0)
myHKLview.SetColumnBinThresholds([3, 3.2, 3.4, 3.8, 4.5, 6, 9, 30])
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

ma1 = miller.array( miller.set(xs, mi), ) [11.205, 6.353, 26.167, 14.94, 2.42, 24.921, 16.185, 11.798, 21.183, 4.98] ),
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

mi3 = flex.miller_index([ (1,-2,3), (0,0,-4), (0, 1, 2), (1, 0, 2), (-1, 1, -2), (-2, 1, 0) , (1, 0, -2), (0, 0, 2) ] )
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
mtz1.add_miller_array(ma4, column_root_label="blip")
mtz1.add_miller_array(ma5, column_root_label="bleep")
mtz1.add_miller_array(mafom, column_root_label="FOM")
mtz1.set_wavelength(1.2)
mtz1.set_name("MyTestData")
mtz1.mtz_object().write("mymtz.mtz")


from crys3d.hklview import cmdlineframes
myHKLview = cmdlineframes.HKLViewFrame(jscriptfname = "myjstr.js", verbose=False)
myHKLview.LoadReflectionsFile("mymtz.mtz")
myHKLview.SetColumn(0)
no
myHKLview.SetRadiiScale(1, nth_power_scale=0.2)
myHKLview.SetColumnBinThresholds([50, 20, 15, 12, 9])
myHKLview.SetOpacity(0, 0.0)
myHKLview.SetOpacity(1, 0.0)
myHKLview.SetOpacity(2, 0.0)


myHKLview.ExpandToP1(True)
myHKLview.ExpandAnomalous(True)
myHKLview.ShowMissing(True)
myHKLview.ShowSystematicAbsences(True)

myHKLview.SetColourColumn(2)
myHKLview.SetRadiusColumn(3)
myHKLview.SetColumn(1)

myHKLview.SetColumn(0)
myHKLview.SetRadiiToSigmas(True)

from crys3d.hklview import cmdlineframes
myHKLview = cmdlineframes.HKLViewFrame(jscriptfname = "myjstr.js", htmlfname = "myhkl.html")
myHKLview.LoadReflectionsFile(r"C:\Users\oeffner\Buser\Tests\MRproblem\MRproblem_15.1.mtz")

myHKLview.SetColumn(2)
myHKLview.SetColoursToPhases(True)


from crys3d.hklview import cmdlineframes
myHKLview = cmdlineframes.HKLViewFrame(jscriptfname = "myjstr.js")
myHKLview.LoadReflectionsFile("3RP2_A.1.mtz")
myHKLview.SetColumn(1,4)
myHKLview.SetRadiiScale(1, nth_power_scale=0.0)
myHKLview.ExpandToP1(True)
myHKLview.ExpandAnomalous(True)
myHKLview.ShowSlice(True, "l", 25)
myHKLview.ShowMissing(True)


myHKLview.LoadReflectionsFile(r"C:\Users\oeffner\Buser\Work\TNCS\4N3E\4n3e_final.mtz")


"""



from cctbx.miller import display2 as display
from crys3d.hklview import jsview_3d as view_3d
from crys3d.hklview.jsview_3d import ArrayInfo
from libtbx.str_utils import format_value
from cctbx.array_family import flex
from libtbx.utils import Sorry, to_str
from libtbx import group_args
import libtbx
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
    self.array_infostrs = []
    self.array_infotpls = []
    self.procarrays = []
    self.merge_answer = [None]
    self.dmin = -1
    self.settings = display.settings()
    self.verbose = True
    if 'verbose' in kwds:
      self.verbose = kwds['verbose']
    kwds['settings'] = self.settings
    kwds['mprint'] = self.mprint
    self.viewer = view_3d.hklview_3d( **kwds )
    self.ResetPhilandViewer()
    self.NewFileLoaded = False
    self.infostr = ""
    self.useSocket=False
    if 'useSocket' in kwds:
      self.useSocket = kwds['useSocket']
      self.context = zmq.Context()
      self.socket = self.context.socket(zmq.PAIR)
      self.socket.connect("tcp://127.0.0.1:7895")
      self.STOP = False
      print("starting socketthread")
      self.msgqueuethrd = threading.Thread(target = self.zmq_listen )
      self.msgqueuethrd.daemon = True
      self.msgqueuethrd.start()


  def zmq_listen(self):
    while not self.STOP:
      philstr = self.socket.recv()
      philstr = str(philstr)
      self.mprint("Received phil string:\n" + philstr)
      new_phil = libtbx.phil.parse(philstr)
      self.update_settings(new_phil)
      #print "in zmq listen"
      time.sleep(1)
    del self.socket


  def ResetPhilandViewer(self, extraphil=None):
    self.master_phil = libtbx.phil.parse( masterphilstr )
    self.currentphil = self.master_phil
    if extraphil:
      self.currentphil = self.currentphil.fetch(source = extraphil)
    self.params = self.currentphil.fetch().extract()
    self.viewer.miller_array = None
    self.viewer.symops = []
    self.viewer.sg = None
    self.viewer.proc_arrays = []
    self.viewer.HKLscenesdict = {}


  def GetNewCurrentPhilFromString(self, philstr, oldcurrentphil):
    user_phil = libtbx.phil.parse(philstr)
    newcurrentphil = oldcurrentphil.fetch(source = user_phil)
    diffphil = oldcurrentphil.fetch_diff(source = user_phil)
    return newcurrentphil, diffphil


  def GetNewCurrentPhilFromPython(self, pyphilobj, oldcurrentphil):
    newcurrentphil = oldcurrentphil.fetch(source = pyphilobj)
    diffphil = oldcurrentphil.fetch_diff(source = pyphilobj)

    oldcolbintrshld = oldcurrentphil.extract().NGL_HKLviewer.column_bin_thresholds
    newcolbintrshld = pyphilobj.extract().NGL_HKLviewer.column_bin_thresholds
    # fetch_diff doesn't seem able to correclty spot changes
    # in the multiple scope phil object "NGL_HKLviewer.column_bin_thresholds"
    # Must do it manually
    if oldcolbintrshld != newcolbintrshld:
      params = newcurrentphil.extract()
      params.NGL_HKLviewer.column_bin_thresholds = newcolbintrshld
      newcurrentphil = self.master_phil.format(python_object = params)
      diffphil = self.master_phil.fetch_diff(source = newcurrentphil)

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
      self.mprint("diff phil:\n" + diff_phil.as_str(), True )

      self.params = self.currentphil.extract()
      phl = self.params.NGL_HKLviewer

      if hasattr(diff, "filename"):
        self.ResetPhilandViewer(diff_phil)
        self.load_reflections_file(phl.filename)

      #if hasattr(diff, "column") or hasattr(diff, "fom_column") \
      if hasattr(diff, "column") \
       or hasattr(diff, "merge_data") or hasattr(diff, "column_bin_thresholds"):
        #if self.set_column(phl.column, phl.fom_column):
        if self.set_column(phl.column):
          self.set_column_bin_thresholds(phl.column_bin_thresholds, phl.bin_array)

      if hasattr(diff, "spacegroup_choice"):
        self.set_spacegroup_choice(phl.spacegroup_choice)

      if hasattr(diff, "using_space_subgroup") and phl.using_space_subgroup==False:
        self.set_default_spacegroup()

      if hasattr(diff, "camera_type"):
        self.set_camera_type(phl.camera_type)

      if hasattr(diff, "shape_primitive"):
        self.set_shape_primitive(phl.shape_primitive)

      if hasattr(diff, "tooltips_in_script"):
        self.viewer.script_has_tooltips = phl.tooltips_in_script

      if hasattr(diff, "viewer"):
        self.viewer.settings = phl.viewer
        self.settings = phl.viewer

      msg = self.viewer.update_settings(diff, phl)
      self.mprint( msg, True)
      for i,e in enumerate(self.spacegroup_choices):
        self.mprint("%d, %s" %(i,e.symbol_and_number()) , True)
      #import code, traceback; code.interact(local=locals(), banner="".join( traceback.format_stack(limit=10) ) )
      self.SendInfo()
      self.NewFileLoaded = False

      if (self.viewer.miller_array is None) :
        self.mprint( "No miller array has been selected", True)
        return False
      return True
    except Exception as e:
      self.mprint(to_str(e) + "\n" + traceback.format_exc(), True)
      return False


  def mprint(self, m, verbose=False):
    if self.verbose or verbose:
      print(m)


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
      raise Sorry("Hendrickson-Lattman coefficients are not supported.")
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
        self.SendInfo()
        while 1:
          philstr = self.socket.recv()
          philstr = str(philstr)
          self.mprint("process_miller_array, Received phil:\n" + philstr, True)
          new_phil = libtbx.phil.parse(philstr)
          #working_phil = self.master_phil.fetch(source = new_phil)
          params = new_phil.extract().NGL_HKLviewer
          if hasattr(params, "merge_data"): # awaiting user to tick a checkbox on the gui
            self.params.NGL_HKLviewer.mergedata = params.merge_data
            break
          time.sleep(1)
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
    if col < 0:
      array_info = procarray_info
    self.merge_answer = [None]
    #import code, traceback; code.interact(local=locals(), banner="".join( traceback.format_stack(limit=10) ) )
    #self.viewer.set_miller_array(col, merge=array_info.merge,
    #   details=array_info.details_str, proc_arrays=self.procarrays)
    #self.viewer.identify_suitable_fomsarrays()
    return array_info


  def set_miller_array(self, col=-1) :
    #print "in set_miller_array"
    #if col >= len(self.valid_arrays):
    if col >= len(self.viewer.hkl_scenes_info ):
      return
    array_info = self.process_all_miller_arrays(col)
    self.viewer.set_miller_array(col, merge=array_info.merge,
       details=array_info.details_str)
    self.viewer.proc_arrays=self.procarrays
    self.viewer.identify_suitable_fomsarrays()
    self.update_space_group_choices()


  def update_space_group_choices(self) :
    if self.viewer.miller_array is None or \
      self.params.NGL_HKLviewer.using_space_subgroup:
      return
    current_miller_array_idx = self.viewer.hkl_scenes_info[self.viewer.iarray][1]
    matching_valid_array = self.valid_arrays[ current_miller_array_idx ]
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


  def set_spacegroup_choice(self, n) :
    if (self.viewer.miller_array is None) :
      raise Sorry("No data loaded!")
    self.current_spacegroup = self.spacegroup_choices[n]
    from cctbx import crystal
    symm = crystal.symmetry(
      space_group_info= self.current_spacegroup,
      unit_cell=self.viewer.miller_array.unit_cell())
    othervalidarrays = []
    for validarray in self.valid_arrays:
      # TODO: check if array is unmerged i.e. not symmetry unique
      #print "Space group casting ", validarray.info().label_string()
      arr = validarray.expand_to_p1().customized_copy(crystal_symmetry=symm)
      arr = arr.merge_equivalents().array().set_info(validarray.info())
      arr = self.detect_Rfree(arr)
      othervalidarrays.append( arr )
    self.mprint( "MERGING 2")
    self.viewer.proc_arrays = othervalidarrays
    self.params.NGL_HKLviewer.using_space_subgroup = True
    #self.viewer.identify_suitable_fomsarrays()
    self.viewer.set_miller_array()
    #self.viewer.DrawNGLJavaScript()


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



  def load_reflections_file(self, file_name, set_array=True, data_only=False):
    file_name = to_str(file_name)
    if (file_name != ""):
      from iotbx.reflection_file_reader import any_reflection_file
      self.viewer.isnewfile = True
      #self.params.NGL_HKLviewer.mergedata = None
      self.viewer.iarray = -1
      self.viewer.icolourcol = -1
      self.viewer.iradiicol = -1
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
        self.mprint(to_str(e), True)
        arrays = []
      valid_arrays = []
      self.array_infostrs = []
      self.array_infotpls = []
      for array in arrays :
        #if array.is_hendrickson_lattman_array() :
        #  continue
        #elif (data_only) :
        if (not array.is_real_array()) and (not array.is_complex_array()) \
         and (not array.is_integer_array()) and (not array.is_bool_array()) :
          self.mprint('Ignoring miller array \"%s\" of %s' \
            %(array.info().label_string(), type(array.data()[0]) ), True )
          continue
        self.array_infostrs.append( ArrayInfo(array, self.mprint).infostr )
        self.array_infotpls.append( ArrayInfo(array, self.mprint).infotpl )
        valid_arrays.append(array)
      self.valid_arrays = valid_arrays
      for i,e in enumerate(self.array_infostrs):
        self.mprint("%d, %s" %(i, e), True)
      self.NewFileLoaded = True
      if (len(valid_arrays) == 0):
        msg = "No arrays of the supported types in this file."
        self.mprint(msg, True)
        self.NewFileLoaded=False
      elif (len(valid_arrays) >= 1):
        if (set_array):
          self.set_miller_array()
      return True


  def LoadReflectionsFile(self, filename):
    #import code, traceback; code.interact(local=locals(), banner="".join( traceback.format_stack(limit=10) ) )
    self.params.NGL_HKLviewer.filename = filename
    self.update_settings()


  def set_camera_type(self, camtype):
    self.viewer.cameratype = self.params.NGL_HKLviewer.cameratype
    self.viewer.DrawNGLJavaScript()


  def SetCameraType(self, camtype):
    self.params.NGL_HKLviewer.cameratype = camtype
    self.update_settings()


  def ExpandToP1(self, val):
    self.params.NGL_HKLviewer.viewer.expand_to_p1 = val
    self.update_settings()


  def ExpandAnomalous(self, val):
    self.params.NGL_HKLviewer.viewer.expand_anomalous = val
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
    self.params.NGL_HKLviewer.viewer.slice_mode = val
    self.params.NGL_HKLviewer.viewer.slice_axis = axis.lower()
    self.params.NGL_HKLviewer.viewer.slice_index = index
    self.update_settings()


  def set_column_bin_thresholds(self, binvals=[], bin_array="Resolution"):
    self.viewer.binarray = bin_array
    if binvals:
      if self.viewer.binarray=="Resolution":
        binvals = list( 1.0/flex.double(binvals) )
      binvals.sort()
    self.viewer.UpdateBinValues( binvals )


  def SetColumnBinThresholds(self, binvals, bin_array="Resolution"):
    self.params.NGL_HKLviewer.column_bin_thresholds = binvals
    self.params.NGL_HKLviewer.bin_array = bin_array
    self.update_settings()


  def SetOpacity(self, bin, alpha):
    self.viewer.SetOpacity(bin, alpha)


  #def set_column(self, column, fom_column=None) :
  def set_column(self, column) :
    self.viewer.binvals = []
    if column is None:
      return False
    self.viewer.icolourcol = column
    self.viewer.iradiicol = column
    """
    if self.valid_arrays[column].is_complex_array():
      if fom_column and self.valid_arrays[fom_column].is_real_array():
        self.viewer.mapcoef_fom_dict[self.valid_arrays[column].info().label_string()] = fom_column
        self.mprint("Using array %d as FOM values for array %d" %(fom_column, column))
      else:
        if self.viewer.mapcoef_fom_dict.get(self.valid_arrays[column].info().label_string()):
          del self.viewer.mapcoef_fom_dict[self.valid_arrays[column].info().label_string()]
        self.mprint("No valid FOM array provided.")
    """
    #self.set_miller_array(self.valid_arrays[column])
    #self.set_miller_array(column)
    self.viewer.set_miller_array(column)
    self.viewer.miller_array = self.viewer.HKLscenes[column].miller_array
    if (self.viewer.miller_array is None):
      raise Sorry("No data loaded!")
    self.mprint( "Miller array %s runs from hkls: %s to %s" \
     %(self.viewer.miller_array.info().label_string(), self.viewer.miller_array.index_span().min(),
        self.viewer.miller_array.index_span().max() ), True )
    self.update_space_group_choices()
    #self.viewer.DrawNGLJavaScript()
    #self.update_settings()
    #import code, traceback; code.interact(local=locals(), banner="".join( traceback.format_stack(limit=10) ) )
    return True


  def SetColumn(self, column):
    self.params.NGL_HKLviewer.column = column
    self.update_settings()


  def SetMergeData(self, val):
    self.params.NGL_HKLviewer.merge_data = val
    self.update_settings()

  def SetColourColumn(self, colourcol):
    self.params.NGL_HKLviewer.viewer.icolourcol = colourcol
    self.update_settings()


  def SetRadiusColumn(self, radiuscol):
    self.params.NGL_HKLviewer.viewer.iradiicol = radiuscol
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


  def GetSpaceGroupChoices(self):
    """
    return array of strings with available subgroups of the space group
    """
    if (self.viewer.miller_array is None) :
      self.mprint( "No miller array has been selected", True)
    if self.spacegroup_choices:
      return [e.symbol_and_number() for e in self.spacegroup_choices]
    return []


  def GetHtmlURL(self):
    return self.viewer.url


  def GetNGLstring(self):
    return self.viewer.NGLscriptstr


  def GetArrayInfotpls(self):
    """
    return array of strings with information on each miller array
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


  def SendInfo(self):
    mydict = { "info": self.infostr,
               "miller_arrays": self.array_infotpls,
               "hklscenes_arrays": self.viewer.hkl_scenes_info,
               "bin_info": self.viewer.binstrs,
               "html_url": self.viewer.url,
               "merge_data": self.params.NGL_HKLviewer.merge_data,
               "spacegroups": [e.symbol_and_number() for e in self.spacegroup_choices],
               "NewFileLoaded": self.NewFileLoaded
            }
    #print str(mydict)
    if self.useSocket:
      self.socket.send( str(mydict).encode("utf-8") )



masterphilstr = """
NGL_HKLviewer {
  filename = None
    .type = path
  column = None
    .type = int
  merge_data = False
    .type = bool
  spacegroup_choice = None
    .type = int
  using_space_subgroup = False
    .type = bool
  column_bin_thresholds = None
    .type = float
    .multiple = True
  bin_array = 'Resolution'
    .type = str
  camera_type = *'orthographic' 'perspective'
    .type = choice
  shape_primitive = *'spheres' 'points'
    .type = choice
  tooltips_in_script = False
    .type = bool
  viewer {
    %s
  }
}

""" %display.philstr



if __name__ == '__main__':
  myHKLview = HKLViewFrame()
