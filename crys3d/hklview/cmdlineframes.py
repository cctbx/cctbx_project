
from __future__ import absolute_import, division, print_function

from iotbx.reflection_file_reader import any_reflection_file
from cctbx.miller import display2 as display
from crys3d.hklview import jsview_3d as view_3d
from crys3d.hklview.jsview_3d import ArrayInfo
from cctbx import miller
from libtbx.math_utils import roundoff
from libtbx.str_utils import format_value
from cctbx.array_family import flex
from libtbx.utils import Sorry, to_str
from scitbx import matrix
from cctbx import sgtbx
from libtbx import group_args
import libtbx
import libtbx.load_env
import traceback
import sys, zmq, threading,  time, cmath, zlib, os.path, math, re


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
    self.origarrays = {}
    self.merge_answer = [None]
    self.dmin = -1
    self.settings = display.settings()
    self.verbose = 0
    if 'verbose' in kwds:
      self.verbose = eval(kwds['verbose'])
    self.guiSocketPort=None
    self.mprint("kwds= " +str(kwds), 1)
    self.mprint("args= " + str(args), 1)
    kwds['settings'] = self.settings
    kwds['mprint'] = self.mprint
    self.infostr = ""
    self.hklfile_history = []
    self.tncsvec = None
    self.uservectors = []
    self.new_miller_array_operations_lst = []
    self.copyrightpaths = [("CCTBX copyright", libtbx.env.under_root(os.path.join("modules","cctbx_project","COPYRIGHT.txt"))),
     ("NGL copyright", libtbx.env.under_root(os.path.join("modules","cctbx_project","crys3d","hklview","License_for_NGL.txt"))),
     ("html2canvas copyright", libtbx.env.under_root(os.path.join("modules","cctbx_project","crys3d","hklview","LICENSE_for_html2canvas.txt")))
    ]
    self.zmqsleeptime = 0.1
    if 'useGuiSocket' in kwds:
      self.guiSocketPort = eval(kwds['useGuiSocket'])
      self.context = zmq.Context()
      self.guisocket = self.context.socket(zmq.PAIR)
      self.guisocket.connect("tcp://127.0.0.1:%s" %self.guiSocketPort )
      self.STOP = False
      self.mprint("starting socketthread", 1)
      # name this thread to ensure any asyncio functions are called only from main thread
      self.msgqueuethrd = threading.Thread(target = self.zmq_listen, name="HKLviewerZmqThread" )
      self.msgqueuethrd.daemon = True
      kwds['send_info_to_gui'] = self.SendInfoToGUI # function also used by hklview_3d
      pyversion = "cctbx.python.version: " + str(sys.version_info[0])
      # tell gui what python version we are
      self.SendInfoToGUI(pyversion )
      self.SendInfoToGUI({"copyrights": self.copyrightpaths } )
    kwds['websockport'] = self.find_free_port()
    kwds['parent'] = self
    self.viewer = view_3d.hklview_3d( **kwds )
    self.ResetPhilandViewer()
    self.idx_data = None
    self.NewFileLoaded = False
    self.loaded_file_name = ""
    self.hklin = None
    if 'hklin' in kwds or 'HKLIN' in kwds:
      self.hklin = kwds.get('hklin', kwds.get('HKLIN') )
      self.LoadReflectionsFile(self.hklin)
    if 'useGuiSocket' in kwds:
      self.msgqueuethrd.start()


  def __exit__(self, exc_type=None, exc_value=0, traceback=None):
    self.viewer.__exit__(exc_type, exc_value, traceback)
    self.mprint("Destroying HKLViewFrame", verbose=0) # this string is expected by HKLviewer.py so don't change
    self.STOP = True
    del self
    sys.exit()


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
    #time.sleep(5)
    while not self.STOP:
      try:
        msgstr = self.guisocket.recv().decode("utf-8")
        self.mprint("Received string:\n" + msgstr, verbose=1)
        msgtype, mstr = eval(msgstr)
        if msgtype=="dict":
          self.viewer.datatypedict = eval(mstr)
        if msgtype=="philstr":
          new_phil = libtbx.phil.parse(mstr)
          self.update_settings(new_phil)
        time.sleep(self.zmqsleeptime)
      except Exception as e:
        self.mprint( str(e) + traceback.format_exc(limit=10), verbose=1)
    self.mprint( "Shutting down zmq_listen() thread", 1)
    self.guiSocketPort=None


  def ResetPhilandViewer(self, extraphil=None):
    self.master_phil = libtbx.phil.parse( masterphilstr )
    self.currentphil = self.master_phil
    if extraphil:
      self.currentphil = self.currentphil.fetch(source = extraphil)
      # Don't retain clip plane values as these are specific to each crystal
      # so use clip plane parameters from the master phil
      default_clipphil = self.master_phil.fetch().extract().clip_plane
      currentparms = self.currentphil.extract()
      currentparms.clip_plane = default_clipphil
      self.currentphil = self.master_phil.format(python_object = currentparms)
    self.params = self.currentphil.fetch().extract()
    self.viewer.viewerparams = self.params.viewer
    self.viewer.params = self.params
    self.params.binner_idx = 0
    self.params.nbins = 1
    self.params.scene_bin_thresholds = ""
    self.params.using_space_subgroup = False
    self.viewer.symops = []
    self.viewer.sg = None
    self.viewer.proc_arrays = []
    self.viewer.HKLscenedict = {}
    self.uservectors = []
    self.viewer.visual_symmxs = []
    self.visual_symHKLs = []
    self.viewer.sceneisdirty = True
    self.viewer.isnewfile = True
    if self.viewer.miller_array:
      self.viewer.params.viewer.scene_id = None
      self.viewer.RemoveStageObjects()
    self.viewer.miller_array = None
    self.viewer.lastviewmtrx = None
    return self.viewer.params


  def GetNewCurrentPhilFromString(self, philstr, oldcurrentphil):
    user_phil = libtbx.phil.parse(philstr)
    newcurrentphil = oldcurrentphil.fetch(source = user_phil)
    diffphil = oldcurrentphil.fetch_diff(source = user_phil)
    return newcurrentphil, diffphil


  def GetNewCurrentPhilFromPython(self, pyphilobj, oldcurrentphil):
    newcurrentphil, unusedphilparms = oldcurrentphil.fetch(source = pyphilobj, track_unused_definitions=True)
    for parm in unusedphilparms:
      self.mprint( "Received unrecognised phil parameter: " + parm.path, verbose=1)
    diffphil = oldcurrentphil.fetch_diff(source = pyphilobj)
    """
    oldcolbintrshld = oldcurrentphil.extract().scene_bin_thresholds
    newcolbintrshld = oldcolbintrshld
    if hasattr(pyphilobj.extract(), "scene_bin_thresholds"):
      newcolbintrshld = pyphilobj.extract().scene_bin_thresholds
    # fetch_diff doesn't seem able to correclty spot changes
    # in the multiple scope phil object "scene_bin_thresholds"
    # Must do it manually
    params = newcurrentphil.extract()
    if oldcolbintrshld != newcolbintrshld: # or old_binopacities != new_binopacities:
      params.scene_bin_thresholds = newcolbintrshld
      newcurrentphil = self.master_phil.format(python_object = params)
      diffphil = self.master_phil.fetch_diff(source = newcurrentphil)
    #import code, traceback; code.interact(local=locals(), banner="".join( traceback.format_stack(limit=10) ) )
    """
    return newcurrentphil, diffphil


  def SetCurrentPhilAsPython(self, pyphil):
    newphil = master_phil.format(python_object= pyphil)
    currphil = master_phil.fetch(source = newphil)


  def update_settings(self, new_phil=None):
    try:
      if not new_phil:
        #self.params = self.viewer.params
        new_phil = self.master_phil.format(python_object = self.params)
      #import code, traceback; code.interact(local=locals(), banner="".join( traceback.format_stack(limit=10) ) )
      self.currentphil, diff_phil = self.GetNewCurrentPhilFromPython(new_phil, self.currentphil)
      #diff = None
      self.params = self.currentphil.extract()
      phl = self.params

      if len(diff_phil.all_definitions()) < 1 and not phl.mouse_moved:
        self.mprint( "Nothing's changed", verbose=1)
        return False

      #diff = diff_phil.extract()
      self.mprint("diff phil:\n" + diff_phil.as_str(), verbose=1 )

      #self.params = self.currentphil.extract()
      #phl = self.params

      if view_3d.has_phil_path(diff_phil, "use_provided_miller_arrays"):
        phl = self.ResetPhilandViewer(self.currentphil)
        if not self.load_miller_arrays():
          return False
        self.viewer.lastscene_id = phl.viewer.scene_id

      if view_3d.has_phil_path(diff_phil, "openfilename"):
        phl = self.ResetPhilandViewer(self.currentphil)
        if not self.load_reflections_file(phl.openfilename):
          return False
        self.viewer.lastscene_id = phl.viewer.scene_id

      if view_3d.has_phil_path(diff_phil, "scene_id", "merge_data", "show_missing", \
         "show_only_missing", "show_systematic_absences", "nbins", "binner_idx",\
         "scene_bin_thresholds"):
        if self.set_scene(phl.viewer.scene_id):
          self.update_space_group_choices()
          self.set_scene_bin_thresholds(strbinvals=phl.scene_bin_thresholds,
                                         binner_idx=phl.binner_idx,
                                         nbins=phl.nbins )
      if phl.spacegroup_choice == None:
        self.mprint("! spacegroup_choice == None")
        #time.sleep(15)

      if view_3d.has_phil_path(diff_phil, "spacegroup_choice"):
        self.set_spacegroup_choice(phl.spacegroup_choice)

      if view_3d.has_phil_path(diff_phil, "tabulate_miller_array_ids"):
        self.tabulate_arrays(phl.tabulate_miller_array_ids)
        #return True

      if view_3d.has_phil_path(diff_phil, "miller_array_operations"):
        self.make_new_miller_array()

      if view_3d.has_phil_path(diff_phil, "using_space_subgroup") and phl.using_space_subgroup==False:
        self.set_default_spacegroup()

      if view_3d.has_phil_path(diff_phil, "shape_primitive"):
        self.set_shape_primitive(phl.shape_primitive)

      if view_3d.has_phil_path(diff_phil, "add_user_vector_hkl_op",
                                         "add_user_vector_abc",
                                         "add_user_vector_hkl"):
        self.add_user_vector()

      if view_3d.has_phil_path(diff_phil, "save_image_name"):
        self.SaveImageName(phl.save_image_name)
        phl.save_image_name = None

      if view_3d.has_phil_path(diff_phil, "action"):
        ret = self.set_action(phl.action)
        phl.action = "is_running" # ensure the same action in succession can be executed
        if not ret:
          return False

      if view_3d.has_phil_path(diff_phil, "savefilename"):
        self.SaveReflectionsFile(phl.savefilename)
      phl.savefilename = None # ensure the same action in succession can be executed

      if view_3d.has_phil_path(diff_phil, "viewer"):
        self.viewer.settings = phl.viewer
        self.settings = phl.viewer

      self.params = self.viewer.update_settings(diff_phil, phl)
      if view_3d.has_phil_path(diff_phil, "scene_id", "spacegroup_choice"):
        self.list_vectors()
      # parameters might have been changed. So update self.currentphil accordingly
      self.currentphil = self.master_phil.format(python_object = self.params)
      self.NewFileLoaded = False
      phl.mouse_moved = False
      self.SendCurrentPhilValues()
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
      if test_flag_value not in array.data():
        return array # for the few cases where a miller array cannot be considered as a valid Rfree array
      array = array.customized_copy(data=(array.data() == test_flag_value))
      array.set_info(info)
      array._data = array.data().as_int()
    return array


  def process_miller_array(self, array) :
    if (array is None) : return
    info = array.info()
    if isinstance(info, str) :
      labels = "TEST DATA"
    else :
      labels = info.label_string()
    if (array.unit_cell() is None) or (array.space_group() is None) :
      raise Sorry("No space group info is present in data")
    details = []
    self.infostr = ""
    array = self.detect_Rfree(array)
    sg = "%s" % array.space_group_info()
    uc = "a=%g b=%g c=%g angles=%g,%g,%g" % array.unit_cell().parameters()
    details_str = ""
    if (len(details) > 0) :
      details_str = "(%s)" % ", ".join(details)
    array_info = group_args(
      labels=labels,
      details_str=details_str,
      merge=self.params.merge_data,
      sg=sg,
      uc=uc)
    return array, array_info


  def process_all_miller_arrays(self, col):
    self.mprint("Processing reflection data...")
    self.procarrays = []
    if self.params.merge_data == False:
      self.settings.expand_to_p1 = False
      self.settings.expand_anomalous = False
    for c,arr in enumerate(self.valid_arrays):
      procarray, procarray_info = self.process_miller_array(arr)
      self.procarrays.append(procarray)
      if c==col:
        array_info = procarray_info
        self.viewer.miller_array = procarray
    if col is None:
      array_info = procarray_info
    return array_info


  def set_miller_array(self, col=None) :
    if col is not None and col >= len(self.viewer.hkl_scenes_info ):
      return
    array_info = self.process_all_miller_arrays(col)
    self.viewer.set_miller_array(col, merge=array_info.merge,
       details=array_info.details_str)
    self.viewer.proc_arrays = self.procarrays
    self.viewer.identify_suitable_fomsarrays()


  def update_space_group_choices(self, col=None) :
    if (self.viewer.miller_array is None and col is None) or \
      self.params.using_space_subgroup:
      return
    if col is None:
      current_miller_array_idx = self.viewer.HKLInfo_from_dict()[1]
    else:
      current_miller_array_idx = col
    matching_valid_array = self.procarrays[ current_miller_array_idx ]
    from cctbx.sgtbx.subgroups import subgroups
    from cctbx import sgtbx
    sg_info = matching_valid_array.space_group_info()
    subgrs = subgroups(sg_info).groups_parent_setting()
    self.spacegroup_choices = []
    for i,subgroup in enumerate(subgrs) :
      subgroup_info = sgtbx.space_group_info(group=subgroup)
      self.spacegroup_choices.append(subgroup_info)
    for i,e in enumerate(self.spacegroup_choices):
      c = None
      if str(sg_info) == str(e):
        self.current_spacegroup = self.spacegroup_choices[i]
        c = i
        break
    if c is None:
      c = 0
      self.spacegroup_choices.insert(c, sg_info)
      self.current_spacegroup = sg_info
    self.params.spacegroup_choice = c
    spglst = [e.symbol_and_number() for e in self.spacegroup_choices] + ["original spacegroup"]
    mydict = { "spacegroups": spglst }
    self.SendInfoToGUI(mydict)


  def set_spacegroup_choice(self, n) :
    if (self.viewer.miller_array is None) :
      raise Sorry("No data loaded!")
    if n == len(self.spacegroup_choices): # selected the unmerged "original spacegroup" in the list
      self.viewer.proc_arrays = self.procarrays
      self.params.using_space_subgroup = False
    else:
      self.current_spacegroup = self.spacegroup_choices[n]
      from cctbx import crystal
      symm = crystal.symmetry(
        space_group_info= self.current_spacegroup,
        unit_cell=self.viewer.miller_array.unit_cell())
      othervalidarrays = []
      for validarray in self.procarrays:
        # TODO: check if array is unmerged i.e. not symmetry unique
        arr = validarray.expand_to_p1().customized_copy(crystal_symmetry=symm)
        arr = arr.merge_equivalents().array().set_info(validarray.info())
        arr = self.detect_Rfree(arr)
        othervalidarrays.append( arr )

      self.mprint( "MERGING 2", verbose=2)
      self.viewer.proc_arrays = othervalidarrays
      self.params.using_space_subgroup = True
    self.viewer.set_miller_array()
    for i,e in enumerate(self.spacegroup_choices):
      self.mprint("%d, %s" %(i,e.symbol_and_number()) , verbose=0)


  def SetSpaceGroupChoice(self, n):
    self.params.spacegroup_choice = n
    self.update_settings()


  def SetDefaultSpaceGroup(self):
    self.params.using_space_subgroup = False
    self.update_settings()


  def set_default_spacegroup(self):
    self.viewer.proc_arrays = self.procarrays
    self.viewer.set_miller_array()
    self.viewer.identify_suitable_fomsarrays()


  def MakeNewMillerArrayFrom(self, operation, label, arrid1, arrid2=None):
    # get list of existing new miller arrays and operations if present
    miller_array_operations_lst = []
    #if self.params.miller_array_operations:
    #  miller_array_operations_lst = eval(self.params.miller_array_operations)
    miller_array_operations_lst = [ ( operation, label, arrid1, arrid2 ) ]
    self.params.miller_array_operations = str( miller_array_operations_lst )
    self.update_settings()


  def make_new_miller_array(self):
    miller_array_operations_lst = eval(self.params.miller_array_operations)
    unique_miller_array_operations_lst = []
    for (operation, label, arrid1, arrid2) in miller_array_operations_lst:
      for arr in self.procarrays:
        if label in arr.info().labels + [ "", None]:
          raise Sorry("Provide an unambiguous label for your new miller array!")
      unique_miller_array_operations_lst.append( (operation, label, arrid1, arrid2) )
    self.params.miller_array_operations = str(unique_miller_array_operations_lst)
    from copy import deepcopy
    millarr1 = deepcopy(self.procarrays[arrid1])
    newarray = None
    if arrid2 != -1:
      millarr2 = deepcopy(self.procarrays[arrid2])
      newarray = self.viewer.OperateOn2MillerArrays(millarr1, millarr2, operation)
    else:
      newarray = self.viewer.OperateOn1MillerArray(millarr1, operation)
    if newarray is not None:
      self.mprint("New dataset has %d reflections." %newarray.size())
      newarray.set_info(millarr1._info )
      newarray._info.labels = [ label ]
      procarray, procarray_info = self.process_miller_array(newarray)
      self.procarrays.append(procarray)
      self.viewer.proc_arrays = self.procarrays
      self.viewer.has_new_miller_array = True
      self.viewer.array_infostrs.append( ArrayInfo(procarray, self.mprint).infostr )
      self.viewer.array_infotpls.append( ArrayInfo(procarray, self.mprint).infotpl )
      #self.viewer.SupersetMillerArrays()

      hkls = self.origarrays["HKLs"]
      nanarr = flex.double(len(hkls), float("nan"))
      m = miller.match_indices(hkls, procarray.indices() )
      indices_of_matched_hkls = m.pairs().column(0)
      for i,e in enumerate(indices_of_matched_hkls):
        nanarr[e] = procarray.data()[i]
      self.origarrays[label] = list(nanarr)
      mydict = { "array_infotpls": self.viewer.array_infotpls, "NewHKLscenes" : True, "NewMillerArray" : True}
      self.SendInfoToGUI(mydict)


  def prepare_dataloading(self):
    self.viewer.isnewfile = True
    #self.params.mergedata = None
    self.params.viewer.scene_id = None
    self.viewer.colour_scene_id = None
    self.viewer.radii_scene_id = None
    self.viewer.match_valarrays = []
    self.viewer.proc_arrays = {}
    self.spacegroup_choices = []
    self.origarrays = {}
    display.reset_settings()
    self.settings = display.settings()
    self.viewer.settings = self.params.viewer
    self.viewer.mapcoef_fom_dict = {}
    self.viewer.sceneid_from_arrayid = []
    self.hklfile_history = []
    self.tncsvec = None
    self.loaded_file_name = ""


  def finish_dataloading(self, arrays):
    valid_arrays = []
    self.viewer.array_infostrs = []
    self.viewer.array_infotpls = []
    spg = arrays[0].space_group()
    for i,array in enumerate(arrays):
      if array.space_group() is None:
        array._space_group_info = spg.info()
        self.mprint("""No space group info present in the %d. miller array.
\nBorrowing space group info from the first miller array""" %i)
      arrayinfo = ArrayInfo(array, self.mprint)
      self.viewer.array_infostrs.append( arrayinfo.infostr )
      self.viewer.array_infotpls.append( arrayinfo.infotpl )
      if i==0:
        mydict = { "spacegroup_info": arrayinfo.spginf, "unitcell_info": arrayinfo.ucellinf }
        self.SendInfoToGUI(mydict)
      valid_arrays.append(array)
    self.valid_arrays = valid_arrays
    self.mprint("%d Miller arrays in this dataset:" %len(arrays))
    for e in self.viewer.array_infostrs:
      self.mprint("%s" %e)
    self.mprint("\n")
    self.NewFileLoaded = True
    if (len(valid_arrays) == 0):
      msg = "No arrays of the supported types present."
      self.mprint(msg)
      self.NewFileLoaded=False
    elif (len(valid_arrays) >= 1):
      self.set_miller_array()
      self.update_space_group_choices(0) # get the default spacegroup choice
      mydict = { "info": self.infostr,
                  "array_infotpls": self.viewer.array_infotpls,
                  "bin_infotpls": self.viewer.bin_infotpls,
                  "html_url": self.viewer.url,
                  "tncsvec": self.tncsvec,
                  "merge_data": self.params.merge_data,
                  "spacegroups": [e.symbol_and_number() for e in self.spacegroup_choices],
                  "NewFileLoaded": self.NewFileLoaded,
                  "file_name": self.params.openfilename
                }
      self.SendInfoToGUI(mydict)
    self.params.openfilename = None


  def load_reflections_file(self, file_name):
    file_name = to_str(file_name)
    ret = False
    if (file_name != ""):
      try :
        self.mprint("Reading file...")
        self.prepare_dataloading()
        hkl_file = any_reflection_file(file_name)
        if hkl_file._file_type == 'cif':
          # use new cif label parser for reflections
          cifreader = hkl_file.file_content()
          arrays = cifreader.as_miller_arrays(merge_equivalents=False, style="new")
          # sanitise labels by removing redundant strings.
          # remove the data name of this cif file from all labels
          dataname = list(hkl_file._file_content.builder._model.keys())
          unwantedstrings = dataname[:]
          # remove "_refln." from all labels
          unwantedstrings.append("_refln.")
          unwantedstrings.append("_refln_")
          for arr in arrays:
            if len(arr.info().labels):
              newlabels = []
              for label in arr.info().labels:
                found = False
                for s in unwantedstrings:
                  if s in label:
                    newlabel = label.replace(s, "")
                    found = True
                    if len(newlabel) > 0:
                      newlabels.append(newlabel)
                    break
                if not found:
                  newlabels.append(label)
                arr.info().labels = newlabels
          self.origarrays = cifreader.as_original_arrays()[dataname[0]]
        else: # some other type of reflection file than cif
          arrays = hkl_file.as_miller_arrays(merge_equivalents=False)
        if hkl_file._file_type == 'ccp4_mtz':
          self.hklfile_history = list(hkl_file._file_content.history())
          self.loaded_file_name = file_name
          for e in self.hklfile_history:
            if "TNCS NMOL" in e and "VECTOR" in e:
              svec = e.split()[-3:]
              t1 = float(svec[0])
              t2 = float(svec[1])
              t3 = float(svec[2])
              if (t1*t1 + t2*t2 + t3*t3) > 0.0:
                self.tncsvec = (t1, t2, t3)
                self.mprint("tNCS vector found in header of mtz file: %s" %str(self.tncsvec) )
          from iotbx import mtz
          mtzobj = mtz.object(file_name)
          nanval = float("nan")
          self.origarrays["HKLs"] = mtzobj.extract_miller_indices()
          for mtzlbl in mtzobj.column_labels():
            col = mtzobj.get_column( mtzlbl )
            newarr = col.extract_values_and_selection_valid().values.deep_copy()
            for i,b in enumerate(col.extract_values_and_selection_valid().selection_valid):
              if not b:
                newarr[i] = nanval
            self.origarrays[mtzlbl] = list(newarr)
        self.finish_dataloading(arrays)
      except Exception as e :
        self.NewFileLoaded=False
        self.mprint("".join(traceback.format_tb(e.__traceback__ )) + e.__repr__())
        arrays = []
      ret = True
    return ret


  def LoadReflectionsFile(self, openfilename):
    self.params.openfilename = openfilename
    self.update_settings()


  def load_miller_arrays(self):
    ret = False
    try:
      self.ResetPhilandViewer(self.currentphil)
      self.prepare_dataloading()
      self.finish_dataloading(self.provided_miller_arrays)
      ret = True
    except Exception as e :
      self.NewFileLoaded=False
      self.mprint("".join(traceback.format_tb(e.__traceback__ )) + e.__repr__())
      arrays = []
    return ret


  def LoadMillerArrays(self, marrays):
    self.provided_miller_arrays = marrays
    self.params.use_provided_miller_arrays = True
    self.update_settings()


  def SaveReflectionsFile(self, savefilename):
    if self.loaded_file_name == savefilename:
      self.mprint("Not overwriting currently loaded file. Choose a different name!")
      return
    self.mprint("Saving file...")
    fileextension = os.path.splitext(savefilename)[1]
    if fileextension == ".mtz":
      mtz1 = self.viewer.proc_arrays[0].as_mtz_dataset(column_root_label= self.viewer.proc_arrays[0].info().labels[0])
      for i,arr in enumerate(self.viewer.proc_arrays):
        if i==0:
          continue
        mtz1.add_miller_array(arr, column_root_label=arr.info().labels[0] )
      try: # python2 or 3
        mtz1.mtz_object().write(savefilename)
      except Exception as e:
        mtz1.mtz_object().write(savefilename.encode("ascii"))
      self.mprint("Miller array(s) saved to: " + savefilename)
    elif fileextension == ".cif":
      import iotbx.cif
      mycif = None
      fname = savefilename
      fnames = []

      def save2cif(filename, mycif):
        with open(filename.encode("ascii"), "w") as f:
          f.write("data_%s\n#\n" %os.path.splitext(os.path.basename(filename))[0])
          print(mycif.cif_block, file= f)

      for i,arr in enumerate(self.viewer.proc_arrays):
        arrtype = None
        colnames = ["_refln.%s" %e for e in arr.info().labels ]
        colname= None
        if self.has_indices_with_multiple_data(arr):
          # if array contains data with more than one data point for the same hkl index iotbx.cif
          # cannot add additional arrays to the cif block so save this array in a separate file
          singlecif = iotbx.cif.miller_arrays_as_cif_block(arr, array_type = arrtype,
                                                       column_name=colname, column_names = colnames )
          fname = os.path.splitext(savefilename)[0] + "_%d"%i + os.path.splitext(savefilename)[1]
          save2cif(fname, singlecif)
          fnames.append(fname)
          continue
        if not mycif:
          mycif = iotbx.cif.miller_arrays_as_cif_block(arr, array_type = arrtype,
                                                       column_name=colname, column_names = colnames )
        else:
          mycif.add_miller_array(arr, column_name= colname, array_type= arrtype,
                                 column_names = colnames)
      if mycif:
        save2cif(savefilename, mycif)
        fnames.append(savefilename)
      self.mprint("Miller array(s) saved to: " + ",\n".join(fnames))
      if len(fnames) > 1:
        self.mprint("Unmerged data put into separate files")
    else:
      self.mprint("Can only save file in MTZ or CIF format. Sorry!")


  def has_indices_with_multiple_data(self, arr):
    return len(set(list(arr.indices()))) < arr.size()


  def tabulate_arrays(self, datalabels):
    if len(self.origarrays) == 0: # if not an mtz file then split columns
      # SupersetMillerArrays may not be necessary if file formats except for cif and mtz can't store multiple data columns
      #self.viewer.SupersetMillerArrays()
      self.origarrays["HKLs"] = self.viewer.proc_arrays[0].indices()
      for arr in self.viewer.proc_arrays:
        if arr.is_complex_array():
          ampls, phases = self.viewer.Complex2AmplitudesPhases(arr.data())
          cmplxlst = [ "%.4f + %.4f * i"%(e.real, e.imag)
                        if not cmath.isnan(e) else display.nanval for e in arr.data() ]
          self.origarrays[arr.info().label_string()] = cmplxlst
          self.origarrays[arr.info().labels[0]] = list(ampls)
          self.origarrays[arr.info().labels[-1]] = list(phases)
        elif arr.is_hendrickson_lattman_array():
          A,B,C,D = arr.data().as_abcd()
          HLlst = [ "%.4f, %.4f, %.4f, %.4f"%(e[0], e[1], e[2], e[3]) for e in arr.data() ]
          self.origarrays[arr.info().label_string()] = HLlst
          self.origarrays[arr.info().labels[0]] = list(A)
          self.origarrays[arr.info().labels[1]] = list(B)
          self.origarrays[arr.info().labels[2]] = list(C)
          self.origarrays[arr.info().labels[3]] = list(D)
        elif arr.sigmas() is not None:
          labels = arr.info().labels
          # Labels could be something like ['I(+)', 'SIGI(+)', 'I(-)', 'SIGI(-)'].
          # So group datalabels and sigmalabels separately assuming that sigma column contain the three letters "sig"
          datalabel = ",".join([ e for e in labels if "sig" not in e.lower()])
          sigmalabel = ",".join([ e for e in labels if "sig" in e.lower()])
          self.origarrays[datalabel] = list(arr.data())
          self.origarrays[sigmalabel] = list(arr.sigmas())
        elif arr.is_integer_array():
          list_with_nans = [ e if not e==display.inanval else display.nanval for e in arr.data() ]
          if self.viewer.array_infotpls[id][0] == 'FreeR_flag': # want True or False back
            list_with_nans = [ 1==e if not cmath.isnan(e) else display.nanval for e in list_with_nans ]
          self.origarrays[arr.info().label_string()] = list_with_nans
        else:
          self.origarrays[arr.info().label_string()] = list(arr.data())

    labels = eval(datalabels)
    indices = self.origarrays["HKLs"]
    dres = self.procarrays[0].unit_cell().d( indices)
    dreslst = [("d_res", roundoff(list(dres)),3)]
    hkls = list(indices)
    hkllst = [ ("H", [e[0] for e in hkls] ), ("K", [e[1] for e in hkls] ), ("L", [e[2] for e in hkls] )]
    datalst = []
    for label in labels:
      datalst.append( (label, list(self.origarrays[label])))
    self.idx_data = hkllst + dreslst + datalst
    self.mprint("Sending table data...", verbose=0)
    mydict = { "tabulate_miller_array": self.idx_data }
    self.params.tabulate_miller_array_ids = "[]" # to allow reopening a closed window again
    self.SendInfoToGUI(mydict)


  def TabulateMillerArray(self, ids):
    self.params.tabulate_miller_array_ids = str(ids)
    self.update_settings()


  def SetCameraType(self, camtype):
    self.params.NGL.camera_type = camtype
    self.update_settings()


  def ExpandToP1(self, val, inbrowser=True):
    self.params.viewer.expand_to_p1 = val
    self.params.viewer.inbrowser = inbrowser
    self.update_settings()


  def ExpandAnomalous(self, val, inbrowser=True):
    self.params.viewer.expand_anomalous = val
    self.params.viewer.inbrowser = inbrowser
    self.update_settings()


  def ShowOnlyMissing(self, val):
    self.params.viewer.show_only_missing = val
    self.update_settings()


  def ShowMissing(self, val):
    self.params.viewer.show_missing = val
    self.update_settings()


  def ShowDataOverSigma(self, val):
    self.params.viewer.show_data_over_sigma = val
    self.update_settings()


  def ShowSystematicAbsences(self, val):
    self.params.viewer.show_systematic_absences = val
    self.update_settings()


  def ShowSlice(self, val, axis="h", index=0):
    axisstr = axis.lower()
    self.params.viewer.slice_mode = val
    self.params.viewer.slice_axis = axisstr
    self.params.viewer.slice_index = index
    self.update_settings()


  def set_scene_bin_thresholds(self, strbinvals = "", binner_idx = 0,  nbins = 6):
    nuniquevalues = -1
    if not strbinvals:
      binvals, nuniquevalues = self.viewer.calc_bin_thresholds(binner_idx, nbins)
    else:
      nan = float("nan")
      binvals = eval(strbinvals)
    if binvals and binner_idx == 0:
      binvals = list( 1.0/flex.double(binvals) )
    self.viewer.UpdateBinValues(binner_idx, binvals, nuniquevalues)


  def SetSceneNbins(self, nbins, binner_idx = 0):
    self.params.nbins = nbins
    self.params.binner_idx = binner_idx
    self.params.NGL.bin_opacities = str([ (1.0, e) for e in range(nbins) ])
    self.update_settings()


  def GetNumberingOfBinners(self):
    return [ (i,e) for i,e in enumerate(self.viewer.bin_labels_type_idxs) ]


  def SetSceneBinThresholds(self, binvals=[]):
    self.params.scene_bin_thresholds = str(binvals)
    self.params.nbins = len(binvals)
    self.update_settings()


  def SetOpacities(self, bin_opacities):
    self.params.NGL.bin_opacities = str(bin_opacities)
    self.update_settings()


  def SetToolTipOpacity(self, val):
    self.params.NGL.tooltip_alpha = val
    self.update_settings()


  def SetShowToolTips(self, val):
    self.params.NGL.show_tooltips = val
    self.update_settings()


  def set_scene(self, scene_id):
    self.viewer.binvals = []
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
    self.params.viewer.scene_id = scene_id
    self.update_settings()


  def SetMergeData(self, val):
    self.params.merge_data = val
    self.update_settings()


  def SetColourScene(self, colourcol):
    self.params.viewer.colour_scene_id = colourcol
    self.update_settings()


  def SetRadiusScene(self, radiuscol):
    self.params.viewer.radii_scene_id = radiuscol
    self.update_settings()


  def SetRadiiScale(self, scale=1.0, nth_power_scale = -1.0):
    """
    Scale radii. Decrease the contrast between large and small radii with nth_root_scale < 1.0
    If nth_power_scale=0.0 then all radii will have the same size regardless of data values.
    If nth_power_scale < 0.0 an automatic power will be computed ensuring the smallest radius
    is 0.1 times the maximum radius
    """
    self.params.viewer.scale = scale
    self.params.viewer.nth_power_scale_radii = nth_power_scale
    self.update_settings()


  def SetColourRadiusToSigmas(self, val):
    self.params.viewer.sigma_color_radius = val
    self.update_settings()


  def SetColourScheme(self, color_scheme, color_powscale=1.0):
    self.params.viewer.color_scheme = color_scheme
    self.params.viewer.color_powscale = color_powscale
    self.update_settings()


  def SetShapePrimitive(self, val):
    self.params.shape_primitive = val
    self.update_settings()


  def set_shape_primitive(self, val):
    if val == "points":
      self.viewer.primitivetype = "PointBuffer"
    else:
      self.viewer.primitivetype = "sphereBuffer"


  def SetAction(self, val):
    self.params.action = val
    self.update_settings()


  def set_action(self, val):
    if val == "reset_view":
      self.viewer.SetAutoView()
    if val == "is_terminating":
      self.__exit__()
      return False
    return True


  def SetFontSize(self, val):
    self.params.NGL.fontsize = val
    self.viewer.SetFontSize(val)


  def list_vectors(self):
    self.viewer.all_vectors = self.viewer.rotation_operators[:]
    if self.tncsvec is not None:
      uc = self.viewer.miller_array.unit_cell()
      # TNCS vector is specified in realspace fractional coordinates. Convert it to cartesian
      cartvec = list( self.tncsvec * matrix.sqr(uc.orthogonalization_matrix()) )
      ln = len(self.viewer.all_vectors)
      self.viewer.all_vectors.append( (ln, "TNCS", 0, cartvec, "", "", str(roundoff(self.tncsvec, 5)) ) )
    self.viewer.all_vectors = self.viewer.all_vectors + self.uservectors
    for (opnr, label, order, cartvec, hkl_op, hkl, abc) in self.viewer.all_vectors:
      # avoid onMessage-DrawVector in HKLJavaScripts.js misinterpreting the commas in strings like "-x,z+y,-y"
      name = label + hkl_op.replace(",", "_")
      self.viewer.RemovePrimitives(name)
    self.SendInfoToGUI( { "all_vectors": self.viewer.all_vectors } )
    return self.viewer.all_vectors


  def add_user_vector(self):
    uc = self.viewer.miller_array.unit_cell()
    ln = len(self.viewer.all_vectors)
    label = self.params.viewer.user_label
    order = 0
    try:
      hklvec = ""
      abcvec = ""
      hklop = ""
      unwantedchars = " |(|)|[|]|{|}"
      # individual characters separated by | substituted with a "" using re.sub()
      if self.params.viewer.add_user_vector_hkl not in [None, "", "()"]:
        hklvec = eval(re.sub(unwantedchars, "", self.params.viewer.add_user_vector_hkl))
        # convert into cartesian space
        cartvec = list( self.viewer.scene.renderscale*(hklvec * matrix.sqr(uc.fractionalization_matrix()).transpose()) )
      elif self.params.viewer.add_user_vector_abc not in [None, "", "()"]:
        abcvec = eval(re.sub(unwantedchars, "", self.params.viewer.add_user_vector_abc))
        # convert into cartesian space
        cartvec = list(abcvec * matrix.sqr(uc.orthogonalization_matrix()))
      elif self.params.viewer.add_user_vector_hkl_op not in [None, ""]:
        hklop = re.sub(unwantedchars, "", self.params.viewer.add_user_vector_hkl_op)
        rt = sgtbx.rt_mx(symbol=hklop, r_den=12, t_den=144)
        self.viewer.symops.append( rt ) #
        (cartvec, a, label, order) = self.viewer.GetVectorAndAngleFromRotationMx( rt.r() )
        if label:
          label = "%s-fold_%s" %(str(int(roundoff(2*math.pi/a, 0))), self.params.viewer.user_label)
          self.mprint("Rotation axis, %s, added" %label)
        if label =="" or order==0:
          self.mprint("Cannot compute a rotation axis from %s" %self.params.viewer.add_user_vector_hkl_op)
          return
      if (self.params.viewer.add_user_vector_hkl in [None, "", "()"] \
       and self.params.viewer.add_user_vector_abc in [None, "", "()"] \
       and self.params.viewer.add_user_vector_hkl_op) in [None, ""]:
        self.mprint("No vector was specified")
      self.uservectors.append( (ln, label, order, cartvec, hklop, str(hklvec), str(abcvec) ))
      self.list_vectors()
    except Exception as e:
      raise Sorry( str(e))
    self.params.viewer.add_user_vector_hkl_op = ""
    self.params.viewer.add_user_vector_hkl = ""
    self.params.viewer.add_user_vector_abc = ""


  def AddUserVector(self, hkl_op="", abc="", hkl="", label=""):
    """
    Vector can be specified as a rotation operator, say "-h-k,k,-l" subject to spacegroup contraints,
    as a fractional vector in real space or as a fractional vector in reciprocal space. If
    specified as a rotation operator the derived vector is the implicit rotation axis.
    """
    self.params.viewer.user_label = label
    self.params.viewer.add_user_vector_hkl_op = str(hkl_op)
    self.params.viewer.add_user_vector_abc = str(abc)
    self.params.viewer.add_user_vector_hkl = str(hkl)
    self.update_settings()


  def ShowRotationAxes(self, val):
    self.params.viewer.show_symmetry_rotation_axes = val
    self.update_settings()


  def ShowVector(self, i, val=True):
    self.params.viewer.show_vector = str([i, val])
    self.update_settings()


  def ShowUnitCell(self, val):
    self.params.show_real_space_unit_cell = val
    self.update_settings()


  def ShowReciprocalUnitCell(self, val):
    self.params.show_reciprocal_unit_cell = val
    self.update_settings()


  def SetClipPlane(self, use=True, hkldist=0.0, clipwidth=2.0):
    if use:
      self.params.clip_plane.hkldist = hkldist
      self.params.clip_plane.clipwidth = clipwidth
      self.params.slice_mode = False
      self.params.inbrowser = True
    else:
      self.params.clip_plane.clipwidth = None
    self.update_settings()


  def SinglePlaneOfReflections(self, use=True, axis="h", slice_index=0 ):
    if use:
      viewer.slice_axis = axis
      viewer.is_parallel = False
      viewer.slice_mode = True
      viewer.inbrowser = False
      viewer.fixorientation = "reflection_slice"
      viewer.slice_index = slice_index
    else:
      viewer.slice_mode = False
      viewer.inbrowser = True
      viewer.fixorientation = "None"
    self.update_settings()


  def OrientVector(self, vecnr, is_parallel, val=True):
    viewer.fixorientation = "None"
    if val:
      viewer.is_parallel = is_parallel
      viewer.fixorientation = "vector"
      viewer.show_vector = '[%d, True]' %vecnr
    self.update_settings()


  def AnimateRotateAroundVector(self, vecnr, speed):
    self.params.clip_plane.animate_rotation_around_vector = str([vecnr, speed])
    self.update_settings()


  def RotateAroundVector(self, vecnr, dgr):
    self.params.clip_plane.angle_around_vector = str([vecnr, dgr])
    self.update_settings()


  def ShowHKL(self, hkl):
    self.params.viewer.show_hkl = str(hkl)
    self.update_settings()


  def SetMouseSpeed(self, trackspeed):
    self.params.NGL.mouse_sensitivity = trackspeed
    self.update_settings()


  def GetMouseSpeed(self):
    self.viewer.GetMouseSpeed()
    return self.params.NGL.mouse_sensitivity


  def GetSpaceGroupChoices(self):
    """
    return array of strings with available subgroups of the space group
    """
    if (self.viewer.miller_array is None) :
      self.mprint( NOREFLDATA)
    if self.spacegroup_choices:
      return [e.symbol_and_number() for e in self.spacegroup_choices]
    return []


  def SaveImageName(self, fname):
    self.viewer.MakeImage(fname)


  def SendCurrentPhilValues(self):
    philstrvalsdict = {}
    for e in self.currentphil.all_definitions():
      philstrvalsdict[e.path] = e.object.extract()
    mydict = { "current_phil_strings": philstrvalsdict }
    self.SendInfoToGUI(mydict)
    if self.viewer.params.viewer.scene_id is not None:
      self.SendInfoToGUI({ "used_nth_power_scale_radii": self.viewer.HKLscene_from_dict().nth_power_scale_radii })


  def GetHtmlURL(self):
    return self.viewer.url


  def GetHtmlstring(self):
    return self.viewer.htmlstr


  def GetArrayInfotpls(self):
    """
    return array of tuples with information on each miller array
    """
    return self.viewer.array_infotpls


  def GetSceneDataLabels(self):
    return [ e[3][0] for e in myHKLview.viewer.hkl_scenes_infos ]


  def GetHklScenesInfos(self):
    """
    return array of strings with information on each processed miller array
    which may have been expanded with anomalous reflections or truncated to non-anomalous reflections
    as to match the currently selected miller array
    """
    return self.viewer.hkl_scenes_infos


  def GetBinInfo(self):
    """
    return array of number of hkls and bin boundaries of the bins the current miller array data has been sorted into.
    Useful when deciding which bin of reflections to make transparent
    """
    return self.viewer.binstrs


  def SendInfoToGUI(self, infodict, binary=True):
    if self.guiSocketPort:
      m = str(infodict).encode("utf-8")
      if not binary:
        self.guisocket.send( m )
      else:
        if type(m) is not bytes:
          m = bytes(m)
        bindict = zlib.compress( m )
        self.guisocket.send( bindict )


masterphilstr = """
  openfilename = None
    .type = path
  use_provided_miller_arrays = False
    .type = bool
  savefilename = None
    .type = path
  save_image_name = None
    .type = path
  merge_data = False
    .type = bool
  miller_array_operations = ''
    .type = str
  spacegroup_choice = 0
    .type = int
  using_space_subgroup = False
    .type = bool
  mouse_moved = False
    .type = bool
  real_space_unit_cell_scale_fraction = None
    .type = float
  reciprocal_unit_cell_scale_fraction = None
    .type = float
  clip_plane {
    angle_around_vector = \"[0,0]\"
      .type = str
    animate_rotation_around_vector = \"[0,0]\"
      .type = str
    hkldist = 0.0
      .type = float
    clipwidth = None
      .type = float
    fractional_vector = reciprocal *realspace
      .type = choice
    bequiet = False
      .type = bool
  }
  scene_bin_thresholds = ''
    .type = str
  binner_idx = 0
    .type = int
  nbins = 1
    .type = int(value_min=1, value_max=40)
  shape_primitive = *'spheres' 'points'
    .type = choice
  viewer {
    scene_id = None
      .type = int
    ncolourlabels = 6
      .type = int
    show_symmetry_rotation_axes = False
      .type = bool
    show_vector = ''
      .type = str
    add_user_vector_hkl_op = ""
      .type = str
    add_user_vector_abc = ""
      .type = str
    add_user_vector_hkl = ""
      .type = str
    user_label = ""
      .type = str
    show_hkl = ""
      .type = str
    is_parallel = False
      .type = bool
    fixorientation = vector reflection_slice *None
      .type = choice
    angle_around_XHKL_vector = 0.0
      .type = float
    angle_around_YHKL_vector = 0.0
      .type = float
    angle_around_ZHKL_vector = 0.0
      .type = float
    %s
  }
  NGL {
    %s
  }
  action = *is_running is_terminating reset_view
    .type = choice
  tabulate_miller_array_ids = "[]"
    .type = str

""" %(display.philstr, view_3d.ngl_philstr)


def run():
  """
  utility function for passing keyword arguments more directly to HKLViewFrame()
  """
  #time.sleep(15)
  # dirty hack for parsing a file path with spaces of a browser if not using default
  args = sys.argv[1:]
  sargs = " ".join(args)
  qchar = "'"
  if sargs.find("'") > -1:
    quote1 = sargs.find(qchar)
    if sargs[ quote1 + 1:].find(qchar) < 0:
      raise Sorry("Missing quote in arguments")
    quote2 = sargs[ quote1 + 1:].find(qchar) + quote1 + 1
    space1 = sargs[ :quote1].rfind(" ")
    arg = sargs[space1 +1: quote2 +1]
    sargs2 = sargs.replace(arg,"")
    args = sargs2.split(" ")
    arg = arg.replace("'","")
    arg = arg.replace('"',"")
    arg = arg.replace('\\', '/') # webbrowser module wants browser paths having unix forward slashes
    args.append(arg)

  kwargs = dict(arg.split('=') for arg in args if '=' in arg)
  #check if any argument is a filename
  for arg in args:
    # if so add it as a keyword argument
    if os.path.isfile(arg) and '=' not in arg:
      kwargs['hklin'] = arg


  myHKLview = HKLViewFrame(**kwargs)


if __name__ == '__main__':
  run()
