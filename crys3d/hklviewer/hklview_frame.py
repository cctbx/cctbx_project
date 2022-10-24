
from __future__ import absolute_import, division, print_function

from iotbx.reflection_file_reader import any_reflection_file
from cctbx.xray import observation_types
from iotbx.gui_tools.reflections import ArrayInfo
from crys3d.hklviewer import display2 as display
from crys3d.hklviewer import jsview_3d
from cctbx import miller
from cctbx import crystal
from libtbx.math_utils import roundoff
from cctbx.array_family import flex
from libtbx.utils import Sorry, to_str
from scitbx import matrix
from cctbx import sgtbx
from libtbx import group_args, version
import libtbx
import libtbx.load_env
import traceback
import sys, zmq, threading,  time, cmath, zlib, os.path, math, re
from pathlib import Path
import importlib.util



NOREFLDATA = "No reflection data has been selected"

class HKLViewFrame() :
  def __init__ (self, *args, **kwds) :
    self.valid_arrays = []
    self.spacegroup_choices = []
    self.procarrays = []
    self.origarrays = {}
    self.ano_spg_tpls =[]
    self.merge_answer = [None]
    self.dmin = -1
    self.verbose = 1
    if 'verbose' in kwds:
      try:
        self.verbose = eval(kwds['verbose'])
      except Exception as e:
        self.verbose = kwds['verbose']
    self.guiSocketPort=None
    kwds['mprint'] = self.mprint
    self.outputmsgtypes = []
    self.userpresetbuttonsfname = os.path.join( Path.home(), ".hkl_viewer_buttons.py")
    self.infostr = ""
    self.allbuttonslist = []
    self.hklfile_history = []
    self.arrayinfos = []
    self.tncsvec = None
    self.aniso1 = None
    self.aniso2 = None
    self.aniso3 = None
    self.uservectors = []
    self.new_miller_array_operations_lst = []
    cfilename = libtbx.env.under_root(os.path.join("modules","cctbx_project","COPYRIGHT.txt"))
    if not os.path.exists( cfilename ):
      cfilename = libtbx.env.under_root(os.path.join("cctbx", "COPYRIGHT.txt")) # conda installations
    self.copyrightpaths = [
     ("CCTBX copyright", cfilename),
     ("NGL copyright", os.path.join(os.path.dirname(os.path.abspath(__file__)), "License_for_NGL.txt")),
     ("html2canvas copyright", os.path.join(os.path.dirname(os.path.abspath(__file__)), "LICENSE_for_html2canvas.txt"))
    ]
    self.zmqsleeptime = 0.1
    buttonsdeflist = []
    if 'useGuiSocket' in kwds:
      self.guiSocketPort = eval(kwds['useGuiSocket'])
      self.context = zmq.Context()
      self.guisocket = self.context.socket(zmq.PAIR)
      self.guisocket.connect("tcp://127.0.0.1:%s" %self.guiSocketPort )
      self.STOP = False
      self.mprint("CCTBX process with pid: %s starting socket thread" %os.getpid(), verbose=1)
      # name this thread to ensure any asyncio functions are called only from main thread
      self.msgqueuethrd = threading.Thread(target = self.zmq_listen, name="HKLviewerZmqThread" )
      self.msgqueuethrd.daemon = True
      kwds['send_info_to_gui'] = self.SendInfoToGUI # function also used by HKLjsview_3d
      pyversion = "cctbx.python.version: " + str(sys.version_info[0])
      # tell gui what python version we are
      self.SendInfoToGUI(pyversion )
      self.SendInfoToGUI({"copyrights": self.copyrightpaths,
                          "cctbxversion": version.get_version()} )
      try:
        from .preset_buttons import cctbx_buttonsdeflist
        buttonsdeflist = cctbx_buttonsdeflist
        try:
          from phasertng.scripts import xtricorder # then we are in phenix and can load phasertng
          from .preset_buttons import phenix_buttonsdeflist
          buttonsdeflist.extend(phenix_buttonsdeflist) # add phenix buttons to Quick View list
        except Exception as e: # otherwise we only have cctbx
          from .preset_buttons import cctbx_buttonsdeflist
          buttonsdeflist = cctbx_buttonsdeflist
      except Exception as e: # don't even have cctbx! Should never get here!
        buttonsdeflist = []
      if "phenix_buttonsdeflist" in dir():
        self.SendInfoToGUI({"AddPhenixButtons": True})
    self.mprint("kwds= " +str(kwds), verbose=1)
    self.mprint("args= " + str(args), verbose=1)
    kwds['websockport'] = self.find_free_port()
    kwds['parent'] = self
    self.viewer = jsview_3d.HKLview_3d( **kwds )
    self.allbuttonslist = buttonsdeflist
    if os.path.exists(self.userpresetbuttonsfname):
      self.mprint("Using user defined preset-buttons from " + self.userpresetbuttonsfname, verbose=1)
    else:
      factorydefault_userbutton_fname = os.path.join(os.path.split(jsview_3d.__file__)[0], "default_user_preset_buttons.py")
      import shutil
      shutil.copyfile(factorydefault_userbutton_fname, self.userpresetbuttonsfname)
      self.mprint("New UserPresetButton file copied to " + self.userpresetbuttonsfname)
      self.mprint("Taylor button definitions to your own needs.")

    spec = importlib.util.spec_from_file_location("UserPresetButtons", self.userpresetbuttonsfname)
    UserPresetButtons_module = importlib.util.module_from_spec(spec)
    sys.modules["UserPresetButtons"] = UserPresetButtons_module
    try:
      spec.loader.exec_module(UserPresetButtons_module)
      self.allbuttonslist = buttonsdeflist + UserPresetButtons_module.buttonsdeflist
    except Exception as e:
      self.mprint( str(e) + traceback.format_exc(limit=10))
    self.ResetPhilandViewer()
    self.firsttime = True
    self.idx_data = None
    self.clipper_crystdict = None
    self.NewFileLoaded = False
    self.loaded_file_name = ""
    self.validated_preset_buttons = False
    self.fileinfo = None
    if 'fileinfo' in kwds:
      self.fileinfo = kwds.get('fileinfo', 1 )
    self.hklin = None
    if 'hklin' in kwds or 'HKLIN' in kwds:
      self.hklin = kwds.get('hklin', kwds.get('HKLIN') )
    self.LoadReflectionsFile(self.hklin)
    if 'useGuiSocket' in kwds:
      self.msgqueuethrd.start()
    self.validate_preset_buttons()
    if 'show_master_phil' in args:
      self.mprint("Default PHIL parameters:\n" + "-"*80 + "\n" + master_phil.as_str(attributes_level=2) + "-"*80)


  def __exit__(self, exc_type=None, exc_value=0, traceback=None):
    self.viewer.__exit__(exc_type, exc_value, traceback)
    self.mprint("Destroying HKLViewFrame", verbose=0) # this string is expected by HKLviewer.py so don't change
    self.STOP = True
    del self
    #sys.exit()


  def mprint(self, msg, verbose=0, end="\n"):
    if self.guiSocketPort:
      if  verbose == 0:
        # say verbose="2threading" then print all messages with verbose=2 or verbose=threading
        self.SendInfoToGUI( { "info": msg + end } )
      if  (isinstance(self.verbose,int) and isinstance(verbose,int) and verbose >= 1 and verbose <= self.verbose) \
       or (isinstance(self.verbose,str) and self.verbose.find(str(verbose))>=0 ):
        # say verbose="2threading" then print all messages with verbose=2 or verbose=threading
        self.SendInfoToGUI( { "alert": msg + end } )
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
    nan = float("nan") # workaround for "evaluating" any NaN values in the messages received
    lastmsgtype = ""
    while not self.STOP:
      try:
        msgstr = self.guisocket.recv().decode("utf-8")
        if msgstr == "":
          continue
        self.mprint("Received message string:\n" + msgstr, verbose=2)
        msgtype, mstr = eval(msgstr)
        if msgtype=="debug_show_phil":
          self.mprint(self.show_current_phil() )
        if msgtype=="datatypedict":
          self.viewer.datatypedict = eval(mstr)
        if msgtype=="clipper_crystdict":
          self.clipper_crystdict = eval(mstr)
          self.convert_clipperdict_to_millerarrays(self.clipper_crystdict)
        if msgtype=="philstr" or msgtype=="preset_philstr":
          self.mprint("Received PHIL string:\n" + mstr, verbose=1)
          new_phil = libtbx.phil.parse(mstr)
          self.update_settings(new_phil, msgtype, lastmsgtype)
        if msgtype=="external_cmd":
          self.external_cmd = mstr
          self.mprint("Received python command string:\n" + mstr, verbose=1)
          self.run_external_cmd()

        lastmsgtype = msgtype
        time.sleep(self.zmqsleeptime)
      except Exception as e:
        self.mprint( str(e) + traceback.format_exc(limit=10), verbose=1)
    self.mprint("Shutting down zmq_listen() thread", verbose=1)
    self.guiSocketPort=None


  def ResetPhilandViewer(self, extraphil=None):
    self.ResetPhil(extraphil)
    self.viewer.symops = []
    self.viewer.sg = None
    self.viewer.proc_arrays = []
    self.viewer.HKLscenedict = {}
    self.uservectors = []
    self.viewer.visual_symmxs = []
    self.viewer.visual_symHKLs = []
    self.viewer.sceneisdirty = True
    self.viewer.isnewfile = True
    self.validated_preset_buttons = False
    if self.viewer.miller_array:
      self.viewer.params.viewer.scene_id = None
      self.viewer.RemoveStageObjects()
    self.viewer.miller_array = None
    self.viewer.lastviewmtrx = None
    return self.viewer.params


  def ResetPhil(self, extraphil=None):
    self.currentphil = master_phil
    if extraphil:
      self.currentphil = self.currentphil.fetch(source = extraphil)
      # Don't retain clip plane values as these are specific to each crystal
      # so use clip plane parameters from the master phil
      default_clipphil = master_phil.fetch().extract().clip_plane
      currentparms = self.currentphil.extract()
      currentparms.clip_plane = default_clipphil
      self.currentphil = master_phil.format(python_object = currentparms)
    self.params = self.currentphil.fetch().extract()
    self.viewer.params = self.params
    self.params.binning.binner_idx = 0
    self.params.binning.nbins = 1
    self.params.binning.scene_bin_thresholds = []
    self.params.using_space_subgroup = False


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
    # bin_opacity, show_vector and user_vector have the "multiple" attribute. In order to retain any existing
    # list elements we copy them from the oldcurrentphil parameters and assign them to the newphil parameters
    # Failure to do this would lead to list of bin_opacity elements not being properly updated when
    # changing one of the bin_opacity elements. Likewise for show_vector and user_vector.
    if jsview_3d.has_phil_path(pyphilobj, "bin_opacity"):
      bin_opacitylst = master_phil.fetch(source = pyphilobj ).extract().binning.bin_opacity
    else:
      bin_opacitylst = master_phil.fetch(source = oldcurrentphil ).extract().binning.bin_opacity
    newpyobj = newcurrentphil.extract()
    newpyobj.binning.bin_opacity = bin_opacitylst
    newcurrentphil = newcurrentphil.format(python_object= newpyobj )

    if jsview_3d.has_phil_path(pyphilobj, "show_vector"):
      show_vectorlst = master_phil.fetch(source = pyphilobj ).extract().viewer.show_vector
    else:
      show_vectorlst = master_phil.fetch(source = oldcurrentphil ).extract().viewer.show_vector
    newpyobj = newcurrentphil.extract()
    newpyobj.viewer.show_vector = show_vectorlst
    newcurrentphil = newcurrentphil.format(python_object= newpyobj )

    if jsview_3d.has_phil_path(pyphilobj, "user_vector"):
      user_vectorlst = master_phil.fetch(source = pyphilobj ).extract().viewer.user_vector \
       + oldcurrentphil.extract().viewer.user_vector
      newpyobj = newcurrentphil.extract()
      newpyobj.viewer.user_vector = user_vectorlst
      newcurrentphil = newcurrentphil.format(python_object= newpyobj )

    return newcurrentphil, diffphil


  def SetCurrentPhilAsPython(self, pyphil):
    newphil = master_phil.format(python_object= pyphil)
    currphil = master_phil.fetch(source = newphil)


  def show_current_phil(self, useful_for_preset_button=True):
    # When useful_for_preset_button=True the output has integer ids for vectors and scene_id
    # replaced with label strings for vectors and data arrays
    diffphil = master_phil.fetch_diff(source = self.currentphil)
    if useful_for_preset_button:
      # Tidy up diffphil by eliminating parameters that are not needed
      # for preset buttons or are being overwritten by user settings (selected_info, colour and radii scheme)
      # First make a copy of all phil parameters some of which may be altered below
      paramscopy = master_phil.format(self.params).copy().extract()

      omitparms = ["viewer.scene_id", "hkls.nth_power_scale_radii", "hkls.scale",
        "hkls.color_scheme", "hkls.color_powscale", "binning.binner_idx"]
        # selected_info and NGL scope are omitted below
      if jsview_3d.has_phil_path(diffphil, "scene_id"):
        # then merge corresponding label and datatype into the diffphil so these can be used for
        # preset button phil strings instead of scene_id which may only apply to the current data file
        label,datatype = self.viewer.get_label_type_from_scene_id(paramscopy.viewer.scene_id)
        paramscopy.viewer.data_array.label = label
        paramscopy.viewer.data_array.datatype = datatype
        paramsobj = master_phil.format(paramscopy)
        diffphil = master_phil.fetch_diff(source =paramsobj)

      if jsview_3d.has_phil_path(diffphil, "angle_around_vector"):
        # convert any integer id into corresponding label for this vector
        # which is needed for userfriendliness when using current phil for preset buttons
        lblvectors = self.viewer.get_vectors_labels_from_ids([paramscopy.viewer.angle_around_vector])
        strvec, angle = lblvectors[0]
        # if angle=0 this amounts to the default of not having rotated at all so omit angle_around_vector altogether
        if angle != 0.0:
          paramscopy.viewer.angle_around_vector = str(lblvectors[0])
        else:
          omitparms = omitparms + ["viewer.angle_around_vector"]
        diffphil = diffphil.format(paramscopy) # adopt new angle_around_vector value

      if jsview_3d.has_phil_path(diffphil, "animate_rotation_around_vector"):
        # convert any integer id into corresponding label for this vector
        # which is needed for userfriendliness when using current phil for preset buttons
        lblvectors = self.viewer.get_vectors_labels_from_ids([paramscopy.viewer.animate_rotation_around_vector])
        strvec, speed = lblvectors[0]
        # if speed <= 0 this amounts to the default of not animating at all so omit animate_rotation_around_vector altogether
        if speed > 0.0:
          paramscopy.viewer.animate_rotation_around_vector = str(lblvectors[0])
        else:
          omitparms = omitparms + ["clip_plane.animate_rotation_around_vector"]
        diffphil = diffphil.format(paramscopy) # adopt new animate_rotation_around_vector value

      if jsview_3d.has_phil_path(diffphil, "show_vector"):
        # convert any integer id into corresponding label for this vector
        # which is needed for userfriendliness when using current phil for preset buttons
        show_vectors_lbl = self.viewer.get_vectors_labels_from_ids(paramscopy.viewer.show_vector)
        showveclst = []
        # viewer.show_vector has multiple=True so use a list to add more vectors to the phil parameter
        for veclbl,isvisible in show_vectors_lbl:
          if isvisible==True: # default is for vectors not to show, so no need to include those
            showveclst.append(str([veclbl,isvisible]) )
        paramscopy.viewer.show_vector = showveclst
        diffphil = diffphil.format(paramscopy) # adopt new show_vector value

      if jsview_3d.has_phil_path(diffphil, "binner_idx"):
        # then merge corresponding binlabel into the diffphil so it can be used for
        # preset button phil strings instead of binner_idx which may only apply to the current data file
        binlabel = self.viewer.get_binlabel_from_binner_idx(paramscopy.binning.binner_idx)
        binlabelphil = libtbx.phil.parse("binning.binlabel = '%s'" %binlabel)
        workingphil = master_phil.fetch(sources=[binlabelphil, diffphil] )
        diffphil = master_phil.fetch_diff(source=workingphil )

      remove_bin_opacities = True
      if jsview_3d.has_phil_path(diffphil, "bin_opacity"):
        bin_opacitieslst = paramscopy.binning.bin_opacity
        for alpha,bin in bin_opacitieslst:
          if alpha < 1.0: # at least 1 bin is not fully opaque.
                          # That's not the default so don't omit opacities
            remove_bin_opacities = False
            break
      if remove_bin_opacities:
        omitparms = omitparms + ["binning.bin_opacity"]
      remainingobjs = []
      for obj in diffphil.objects:
        vobjs = []
        if hasattr(obj, "objects"):
          for vobj in obj.objects:
            if vobj.full_path() not in omitparms:
              vobjs.append(vobj)
          obj.objects = vobjs
          # The miller table column layout irrelevant for preset buttons so skip phil values governing it.
          # Applies as well to NGL phil scope
          if len(obj.objects) > 0 and obj.full_path() !=  "selected_info" and obj.full_path() !=  "NGL":
            remainingobjs.append(obj)
        else:
          remainingobjs.append(obj)
      diffphil.objects = remainingobjs
    return "\nCurrent non-default phil parameters:\n\n" + diffphil.as_str()


  def update_settings(self, new_phil=None, msgtype="philstr", lastmsgtype="philstr"):
    try:
      oldsceneid = self.params.viewer.scene_id
      currentNGLscope = None
      currentSelectInfoscope = None
      if msgtype=="preset_philstr":
        currentNGLscope = self.currentphil.extract().NGL
        currentSelectInfoscope = self.currentphil.extract().selected_info
        self.ResetPhil()
        self.viewer.sceneisdirty = True
        self.viewer.executing_preset_btn = True
      # selecting a new scene_id resets phil parameters if the previous phil was from a preset button
      if lastmsgtype=="preset_philstr" and jsview_3d.has_phil_path(new_phil, "scene_id"):
        currentNGLscope = self.currentphil.extract().NGL
        currentSelectInfoscope = self.currentphil.extract().selected_info

      if not new_phil:
        new_phil = master_phil.format(python_object = self.params)
      self.currentphil, diff_phil = self.GetNewCurrentPhilFromPython(new_phil, self.currentphil)

      self.params = self.currentphil.extract()
      phl = self.params
      if msgtype=="preset_philstr":  # override default NGL and selected_info scopes with user settings
        self.params.NGL = currentNGLscope
        self.params.selected_info = currentSelectInfoscope
      self.viewer.params = phl
      # once a preset phil setting has been enabled allow changing a phil parameter
      # without having to change scene_id
      if (msgtype=="philstr") and (lastmsgtype=="preset_philstr") and (oldsceneid is not None) and \
         jsview_3d.has_phil_path(diff_phil, "scene_id") == False:
        self.params.viewer.scene_id = oldsceneid
        self.viewer.sceneisdirty = True

      if len(diff_phil.all_definitions()) < 1 and not self.viewer.mouse_moved:
        self.mprint( "No change in PHIL parameters\n", verbose=1)
        return

      self.mprint("diff phil:\n" + diff_phil.as_str(), verbose=1 )

      if jsview_3d.has_phil_path(diff_phil, "miller_array_operation"):
        phl.viewer.scene_id = self.make_new_miller_array( msgtype=="preset_philstr" )
        self.set_scene(phl.viewer.scene_id)
        phl.hkls.sigma_color_radius = False

      # preset phil usually comes with data_array.label, data_array.phasertng_tag or data_array.datatype.
      # Scene_id is then inferred from data_array and used throughout
      if jsview_3d.has_phil_path(diff_phil, "data_array"):
        if jsview_3d.has_phil_path(diff_phil, "phasertng_tag"):
          phl.viewer.data_array.label = self.get_label_from_phasertng_tag(phl.viewer.data_array.phasertng_tag)
        phl.viewer.scene_id = self.viewer.get_scene_id_from_label_or_type(phl.viewer.data_array.label,
                                                                          phl.viewer.data_array.datatype)
      if jsview_3d.has_phil_path(diff_phil, "binlabel"):
        phl.binning.binner_idx = self.viewer.get_binner_idx_from_label(phl.binning.binlabel)

      if jsview_3d.has_phil_path(diff_phil, "scene_id"):
        phl.viewer.data_array.label = None
        phl.viewer.data_array.datatype = None

      if jsview_3d.has_phil_path(diff_phil, "use_provided_miller_arrays"):
        if not self.load_miller_arrays():
          return
        self.viewer.lastscene_id = phl.viewer.scene_id
        phl.use_provided_miller_arrays = False # ensure we can do this again

      if jsview_3d.has_phil_path(diff_phil, "openfilename"):
        fname = phl.openfilename
        currentNGLscope = self.currentphil.extract().NGL
        currentSelectInfoscope = self.currentphil.extract().selected_info
        phl = self.ResetPhilandViewer()
        phl.openfilename = fname # as openfilename was reset above
        if not self.load_reflections_file(fname):
          return
        self.params.NGL = currentNGLscope # override default NGL and selected_info scopes with user settings
        self.params.selected_info = currentSelectInfoscope
        self.viewer.lastscene_id = phl.viewer.scene_id
        self.validated_preset_buttons = False

      if jsview_3d.has_phil_path(diff_phil, "scene_id",
                                            "show_missing",
                                            "show_only_missing",
                                            "show_systematic_absences",
                                            "binning",
                                            "data_array"):
        if self.set_scene(phl.viewer.scene_id):
          self.update_space_group_choices()
          self.set_scene_bin_thresholds(phl.binning.scene_bin_thresholds,
                                         binner_idx=phl.binning.binner_idx,
                                         nbins=phl.binning.nbins )

      if jsview_3d.has_phil_path(diff_phil, "spacegroup_choice"):
        self.set_spacegroup_choice(phl.spacegroup_choice)

      if jsview_3d.has_phil_path(diff_phil, "tabulate_miller_array_ids"):
        self.tabulate_arrays(phl.tabulate_miller_array_ids)

      if jsview_3d.has_phil_path(diff_phil, "using_space_subgroup") and phl.using_space_subgroup==False:
        self.set_default_spacegroup()

      if jsview_3d.has_phil_path(diff_phil, "user_vector"):
        self.add_user_vector(self.params.viewer.user_vector)
        self.validated_preset_buttons = False

      make_new_info_tuples=False
      if jsview_3d.has_phil_path(diff_phil, "miller_array_operation"):
        miller_array_operations_lst = eval(self.params.miller_array_operation)
        (operation, millarroplabel, [labl1, type1], [labl2, type2]) = miller_array_operations_lst
        if millarroplabel not in [arr.info().label_string() for arr in self.procarrays]:
          make_new_info_tuples=True

      if jsview_3d.has_phil_path(diff_phil, "selected_info", "openfilename") or make_new_info_tuples:
        self.viewer.array_info_format_tpl = []
        for array in self.procarrays:
          if type(array.data()) == flex.std_string: # in case of status array from a cif file
            uniquestrings = list(set(array.data()))
            info = array.info()
            array = array.customized_copy(data=flex.int([uniquestrings.index(d) for d in array.data()]))
            array.set_info(info)
          if array.space_group() is None:
            array._unit_cell = uc
            array._space_group_info = spg.info()
            self.mprint("""No unit cell or space group info present in the %d. miller array.
    Borrowing them from the first miller array""" %i)
          wrap_labels = 25
          arrayinfo = ArrayInfo(array,wrap_labels)
          info_fmt, dummy, dummy2 = arrayinfo.get_selected_info_columns_from_phil(self.params )
          self.viewer.array_info_format_tpl.append( info_fmt )
        self.SendInfoToGUI({"array_infotpls": self.viewer.array_info_format_tpl})

        colnames_select_lst = []
        for philname,selected in list(self.params.selected_info.__dict__.items()):
          if not philname.startswith("__"):
            colnames_select_lst.append((philname, arrayinfo.caption_dict[philname], selected))
        self.SendInfoToGUI({ "colnames_select_lst": colnames_select_lst })

      if jsview_3d.has_phil_path(diff_phil, "save_image_name"):
        self.SaveImageName(phl.save_image_name)
        phl.save_image_name = None

      if jsview_3d.has_phil_path(diff_phil, "action"):
        ret = self.set_action(phl.action)
        phl.action = "is_running" # ensure the same action in succession can be executed
        if not ret:
          return

      if jsview_3d.has_phil_path(diff_phil, "savefilename"):
        self.SaveReflectionsFile(phl.savefilename)
      phl.savefilename = None # ensure the same action in succession can be executed

      if jsview_3d.has_phil_path(diff_phil, "hkls"):
        self.HKLsettings = phl.hkls

      #if jsview_3d.has_phil_path(diff_phil, "openfilename", "scene_id", "spacegroup_choice", "data_array"):
      if jsview_3d.has_phil_path(diff_phil, "openfilename"):
        self.list_vectors()
        self.validated_preset_buttons = False

      self.params = self.viewer.update_settings(diff_phil, phl)
      # parameters might have been changed. So update self.currentphil accordingly

      self.SendCurrentPhilValues()
      self.NewFileLoaded = False
      self.viewer.mouse_moved = False
      self.validate_preset_buttons()
      if (self.viewer.miller_array is None) :
        self.mprint( NOREFLDATA, verbose=1)
      self.mprint( "Ready")
    except Exception as e:
      self.mprint(to_str(e) + "\n" + traceback.format_exc())


  def SendCurrentPhilValues(self):
    self.currentphil = master_phil.format(python_object = self.params)
    philstrvalsdict = {}
    lst = []
    for e in self.currentphil.all_definitions():
      # deal with multiple definitions of a phil parameter by appending them to a list and
      # then assigning that list to the dictionary value with the key e.path. This assumes
      # that e.object is a phil parameter and not a phil scope.
      # user_vector is a multiple scope and cannot be cast into a dictionary. Instead it is
      # sent as part of self.viewer.all_vectors whenever self.list_vectors() is called
      if e.object.multiple == True:
        lst.append(e.object.extract())
        philstrvalsdict[e.path] = lst
      else:
        philstrvalsdict[e.path] = e.object.extract()
        lst = []
    mydict = { "current_phil_strings": philstrvalsdict }
    self.SendInfoToGUI(mydict)
    if self.viewer.params.viewer.scene_id is not None:
      self.SendInfoToGUI({ "used_nth_power_scale_radii": self.viewer.HKLscene_from_dict().nth_power_scale_radii })


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
      self.mprint(array.info().label_string() +  " looks like R-free flags. Mapping to zeros and ones", verbose=1)
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
    #self.mprint("Processing reflection data")
    self.procarrays = []
    if self.params.merge_data == False:
      self.hkls.expand_to_p1 = False
      self.hkls.expand_anomalous = False
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
    self.viewer.proc_arrays = self.procarrays
    self.viewer.identify_suitable_fomsarrays()
    self.viewer.set_miller_array(col, merge=array_info.merge,
       details=array_info.details_str)


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
    #self.params.spacegroup_choice = c
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
    miller_array_operations_lst = [ ( operation, label, arrid1, arrid2 ) ]
    self.params.miller_array_operations = str( miller_array_operations_lst )
    self.update_settings()


  def make_new_miller_array(self, is_preset_philstr=False):
    miller_array_operations_lst = eval(self.params.miller_array_operation)
    (operation, label, [labl1, type1], [labl2, type2]) = miller_array_operations_lst

    arrid1 = self.viewer.get_scene_id_from_label_or_type(labl1, type1)
    arrid2 = -1
    if labl2 != "":
      arrid2 = self.viewer.get_scene_id_from_label_or_type(labl2, type2)

    for arr in self.procarrays:
      if label in arr.info().labels + [ "", None]:
        if is_preset_philstr: # miller_array created by a preset  button. Just return the scene_id
          return self.viewer.get_scene_id_from_label_or_type(label)
        raise Sorry("Provide a label for the new miller array that isn't already used.")
    from copy import deepcopy
    millarr1 = deepcopy(self.procarrays[arrid1])
    newarray = None
    try:
      if arrid2 != -1:
        millarr2 = deepcopy(self.procarrays[arrid2])
        self.mprint("Creating %s data with array1 as %s and array2 as %s through the operation:\n\n%s" \
                     %(label, millarr1.info().label_string(), millarr2.info().label_string(), operation))
        newarray = self.viewer.OperateOn2MillerArrays(millarr1, millarr2, operation)
      else:
        self.mprint("Creating %s data with array1 as %s through the operation:\n\n%s" \
                     %(label, millarr1.info().label_string(), operation))
        newarray = self.viewer.OperateOn1MillerArray(millarr1, operation)
    except Exception as e:
      self.mprint( str(e) + traceback.format_exc(limit=10), verbose=0)

    if newarray is None:
       # allow user to quickly amend his broken python code without having to enter a new column label
      self.params.miller_array_operation = "" # do this by resetting phil parameter to the master default value
      # and update_settings() won't bail out with a "No change in PHIL parameters" message
    else:
      self.mprint("New dataset has %d reflections." %newarray.size())
      newarray.set_info(millarr1._info )
      newarray._info.labels = [ label ]
      if isinstance( newarray.sigmas(), flex.double):
        newarray._info.labels = [ label, "Sig" +label ]
      procarray, procarray_info = self.process_miller_array(newarray)
      self.procarrays.append(procarray)
      self.viewer.proc_arrays = self.procarrays
      self.viewer.has_new_miller_array = True

      wrap_labels = 25
      arrayinfo = ArrayInfo(procarray,wrap_labels)
      info_fmt, headerstr, infostr = arrayinfo.get_selected_info_columns_from_phil(self.params )
      self.viewer.array_info_format_tpl.append( info_fmt )
      # isanomalous and spacegroup might not have been selected for displaying so send them separatately to GUI
      self.ano_spg_tpls.append((arrayinfo.isanomalous, arrayinfo.spginf) )
      # Storing this new miller_array in the origarrays dictionary allows making a table of the data later.
      # First create a superset of HKLs existing miller arrays and the new procarray.
      hkls = self.origarrays["HKLs"]
      m = miller.match_indices(procarray.indices(), hkls )
      # get subset of indices in hkls matching procarray.indices()
      indices_of_matched_hkls = m.pairs().column(1)
      # pad hkls with the indices only present in procarray.indices()
      hkls.extend(  procarray.indices().select( m.singles(0)) )
      # hkls is now a superset of indices.
      # Make temporary data array the size of hkls. This will be filled with datavalues
      # from procarray matching the order of indices in hkls
      datarr = flex.double(len(hkls), float("nan"))
      # assign data values corresponding to matching indices to datarr
      m = miller.match_indices(procarray.indices(), hkls )
      # get single indices in hkls matching procarray.indices()
      indices_of_matched_hkls = m.pairs().column(1)
      for i,e in enumerate(indices_of_matched_hkls):
        datarr[e] = procarray.data()[i]
      # datarr is now a copy of data values in procarray but ordered to match the indices in hkls
      # join datarr to dictionary so it can be tabulated together with other data sets
      self.origarrays[newarray._info.labels[0]] = list(datarr)
      # If we have Sigmas then also store values and label for these in origarrays
      if isinstance( newarray.sigmas(), flex.double):
        sigarr = flex.double(len(hkls), float("nan"))
        for i,e in enumerate(indices_of_matched_hkls):
          sigarr[e] = procarray.sigmas()[i]
        self.origarrays[newarray._info.labels[1]] = list(sigarr)

      self.arrayinfos.append(arrayinfo)
      self.viewer.get_labels_of_data_for_binning(self.arrayinfos)
      mydict = { "array_infotpls": self.viewer.array_info_format_tpl,
                "ano_spg_tpls": self.ano_spg_tpls,
                "NewHKLscenes" : True,
                "NewMillerArray" : True
                }
      self.SendInfoToGUI(mydict)
      self.validated_preset_buttons = False
      self.viewer.include_tooltip_lst = [True] * len(self.viewer.proc_arrays)
      self.SendInfoToGUI({ "include_tooltip_lst": self.viewer.include_tooltip_lst })
    return len(self.viewer.hkl_scenes_infos)-1 # return scene_id of this new miller_array


  def run_external_cmd(self):
    # Run some python script like xtricorder with the exec function. Script can manipulate HKLViewFrame
    # by accessing functions and attributes on 'self' that is exported as a local variable.
    # Get logfile name and tabname assigned
    # by the script and send these to the HKLviewer GUI. Also expecting retval and errormsg to be defined
    # in the script
    try:
      ldic= {'retval': None, 'errormsg': None, 'self': self, 'master_phil': master_phil }
      exec(self.external_cmd, globals(), ldic)
      retval = ldic.get("retval", None)
      errormsg = ldic.get("errormsg", None)
      if retval != 0:
        raise Sorry(errormsg)
      tabname = ldic.get("tabname", None)
      logfname = ldic.get("logfname", None)
      self.SendInfoToGUI( {"show_log_file_from_external_cmd": [tabname, logfname ]  } )
      self.validated_preset_buttons = False
      self.validate_preset_buttons()
    except Exception as e:
      self.SendInfoToGUI( {"show_log_file_from_external_cmd": -42 } )
      raise Sorry(str(e))


  def prepare_dataloading(self):
    self.viewer.isnewfile = True
    self.params.viewer.scene_id = None
    self.viewer.match_valarrays = []
    self.viewer.proc_arrays = {}
    self.spacegroup_choices = []
    self.origarrays = {}
    display.reset_settings()
    self.hkls = display.settings()
    self.viewer.mapcoef_fom_dict = {}
    self.viewer.sceneid_from_arrayid = []
    self.hklfile_history = []
    self.tncsvec = None
    self.aniso1 = None
    self.aniso2 = None
    self.aniso3 = None
    self.loaded_file_name = ""


  def finish_dataloading(self, arrays):
    valid_arrays = []
    self.viewer.array_info_format_tpl = []
    spg = arrays[0].space_group()
    uc = arrays[0].unit_cell()
    self.ano_spg_tpls =[]
    self.mprint("%d Miller arrays in this dataset:" %len(arrays))
    spgset = set([])
    self.arrayinfos = []
    previous_ucell = None
    for i,array in enumerate(arrays):
      if type(array.data()) == flex.std_string: # in case of status array from a cif file
        uniquestrings = list(set(array.data()))
        info = array.info()
        array = array.customized_copy(data=flex.int([uniquestrings.index(d) for d in array.data()]))
        array.set_info(info)
      if i>0:
        if arrays[i-1].unit_cell() is not None:
          previous_ucell = arrays[i-1].unit_cell()
        if arrays[i-1].space_group() is not None:
          previous_spg = arrays[i-1].space_group()

      # A cif file might lack unit cell or space group for all the crystals in the file
      if array.unit_cell() is None:
        if previous_ucell is None:
          raise Sorry("No unit cell found in the first miller array.")
        array._unit_cell = previous_ucell
        self.mprint("""No unit cell present in the %d. miller array. Borrowing from previous miller array""" %i)
      if array.space_group() is None:
        symm_new = crystal.symmetry( unit_cell = previous_ucell,
                                    space_group_info = previous_spg.info()
                                    )
        info = array.info()
        array = array.customized_copy(crystal_symmetry = symm_new)
        array.set_info(info)
        self.mprint("""No space group present in the %d. miller array. Borrowing from previous miller array""" %i)
      if array.space_group() is None:
        raise Sorry("No space group definition found in the first miller array.")

      wrap_labels = 25
      arrayinfo = ArrayInfo(array,wrap_labels)
      info_fmt, headerstr, infostr = arrayinfo.get_selected_info_columns_from_phil(self.params )
      if i==0: # print formatted table of array info here
        self.mprint(headerstr)
      self.mprint(infostr)
      self.viewer.array_info_format_tpl.append( info_fmt )
      # isanomalous and spacegroup might not have been selected for displaying so send them separatately to GUI
      self.ano_spg_tpls.append((arrayinfo.isanomalous, arrayinfo.spginf) )
      spgset.add(arrayinfo.ucellinf)
      if i==0:
        # convert philstring of selected_info into a list so GUI can make a selection settings dialog
        # for what columns to show in the millertable
        colnames_select_lst = []
        for philname,selected in list(self.params.selected_info.__dict__.items()):
          if not philname.startswith("__"):
            colnames_select_lst.append((philname, arrayinfo.caption_dict[philname], selected))
        self.SendInfoToGUI({ "colnames_select_lst": colnames_select_lst })
      valid_arrays.append(array)
      self.arrayinfos.append(arrayinfo)
    self.valid_arrays = valid_arrays
    self.SendInfoToGUI({"spacegroup_info": arrayinfo.spginf, "unitcell_info": list(spgset) })

    if self.fileinfo:
      return
    if (len(valid_arrays) == 0):
      msg = "No arrays of the supported types present."
      self.mprint(msg)
      self.NewFileLoaded=False
    elif (len(valid_arrays) >= 1):
      self.set_miller_array()
      self.viewer.get_labels_of_data_for_binning(self.arrayinfos)
      self.update_space_group_choices(0) # get the default spacegroup choice
      mydict = { "info": self.infostr,
                  "bin_infotpls": self.viewer.bin_infotpls,
                  "ano_spg_tpls": self.ano_spg_tpls,
                  "html_url": self.viewer.url,
                  "tncsvec": self.tncsvec,
                  "merge_data": self.params.merge_data,
                  "spacegroups": [e.symbol_and_number() for e in self.spacegroup_choices],
                  "NewFileLoaded": self.NewFileLoaded,
                  "file_name": self.params.openfilename
                }
      self.SendInfoToGUI(mydict)
    self.viewer.include_tooltip_lst = [True] * len(self.viewer.proc_arrays)
    self.SendInfoToGUI({ "include_tooltip_lst": self.viewer.include_tooltip_lst })


  def load_reflections_file(self, file_name):
    file_name = to_str(file_name)
    ret = False
    self.NewFileLoaded=True
    if (file_name != ""):
      try :
        self.mprint("\nReading file %s..." %file_name)
        self.prepare_dataloading()
        hkl_file = any_reflection_file(file_name)
        arrays = hkl_file.as_miller_arrays(merge_equivalents=False, reconstruct_amplitudes=False)
        self.origarrays = {}
        if hkl_file._file_type == 'cif':
          # use new cif label parser for reflections
          cifreader = hkl_file.file_content()
          cifarrays = cifreader.as_miller_arrays(merge_equivalents=False)
          arrays = [] # overwrite with simplified label strings
          for arr in cifarrays: # avoid these un-displayable arrays
            if arr.info().labels[-1] not in ['_refln.crystal_id',
                      'HKLs','_refln.wavelength_id', '_refln.scale_group_code']:
              arrays.append(arr)
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
          ciforigarrays = cifreader.as_original_arrays()[dataname[0]]
          for key in ciforigarrays:
            if key not in ['_refln.crystal_id', # avoid these un-displayable arrays
                      '_refln.wavelength_id', '_refln.scale_group_code']:
              self.origarrays[key] = ciforigarrays[key]
          # replace ? with nan in self.origarrays to allow sorting tables of data in HKLviewer
          for labl in self.origarrays.keys():
            origarray = self.origarrays[labl]
            for i,e in enumerate(self.origarrays[labl]):
              if e=="?":
                origarray[i] = "nan"
            try:
              self.origarrays[labl] = flex.double(origarray)
            except Exception as e:
              self.origarrays[labl] = origarray
        if hkl_file._file_type == 'ccp4_mtz':
          self.hklfile_history = list(hkl_file._file_content.history())
          self.loaded_file_name = file_name
          for e in self.hklfile_history:
            if "TNCS" in e and "VECTOR" in e:
              svec = e.split()[-3:]
              t1 = float(svec[0])
              t2 = float(svec[1])
              t3 = float(svec[2])
              if (t1*t1 + t2*t2 + t3*t3) > 0.0:
                self.tncsvec = (t1, t2, t3)
                self.mprint("tNCS vector found in header of mtz file: %s" %str(self.tncsvec) )
            if "PHASER A" in e and len(e.split()) == 11:
              self.aniso1 = [ eval(f) for f in e.split()[2:5] ]
              self.aniso2 = [ eval(f) for f in e.split()[5:8] ]
              self.aniso3 = [ eval(f) for f in e.split()[8:11] ]
              self.mprint("Anisotropic principal axes found in header of mtz file: %s, %s, %s" \
                %(str(self.aniso1),str(self.aniso2),str(self.aniso3) ))
          from iotbx import mtz
          mtzobj = mtz.object(file_name)
          nanval = float("nan")
          # deep copy to avoid out of range errors elsewhere if user merges reflections and
          # we need to extend list of reflections
          self.origarrays["HKLs"] = mtzobj.extract_miller_indices()[:]
          for mtzlbl in mtzobj.column_labels():
            col = mtzobj.get_column( mtzlbl )
            newarr = col.extract_values_and_selection_valid().values.deep_copy()
            for i,b in enumerate(col.extract_values_and_selection_valid().selection_valid):
              if not b:
                newarr[i] = nanval
            self.origarrays[mtzlbl] = list(newarr)

        if len(self.origarrays.items()) == 0:
          self.origarrays["HKLs"] = arrays[0].indices()[:]
          for arr in arrays:
            if (arr.is_complex_array() or arr.is_hendrickson_lattman_array())==False:
              if arr.sigmas() == None:
                self.origarrays[arr.info().label_string()] = arr.data()
              else:
                self.origarrays[arr.info().labels[0]] = arr.data()
                self.origarrays[arr.info().labels[1]] = arr.sigmas()
        self.finish_dataloading(arrays)
        self.SendInfoToGUI({"NewFileLoaded": self.NewFileLoaded})
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
      if self.firsttime:
        self.ResetPhilandViewer(self.currentphil)
        self.prepare_dataloading()
        self.firsttime = False
      self.finish_dataloading(self.provided_miller_arrays)
      self.viewer.sceneisdirty = True
      self.viewer.has_new_miller_array = True
      ret = True
    except Exception as e :
      self.NewFileLoaded=False
      self.mprint("".join(traceback.format_tb(e.__traceback__ )) + e.__repr__())
      arrays = []
    return ret


  def LoadMillerArrays(self, marrays):
    self.provided_miller_arrays = marrays
    self.update_settings()
    self.params.use_provided_miller_arrays = True


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


  def get_label_from_phasertng_tag(self, tngcolumn_tags):
    tngcols = []
    # Say tngcolumn_tags = "INAT,SIGINAT" and mtz history looks like:
    # PHASER LABIN INAT/I<<FW/J SIGINAT/SIGI<<FW/Q IPOS/I(+)<<I/K
    # PHASER LABIN SIGIPOS/SIGI(+)<<SIGI/M INEG/I(-)<<I/K SIGINEG/SIGI(-)<<SIGI/M
    # then m would look like [('INAT', 'SIGINAT')]. Concatenate the strings before returning
    label = ""
    for e in self.hklfile_history:
      for tngcolumn_tag in tngcolumn_tags.split(","):
        m =  re.findall(tngcolumn_tag + '/(\S*)/', e, re.VERBOSE)
        if len(m) > 0:
          tngcols.append(m[0])
    if len(tngcols):
      label = ",".join(tngcols)
    return label


  def validate_preset_buttons(self):
    if not self.validated_preset_buttons:
      activebtns = []
      # look for strings like data_array.label="F,SIGFP" and see if that data column exists in the file
      uniquebtnids = set([])
      if self.viewer.miller_array is not None:
        ma = self.viewer.miller_array
      elif len(self.procarrays):
        ma = self.procarrays[0]
      else:
        ma = -1

      self.mprint("Preset buttons:", verbose=1)
      for ibtn,(btn_id, btnlabel, philstr) in enumerate(self.allbuttonslist):
        button_fate_decided = False
        if btn_id not in uniquebtnids:
          uniquebtnids.add(btn_id)
        else:
          raise Sorry("Button ID, %s, has already been used by another button." %btn_id)

        btnphil = libtbx.phil.parse(philstr)
        philstr_label = None
        philstr_type = None
        phasertng_tag = None
        philstr_showvectors = []
        philstr_user_vectors_labels = []
        millaroperationstr = None
        if jsview_3d.has_phil_path(btnphil, "data_array", "show_vector", "miller_array_operation"):
          btnphilobj, unusedparms = master_phil.fetch(btnphil, track_unused_definitions=True)
          for parm in unusedparms:
            self.mprint( "Preset button, %s, has unrecognised phil parameter:\n   %s" %(btn_id, parm.path), verbose=1)
          btnphilextract = btnphilobj.extract()
          if btnphilextract.viewer.data_array.label is not None:
            philstr_label = btnphilextract.viewer.data_array.label
          if btnphilextract.viewer.data_array.datatype is not None:
            philstr_type = btnphilextract.viewer.data_array.datatype
          if btnphilextract.viewer.data_array.phasertng_tag is not None:
            phasertng_tag = btnphilextract.viewer.data_array.phasertng_tag
            # find the miller array used by phasertng as specified in the mtz history header
            philstr_label = [ self.get_label_from_phasertng_tag(",".join(phasertng_tag)) ]
          if len(btnphilextract.viewer.show_vector) > 0:
            philstr_showvectors = btnphilextract.viewer.show_vector
          if len(btnphilextract.viewer.user_vector) > 0:
            if ma==-1:
              continue
            for uvec in btnphilextract.viewer.user_vector:
              if uvec.hkl_op != "": # then verify this operation is commensurate with the spacegroup
                rt = sgtbx.rt_mx(symbol=uvec.hkl_op, r_den=12, t_den=144)
                (cartvec, a, rotlabel, order) = self.viewer.GetVectorAndAngleFromRotationMx( rt.r(), ma)
                if rotlabel =="" or order==0:
                  self.mprint("\"%s\" is disabled because HKL operation, \"%s\", is not a rotation in space group %s" \
                   %(btnlabel, uvec.hkl_op, ma.space_group().info().symbol_and_number()), verbose=1)
                  activebtns.append((self.allbuttonslist[ibtn],False, "", None))
                  button_fate_decided = True
                  break
              else:
                philstr_user_vectors_labels.append( uvec.label)
          if button_fate_decided:
            continue

          if btnphilextract.miller_array_operation != "":
# The miller array operation part in the philstr could look like:
# miller_array_operation = "('newarray._data = array1.data()/array1.sigmas()\\nnewarray._sigmas = None', 'IoverSigI', ['I<<FSQ,SIGI<<FSQ', 'Intensity'], ['', ''])"
# We want to capture 'I<<FSQ,SIGI<<FSQ' and 'Intensity' strings which will be in arr1label and arr1type
            millaroperationstr, millarrlabel, (arr1label, arr1type), (arr2label, arr2type) = \
                                                     eval( btnphilextract.miller_array_operation)
        nvectorsfound = len(philstr_showvectors)
        veclabels = []
        philveclabel = ""
        if philstr_showvectors:
          nvectorsfound = 0
          for iphilvec,philstrvec in enumerate(philstr_showvectors):
            philveclabel, philshowvec = eval(philstrvec)
            # see if any of the user_vectors_labels is a substring of the label for the vectors to display
            #if True in [ lbl in philveclabel for lbl in philstr_user_vectors_labels ] :
            #  nvectorsfound = len(philstr_showvectors)
            #  continue # button phil defines a user vector matching the show vector
            for opnr, veclabel, order, cartvec, hklop, hkl, abc, length in self.viewer.all_vectors:
              # allow label to be just a substring of veclabel
              philstr_userlbl = ""
              for lbl in philstr_user_vectors_labels:
                if lbl == veclabel: # button phil defines a user vector matching the show vector
                  philstr_userlbl = lbl
                  break
              if philshowvec and philveclabel in veclabel:
                nvectorsfound +=1
                if len(philveclabel) < len(veclabel):
                  if philstr_userlbl:
                    #veclabels += "," + philstr_userlbl
                    veclabels.append(philstr_userlbl)
                  else:
                    #veclabels += "," + veclabel
                    veclabels.append(veclabel)
            if (iphilvec+1) > nvectorsfound:
              self.mprint("\"%s\" is disabled until a vector, \"%s\", has been " \
                   "found in a dataset or by manually adding this vector." %(btnlabel, philveclabel), verbose=1)
        miller_array_operation_can_be_done = False
        if millaroperationstr:
          for inflst, pidx, fidx, datalabel, datatype, hassigmas, sceneid in self.viewer.hkl_scenes_infos:
            if datalabel == arr1label or datatype == arr1type:
              miller_array_operation_can_be_done = True
              break
          if miller_array_operation_can_be_done:
            self.mprint("\"%s\" declared using %s and %s is assigned to data %s of type %s." \
                          %(btnlabel, arr1label, arr1type, datalabel, datatype), verbose=1)
            activebtns.append((self.allbuttonslist[ibtn],True, datalabel, None))
          else:
            self.mprint("\"%s\" declared using %s and %s is not assigned to any dataset." \
                            %(btnlabel, arr1label, arr1type), verbose=1)
            #activebtns.append((self.allbuttonslist[ibtn], False, "", None))
        if philstr_label is not None and millaroperationstr is None:
          labeltypefound = False
          for inflst, pidx, fidx, datalabel, datatype, hassigmas, sceneid in self.viewer.hkl_scenes_infos:
            if datalabel == philstr_label:
              labeltypefound = True
              break
          if not labeltypefound:
            for inflst, pidx, fidx, datalabel, datatype, hassigmas, sceneid in self.viewer.hkl_scenes_infos:
              if philstr_type is not None and philstr_type == datatype:
                labeltypefound = True
                break
          if labeltypefound and nvectorsfound >= len(philstr_showvectors):
            self.mprint("\"%s\" assigned to dataset %s of type %s." \
                          %(btnlabel + str(veclabels), datalabel, datatype), verbose=1)
            activebtns.append((self.allbuttonslist[ibtn], True, datalabel, (philveclabel, veclabels) ))
          else:
            self.mprint("\"%s\" expecting dataset of type \"%s\" has not been assigned to any dataset." \
                              %(btnlabel, philstr_type), verbose=1)
            #activebtns.append((self.allbuttonslist[ibtn], False, "", None))

      self.SendInfoToGUI({"enable_disable_preset_buttons": str(activebtns)})
    self.validated_preset_buttons = True


  def convert_clipperdict_to_millerarrays(self, crystdict):
    """
    Called in zmq_listen() when Chimerax with Isolde sends clipper arrays to
    our zmqsocket from HKLviewer.ProcessMessages()
    """
    xs = crystal.symmetry(unit_cell=crystdict["unit_cell"], space_group_symbol= crystdict["spg_number"] )
    mi = flex.miller_index(crystdict["HKL"])
    lst = list(crystdict.keys())
    lst.remove('HKL')
    lst.remove('unit_cell')
    lst.remove('spg_number')
    clipperlabel = lst[0]
    Flabl, Siglabl = clipperlabel.split(", ")
    data = flex.double(crystdict[clipperlabel][0])
    sigmas = flex.double(crystdict[clipperlabel][1])
    marray = miller.array( miller.set(xs, mi, anomalous_flag=False),
                         data, sigmas).set_observation_type( observation_types.amplitude())
    marray.set_info(miller.array_info(source="Isolde", labels=[Flabl, Siglabl]))

    fcamplitudes = flex.double(crystdict["FCALC,PHFCALC"][0])
    fcphases = flex.double(crystdict["FCALC,PHFCALC"][1])
    mapcoeffarray = miller.array( miller.set(xs, mi, anomalous_flag=False), fcamplitudes)
    mapcoeffarray = mapcoeffarray.phase_transfer(fcphases, deg=False)
    mapcoeffarray.set_info(miller.array_info(source="Isolde", labels=["FCALC", "PHFCALC"]))

    wamplitudes = flex.double(crystdict["2FOFC,PH2FOFC"][0])
    wphases = flex.double(crystdict["2FOFC,PH2FOFC"][1])
    wmapcoeffarray = miller.array( miller.set(xs, mi, anomalous_flag=False), wamplitudes)
    wmapcoeffarray = wmapcoeffarray.phase_transfer(wphases, deg=False)
    wmapcoeffarray.set_info(miller.array_info(source="Isolde", labels=["2FOFC", "PH2FOFC"]))

    self.LoadMillerArrays([marray, mapcoeffarray, wmapcoeffarray])


  def has_indices_with_multiple_data(self, arr):
    return len(set(list(arr.indices()))) < arr.size()


  def tabulate_arrays(self, datalabels):
    if len(self.origarrays) == 0: # if not an mtz file then split columns
      self.origarrays["HKLs"] = self.viewer.proc_arrays[0].indices()[:]
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
          if self.viewer.array_info_format_tpl[id][0] == 'FreeR_flag': # want True or False back
            list_with_nans = [ 1==e if not cmath.isnan(e) else display.nanval for e in list_with_nans ]
          self.origarrays[arr.info().label_string()] = list_with_nans
        else:
          self.origarrays[arr.info().label_string()] = list(arr.data())

    indices = self.origarrays["HKLs"]
    dres = self.procarrays[0].unit_cell().d( indices)
    dreslst = [("d_res", roundoff(list(dres)),3)]
    hkls = list(indices)
    hkllst = [ ("H", [e[0] for e in hkls] ), ("K", [e[1] for e in hkls] ), ("L", [e[2] for e in hkls] )]
    datalst = []
    labellists = eval(datalabels)
    for labels in labellists:
      crystlbl = ""; wavelbl = ""; scalelbl =""
      for i,label in enumerate(labels):
        if "crystal_id" in label:
          crystlbl = "," + label
        if "wavelength_id" in label:
          wavelbl = "," + label
        if "scale_group_code" in label:
          scalelbl = "," + label
      for label in labels:
        if "crystal_id" in label or "wavelength_id" in label or "scale_group_code" in label:
          continue
        fulllabel = label + crystlbl + wavelbl + scalelbl
        pydatlst = list(self.origarrays[fulllabel])
        # If a merged array has been created by the user pydatlst could be shorter than len(hkls).
        # pydatlst must have the same size as hkls when received by helpers.MillerArrayTableModel()
        # If it is not, then pad nan values at the end to make up for it.
        # Otherwise the tabulated reflection data won't display correctly
        if len(pydatlst) < len(hkls):
          pydatlst.extend( [float("nan")]*(len(hkls)-len(pydatlst) ))
        datalst.append( (label, pydatlst))
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


  def ExpandToP1(self, val):
    self.params.hkls.expand_to_p1 = val
    self.update_settings()


  def ExpandAnomalous(self, va):
    self.params.hkls.expand_anomalous = val
    self.update_settings()


  def ShowOnlyMissing(self, val):
    self.params.hkls.show_only_missing = val
    self.update_settings()


  def ShowMissing(self, val):
    self.params.hkls.show_missing = val
    self.update_settings()


  def ShowDataOverSigma(self, val):
    self.params.hkls.show_data_over_sigma = val
    self.update_settings()


  def ShowSystematicAbsences(self, val):
    self.params.hkls.show_systematic_absences = val
    self.update_settings()


  def ShowSlice(self, val, axis="h", index=0):
    axisstr = axis.lower()
    self.params.hkls.slice_mode = val
    self.params.hkls.slice_axis = axisstr
    self.params.hkls.slice_index = index
    self.update_settings()


  def set_scene_bin_thresholds(self, thresholds = None, binner_idx = 0,  nbins = 6):
    nuniquevalues = -1
    if not thresholds:
      binvals, nuniquevalues = self.viewer.calc_bin_thresholds(binner_idx, nbins)
    else:
      binvals = thresholds[:]
    if binvals and binner_idx == 0: # binner_idx=0 is for binning against resolution
      binvals = list( 1.0/flex.double(binvals) )
    self.viewer.UpdateBinValues(binner_idx, binvals, nuniquevalues)


  def SetSceneNbins(self, nbins, binner_idx = 0):
    self.params.nbins = nbins
    self.params.binning.binner_idx = binner_idx
    self.params.binning.bin_opacity = [ [1.0, e] for e in range(nbins) ]
    self.update_settings()


  def GetNumberingOfBinners(self):
    return [ (i,e) for i,e in enumerate(self.viewer.bin_labels_type_idxs) ]


  def SetSceneBinThresholds(self, thresholds=None):
    if thresholds==None:
      self.params.scene_bin_thresholds = []
    else:
      self.params.scene_bin_thresholds = thresholds[:]
    self.params.nbins = len(binvals)
    self.update_settings()


  def SetOpacities(self, bin_opacities):
    #self.params.bin_opacities = str(bin_opacities)
    self.params.bin_opacities = bin_opacity
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


  def SetRadiiScale(self, scale=1.0, nth_power_scale = float("nan")):
    """
    Scale radii. Decrease the contrast between large and small radii with nth_root_scale < 1.0
    If nth_power_scale=0 then all radii will have the same size regardless of data values.
    If nth_power_scale=NaN an automatic power will be computed ensuring the smallest radius
    is 0.1 times the maximum radius
    """
    self.params.hkls.scale = scale
    self.params.hkls.nth_power_scale_radii = nth_power_scale
    self.update_settings()


  def SetColourRadiusToSigmas(self, val):
    self.params.hkls.sigma_color_radius = val
    self.update_settings()


  def SetColourScheme(self, color_scheme, color_powscale=1.0):
    self.params.hkls.color_scheme = color_scheme
    self.params.hkls.color_powscale = color_powscale
    self.update_settings()


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
    self.viewer.calc_rotation_axes()
    self.viewer.all_vectors = self.viewer.rotation_operators[:]
    if self.viewer.miller_array is not None:
      uc = self.viewer.miller_array.unit_cell()
    else: # a fallback
      uc = self.procarrays[0].unit_cell()

    tncsvec = []
    if self.tncsvec is not None:
      # TNCS vector is specified in realspace fractional coordinates. Convert it to cartesian
      cartvec = list( self.tncsvec * matrix.sqr(uc.orthogonalization_matrix()) )
      ln = len(self.viewer.all_vectors)
      # Use half the length of the tncs vector to allow stepping through alternating weak and strong layers
      # of reflections in the GUI when orienting clip plane perpendicular to the tncs vector
      veclength = self.viewer.renderscale*0.5/math.sqrt( cartvec[0]*cartvec[0] + cartvec[1]*cartvec[1] + cartvec[2]*cartvec[2] )
      tncsvec = [("TNCS", 0, cartvec, "", "", str(roundoff(self.tncsvec, 5)), veclength )]

    anisovectors = []
    if self.aniso1 is not None:
      # anisotropic principal axes vector are specified in cartesian coordinates.
      cartvec = [self.aniso1[0]*self.viewer.renderscale, self.aniso1[1]*self.viewer.renderscale, self.aniso1[2]*self.viewer.renderscale]
      # convert vector to fractional coordinates for the table in the GUI
      aniso1frac = list(cartvec * matrix.sqr(uc.fractionalization_matrix()))
      veclength = self.viewer.renderscale/math.sqrt( cartvec[0]*cartvec[0] + cartvec[1]*cartvec[1] + cartvec[2]*cartvec[2] )
      anisovectors = [("Anisotropy1", 0, cartvec, "", "", str(roundoff(aniso1frac, 5)), veclength )]

      cartvec = [self.aniso2[0]*self.viewer.renderscale, self.aniso2[1]*self.viewer.renderscale, self.aniso2[2]*self.viewer.renderscale]
      # convert vector to fractional coordinates for the table in the GUI
      aniso2frac = list(cartvec * matrix.sqr(uc.fractionalization_matrix()))
      veclength = self.viewer.renderscale/math.sqrt( cartvec[0]*cartvec[0] + cartvec[1]*cartvec[1] + cartvec[2]*cartvec[2] )
      anisovectors.append(("Anisotropy2", 0, cartvec, "", "", str(roundoff(aniso2frac, 5)), veclength ) )

      cartvec = [self.aniso3[0]*self.viewer.renderscale, self.aniso3[1]*self.viewer.renderscale, self.aniso3[2]*self.viewer.renderscale]
      # convert vector to fractional coordinates for the table in the GUI
      aniso3frac = list(cartvec * matrix.sqr(uc.fractionalization_matrix()))
      veclength = self.viewer.renderscale/math.sqrt( cartvec[0]*cartvec[0] + cartvec[1]*cartvec[1] + cartvec[2]*cartvec[2] )
      anisovectors.append( ("Anisotropy3", 0, cartvec, "", "", str(roundoff(aniso3frac, 5)), veclength ) )

    ln = len(self.viewer.all_vectors)
    Hcartvec = list( self.viewer.renderscale*( (1,0,0)*matrix.sqr(uc.fractionalization_matrix()).transpose()) )
    Kcartvec = list( self.viewer.renderscale*( (0,1,0)*matrix.sqr(uc.fractionalization_matrix()).transpose()) )
    Lcartvec = list( self.viewer.renderscale*( (0,0,1)*matrix.sqr(uc.fractionalization_matrix()).transpose()) )
    Hlength = math.sqrt( Hcartvec[0]*Hcartvec[0] + Hcartvec[1]*Hcartvec[1] + Hcartvec[2]*Hcartvec[2] )
    Klength = math.sqrt( Kcartvec[0]*Kcartvec[0] + Kcartvec[1]*Kcartvec[1] + Kcartvec[2]*Kcartvec[2] )
    Llength = math.sqrt( Lcartvec[0]*Lcartvec[0] + Lcartvec[1]*Lcartvec[1] + Lcartvec[2]*Lcartvec[2] )
    hklunit_vectors = [ ("H-axis (1,0,0)", 0, Hcartvec, "", "(1,0,0)", "", Hlength ),
                        ("K-axis (0,1,0)", 0, Kcartvec, "", "(0,1,0)", "", Klength ),
                        ("L-axis (0,0,1)", 0, Lcartvec, "", "(0,0,1)", "", Llength )]

    all_vecs = hklunit_vectors + tncsvec + anisovectors + self.viewer.rotation_operators[:] + self.uservectors
    self.viewer.all_vectors = []
    for opnr,(label, order, cartvec, hkl_op, hkl, abc, length) in enumerate(all_vecs):
      self.viewer.all_vectors.append( (opnr, label, order, cartvec, hkl_op, hkl, abc, length) )

    for (opnr, label, order, cartvec, hkl_op, hkl, abc, length) in self.viewer.all_vectors:
      # avoid onMessage-DrawVector in HKLJavaScripts.js misinterpreting the commas in strings like "-x,z+y,-y"
      name = label + hkl_op.replace(",", "_")
      self.viewer.RemovePrimitives(name)
    self.SendInfoToGUI( { "all_vectors": self.viewer.all_vectors } )
    return self.viewer.all_vectors


  def add_user_vector(self, philuser_vectors, rectify_improper_rotation=False):
    uc = self.viewer.miller_array.unit_cell()
    try:
      for phil_uvec in philuser_vectors:
        label = phil_uvec.label
        userveclabels = [ e[0] for e in self.uservectors ]
        if label in userveclabels:
          continue
        order = 0
        hklvec = ""
        abcvec = ""
        hklop = ""
        unwantedchars = " |(|)|[|]|{|}"
        # individual characters separated by | substituted with a "" using re.sub()
        if phil_uvec.hkl not in [None, "", "()"]:
          hklvec = eval(re.sub(unwantedchars, "", phil_uvec.hkl))
          # convert into cartesian space
          cartvec = list( self.viewer.renderscale*(hklvec * matrix.sqr(uc.fractionalization_matrix()).transpose()) )
          veclength = math.sqrt( cartvec[0]*cartvec[0] + cartvec[1]*cartvec[1] + cartvec[2]*cartvec[2] )
        elif phil_uvec.abc not in [None, "", "()"]:
          abcvec = eval(re.sub(unwantedchars, "", phil_uvec.abc))
          # convert into cartesian space
          cartvec = list(abcvec * matrix.sqr(uc.orthogonalization_matrix()))
          # length unit used by a realspace vector in reciprocal space is the inverse of its realspace length
          veclength = self.viewer.renderscale/math.sqrt( cartvec[0]*cartvec[0] + cartvec[1]*cartvec[1] + cartvec[2]*cartvec[2] )
        elif phil_uvec.hkl_op not in [None, ""]:
          hklop = re.sub(unwantedchars, "", phil_uvec.hkl_op)
          rt = sgtbx.rt_mx(symbol=hklop, r_den=12, t_den=144)
          self.viewer.symops.append( rt ) #
          (cartvec, a, rotlabel, order) = self.viewer.GetVectorAndAngleFromRotationMx( rt.r(),
                                                    rectify_improper_rotation=rectify_improper_rotation )
          veclength = math.sqrt( cartvec[0]*cartvec[0] + cartvec[1]*cartvec[1] + cartvec[2]*cartvec[2] )
          if rotlabel:
            label = "%s-fold_%s" %(str(int(roundoff(2*math.pi/a, 0))), label)
            if label in userveclabels:
              continue # this vector label is already there. Don't add it
            self.mprint("Rotation axis, %s, added" %rotlabel)
          if rotlabel =="" or order==0:
            self.mprint("Cannot compute a rotation axis from %s" %phil_uvec.hkl_op)
            return
        vecundefined = (phil_uvec.hkl in [None, "", "()"] \
         and phil_uvec.abc in [None, "", "()"] \
         and phil_uvec.hkl_op in [None, ""])
        if not vecundefined:
          if phil_uvec.label in [None, ""]:
            raise Sorry("Specify user_vector properly!")
          self.uservectors.append( (label, order, cartvec, hklop, str(hklvec), str(abcvec), veclength ))
      self.list_vectors()
    except Exception as e:
      raise Sorry( str(e))


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


  def ShowVector(self, val, b=True):
    self.params.viewer.show_vector = [str([val, b])]
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
      self.params.clip_plane.clip_width = clipwidth
      self.params.hkls.slice_mode = False
    else:
      self.params.clip_plane.clip_width = None
    self.update_settings()


  def AnimateRotateAroundVector(self, vecnr, speed):
    self.params.viewer.animate_rotation_around_vector = str([vecnr, speed])
    self.update_settings()


  def RotateAroundVector(self, vecnr, dgr):
    self.params.viewer.angle_around_vector = str([vecnr, dgr])
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


  def GetHtmlURL(self):
    return self.viewer.url


  def GetHtmlstring(self):
    return self.viewer.htmlstr


  def GetArrayInfotpls(self):
    """
    return array of tuples with information on each miller array
    """
    return self.viewer.array_info_format_tpl


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


master_phil_str = """
  openfilename = None
    .type = path
    .help = "Name of file with one or more datasets"
  use_provided_miller_arrays = False
    .type = bool
    .help = "internal flag"
  savefilename = None
    .type = path
    .help = "Name of file where the user wants to save datasets. Optionally used after making new datasets from existing ones"
  save_image_name = None
    .type = path
    .help = "Name of image file (PNG format) where the current displayed reflections will be saved to at the users request"
  merge_data = False
    .type = bool
    .help = "internal flag"
  miller_array_operation = ''
    .type = str
    .help = "Python syntax string defining a new cctbx.miller_array object from one or two exisitng miller arrays " \
            "in the loaded data file. The CCTBX API is used for this"
  spacegroup_choice = None
    .type = int
  using_space_subgroup = False
    .type = bool
    .help = "internal flag"
  real_space_unit_cell_scale_fraction = None
    .type = float(value_min=0.0, value_max=1.0)
    .help = "Parameter specifying the scale at which to display the unit cell compared to reciprocal space. 0 means " \
            "true to scale of the size of the reciprocal lattice. 1 means close to the radius of the displayed sphere of reflections."
  reciprocal_unit_cell_scale_fraction = None
    .type = float(value_min=0.0, value_max=1.0)
    .help = "Parameter specifying the scale at which to display the reciprocal unit cell compared to reciprocal space. " \
            "0 means true to scale of the size of the reciprocal lattice. 1 means close to the radius of the displayed sphere of reflections."
  clip_plane
    .help = "Optionally imposed clip plane which is always parallel to the screen."
  {
    hkldist = 0.0
      .type = float
      .help = "Distance from origin of the center of the clip plane. " \
              "If the clip plane is normal to a reciprocal lattice vector or the associated real space " \
              "vector the unit of hkldist is the length of the reciprocal lattice vector. " \
              "If the clip plane is normal to a real space vector the unit of hkldist is the " \
              "inverse of its real space length."
    normal_vector = "-1"
      .type = str
    is_assoc_real_space_vector = False
      .type = bool
      .help = "Indicate if using associated real space vector to the selected reciprocal space vector"
    normal_vector_length_scale = -1
      .type = float
      .help = "If value is negative the unit length of hkldist is used as the scale."
    clip_width = None
      .type = float
      .help = "If value is not None then we are clipping. If auto_clip_width is True this value is ignored."
    auto_clip_width = True
      .type = bool
      .help = "If true compute appropriate clip plane width. Otherwise use clip_width value"
    fractional_vector = reciprocal *realspace
      .type = choice
  }

  %s

  binning {
    scene_bin_thresholds = []
      .type = floats
    binner_idx = 0
      .type = int
      .help = "Index in list of binners, say ['Resolution', 'Singletons', 'I,SIGI', 'Sigmas of I,SIGI',..] "
    binlabel = None
      .type = str
      .help = "Element in list of binners, say ['Resolution', 'Singletons', 'I,SIGI', 'Sigmas of I,SIGI',..] "
    bin_opacity = None
      .type = floats(size=2)
      .multiple = True
      .help = "A list of tuples (alpha, idx) with as many or more elements as the current number of binners. List is cast to a string"
    nbins = 1
      .type = int(value_min=1, value_max=40)
  }

  viewer {
    data_array {
      label = None
        .type = str
        .help = "If provided this assigns scene_id with a value corresponding to the numbering " \
                   "order the miller array with this label is found in the reflection data file."
      phasertng_tag = None
        .type = str
        .help = "If provided this assigns scene_id with a value corresponding to the numbering " \
                   "order the miller array with a label found in the parsed history of the MTZ header."
      datatype = None
        .type = str
        .help = "In case label is not found this assigns scene_id with a value corresponding to " \
                   "the first miller array of this data type found in the reflection data file."
    }
    scene_id = None
      .type = int
    ncolourlabels = 6
      .type = int
      .help = "internal"
    #show_symmetry_rotation_axes = False
    #  .type = bool
    show_vector = ''
      .type = str
      .multiple = True
      .help = "Vectors to display. Each show_vector is a stringified python list consisting of the name " \
              "of the vector as the first element and a boolean value as the second element indicating " \
              "visibility of the vector, say show_vector = "['4-fold#2', True]""
    show_all_vectors = 0
      .type = int(value_min=-1, value_max=1)
    user_vector
      .multiple = True
      .help = "Vectors the user add in addition to existing vectors (rotations, TNCS, anisotropy principal axes). " \
              "A vector has to be entered either as a rotation, a real space or a reciprocal space vector. " \
              "The label is required but only one of hkl_op, abc or hkl must be specified"
    {
      label = ""
        .type = str
      hkl_op = ""
        .type = str
        .help = "Rotation operation specified with h,k and l"
      abc = ""
        .type = str
        .help = "Real space vector in real space fractional coordinates"
      hkl = ""
        .type = str
        .help = "Reciprocal space vector in reciprocal space fractional coordinates"
    }
    show_hkl = ""
      .type = str
      .help = "Highlight a reflection with a red meshed wire net surrounding it."
    is_parallel = False
      .type = bool
      .help = "Specifies if reciprocal space is rotated to have a selected and displayed vector being parallel " \
              "or perpendicular to the screen."
    fixorientation = vector *None
      .type = choice
      .help = "fixes orientation of reciprocal space to be aligned with a vector so only mouse zoom will work"
    angle_around_XHKL_vector = 0.0
      .type = float
    angle_around_YHKL_vector = 0.0
      .type = float
    angle_around_ZHKL_vector = 0.0
      .type = float
    angle_around_vector = \"[0,0]\"
      .type = str
      .help = "Rotation with a specified angle of all reflections around a specified vector, " \
              "say angle_around_vector = \"['2-fold#5', 13.0]\""
    animate_rotation_around_vector = \"[0,0]\"
      .type = str
      .help = "Continuous rotation of all reflections around a specified vector at a certain speed, " \
              "say animate_rotation_around_vector = \"['2-fold#5', 10.0]\""
  }
  hkls {
    %s
  }
  NGL {
    %s
  }
  action = *is_running is_terminating reset_view
    .type = choice
  tabulate_miller_array_ids = "[]"
    .type = str
  tooltip_data = "[]"
    .type = str

""" %(ArrayInfo.arrayinfo_phil_str, display.philstr, jsview_3d.ngl_philstr)

master_phil = libtbx.phil.parse( master_phil_str )
