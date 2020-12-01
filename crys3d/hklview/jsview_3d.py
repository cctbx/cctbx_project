
from __future__ import absolute_import, division, print_function
from libtbx.math_utils import roundoff
import traceback
from cctbx.miller import display2 as display
from cctbx.array_family import flex
from cctbx import miller
from scitbx import graphics_utils
from scitbx import matrix
import scitbx.math
from libtbx.utils import Sorry, to_str
import threading, math, sys, cmath
if sys.version_info[0] > 2: # using websockets which is superior to websocket_server
  from crys3d.hklview.WebBrowserMessengerPy3 import WBmessenger
else: # using websocket_server
  from crys3d.hklview.WebBrowserMessengerPy2 import WBmessenger

import os.path, time, copy
import libtbx
import webbrowser, tempfile
from six.moves import range


def has_phil_path(philobj, *paths): # variable number of arguments
  for path in paths:
    if len([ e.path for e in philobj.all_definitions() if path in e.path.split(".") ]):
      return True
  return False


class ArrayInfo:
  def __init__(self, millarr, mprint=sys.stdout.write, fomlabel=None):
    from iotbx.gui_tools.reflections import get_array_description
    if (millarr.unit_cell() is None) or (millarr.space_group() is None) :
      raise Sorry("No space group info is present in data")
    data = millarr.data()
    self.datatype = ""
    self.maxdata = self.mindata = self.maxsigmas = self.minsigmas = None
    self.minmaxdata = (None, None)
    self.minmaxsigs = (None, None)
    arrsize = data.size()
    if not isinstance(data, flex.std_string):
      if (isinstance(data, flex.hendrickson_lattman)):
        data = graphics_utils.NoNansHL( data )
        # for now display HL coefficients as a simple sum
        self.maxdata = max([e[0]+e[1]+e[2]+e[3] for e in data ])
        self.mindata = min([e[0]+e[1]+e[2]+e[3] for e in data ])
        arrsize = len([42 for e in millarr.data() if not math.isnan(e[0]+e[1]+e[2]+e[3])])
        self.datatype = "ishendricksonlattman"
      else:
        arrsize = len([42 for e in millarr.data() if not math.isnan(abs(e))])
        if (isinstance(data, flex.int)):
          data = flex.double([e for e in data if e!= display.inanval])
          self.datatype = "isint"
        if millarr.is_complex_array():
          data = flex.abs(millarr.data())
          self.datatype = "iscomplex"
        #data = [e for e in data if not math.isnan(e)]
        data = graphics_utils.NoNansArray( data, data[0] ) # assuming data[0] isn't NaN
        self.maxdata = flex.max( data )
        self.mindata = flex.min( data )

      if millarr.sigmas() is not None:
        data = millarr.sigmas()
        self.datatype = "hassigmas"
        #data = [e for e in data if not math.isnan(e)]
        data = graphics_utils.NoNansArray( data, data[0] ) # assuming data[0] isn't NaN
        self.maxsigmas = flex.max( data )
        self.minsigmas = flex.min( data )
      self.minmaxdata = (roundoff(self.mindata), roundoff(self.maxdata))
      self.minmaxsigs = (roundoff(self.minsigmas), roundoff(self.maxsigmas))

    self.labels = self.desc = ""
    #import code, traceback; code.interact(local=locals(), banner="".join( traceback.format_stack(limit=10) ) )
    if millarr.info():
      #self.labels = millarr.info().label_string()
      self.labels = millarr.info().labels
      if fomlabel:
        self.labels = [ millarr.info().label_string() + " + " + fomlabel ]
        self.datatype = "iscomplex_fom"
      self.desc = get_array_description(millarr)
    self.span = ("?" , "?")
    self.spginf = millarr.space_group_info().symbol_and_number()
    dmin = 0.0
    dmax = 0.0
    try:
      self.span = ( millarr.index_span().min(), millarr.index_span().max())
      dmin = millarr.d_max_min()[1]
      dmax = millarr.d_max_min()[0]
    except Exception as e:
      mprint(to_str(e))
    issymunique = millarr.is_unique_set_under_symmetry()
    isanomalous = millarr.anomalous_flag()

    self.infotpl = ( ",".join(self.labels), self.desc, self.spginf, arrsize, self.span,
     self.minmaxdata, self.minmaxsigs, (roundoff(dmin), roundoff(dmax)), issymunique, isanomalous )
    self.infostr = "%s (%s), space group: %s, %s HKLs: %s, MinMax: %s, MinMaxSigs: %s, d_minmax: %s, SymUnique: %d, Anomalous: %d" %self.infotpl


def MakeHKLscene( proc_array, pidx, setts, mapcoef_fom_dict, merge, mprint=sys.stdout.write):
  scenemaxdata =[]
  scenemindata =[]
  scenemaxsigmas = []
  sceneminsigmas = []
  scenearrayinfos = []
  hklscenes = []
  fomsarrays_idx = [(None, None)]
  #mprint("in MakeHKLscene", verbose=True)
  #import code, traceback; code.interact(local=locals(), banner="".join( traceback.format_stack(limit=10) ) )
  if proc_array.is_complex_array():
    fomsarrays_idx.extend( mapcoef_fom_dict.get(proc_array.info().label_string()) )
  settings = setts
  if (settings.expand_anomalous or settings.expand_to_p1) \
      and not proc_array.is_unique_set_under_symmetry() and not merge:
    #settings = copy.deepcopy(settings)
    settings.expand_anomalous = False
    settings.expand_to_p1 = False
    mprint("The " + proc_array.info().label_string() + \
         " array is not symmetry unique and therefore won't be expanded")
  if (settings.inbrowser==True):
    settings.expand_anomalous = False
    settings.expand_to_p1 = False

  for (fomsarray, fidx) in fomsarrays_idx:
    hklscene = display.scene(miller_array=proc_array, merge=merge,
      settings=settings, foms_array=fomsarray, fullprocessarray=True )
    if not hklscene.SceneCreated:
      mprint("The " + proc_array.info().label_string() + " array was not processed")
      #return False
      continue
    #import code, traceback; code.interact(local=locals(), banner="".join( traceback.format_stack(limit=10) ) )
    # cast any NAN values to 1 of the colours and radii to 0.2 before writing javascript
    if hklscene.SceneCreated:
      hklscenes.append( hklscene)
      #b = flex.bool([math.isnan(e[0] + e[1] + e[2]) for e in hklscene.colors])
      #hklscene.colors = hklscene.colors.set_selected(b, (1.0, 1.0, 1.0))
      hklscene.colors = graphics_utils.NoNansvec3( hklscene.colors, 1.0, 1.0, 1.0)
      #b = flex.bool([math.isnan(e) for e in hklscene.radii])
      #hklscene.radii = hklscene.radii.set_selected(b, 0.2)
      hklscene.radii = graphics_utils.NoNansArray( hklscene.radii, 0.2)
      fomslabel = None
      if fomsarray:
        fomslabel = fomsarray.info().label_string()
      ainf = ArrayInfo(hklscene.work_array, fomlabel=fomslabel)
      scenemaxdata.append( ainf.maxdata )
      scenemindata.append( ainf.mindata )
      scenemaxsigmas.append(ainf.maxsigmas)
      sceneminsigmas.append(ainf.minsigmas)
      scenearrayinfos.append([ainf.infostr, pidx, fidx, ainf.labels, ainf.datatype])
      #self.mprint("%d, %s" %(i, infostr) )
      #i +=1
  return (hklscenes, scenemaxdata, scenemindata, scenemaxsigmas, sceneminsigmas, scenearrayinfos)


def MakeTtips(hklscene, j):
  tooltipstringsdict = {}
  colstraliases = ""
  if hklscene.isUsingFOMs():
    return tooltipstringsdict, colstraliases # already have tooltips for the scene without the associated fom
  colstraliases += "\n  var st%d = '\\n%s: '" %(j, hklscene.work_array.info().label_string() )
  ocolstr = hklscene.work_array.info().label_string()
  if hklscene.work_array.is_complex_array():
    ampl = flex.abs(hklscene.data)
    phases = flex.arg(hklscene.data) * 180.0/math.pi
    # purge nan values from array to avoid crash in fmod_positive()
    #b = flex.bool([bool(math.isnan(e)) for e in phases])
    # replace the nan values with an arbitrary float value
    #phases = phases.set_selected(b, 42.4242)
    phases = graphics_utils.NoNansArray( phases, 42.4242)
    # Cast negative degrees to equivalent positive degrees
    phases = flex.fmod_positive(phases, 360.0)
  sigmas = hklscene.sigmas
  for i,datval in enumerate(hklscene.data):
    od =""
    if hklscene.work_array.is_complex_array():
      od = str(roundoff(ampl[i], 2)) + ", " + str(roundoff(phases[i], 1)) + \
        "\'+DGR+\'"
    elif sigmas is not None:
      od = str(roundoff(datval, 2)) + ", " + str(roundoff(sigmas[i], 2))
    else:
      od = str(roundoff(datval, 2))
    if not (math.isnan( abs(datval) ) or datval == display.inanval):
      hkl = hklscene.indices[i]
      if not hkl in tooltipstringsdict:
        spbufttip = '\'+hk+\'%s, %s, %s' %(hkl[0], hkl[1], hkl[2])
        spbufttip += '\ndres: %s ' %str(roundoff(hklscene.dres[i], 2) )
        spbufttip += '\'+AA+\'' # javascript alias for angstrom
        tooltipstringsdict[hkl] = spbufttip
      # st1, st2,... are javascript aliases for miller array labelstrings as declared in colstraliases
      tooltipstringsdict[hkl] += '\'+st%d+\'%s' %(j, od)
  return tooltipstringsdict, colstraliases


class hklview_3d:
  def __init__ (self, *args, **kwds) :
    self.settings = kwds.get("settings")
    self.ngl_settings = None #NGLsettings()
    self.viewerparams = kwds.get("settings")
    self.diff_phil = None
    self.params = None
    self.miller_array = None
    self.symops = []
    self.sg = None
    self.tooltipstrings = []
    self.tooltipstringsdict = {}
    self.d_min = None
    self.scene = None
    self.lastscene_id = None
    self.merge = False
    self.primitivetype = "SphereBuffer"
    self.url = ""
    self.bin_labels_type_idxs = []
    self.colour_scene_id = None
    self.radii_scene_id = None
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
    self.boundingX = None
    self.boundingY = None
    self.boundingZ = None
    self.OrigClipNear = None
    self.OrigClipFar = None
    self.cameratranslation = ( 0,0,0 )
    #self.angle_x_svec = 0.0
    #self.angle_y_svec = 0.0
    self.angle_z_svec = 0.0
    #self.angle_z_yzvec = 0.0
    #self.angle_y_yzvec = 0.0
    #self.angle_y_xyvec = 0.0
    self.angle_x_xyvec = 0.0
    self.vecrotmx = None
    self.unit_h_axis = None
    self.unit_k_axis = None
    self.unit_l_axis = None
    self.normal_hk = None
    self.normal_kl = None
    self.normal_lh = None
    self.isnewfile = False
    self.has_new_miller_array = False
    self.sleeptime = 0.01 # 0.025
    self.colstraliases = ""
    self.binvals = []
    self.binvalsboundaries = []
    self.oldnbinvalsboundaries = None
    self.proc_arrays = []
    self.HKLscene = []
    self.HKLscenes = []
    self.HKLscenedict = {}
    self.HKLscenesdict = {}
    self.HKLscenesMaxdata = []
    self.HKLscenesMindata = []
    self.HKLscenesMaxsigmas = []
    self.HKLscenesMinsigmas = []
    self.bindata = None
    self.reciproc_scale = 1.0
    self.realspace_scale = 1.0
    self.sceneisdirty = True
    self.imagename = None
    self.imgdatastr = ""
    self.hkl_scenes_info = []
    self.hkl_scenes_infos = []
    self.match_valarrays = []
    self.array_infostrs = []
    self.array_infotpls = []
    self.binstrs = []
    self.nuniqueval = 0
    self.bin_infotpls = []
    self.mapcoef_fom_dict = {}
    self.sceneid_from_arrayid = []
    self.parent = None
    if 'parent' in kwds:
      self.parent = kwds['parent']
    self.verbose = 0
    if 'verbose' in kwds:
      self.verbose = eval(kwds['verbose'])
    self.debug = None
    if 'debug' in kwds:
      self.debug = kwds['debug']
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
    self.send_info_to_gui = None
    if 'send_info_to_gui' in kwds:
      self.send_info_to_gui = kwds['send_info_to_gui']
    self.mprint('Output will be written to \"%s\"'  %self.hklfname)
    self.hklhtml = r"""
<!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 4.0 Transitional//EN">
<html><head><meta charset="utf-8" /></head>
<body>
<script>var websocket_portnumber = %s; </script>
<script src="%s" type="text/javascript"></script>
<script src="%s" type="text/javascript"></script>
<div id="viewport" style="width:100%%; height:100%%;"></div>
</body></html>

    """
    NGLlibpath = libtbx.env.under_root(os.path.join("modules","cctbx_project","crys3d","hklview","ngl.js") )
    HKLjscriptpath = libtbx.env.under_root(os.path.join("modules","cctbx_project","crys3d","hklview","HKLJavaScripts.js") )
    self.htmlstr = self.hklhtml %(self.websockport, NGLlibpath, os.path.abspath( HKLjscriptpath))
    self.colourgradientvalues = []
    self.UseOSBrowser = ""
    ldic=locals()
    if 'UseOSBrowser' in kwds:
      exec("UseOSBrowser = kwds['UseOSBrowser']", globals(), ldic)
      self.UseOSBrowser = ldic["UseOSBrowser"]
    self.viewmtrx = None
    self.lastviewmtrx = None
    self.currentRotmx = matrix.identity(3)
    self.HKLsceneKey = ( 0, False, self.viewerparams.expand_anomalous, self.viewerparams.expand_to_p1  )
    self.handshakewait = 5
    if 'handshakewait' in kwds:
      self.handshakewait = eval(kwds['handshakewait'])
    self.lastmsg = "" # "Ready"
    self.boundingbox_msg_sem = threading.Semaphore()
    self.clipplane_msg_sem = threading.Semaphore()
    self.mousespeed_msg_sem = threading.Semaphore()
    self.WBmessenger = WBmessenger(self)
    self.AddToBrowserMsgQueue = self.WBmessenger.AddToBrowserMsgQueue
    self.WBmessenger.StartWebsocket()
    self.javascriptcleaned = False


  def __exit__(self, exc_type, exc_value, traceback):
    # not called unless instantiated with a "with hklview_3d ... " statement
    self.JavaScriptCleanUp()
    nwait = 0
    while not self.WBmessenger.isterminating and nwait < 5:
      time.sleep(self.sleeptime)
      nwait += self.sleeptime
    if os.path.isfile(self.hklfname):
      os.remove(self.hklfname)
    self.mprint("Destroying hklview_3d", 1)


  def SendInfoToGUI(self, mydict):
    if self.send_info_to_gui:
      self.send_info_to_gui( mydict )


  def update_settings(self, diff_phil, curphilparam) :
    self.ngl_settings = curphilparam.NGL
    self.viewerparams = curphilparam.viewer
    self.params = curphilparam
    self.diff_phil = diff_phil

    if has_phil_path(diff_phil,
                     "openfilename",
                     "use_provided_miller_arrays",
                     "spacegroup_choice",
                     "using_space_subgroup",
                     "merge_data",
                     "camera_type",
                     "miller_array_operations",
                     ) \
     or has_phil_path(diff_phil, "viewer") \
     and has_phil_path(diff_phil,
                       "show_missing",
                       "show_only_missing",
                       "show_systematic_absences",
                       "slice_axis",
                       "slice_mode",
                       "slice_index",
                       "sigma_color",
                       "sigma_radius",
                       "scene_id",
                       "color_scheme",
                       "color_powscale",
                       "scale",
                       "nth_power_scale_radii"
                       ) \
     or self.viewerparams.inbrowser==False and \
      ( has_phil_path(diff_phil,
                     "expand_anomalous",
                     "expand_to_p1",
                     "show_anomalous_pairs")
       ):
        self.sceneisdirty = True
        if has_phil_path(diff_phil, "scene_id",
                         "scale",
                         "nth_power_scale_radii"
            ):
          self.ConstructReciprocalSpace(curphilparam, scene_id=self.viewerparams.scene_id )
        else:
          self.ConstructReciprocalSpace(curphilparam )
    msg = ""
    if self.viewerparams.scene_id is not None and \
      ( has_phil_path(diff_phil,
                      "show_missing",
                      "show_only_missing",
                      "show_systematic_absences",
                      "binner_idx",
                      "nbins",
                      "fontsize",
                      )
       ) and not has_phil_path(diff_phil, "scene_bin_thresholds") :
      self.binvals, self.nuniqueval = self.calc_bin_thresholds(curphilparam.binner_idx, 
                                                               curphilparam.nbins)
      self.sceneisdirty = True
    if has_phil_path(diff_phil, "scene_bin_thresholds"):
      self.sceneisdirty = True

    if has_phil_path(diff_phil, "camera_type"):
      self.set_camera_type()

    if has_phil_path(diff_phil, "show_tooltips"):
      self.set_show_tooltips()

    if has_phil_path(diff_phil, "fixorientation"):
      self.fix_orientation()

    if has_phil_path(diff_phil, "expand_to_p1", "expand_anomalous") \
     and self.viewerparams.inbrowser and not self.viewerparams.slice_mode:
      self.ExpandInBrowser()

    if has_phil_path(diff_phil, "miller_array_operations"):
      self.viewerparams.scene_id = len(self.HKLscenedict)-1
      self.set_scene(self.viewerparams.scene_id)
      self.params.miller_array_operations = ""

    if has_phil_path(diff_phil, 
                      "scene_bin_thresholds", # TODO: group bin phil parameters together in subscope
                      "bin_opacities",
                      "binner_idx",
                      "nbins",
                      "fontsize",
                      "mouse_sensitivity",
                      "real_space_unit_cell_scale_fraction",
                      "reciprocal_unit_cell_scale_fraction",
                      "clip_plane",
                      "viewer") and self.viewerparams.scene_id is not None:
       # any change to parameters in the master phil in display2.py
      self.scene = self.HKLscene_from_dict(self.viewerparams.scene_id)
      self.DrawNGLJavaScript()
      self.mprint( "Rendered %d reflections" % self.scene.points.size(), verbose=1)
      self.set_volatile_params()

    if self.viewerparams.scene_id is None:
      self.DrawNGLJavaScript(blankscene=True)
    return curphilparam


  def set_volatile_params(self):
    if self.viewerparams.scene_id is not None:
      if has_phil_path(self.diff_phil, "angle_around_vector"): # no need to redraw any clip plane
        return
      self.fix_orientation()
      self.SetMouseSpeed(self.ngl_settings.mouse_sensitivity)
      R = flex.vec3_double( [(0,0,0)])
      hkldist = -1
      clipwidth = None
      isreciprocal = True
      if self.viewerparams.slice_mode: # explicit slicing
        if self.viewerparams.slice_axis=="h": hkl = [1,0,0]
        if self.viewerparams.slice_axis=="k": hkl = [0,1,0]
        if self.viewerparams.slice_axis=="l": hkl = [0,0,1]
        R = hkl[0] * self.normal_kl + hkl[1] * self.normal_lh - hkl[2] * self.normal_hk
        clipwidth = 200
      if self.params.clip_plane.clipwidth and not self.viewerparams.slice_mode:
        clipwidth = self.params.clip_plane.clipwidth
        hkldist = self.params.clip_plane.hkldist
        R = flex.vec3_double( [(self.params.clip_plane.h, self.params.clip_plane.k, self.params.clip_plane.l)])
        if self.params.clip_plane.fractional_vector == "realspace" or self.params.clip_plane.fractional_vector == "tncs":
          isreciprocal = False
      self.clip_plane_vector(R[0][0], R[0][1], R[0][2], hkldist,
          clipwidth, self.ngl_settings.fixorientation, self.params.clip_plane.is_parallel,
          isreciprocal)
      if self.viewerparams.inbrowser and not self.viewerparams.slice_mode:
        self.ExpandInBrowser()
      self.SetOpacities(self.ngl_settings.bin_opacities )
      if self.params.real_space_unit_cell_scale_fraction is None:
        scale = None
      else:
        scale = (self.realspace_scale - 1.0)*self.params.real_space_unit_cell_scale_fraction + 1.0
      self.DrawUnitCell(scale )
      if self.params.reciprocal_unit_cell_scale_fraction is None:
        scale = None
      else:
        scale = (self.reciproc_scale - 1.0)*self.params.reciprocal_unit_cell_scale_fraction + 1.0
      self.DrawReciprocalUnitCell(scale )
      self.set_tooltip_opacity()
      self.set_show_tooltips()


  def set_scene(self, scene_id):
    self.binvals = []
    if scene_id is None:
      return False
    self.colour_scene_id = scene_id
    self.radii_scene_id = scene_id
    self.set_miller_array(scene_id)
    if (self.miller_array is None):
      raise Sorry("No data loaded!")
    self.mprint( "Miller array %s runs from hkls: %s to %s" \
     %(self.miller_array.info().label_string(), self.miller_array.index_span().min(),
        self.miller_array.index_span().max() ) )
    self.mprint("Spacegroup: %s" %self.miller_array.space_group().info().symbol_and_number())
    return True


  def set_miller_array(self, scene_id=None, merge=None, details=""):
    if scene_id is not None:
      self.viewerparams.scene_id = scene_id
    if self.viewerparams and self.viewerparams.scene_id is not None and self.viewerparams.scene_id >= 0 and self.HKLscene:
      self.miller_array = self.HKLscene_from_dict(self.viewerparams.scene_id).miller_array
      self.scene = self.HKLscene_from_dict(self.viewerparams.scene_id)
    self.merge = merge
    if (self.miller_array is None):
      return
    self.identify_suitable_fomsarrays()
    self.GetUnitcellScales()
    self.d_min = self.miller_array.d_min()
    array_info = self.miller_array.info()
    self.sg = self.miller_array.space_group()
    self.symops = self.sg.all_ops()
    if len(self.binvals) == 0:
      self.binvals = [ 1.0/self.miller_array.d_max_min()[0], 1.0/self.miller_array.d_max_min()[1]  ]
    #import code, traceback; code.interact(local=locals(), banner="".join( traceback.format_stack(limit=10) ) )
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


  def GetTooltipOnTheFly(self, id, sym_id, anomalous=False):
    hkl = self.scene.indices[id]
    hklvec = flex.vec3_double( [(hkl[0], hkl[1], hkl[2])])
    rotmx=None
    if sym_id >= 0 and sym_id < len(self.symops):
      rotmx = self.symops[sym_id].r()

    Rhkl = hklvec[0]
    if rotmx:
      Rhkl = hklvec[0] * rotmx
    rothkl = Rhkl
    #import code, traceback; code.interact(local=locals(), banner="".join( traceback.format_stack(limit=10) ) )
    if anomalous:
      rothkl =  (-Rhkl[0], -Rhkl[1], -Rhkl[2])
    spbufttip = '\'H,K,L: %d, %d, %d' %(rothkl[0], rothkl[1], rothkl[2])
    # resolution and angstrom character
    spbufttip += '\\ndres: %s \'+ String.fromCharCode(197) +\'' \
      %str(roundoff(self.miller_array.unit_cell().d(hkl), 2) )
    for hklscene in self.HKLscenes:
      if hklscene.isUsingFOMs():
        continue # already have tooltips for the scene without the associated fom
      sigvals = None
      if hklscene.work_array.sigmas() is not None:
        sigvals = list( hklscene.work_array.select(hklscene.work_array.indices() == hkl).sigmas() )
      datval = None
      if hkl in hklscene.work_array.indices():
        datvals = list( hklscene.work_array.select(hklscene.work_array.indices() == hkl).data() )
      else:
        if id >= hklscene.data.size():
          continue
        datvals = [ hklscene.data[id] ]
      for i,datval in enumerate(datvals):
        if isinstance(datval, tuple) and math.isnan(datval[0] + datval[1] + datval[2] + datval[3]):
          continue
        if not isinstance(datval, tuple) and (math.isnan( abs(datval) ) or datval == display.inanval):
          continue
        if hklscene.work_array.is_complex_array():
          ampl = abs(datval)
          phase = cmath.phase(datval) * 180.0/math.pi
          # purge nan values from array to avoid crash in fmod_positive()
          # and replace the nan values with an arbitrary float value
          if math.isnan(phase):
            phase = 42.4242
          # Cast negative degrees to equivalent positive degrees
          phase = phase % 360.0
        spbufttip +="\\n" + hklscene.work_array.info().label_string() + ': '
        if hklscene.work_array.is_complex_array():
          spbufttip += str(roundoff(ampl, 2)) + ", " + str(roundoff(phase, 1)) + \
            "\'+ String.fromCharCode(176) +\'" # degree character
        elif sigvals:
          sigma = sigvals[i]
          spbufttip += str(roundoff(datval, 2)) + ", " + str(roundoff(sigma, 2))
        else:
          spbufttip += str(roundoff(datval, 2))
    spbufttip += '\\n\\n%d,%d,%d' %(id, sym_id, anomalous) # compared by the javascript
    spbufttip += '\''
    return spbufttip


  def get_col_fomcol(self, idx):
    if len(self.HKLInfo_from_dict()) == 0:
      return -1, -1
    return self.HKLInfo_from_dict(idx)[6], self.HKLInfo_from_dict(idx)[7]


  def SupersetMillerArrays(self):
    self.match_valarrays = []
    # First loop over all miller arrays to make a superset of hkls of all
    # miller arrays. Then loop over all miller arrays and extend them with NaNs
    # as to contain the same hkls as the superset
    self.mprint("Gathering superset of miller indices...")
    superset_array = self.proc_arrays[0].deep_copy()
    set_of_indices = set([])
    for i,procarray in enumerate(self.proc_arrays):
      set_of_indices |= set( list(procarray.indices()) )
    self.mprint("Extending miller arrays to match superset of miller indices...")
    indiceslst = flex.miller_index( list( set_of_indices ) )
    for i,procarray in enumerate(self.proc_arrays):
      # first match indices in currently selected miller array with indices in the other miller arrays
      matchindices = miller.match_indices(indiceslst, procarray.indices() )
      #matchindices = miller.match_indices( procarray.indices(), indiceslst )
      valarray = procarray.select( matchindices.pairs().column(1) )
      #if valarray.anomalous_flag() != superset_array.anomalous_flag():
      #  superset_array._anomalous_flag = valarray._anomalous_flag
      #missing = indiceslst.lone_set( valarray.indices() )

      missing = indiceslst.select( miller.match_indices(valarray.indices(), indiceslst ).singles(1))

      # insert NAN values for reflections in self.miller_array not found in procarray
      valarray = display.ExtendMillerArray(valarray, missing.size(), missing )
      #match_valarray = valarray
      match_valindices = miller.match_indices(superset_array.indices(), valarray.indices() )
      match_valarray = valarray.select( match_valindices.pairs().column(1) )
      match_valarray.sort(by_value="packed_indices")
      match_valarray.set_info(procarray.info() )
      self.match_valarrays.append( match_valarray )
    self.mprint("Done making superset")

  """
  def SupersetMillerArrays(self):
    self.match_valarrays = []
    # First loop over all miller arrays to make a superset of hkls of all
    # miller arrays. Then loop over all miller arrays and extend them with NaNs
    # as to contain the same hkls as the superset
    self.mprint("Gathering superset of miller indices...")
    superset_array = self.proc_arrays[0].deep_copy()
    for i,procarray in enumerate(self.proc_arrays):
      if i==0:
        continue
      # first match indices in currently selected miller array with indices in the other miller arrays
      matchindices = miller.match_indices(superset_array.indices(), procarray.indices() )
      valarray = procarray.select( matchindices.pairs().column(1) )
      if valarray.anomalous_flag() != superset_array.anomalous_flag():
        superset_array._anomalous_flag = valarray._anomalous_flag
      missing = procarray.lone_set( superset_array )
      superset_array = display.ExtendMillerArray(superset_array, missing.size(), missing.indices())
    self.mprint("Extending miller arrays to match superset of miller indices...")
    for i,procarray in enumerate(self.proc_arrays):
      # first match indices in currently selected miller array with indices in the other miller arrays
      matchindices = miller.match_indices(superset_array.indices(), procarray.indices() )
      valarray = procarray.select( matchindices.pairs().column(1) )
      if valarray.anomalous_flag() != superset_array.anomalous_flag():
        superset_array._anomalous_flag = valarray._anomalous_flag
      missing = superset_array.lone_set( valarray )
      # insert NAN values for reflections in self.miller_array not found in procarray
      valarray = display.ExtendMillerArray(valarray, missing.size(), missing.indices())
      match_valindices = miller.match_indices(superset_array.indices(), valarray.indices() )
      match_valarray = valarray.select( match_valindices.pairs().column(1) )
      match_valarray.sort(by_value="packed_indices")
      match_valarray.set_info(procarray.info() )
      self.match_valarrays.append( match_valarray )
    self.mprint("Done making superset")
  """

  def ConstructReciprocalSpace(self, curphilparam, scene_id=None):
    sceneid = scene_id
    if sceneid is None:
      sceneid = self.viewerparams.scene_id

    self.HKLsceneKey = (curphilparam.spacegroup_choice,
                         curphilparam.using_space_subgroup,
                         curphilparam.merge_data,
                         self.viewerparams.expand_anomalous or self.viewerparams.inbrowser,
                         self.viewerparams.expand_to_p1 or self.viewerparams.inbrowser,
                         self.viewerparams.inbrowser,
                         self.viewerparams.slice_axis,
                         self.viewerparams.slice_mode,
                         self.viewerparams.slice_index,
                         self.viewerparams.show_missing,
                         self.viewerparams.show_only_missing,
                         self.viewerparams.show_systematic_absences,
                         self.viewerparams.sigma_radius,
                         self.viewerparams.sigma_color,
                         self.viewerparams.color_scheme,
                         self.viewerparams.color_powscale,
                         self.ngl_settings.fontsize,
                         sceneid,
                         self.viewerparams.scale,
                         self.viewerparams.nth_power_scale_radii
                         )
    if self.HKLsceneKey in self.HKLscenedict and not self.has_new_miller_array:
      self.HKLscene = self.HKLscenedict.get(self.HKLsceneKey, False)
      if self.HKLscene:
        self.mprint("Using cached HKL scene", verbose=1)
        return True
    if self.has_new_miller_array:
      self.identify_suitable_fomsarrays()
    self.mprint("Constructing HKL scenes", verbose=0)
    assert(self.proc_arrays)
    if scene_id is None:
      hkl_scenes_infos = []
      self.HKLscenes = []
      sceneid = 0
      for (idx, arr) in enumerate(self.proc_arrays):
        (hklscenes, scenemaxdata,
          scenemindata, scenemaxsigmas,
            sceneminsigmas, scenearrayinfos
         ) = MakeHKLscene( arr.deep_copy(), idx, copy.deepcopy(self.viewerparams), self.mapcoef_fom_dict, None, self.mprint )

        for i,inf in enumerate(scenearrayinfos):
          self.mprint("%d, %s" %(idx+i+1, inf[0]), verbose=0)
          self.HKLsceneKey = (curphilparam.spacegroup_choice,
                                curphilparam.using_space_subgroup,
                                curphilparam.merge_data,
                                self.viewerparams.expand_anomalous or self.viewerparams.inbrowser,
                                self.viewerparams.expand_to_p1 or self.viewerparams.inbrowser,
                                self.viewerparams.inbrowser,
                                self.viewerparams.slice_axis,
                                self.viewerparams.slice_mode,
                                self.viewerparams.slice_index,
                                self.viewerparams.show_missing,
                                self.viewerparams.show_only_missing,
                                self.viewerparams.show_systematic_absences,
                                self.viewerparams.sigma_radius,
                                self.viewerparams.sigma_color,
                                self.viewerparams.color_scheme,
                                self.viewerparams.color_powscale,
                                self.ngl_settings.fontsize,
                                sceneid,
                                self.viewerparams.scale,
                                self.viewerparams.nth_power_scale_radii
                                )
          self.HKLscenedict[self.HKLsceneKey] = ( hklscenes[i], scenemaxdata[i],
          scenemindata[i], scenemaxsigmas[i], sceneminsigmas[i], inf )
          hkl_scenes_infos.append(inf)
          self.HKLscenes.append(hklscenes[i])
          sceneid += 1
      self.hkl_scenes_infos = hkl_scenes_infos
      if self.viewerparams.scene_id is not None:
        self.HKLsceneKey = (curphilparam.spacegroup_choice,
                              curphilparam.using_space_subgroup,
                              curphilparam.merge_data,
                              self.viewerparams.expand_anomalous or self.viewerparams.inbrowser,
                              self.viewerparams.expand_to_p1 or self.viewerparams.inbrowser,
                              self.viewerparams.inbrowser,
                              self.viewerparams.slice_axis,
                              self.viewerparams.slice_mode,
                              self.viewerparams.slice_index,
                              self.viewerparams.show_missing,
                              self.viewerparams.show_only_missing,
                              self.viewerparams.show_systematic_absences,
                              self.viewerparams.sigma_radius,
                              self.viewerparams.sigma_color,
                              self.viewerparams.color_scheme,
                              self.viewerparams.color_powscale,
                              self.ngl_settings.fontsize,
                              self.viewerparams.scene_id,
                              self.viewerparams.scale,
                              self.viewerparams.nth_power_scale_radii
                              )
      scenearraylabeltypes = [ (e[3], e[4], e[1], sceneid) for sceneid,e in enumerate(hkl_scenes_infos) ]
      self.SendInfoToGUI({ "scene_array_label_types": scenearraylabeltypes, "NewHKLscenes" : True })

      self.bin_labels_type_idxs = []
      self.bin_labels_type_idxs.append(("Resolution",  "", -1 ))
      self.bin_labels_type_idxs.append(("Singletons", "", -1 ))
      for labels,labeltype,idx,sceneid in scenearraylabeltypes:
        label = ",".join(labels)
        if labeltype not in  ["iscomplex", "iscomplex_fom", "ishendricksonlattman"]:
          self.bin_labels_type_idxs.append((label, labeltype, sceneid))
        if labeltype == "hassigmas":
          self.bin_labels_type_idxs.append(("Sigmas of " + label, labeltype, sceneid))
        if labeltype == "iscomplex":
          self.bin_labels_type_idxs.append(("Phases of " + label, labeltype, sceneid))
          self.bin_labels_type_idxs.append(("Amplitudes of " + label, labeltype, sceneid))
      self.SendInfoToGUI({ "bin_labels_type_idxs": self.bin_labels_type_idxs})
      self.get_labels_of_data_for_binning()
    else:
      idx = self.scene_id_to_array_id(scene_id)
      (hklscenes, scenemaxdata,
        scenemindata, scenemaxsigmas,
          sceneminsigmas, scenearrayinfos
      ) = MakeHKLscene( self.proc_arrays[idx].deep_copy(), idx, copy.deepcopy(self.viewerparams), self.mapcoef_fom_dict, None, self.mprint )
      for i,inf in enumerate(scenearrayinfos):
        self.mprint("%d, %s" %(idx+i+1, inf[0]), verbose=0)
        self.HKLsceneKey = (curphilparam.spacegroup_choice,
                              curphilparam.using_space_subgroup,
                              curphilparam.merge_data,
                              self.viewerparams.expand_anomalous or self.viewerparams.inbrowser,
                              self.viewerparams.expand_to_p1 or self.viewerparams.inbrowser,
                              self.viewerparams.inbrowser,
                              self.viewerparams.slice_axis,
                              self.viewerparams.slice_mode,
                              self.viewerparams.slice_index,
                              self.viewerparams.show_missing,
                              self.viewerparams.show_only_missing,
                              self.viewerparams.show_systematic_absences,
                              self.viewerparams.sigma_radius,
                              self.viewerparams.sigma_color,
                              self.viewerparams.color_scheme,
                              self.viewerparams.color_powscale,
                              self.ngl_settings.fontsize,
                              sceneid,
                              self.viewerparams.scale,
                              self.viewerparams.nth_power_scale_radii
                              )
        self.HKLscenedict[self.HKLsceneKey] =  ( hklscenes[i], scenemaxdata[i],
          scenemindata[i], scenemaxsigmas[i], sceneminsigmas[i], inf )
        sceneid += 1
    (
      self.HKLscene,
      self.HKLscenesMaxdata,
      self.HKLscenesMindata,
      self.HKLscenesMaxsigmas,
      self.HKLscenesMinsigmas,
      self.hkl_scenes_info
    ) =  self.HKLscenedict[self.HKLsceneKey]
    self.sceneisdirty = True
    self.has_new_miller_array = False
    return True


  def Sceneid_to_SceneKey(self, sceneid):
    return (self.params.spacegroup_choice,
                      self.params.using_space_subgroup,
                      self.params.merge_data,
                      self.viewerparams.expand_anomalous or self.viewerparams.inbrowser,
                      self.viewerparams.expand_to_p1 or self.viewerparams.inbrowser,
                      self.viewerparams.inbrowser,
                      self.viewerparams.slice_axis,
                      self.viewerparams.slice_mode,
                      self.viewerparams.slice_index,
                      self.viewerparams.show_missing,
                      self.viewerparams.show_only_missing,
                      self.viewerparams.show_systematic_absences,
                      self.viewerparams.sigma_radius,
                      self.viewerparams.sigma_color,
                      self.viewerparams.color_scheme,
                      self.viewerparams.color_powscale,
                      self.ngl_settings.fontsize,
                      sceneid,
                      self.viewerparams.scale,
                      self.viewerparams.nth_power_scale_radii
                      )


  def HKLscene_from_dict(self, sceneid=None):
    if sceneid is None:
      sceneid = self.viewerparams.scene_id
    HKLsceneKey = self.Sceneid_to_SceneKey(sceneid)
    if not self.HKLscenedict.get(HKLsceneKey, False):
      self.ConstructReciprocalSpace(self.params, scene_id=sceneid)
    return self.HKLscenedict[HKLsceneKey][0]


  def HKLMaxData_from_dict(self, sceneid=None):
    if sceneid is None:
      sceneid = self.viewerparams.scene_id
    HKLsceneKey = self.Sceneid_to_SceneKey(sceneid)
    return self.HKLscenedict[HKLsceneKey][1]


  def HKLMinData_from_dict(self, sceneid=None):
    if sceneid is None:
      sceneid = self.viewerparams.scene_id
    HKLsceneKey = self.Sceneid_to_SceneKey(sceneid)
    return self.HKLscenedict[HKLsceneKey][2]


  def HKLMaxSigmas_from_dict(self, sceneid=None):
    if sceneid is None:
      sceneid = self.viewerparams.scene_id
    HKLsceneKey = self.Sceneid_to_SceneKey(sceneid)
    return self.HKLscenedict[HKLsceneKey][3]


  def HKLMinSigmas_from_dict(self, sceneid=None):
    if sceneid is None:
      sceneid = self.viewerparams.scene_id
    HKLsceneKey = self.Sceneid_to_SceneKey(sceneid)
    return self.HKLscenedict[HKLsceneKey][4]


  def HKLInfo_from_dict(self, sceneid=None):
    if sceneid is None:
      sceneid = self.viewerparams.scene_id
    HKLsceneKey = self.Sceneid_to_SceneKey(sceneid)
    return self.HKLscenedict[HKLsceneKey][5]


  def identify_suitable_fomsarrays(self):
    self.mprint("Matching complex arrays to suitable FOM arrays")
    self.mapcoef_fom_dict = {}
    self.sceneid_from_arrayid = []
    for k,proc_array in enumerate(self.proc_arrays):
      #if not proc_array.is_complex_array() or not proc_array.is_real_array():
      #  continue
      fom_arrays_idx = []
      array_scene_ids = [(k,k)]
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


  def scene_id_to_array_id(self, scene_id):
    for i,array_scene_id in enumerate(self.sceneid_from_arrayid):
      if scene_id == i:
        return array_scene_id[0]
    return None


  def calc_bin_thresholds(self, binner_idx, nbins):
    # make default bin thresholds if scene_bin_thresholds is not set
    #self.bin_labels_type_idx = self.bin_labels_type_idxs[binner_idx]
    binscenelabel = self.bin_labels_type_idxs[binner_idx][0]
    self.mprint("Using %s for binning" %binscenelabel)
    if binscenelabel=="Resolution":
      warray = self.HKLscene_from_dict(int(self.viewerparams.scene_id)).work_array
      dres = self.HKLscene_from_dict(int(self.viewerparams.scene_id)).dres
      uc = warray.unit_cell()
      indices = self.HKLscene_from_dict(int(self.viewerparams.scene_id)).indices
      if flex.max(dres) == flex.min(dres): # say if only one reflection
        binvals = [dres[0]-0.1, flex.min(dres)+0.1]
        nuniquevalues = 2
      else: # use generic binning function from cctbx
        binning = miller.binning( uc, nbins, indices, flex.max(dres), flex.min(dres) )
        binvals = [ binning.bin_d_range(n)[0] for n in binning.range_all() ]
        binvals = [ e for e in binvals if e != -1.0] # delete dummy limit
        binvals = list( 1.0/flex.double(binvals) )
        nuniquevalues = len(set(list(dres)))
    elif binscenelabel=="Singletons":
        binvals = [ -1.5, -0.5, 0.5, 1.5 ]
        nuniquevalues = len(binvals)
    else:
      bindata, dummy = self.get_matched_binarray(binner_idx)
      selection = flex.sort_permutation( bindata )
      bindata_sorted = bindata.select(selection)
      # get binvals by dividing bindata_sorted with nbins
      binvals = [bindata_sorted[0]] * nbins #
      for i,e in enumerate(bindata_sorted):
        idiv = int(nbins*float(i)/len(bindata_sorted))
        binvals[idiv] = e
      nuniquevalues = len(set(list(bindata)))
    binvals.sort()
    self.mprint("Bin thresholds are:\n" + str(binvals))
    return binvals, nuniquevalues


  def UpdateBinValues(self, binner_idx, binvals = [], nuniquevalues = -1):
    if binvals:
      binvals.sort()
      self.binvals = binvals
    else: # ensure default resolution interval includes all data by avoiding rounding errors
      self.binvals = [ 1.0/(self.miller_array.d_max_min()[0]*1.001),
                       1.0/(self.miller_array.d_max_min()[1]*0.999) ]
    if nuniquevalues == -1:
      if binner_idx==0:
        nuniquevalues = len(set(list( self.HKLscene_from_dict(int(self.viewerparams.scene_id)).dres )))
      else:
        bindata, dummy = self.get_matched_binarray(binner_idx)
        nuniquevalues = len(set(list(bindata)))
    self.nuniqueval = nuniquevalues


  def get_matched_binarray(self, binner_idx):
    binscenelabel, datatype, sceneid = self.bin_labels_type_idxs[binner_idx]
    label = self.HKLscene_from_dict(sceneid).work_array.info().label_string()
    if datatype == "hassigmas" and binscenelabel == "Sigmas of " + label:
      bindata = self.HKLscene_from_dict(sceneid).sigmas.deep_copy()
      binvalsboundaries = [ self.HKLMinSigmas_from_dict(sceneid) - 0.1 , self.HKLMaxSigmas_from_dict(sceneid) + 0.1 ]
    elif datatype == "iscomplex" and "Phases of " + label in binscenelabel:
      bindata = self.HKLscene_from_dict(sceneid).phases.deep_copy()
      # preselect centric reflections, i.e. those with phi = 0 or 180
      binvalsboundaries = [-0.01, 0.01, 179.99, 180.01, 359.99, 360]
    elif datatype == "iscomplex" and "Amplitudes of " + label in binscenelabel:
      bindata = self.HKLscene_from_dict(sceneid).ampl.deep_copy()
      binvalsboundaries = [ self.HKLMinData_from_dict(sceneid) - 0.1 , self.HKLMaxData_from_dict(sceneid) + 0.1 ]
    else:
      bindata = self.HKLscene_from_dict(sceneid).data.deep_copy()
      binvalsboundaries = [ self.HKLMinData_from_dict(sceneid) - 0.1 , self.HKLMaxData_from_dict(sceneid) + 0.1 ]
    return bindata, binvalsboundaries


  def MatchBinArrayToSceneArray(self):
    # match bindata with data or sigmas
    if self.bin_labels_type_idxs[self.params.binner_idx][0] == "Resolution":
      return 1.0/self.scene.dres
    binarraydata, dummy = self.get_matched_binarray(self.params.binner_idx)
    scenearraydata = self.HKLscene_from_dict(self.viewerparams.scene_id).data
    ibinarray = self.bin_labels_type_idxs[self.params.binner_idx][2]
    matchindices = miller.match_indices(self.HKLscene_from_dict(self.viewerparams.scene_id).indices,
                                        self.HKLscene_from_dict(ibinarray).indices )
    matched_binarray = binarraydata.select( matchindices.pairs().column(1) )
    #valarray.sort(by_value="packed_indices")
    #import code, traceback; code.interact(local=locals(), banner="".join( traceback.format_stack(limit=10) ) )
    #missing = scenearraydata.lone_set( valarray )
    # insert NAN values for reflections in self.miller_array not found in binarray
    #valarray = display.ExtendMillerArray(valarray, missing.size(), missing.indices() )
    #match_valindices = miller.match_indices(scenearray.indices(), valarray.indices() )
    #match_valarray = valarray.select( match_valindices.pairs().column(1) )
    #match_valarray.sort(by_value="packed_indices")
    #match_valarray.set_info(binarraydata.info() )
    # patch the bin array so its sequence matches the scene array
    patched_binarraydata = []
    c = 0
    for b in matchindices.pair_selection(0):
      if b:
        patched_binarraydata.append(matched_binarray[c])
        c +=1
      else:
        patched_binarraydata.append(float("nan"))
    return flex.double(patched_binarraydata)


  def OperateOn1MillerArray(self, millarr, operation):
    # lets user specify a one line python expression operating on data, sigmas
    newarray = millarr.deep_copy()
    data = newarray.data()
    sigmas = newarray.sigmas()
    dres = newarray.unit_cell().d( newarray.indices() )
    self.mprint("Creating new miller array through the operation: %s" %operation)
    try:
      newdata = None
      newsigmas = None
      ldic= {'data': data, 'sigmas': sigmas, 'dres': dres }
      exec(operation, globals(), ldic)
      newdata = ldic.get("newdata", None)
      newarray._data = newdata
      newsigmas = ldic.get("newsigmas", None)
      newarray._sigmas = newsigmas
      return newarray
    except Exception as e:
      self.mprint( str(e), verbose=0)
      return None


  def OperateOn2MillerArrays(self, millarr1, millarr2, operation):
    # lets user specify a one line python expression operating on data1 and data2
    matchindices = miller.match_indices(millarr1.indices(), millarr2.indices() )
    matcharr1 = millarr1.select( matchindices.pairs().column(0) ).deep_copy()
    matcharr2 = millarr2.select( matchindices.pairs().column(1) ).deep_copy()
    data1 = matcharr1.data()
    data2 = matcharr2.data()
    sigmas1 = matcharr1.sigmas()
    sigmas2 = matcharr2.sigmas()
    dres = matcharr1.unit_cell().d( matcharr1.indices() )
    newarray = matcharr2.deep_copy()
    newarray._sigmas = None
    self.mprint("Creating new miller array through the operation: %s" %operation)
    try:
      newdata = None
      newsigmas = None
      ldic= {'data1': data1, 'sigmas1': sigmas1, 'data2': data2, 'sigmas2': sigmas2, 'dres': dres }
      exec(operation, globals(), ldic)
      newdata = ldic.get("newdata", None)
      newarray._data = newdata
      newsigmas = ldic.get("newsigmas", None)
      newarray._sigmas = newsigmas
      return newarray
    except Exception as e:
      self.mprint( str(e), verbose=0)
      return None


  def DrawNGLJavaScript(self, blankscene=False):
    if not self.scene or not self.sceneisdirty:
      return
    if self.miller_array is None :
      self.mprint( "Select a data set to display reflections" )
      blankscene = True
      #return
    else:
      self.mprint("Rendering reflections...")

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
    l1 = self.scene.renderscale * maxnorm * 1.1
    l2= self.scene.renderscale * maxnorm * 1.15
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
      # Make colour gradient array used for drawing a bar of colours next to associated values on the rendered html
      mincolourscalar = self.HKLMinData_from_dict(self.colour_scene_id)
      maxcolourscalar = self.HKLMaxData_from_dict(self.colour_scene_id)
      if self.viewerparams.sigma_color:
        mincolourscalar = self.HKLMinSigmas_from_dict(self.colour_scene_id)
        maxcolourscalar = self.HKLMaxSigmas_from_dict(self.colour_scene_id)
      span = maxcolourscalar - mincolourscalar
      ln = 60
      incr = span/ln
      colourgradarrays = []
      val = mincolourscalar
      colourscalararray = flex.double()
      colourscalararray.append( val )
      for j,sc in enumerate(range(ln)):
        val += incr
        colourscalararray.append( val )
      if self.HKLscene_from_dict(self.colour_scene_id).miller_array.is_complex_array():
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
        COL = display.MplColorHelper(self.viewerparams.color_scheme, 0, 360)
        rgbcolarray = flex.vec3_double( [ COL.get_rgb(d)[0:3] for d in colourscalararray ] )

        if self.HKLscene_from_dict(self.colour_scene_id).isUsingFOMs():
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
        COL = display.MplColorHelper(self.viewerparams.color_scheme, mincolourscalar, maxcolourscalar)
        rgbcolarray = flex.vec3_double( [ COL.get_rgb(d)[0:3] for d in colourscalararray ])

        arr = graphics_utils.map_to_rgb_colourmap(
            data_for_colors= colourscalararray,
            colormap = rgbcolarray,
            selection=flex.bool(colourscalararray.size(), True),
            powscale = self.viewerparams.color_powscale
          )

        colourgradarrays.append(arr*256)
      colors = self.HKLscene_from_dict(self.colour_scene_id).colors
      radii = self.HKLscene_from_dict(self.radii_scene_id).radii
      self.meanradius = flex.mean(radii)

    bin_labels_type_idx = self.bin_labels_type_idxs[self.params.binner_idx]
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
    elif bin_labels_type_idx[0] =="Singletons":
      colstr = "Singleton"
    else:
      if not blankscene:
        colstr = self.HKLscene_from_dict(bin_labels_type_idx[2]).work_array.info().label_string()
    data = self.scene.data
    if not blankscene:
      colourlabel = self.HKLscene_from_dict(self.colour_scene_id).colourlabel
      fomlabel = self.HKLscene_from_dict(self.colour_scene_id).fomlabel
    #import code, traceback; code.interact(local=locals(), banner="".join( traceback.format_stack(limit=10) ) )
    assert (colors.size() == radii.size() == nrefls)
    self.colours = []
    self.positions = []
    self.radii2 = []
    self.spbufttips = []

    self.binvalsboundaries = []
    if not blankscene:
      if bin_labels_type_idx[0] =="Resolution":
        self.binvalsboundaries = self.binvals
        self.bindata = 1.0/self.scene.dres
      elif bin_labels_type_idx[0] =="Singletons":
        self.binvalsboundaries = self.binvals
        self.bindata = self.scene.singletonsiness
      else:
        dummy, self.binvalsboundaries = self.get_matched_binarray(self.params.binner_idx)
        self.binvalsboundaries.extend( self.binvals )
        self.binvalsboundaries.sort()
        if self.binvalsboundaries[0] < 0.0:
          self.binvalsboundaries.append(0.0)
          self.binvalsboundaries.sort()
        self.bindata = self.MatchBinArrayToSceneArray()

    self.nbinvalsboundaries = len(self.binvalsboundaries)
    # avoid resetting opacities of bins unless we change the number of bins
    if self.oldnbinvalsboundaries != self.nbinvalsboundaries:
      self.ngl_settings.bin_opacities = str([ (1.0, e) for e in range(self.nbinvalsboundaries + 1) ])
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
          return nbinvalsboundaries
        if (ibin+1) == nbinvalsboundaries:
          return ibin
        if d > binval and d <= binvalsboundaries[ibin+1]:
          return ibin
      raise Sorry("data2bin: Should never get here")

    def getprecision(v1,v2):
      diff = abs(v1-v2); precision = 1; e = 1
      while diff*e < 1.0:
        e *= 10
        precision += 1
      return precision

    if nrefls > 0 and self.bindata.size() != points.size():
      raise Sorry("Not the same number of reflections in bin-data and displayed data")

    start_time = time.time()
    for i, hklstars in enumerate(points):
      # bin currently displayed data according to the values of another miller array
      ibin = data2bin( self.bindata[i], self.binvalsboundaries, self.nbinvalsboundaries )
      self.positions[ibin].extend( graphics_utils.flt_roundoffvec3(hklstars, 2) )
      self.colours[ibin].extend( graphics_utils.flt_roundoffvec3(colors[i], 2) )
      self.radii2[ibin].append( graphics_utils.flt_roundoff(radii[i], 2) )
      self.spbufttips[ibin].append( i )

    elapsed_time = time.time() - start_time
    self.mprint("elapsed time: %s" %elapsed_time, verbose=2)

    spherebufferstr = self.colstraliases
    cntbin = 0
    self.binstrs = []
    self.bin_infotpls = []
    if self.nuniqueval < self.params.nbins:
      self.mprint("%d bins was requested but %s data has only %d unique value(s)!" %(self.params.nbins, colstr, self.nuniqueval), 0)
    for ibin in range(self.nbinvalsboundaries+1):
      mstr =""
      nreflsinbin = len(self.radii2[ibin])
      if nreflsinbin == 0:
        continue
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
      self.bin_infotpls.append( roundoff((nreflsinbin, bin1, bin2 ), precision) )
      self.binstrs.append(mstr)
      self.mprint(mstr, verbose=0)
      cntbin += 1

    if self.ngl_settings.bin_opacities == "":
      self.ngl_settings.bin_opacities = str([ (1.0, e) for e in range(cntbin) ])

    self.SendInfoToGUI( { "bin_opacities": self.ngl_settings.bin_opacities,
                          "bin_infotpls": self.bin_infotpls,
                          "binner_idx": self.params.binner_idx,
                          "tooltip_opacity": self.ngl_settings.tooltip_alpha
                         } )
    colourgradstr = []
    if not blankscene:
      colourgradstrs = []
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
          if j%10 == 0 or j==len(self.colourgradientvalues)-1 :
            vstr = roundoff(val[0], 2)
          colourgradstr.append([vstr, rgb[0], rgb[1], rgb[2] ])
        colourgradstrs.append(colourgradstr)  

    if not self.WBmessenger.browserisopen:
      self.ReloadNGL()
    if not blankscene:
      self.RemoveStageObjects()
      for ibin in range(self.nbinvalsboundaries+1):
        nreflsinbin = len(self.radii2[ibin])
        if nreflsinbin == 0:
          continue
        self.SetBrowserDebug(str(self.verbose>=2).lower())
        self.SetFontSize(self.ngl_settings.fontsize)
        self.DefineHKL_Axes(str(Hstararrowstart), str(Hstararrowend),
          str(Kstararrowstart), str(Kstararrowend), 
          str(Lstararrowstart), str(Lstararrowend), 
          Hstararrowtxt, Kstararrowtxt, Lstararrowtxt )
        self.SendCoordinates2Browser(self.positions[ibin], self.colours[ibin], 
                                     self.radii2[ibin], self.spbufttips[ibin] )
      self.RenderStageObjects()
      self.MakeColourChart(10, 10, colourlabel, fomlabel, colourgradstrs)

      if self.WaitforHandshake():
        nwait = 0
        while self.viewmtrx is None and nwait < self.handshakewait:
          time.sleep(self.sleeptime)
          nwait += self.sleeptime
      self.GetClipPlaneDistances()
      self.GetBoundingBox()
      self.OrigClipFar = self.clipFar
      self.OrigClipNear = self.clipNear
      self.SetMouseSpeed( self.ngl_settings.mouse_sensitivity )
      self.isnewfile = False

    self.sceneisdirty = False
    self.lastscene_id = self.viewerparams.scene_id


  def ProcessMessage(self, message):
    try:
      if sys.version_info[0] > 2:
        ustr = str
      else:
        ustr = unicode
      if isinstance(message, bytes) and isinstance(self.lastmsg, ustr) and "Imageblob" in self.lastmsg:
        self.mprint( "Saving image to file", verbose=1)
        with open( self.imagename, "wb") as imgfile:
          imgfile.write( message)

      if isinstance(message, ustr) and message != "":
        if "Orientation" in message:
          self.ProcessOrientationMessage(message)
        elif 'Received message:' in message:
          self.mprint( message, verbose=2)
        elif 'Browser: Got' in message:
          self.mprint( message, verbose=2)
        elif "websocket" in message:
          self.mprint( message, verbose=1)
        elif "Refreshing" in message or "disconnecting" in message:
          self.mprint( message, verbose=1)
          time.sleep(self.sleeptime)
        elif "AutoViewSet" in message:
          self.set_volatile_params()
        elif "JavaScriptCleanUpDone:" in message:
          self.mprint( message, verbose=1)
          time.sleep(0.5) # time for browser to clean up
          if not self.isnewfile:
            self.WBmessenger.StopWebsocket()
        elif "JavaScriptError:" in message:
          self.mprint( message, verbose=0)
        elif "Expand" in message:
          self.mprint( message, verbose=2)
        elif "Connection lost" in message:
          self.mprint( message, verbose=1)
        elif "Warning!: Web browser closed unexpectedly" in message:
          self.mprint( message, verbose=1)
        elif "Imageblob" in message:
          self.mprint( "Image to be received", verbose=1)
        elif "ImageWritten" in message:
          self.mprint( "Image saved to file", verbose=0)
        elif "ReturnClipPlaneDistances:" in message:
          datastr = message[ message.find("\n") + 1: ]
          lst = datastr.split(",")
          flst = [float(e) for e in lst]
          self.clipNear = flst[0]
          self.clipFar = flst[1]
          self.cameraPosZ = flst[2]
          self.clipplane_msg_sem.release()
        elif "ReturnBoundingBox:" in message:
          datastr = message[ message.find("\n") + 1: ]
          lst = datastr.split(",")
          flst = [float(e) for e in lst]
          self.boundingX = flst[0]
          self.boundingY = flst[1]
          self.boundingZ = flst[2]
          self.boundingbox_msg_sem.release()
        elif "ReturnMouseSpeed" in message:
          datastr = message[ message.find("\n") + 1: ]
          lst = datastr.split(",")
          flst = [float(e) for e in lst]
          if flst[0] is not None and not cmath.isnan(flst[0]):
            self.ngl_settings.mouse_sensitivity = flst[0]
          self.mousespeed_msg_sem.release()
        elif "tooltip_id:" in message:
          ttipids = message.split("tooltip_id:")[1]
          hklid = eval(message.split("tooltip_id:")[1])[0]
          sym_id = eval(message.split("tooltip_id:")[1])[1]
          is_friedel_mate = eval(message.split("tooltip_id:")[1])[2]
          rotmx = None
          hkls = self.scene.indices
          if not is_friedel_mate:
            ttip = self.GetTooltipOnTheFly(hklid, sym_id)
          else:
            hklid = hklid % len(hkls)
            ttip = self.GetTooltipOnTheFly(hklid, sym_id, anomalous=True)
          self.AddToBrowserMsgQueue("ShowThisTooltip", ttip)
        elif "doubleclick colour chart" in message:
          self.onDoubleClickColourChart()
        elif "SelectedBrowserDataColumnComboBox" in message:
          sceneid = int(message.split(":")[1])
          self.parent.SetScene(sceneid)
        else:
          if "Ready " in message:
            self.mprint( message, verbose=5)
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
    self.mprint("translation: %s" %str(roundoff(cameratranslation)), verbose=3)
    alllst = roundoff(flst)
    self.mprint("""OrientationMatrix matrix:
  %s,  %s,  %s,  %s
  %s,  %s,  %s,  %s
  %s,  %s,  %s,  %s
  %s,  %s,  %s,  %s
Distance: %s
    """ %tuple(alllst), verbose=4)
    rotdet = ScaleRotMx.determinant()
    if rotdet <= 0.0:
      self.mprint("Negative orientation matrix determinant!!", verbose=1)
      self.SetAutoView() # return old values as a fall back even if they're out of date
      return self.cameraPosZ, self.currentRotmx, self.cameratranslation
    else:
      cameradist = math.pow(rotdet, 1.0/3.0)
    self.mprint("Scale distance: %s" %roundoff(cameradist), verbose=3)
    currentRotmx = matrix.identity(3)
    if cameradist > 0.0:
      currentRotmx = ScaleRotMx/cameradist
      cameraPosZ = cameradist
    return cameraPosZ, currentRotmx, cameratranslation


  def ProcessOrientationMessage(self, message):
    if self.viewerparams.scene_id is None or self.miller_array is None:
      return
    if message.find("NaN")>=0 or message.find("undefined")>=0:
      return
    if "OrientationBeforeReload:" in message:
      if not self.isnewfile:
        self.viewmtrx = message[ message.find("\n") + 1: ]
        self.lastviewmtrx = self.viewmtrx
      self.isnewfile = False
    self.viewmtrx = message[ message.find("\n") + 1: ]
    self.cameraPosZ, self.currentRotmx, self.cameratranslation = self.GetCameraPosRotTrans( self.viewmtrx)
    rotlst = roundoff(self.currentRotmx.elems)
    self.mprint("""Rotation matrix:
  %s,  %s,  %s
  %s,  %s,  %s
  %s,  %s,  %s
    """ %rotlst, verbose=3)
    if "MouseMovedOrientation:" in message:
      self.params.mouse_moved = True
    if self.currentRotmx.is_r3_rotation_matrix():
      # Round off matrix elements to avoid machine imprecision errors that might cast
      # any matrix element into a number strictly larger than 1 which would
      # crash r3_rotation_matrix_as_x_y_z_angles()
      self.currentRotmx = matrix.sqr(roundoff(self.currentRotmx.elems, 9) )
      angles = self.currentRotmx.r3_rotation_matrix_as_x_y_z_angles(deg=True)
      self.mprint("angles: %s" %str(roundoff(angles)), verbose=3)
      z_vec = flex.vec3_double( [(0,0,1)])
      self.rot_zvec = z_vec * self.currentRotmx
      self.mprint("Rotated cartesian Z direction : %s" %str(roundoff(self.rot_zvec[0])), verbose=3)
      rfracmx = matrix.sqr( self.miller_array.unit_cell().reciprocal().fractionalization_matrix() )
      self.rot_recip_zvec = self.rot_zvec * rfracmx
      self.rot_recip_zvec = (1.0/self.rot_recip_zvec.norm()) * self.rot_recip_zvec
      self.mprint("Rotated reciprocal L direction : %s" %str(roundoff(self.rot_recip_zvec[0])), verbose=3)


  def WaitforHandshake(self, sec=5):
    nwait = 0
    while not self.WBmessenger.browserisopen:
      time.sleep(self.sleeptime)
      nwait += self.sleeptime
      if nwait > sec:
        return False
    return True


  def OpenBrowser(self):
    if self.viewerparams.scene_id is not None and not self.WBmessenger.websockclient \
       and not self.WBmessenger.browserisopen or self.isnewfile:
      with open(self.hklfname, "w") as f:
        f.write( self.htmlstr )
      self.url = "file:///" + os.path.abspath( self.hklfname )
      self.url = self.url.replace("\\", "/")
      self.mprint( "Writing %s and connecting to its websocket client" %self.hklfname, verbose=1)
      if self.UseOSBrowser=="default":
        if not webbrowser.open(self.url, new=0):
          self.mprint("Could not open the default web browser")
          return False
      if self.UseOSBrowser != "default" and self.UseOSBrowser != "":
        browserpath = self.UseOSBrowser + " %s"
        if not webbrowser.get(browserpath).open(self.url, new=0):
          self.mprint("Could not open web browser, %s" %self.UseOSBrowser)
          return False
      self.SendInfoToGUI({ "html_url": self.url } )
      self.WBmessenger.browserisopen = True
      #self.isnewfile = False
      return True
    return False


  def set_show_tooltips(self):
    msg = "%s" %self.ngl_settings.show_tooltips
    self.AddToBrowserMsgQueue("DisplayTooltips", msg)


  def set_tooltip_opacity(self):
    msg = "%f" %self.ngl_settings.tooltip_alpha
    self.AddToBrowserMsgQueue("TooltipOpacity", msg)


  def SetOpacities(self, bin_opacities_str):
    retstr = ""
    if self.miller_array and bin_opacities_str:
      self.ngl_settings.bin_opacities = bin_opacities_str
      bin_opacitieslst = eval(self.ngl_settings.bin_opacities)
      for binopacity in bin_opacitieslst:
        alpha = binopacity[0] # float(binopacity.split(",")[0])
        bin = binopacity[1] # int(binopacity.split(",")[1])
        retstr += self.set_opacity(bin, alpha)
      self.SendInfoToGUI( { "bin_opacities": self.ngl_settings.bin_opacities } )
    self.mprint( retstr, verbose=1)


  def set_opacity(self, bin, alpha):
    if bin > self.nbinvalsboundaries-1:
      return "There are only %d bins present\n" %self.nbinvalsboundaries
    msg = "%d, %f" %(bin, alpha)
    self.AddToBrowserMsgQueue("alpha", msg)
    return "Opacity %s set on bin[%s]\n" %(alpha, bin)


  def RedrawNGL(self):
    self.AddToBrowserMsgQueue("Redraw")


  def ReloadNGL(self): # expensive as javascript may be several Mbytes large
    self.mprint("Rendering JavaScript...", verbose=1)
    if not self.OpenBrowser():
      self.AddToBrowserMsgQueue("Reload")


  def JavaScriptCleanUp(self, ):
    self.AddToBrowserMsgQueue("JavaScriptCleanUp")


  def ExpandInBrowser(self):
    if self.sceneisdirty:
      self.mprint( "Not expanding in browser", verbose=1)
      return
    uc = self.miller_array.unit_cell()
    OrtMx = matrix.sqr( uc.orthogonalization_matrix())
    InvMx = OrtMx.inverse()
    msgtype = "Expand"
    msg = ""
    unique_rot_ops = []
    if self.viewerparams.expand_to_p1:
      msgtype += "P1"
      unique_rot_ops = self.symops[ 0 : self.sg.order_p() ] # avoid duplicate rotation matrices
      retmsg = "Expanding to P1 in browser"
      if not self.miller_array.is_unique_set_under_symmetry():
        retmsg += "\nNot all reflections are in the same asymmetric unit in reciprocal space.\n"
        retmsg += "Some reflections might be displayed on top of one another.\n"
      self.mprint( retmsg, verbose=1)
    else:
      unique_rot_ops = [ self.symops[0] ] # No P1 expansion. So only submit the identity matrix
    if self.viewerparams.expand_anomalous and not self.miller_array.anomalous_flag():
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
    self.AddToBrowserMsgQueue(msgtype, msg)
    self.GetBoundingBox() # bounding box changes when the extent of the displayed lattice changes


  def AddVector(self, s1, s2, s3, t1, t2, t3, isreciprocal=True, label="",
                  r=0, g=0, b=0, name="", radius = 0.15):
    """
    Place vector from {s1, s2, s3] to [t1, t2, t3] with colour r,g,b and label
    If name=="" creation is deferred until AddVector is eventually called with name != ""
    These vectors are then joined in the same NGL representation
    """
    uc = self.miller_array.unit_cell()
    vec1 = (s1*self.scene.renderscale, s2*self.scene.renderscale, s3*self.scene.renderscale)
    vec2 = (t1*self.scene.renderscale, t2*self.scene.renderscale, t3*self.scene.renderscale)
    #svec = list(vec)
    if isreciprocal:
      # uc.reciprocal_space_vector() only takes integer miller indices so compute the cartesian coordinates
      # for floating valued miller indices with the transpose of the fractionalization matrix
      vec1 = list( vec1 * matrix.sqr(uc.fractionalization_matrix()).transpose() )
      vec2 = list( vec2 * matrix.sqr(uc.fractionalization_matrix()).transpose() )
      svec1 = [ vec1[0], vec1[1], vec1[2] ]
      svec2 = [ vec2[0], vec2[1], vec2[2] ]
    else:
      vec1 = list( vec1 * matrix.sqr(uc.orthogonalization_matrix()) )
      vec2 = list( vec2 * matrix.sqr(uc.orthogonalization_matrix()) )
      vscale =  1.0/self.scene.renderscale
      # TODO: find suitable scale factor for displaying real space vector together with reciprocal vectors
      svec1 = [ vscale*vec1[0], vscale*vec1[1], vscale*vec1[2] ]
      svec2 = [ vscale*vec2[0], vscale*vec2[1], vscale*vec2[2] ]
    self.mprint("cartesian vector is: %s to %s" %(str(roundoff(svec1)), str(roundoff(svec2))), verbose=2)
    svec = [svec2[0]-svec1[0], svec2[1]-svec1[1], svec2[2]-svec1[2] ]
    xyvec = svec[:] # deep copying
    xyvec[2] = 0.0 # projection vector of svec in the xy plane
    xyvecnorm = math.sqrt( xyvec[0]*xyvec[0] + xyvec[1]*xyvec[1] )
    if xyvecnorm > 0.0:
      angle_x_xyvec = math.acos( xyvec[0]/xyvecnorm )*180.0/math.pi
      angle_y_xyvec = math.acos( xyvec[1]/xyvecnorm )*180.0/math.pi
    else:
      angle_x_xyvec = 90.0
      angle_y_xyvec = 90.0
    yzvec = svec[:]
    yzvec[0] = 0.0 # projection vector of svec in the yz plane
    yzvecnorm = math.sqrt( yzvec[1]*yzvec[1] + yzvec[2]*yzvec[2] )
    if yzvecnorm > 0.0:
      angle_y_yzvec = math.acos( yzvec[1]/yzvecnorm )*180.0/math.pi
      angle_z_yzvec = math.acos( yzvec[2]/yzvecnorm )*180.0/math.pi
    else:
      angle_y_yzvec = 90.0
      angle_z_yzvec = 90.0
    svecnorm = math.sqrt( svec[0]*svec[0] + svec[1]*svec[1] + svec[2]*svec[2] )
    angle_x_svec = math.acos( svec[0]/svecnorm )*180.0/math.pi
    angle_y_svec = math.acos( svec[1]/svecnorm )*180.0/math.pi
    angle_z_svec = math.acos( svec[2]/svecnorm )*180.0/math.pi
    if angle_y_svec > 90.0:
      angle_x_xyvec = -angle_x_xyvec
    self.mprint("angles in xy plane to x,y axis are: %s, %s" %(angle_x_xyvec, angle_y_xyvec), verbose=2)
    self.mprint("angles in yz plane to y,z axis are: %s, %s" %(angle_y_yzvec, angle_z_yzvec), verbose=2)
    self.mprint("angles to x,y,z axis are: %s, %s, %s" %(angle_x_svec, angle_y_svec, angle_z_svec ), verbose=2)
    self.mprint("deferred rendering vector from (%s, %s, %s) to (%s, %s, %s)" %(s1, s2, s3, t1, t2, t3), verbose=2)
    self.AddToBrowserMsgQueue("AddVector", "%s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s" \
         %tuple(svec1 + svec2 + [r, g, b, label, name, radius]) )
    return angle_x_xyvec, angle_z_svec


  def PointVectorPerpendicularToClipPlane(self):
    rotmx = self.Euler2RotMatrix(( self.angle_x_xyvec, self.angle_z_svec, 0.0 ))
    if rotmx.determinant() < 0.99999:
      self.mprint("Rotation matrix determinant is less than 1")
      return rotmx
    self.currentRotmx = rotmx
    self.RotateMxStage(rotmx)
    return rotmx


  def PointVectorParallelToClipPlane(self):
    rotmx = self.Euler2RotMatrix(( self.angle_x_xyvec, self.angle_z_svec+90.0, 90.0 ))
    if rotmx.determinant() < 0.99999:
      self.mprint("Rotation matrix determinant is less than 1")
      return rotmx
    self.currentRotmx = rotmx
    self.RotateMxStage(rotmx)
    return rotmx


  def RotateAroundFracVector(self, phi, r1,r2,r3, prevrotmx = matrix.identity(3), isreciprocal=False, quietbrowser=True):
    if isreciprocal:
    # Assuming vector is in reciprocal space coordinates turn it into cartesian
      cartvec = list( (r1,r2,r3) * matrix.sqr(self.miller_array.unit_cell().fractionalization_matrix()).transpose() )
    else:
      # Assuming vector is in real space fractional coordinates turn it into cartesian
      cartvec = list( (r1,r2,r3) * matrix.sqr(self.miller_array.unit_cell().orthogonalization_matrix()) )
    #  Rodrigues rotation formula for rotation by phi angle around a vector going through origo
    #  See http://mathworld.wolfram.com/RodriguesRotationFormula.html
    # \mathbf I+\left(\sin\,\varphi\right)\mathbf W+\left(2\sin^2\frac{\varphi}{2}\right)\mathbf W^2
    normR = math.sqrt(cartvec[0]*cartvec[0] + cartvec[1]*cartvec[1] + cartvec[2]*cartvec[2] )
    ux = cartvec[0]/normR
    uy = cartvec[1]/normR
    uz = cartvec[2]/normR
    W = matrix.sqr([0, -uz, uy, uz, 0, -ux, -uy, ux, 0])
    I = matrix.identity(3)
    sin2phi2 = math.sin(phi/2)
    sin2phi2 *= sin2phi2
    RotMx = I + math.sin(phi)*W + 2* sin2phi2 * W*W
    self.currentRotmx = RotMx * prevrotmx # impose any other rotation already performed
    self.RotateMxStage(self.currentRotmx, quietbrowser)
    return self.currentRotmx, [ux, uy, uz]


  def SpinAnimate(self, r1, r2, r3):
    self.AddToBrowserMsgQueue("SpinAnimate", "%s, %s, %s" %(r1, r2, r3) )


  def DrawUnitCell(self, scale=1):
    if scale is None:
      self.RemoveVectors("unitcell")
      self.mprint( "Removing real space unit cell", verbose=1)
      return
    uc = self.miller_array.unit_cell()
    rad = 0.2 # scale # * 0.05 #  1000/ uc.volume()
    self.AddVector(0,0,0, scale,0,0, False, label="a", r=0.5, g=0.8, b=0.8, radius=rad)
    self.AddVector(0,0,0, 0,scale,0, False, label="b", r=0.8, g=0.5, b=0.8, radius=rad)
    self.AddVector(0,0,0, 0,0,scale, False, label="c", r=0.8, g=0.8, b=0.5, radius=rad)
    self.AddVector(scale,0,0, scale,scale,0, False, r=0.8, g=0.5, b=0.8, radius=rad)
    self.AddVector(0,scale,0, scale,scale,0, False, r=0.5, g=0.8, b=0.8, radius=rad)
    self.AddVector(0,0,scale, scale,0,scale, False, r=0.5, g=0.8, b=0.8, radius=rad)
    self.AddVector(0,0,scale, 0,scale,scale, False, r=0.8, g=0.5, b=0.8, radius=rad)
    self.AddVector(0,scale,scale, scale,scale,scale, False, r=0.5, g=0.8, b=0.8, radius=rad)
    self.AddVector(scale,0,scale, scale,scale,scale, False, r=0.8, g=0.5, b=0.8, radius=rad)
    self.AddVector(scale,0,0, scale,0,scale, False, r=0.8, g=0.8, b=0.5, radius=rad)
    self.AddVector(0,scale,0, 0,scale,scale, False, r=0.8, g=0.8, b=0.5, radius=rad)
    self.AddVector(scale,scale,0, scale,scale,scale, False, r=0.8, g=0.8, b=0.5, radius=rad, name="unitcell")
    self.mprint( "Adding real space unit cell", verbose=1)



  def DrawReciprocalUnitCell(self, scale=1):
    if scale is None:
      self.RemoveVectors("reciprocal_unitcell")
      self.mprint( "Removing reciprocal unit cell", verbose=1)
      return
    rad = 0.2 # 0.05 * scale
    self.AddVector(0,0,0, scale,0,0, label="a*", r=0.5, g=0.3, b=0.3, radius=rad)
    self.AddVector(0,0,0, 0,scale,0, label="b*", r=0.3, g=0.5, b=0.3, radius=rad)
    self.AddVector(0,0,0, 0,0,scale, label="c*", r=0.3, g=0.3, b=0.5, radius=rad)
    self.AddVector(scale,0,0, scale,scale,0, r=0.3, g=0.5, b=0.3, radius=rad)
    self.AddVector(0,scale,0, scale,scale,0, r=0.5, g=0.3, b=0.3, radius=rad)
    self.AddVector(0,0,scale, scale,0,scale, r=0.5, g=0.3, b=0.3, radius=rad)
    self.AddVector(0,0,scale, 0,scale,scale, r=0.3, g=0.5, b=0.3, radius=rad)
    self.AddVector(0,scale,scale, scale,scale,scale, r=0.5, g=0.3, b=0.3, radius=rad)
    self.AddVector(scale,0,scale, scale,scale,scale, r=0.3, g=0.5, b=0.3, radius=rad)
    self.AddVector(scale,0,0, scale,0,scale, r=0.3, g=0.3, b=0.5, radius=rad)
    self.AddVector(0,scale,0, 0,scale,scale, r=0.3, g=0.3, b=0.5, radius=rad)
    self.AddVector(scale,scale,0, scale,scale,scale, r=0.3, g=0.3, b=0.5, radius=rad, name="reciprocal_unitcell")
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
    self.realspace_scale = self.scene.renderscale * reciprocspan_length / bodydiagonal_length


  def fix_orientation(self):
    if self.ngl_settings.fixorientation:
      self.DisableMouseRotation()
    else:
      self.EnableMouseRotation()


  def clip_plane_vector(self, a, b, c, hkldist=0.0,
             clipwidth=None, fixorientation=True, is_parallel=False, isreciprocal=False):
    # create clip plane oriented parallel or perpendicular to abc vector
    if a==0.0 and b==0.0 and c==0.0 or clipwidth is None:
      self.RemoveVectorsNoClipPlane()
      return
    self.mprint("Applying clip plane to reflections", verbose=1)
    self.RemoveVectors("clip_vector")
    self.angle_x_xyvec, self.angle_z_svec = self.AddVector(0, 0, 0,
                            a, b, c, isreciprocal=isreciprocal, name="clip_vector")
    if fixorientation:
      self.DisableMouseRotation()
    else:
      self.EnableMouseRotation()
    if is_parallel:
      self.vecrotmx = self.PointVectorParallelToClipPlane()
    else:
      self.vecrotmx = self.PointVectorPerpendicularToClipPlane()
    if self.cameraPosZ is None and self.viewmtrx is not None:
      self.cameraPosZ, self.currentRotmx, self.cameratranslation = self.GetCameraPosRotTrans( self.viewmtrx)
    halfdist = self.cameraPosZ  + hkldist # self.viewer.boundingZ*0.5
    if clipwidth == 0.0:
      clipwidth = self.meanradius
    clipNear = halfdist - clipwidth # 50/self.viewer.boundingZ
    clipFar = halfdist + clipwidth  #50/self.viewer.boundingZ
    self.SetClipPlaneDistances(clipNear, clipFar, -self.cameraPosZ)
    #if hkldist < 0.0:
    #  self.TranslateHKLpoints(a, b, c, hkldist)
    scale = max(self.miller_array.index_span().max())/10


  def RemoveVectorsNoClipPlane(self):
    self.EnableMouseRotation()
    self.RemoveVectors()
    self.SetClipPlaneDistances(0, 0)
    self.TranslateHKLpoints(0, 0, 0, 0.0)


  def set_camera_type(self):
    self.AddToBrowserMsgQueue("SetCameraType", self.ngl_settings.camera_type)


  def get_labels_of_data_for_binning(self):
    self.mprint("Data can be binned according to:")
    for i,e in enumerate(self.bin_labels_type_idxs):
      self.mprint("%d, %s" %(i, e[0]))


  def SetFontSize(self, fontsize):
    msg = str(fontsize)
    self.AddToBrowserMsgQueue("SetFontSize", msg)


  def SetBrowserDebug(self, isdebug):
    msg = str(isdebug)
    self.AddToBrowserMsgQueue("SetBrowserDebug", msg)


  def SetMouseSpeed(self, trackspeed):
    msg = str(trackspeed)
    self.AddToBrowserMsgQueue("SetMouseSpeed", msg)
    #self.GetMouseSpeed() # TODO: fix wait time


  def GetMouseSpeed(self):
    self.ngl_settings.mouse_sensitivity = None
    self.mousespeed_msg_sem.acquire()
    self.AddToBrowserMsgQueue("GetMouseSpeed", "")
    if self.WaitforHandshake():
      nwait = 0
      if not self.mousespeed_msg_sem.acquire(blocking=False) and nwait < 5:
        nwait += self.sleeptime
        self.mprint("mousespeed_msg_sem, wait= %s" %nwait, verbose=2)
    self.mousespeed_msg_sem.release()


  def SetClipPlaneDistances(self, near, far, cameraPosZ=None):
    if cameraPosZ is None:
      cameraPosZ = self.cameraPosZ
    msg = str(near) + ", " + str(far) + ", " + str(cameraPosZ)
    self.AddToBrowserMsgQueue("SetClipPlaneDistances", msg)


  def GetClipPlaneDistances(self):
    self.clipNear = None
    self.clipFar = None
    self.cameraPosZ = None
    self.clipplane_msg_sem.acquire()
    self.AddToBrowserMsgQueue("GetClipPlaneDistances", "") # 
    if self.WaitforHandshake():
      nwait = 0
      if not self.clipplane_msg_sem.acquire(blocking=False) and nwait < 5:
        nwait += self.sleeptime
        self.mprint("clipplane_msg_sem, wait= %s" %nwait, verbose=2)
      self.mprint("clipnear, clipfar, cameraPosZ: %s, %s %s" \
                 %(self.clipNear, self.clipFar, self.cameraPosZ), 2)
    self.clipplane_msg_sem.release()
    return (self.clipNear, self.clipFar, self.cameraPosZ)


  def GetBoundingBox(self):
    self.boundingX = 0.0
    self.boundingY = 0.0
    self.boundingZ = 0.0
    self.boundingbox_msg_sem.acquire()
    self.AddToBrowserMsgQueue("GetBoundingBox", "")
    if self.WaitforHandshake():
      nwait = 0
      if not self.boundingbox_msg_sem.acquire(blocking=False) and nwait < 5:
        nwait += self.sleeptime
        self.mprint("boundingbox_msg_sem, wait= %s" %nwait, verbose=2)
      self.mprint("boundingXYZ: %s, %s %s" \
         %(self.boundingX, self.boundingY, self.boundingZ), verbose=2)
    self.boundingbox_msg_sem.release()
    return (self.boundingX, self.boundingY, self.boundingZ)


  def RemoveVectors(self, reprname=""):
    self.AddToBrowserMsgQueue("RemoveVectors", reprname )


  def SetAutoView(self):
    rotmx = self.Euler2RotMatrix( ( 0.0, 0.0, 0.0 ) )
    self.currentRotmx = rotmx
    self.RotateMxStage(rotmx)
    self.AddToBrowserMsgQueue("SetAutoView" )


  def TestNewFunction(self):
    self.AddToBrowserMsgQueue("Testing")


  def MakeImage(self, filename):
    self.imagename = filename
    self.AddToBrowserMsgQueue("MakeImage", "HKLviewer.png,"+ str(sys.version_info[0]) )


  def DisableMouseRotation(self): # disable rotating with the mouse
    self.AddToBrowserMsgQueue("DisableMouseRotation")


  def EnableMouseRotation(self): # enable rotating with the mouse
    self.AddToBrowserMsgQueue("EnableMouseRotation")


  def ReOrientStage(self):
    if self.viewmtrx:
      self.AddToBrowserMsgQueue("SetAutoView", self.viewmtrx)


  def Euler2RotMatrix(self, eulerangles):
    eulerangles1 = eulerangles
    radangles = [e*math.pi/180.0 for e in eulerangles1]
    RotMx = scitbx.math.euler_angles_as_matrix(radangles)
    return RotMx


  def RotateMxStage(self, rotmx, quietbrowser=True):
    if self.cameraPosZ is None:
      return
    scaleRot = rotmx * self.cameraPosZ
    ortrot = scaleRot.as_mat3()
    str_rot = str(ortrot)
    str_rot = str_rot.replace("(", "")
    str_rot = str_rot.replace(")", "")
    msg = str_rot + ", quiet\n"
    if not quietbrowser:
      msg = str_rot + ", verbose\n"
    self.AddToBrowserMsgQueue("RotateStage", msg)


  def TranslateHKLpoints(self, h, k, l, mag):
    # cast this reciprocal vector into cartesian before messaging NGL to translate our HKL points
    #vec = self.miller_array.unit_cell().reciprocal_space_vector((h, k, l))
    hkl_vec = flex.vec3_double( [(h,k,l)])
    rfracmx = matrix.sqr( self.miller_array.unit_cell().reciprocal().orthogonalization_matrix() )
    cartvec = hkl_vec * rfracmx
    if cartvec.norm()==0.0 or mag==0.0:
      svec = (0, 0, 0)
    else:
      #cartvec = (mag/cartvec.norm()) * cartvec
      cartvec = (-mag*self.scene.renderscale/hkl_vec.norm()) * cartvec
      #svec = [cartvec[0][0]*self.scene.renderscale, cartvec[0][1]*self.scene.renderscale, cartvec[0][2]*self.scene.renderscale ]
      svec = cartvec[0]
    #import code, traceback; code.interact(local=locals(), banner="".join( traceback.format_stack(limit=10) ) )
    self.mprint("cartesian translation vector is: " + str(roundoff(svec)), verbose=1)
    str_vec = str(svec)
    str_vec = str_vec.replace("(", "")
    str_vec = str_vec.replace(")", "")
    msg = str_vec + "\n"
    self.AddToBrowserMsgQueue("TranslateHKLpoints", msg)


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
    strdata = ""
    strdata += "%s\n\n" %roundoff(positions, 2)
    strdata += "%s\n\n" %roundoff(colours, 2)
    strdata += "%s\n\n" %roundoff(radii, 2)
    strdata += "%s" %ttipids
    self.AddToBrowserMsgQueue("AddSpheresBin2ShapeBuffer", strdata)


  def RenderStageObjects(self):
    self.AddToBrowserMsgQueue("RenderStageObjects")


  def MakeColourChart(self, ctop, cleft, label, fomlabel, colourgradarray):
    msg = "%s\n\n%s\n\n%s\n\n%s\n\n%s" %(ctop, cleft, label, fomlabel, str(colourgradarray) )
    self.AddToBrowserMsgQueue("MakeColourChart", msg )


  def onDoubleClickColourChart(self):
    # if running the GUI show the colour chart selection dialog
    self.SendInfoToGUI( { "ColourChart": self.viewerparams.color_scheme, 
                          "ColourPowerScale": self.viewerparams.color_powscale } )

  
  def MakeBrowserDataColumnComboBox(self):
    datcolstr =""
    for i,lbl in enumerate(self.hkl_scenes_infos):
      datcolstr = datcolstr + ",".join(lbl[3]) + "\n" + str(i)
      if i < len(self.hkl_scenes_infos)-1:
        datcolstr = datcolstr + "\n\n"

    self.AddToBrowserMsgQueue("MakeBrowserDataColumnComboBox", datcolstr)




ngl_philstr = """
  mouse_sensitivity = 0.04
    .type = float
  bin_opacities = ""
    .type = str
  tooltip_alpha = 0.70
    .type = float
  fontsize = 9
    .type = int
  show_tooltips = none *click hover
    .type = choice
  fixorientation = False
    .type = bool
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
