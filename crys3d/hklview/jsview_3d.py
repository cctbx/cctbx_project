
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
from websocket_server import WebsocketServer
import threading, math, sys, cmath
from time import sleep
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
    if (isinstance(data, flex.int)):
      data = flex.double([e for e in data if e!= display.inanval])
    if millarr.is_complex_array():
      data = flex.abs(millarr.data())
    #data = [e for e in data if not math.isnan(e)]
    data = graphics_utils.NoNansArray( data, data[0] ) # assuming data[0] isn't NaN
    self.maxdata = flex.max( data )
    self.mindata = flex.min( data )
    self.maxsigmas = self.minsigmas = None
    if millarr.sigmas() is not None:
      data = millarr.sigmas()
      #data = [e for e in data if not math.isnan(e)]
      data = graphics_utils.NoNansArray( data, data[0] ) # assuming data[0] isn't NaN
      self.maxsigmas = flex.max( data )
      self.minsigmas = flex.min( data )
    self.minmaxdata = (roundoff(self.mindata), roundoff(self.maxdata))
    self.minmaxsigs = (roundoff(self.minsigmas), roundoff(self.maxsigmas))
    self.labels = self.desc = ""
    #import code, traceback; code.interact(local=locals(), banner="".join( traceback.format_stack(limit=10) ) )
    if millarr.info():
      self.labels = millarr.info().label_string()
      if fomlabel:
        self.labels = millarr.info().label_string() + " + " + fomlabel
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
    self.infotpl = ( self.labels, self.desc, self.spginf, millarr.indices().size(), self.span,
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
      scenearrayinfos.append((ainf.infostr, pidx, fidx, ainf.labels))
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
    self.NGLscriptstr = ""
    self.camera_type = "orthographic"
    self.primitivetype = "SphereBuffer"
    self.url = ""
    self.binscenelabel = "Resolution"
    self.colour_scene_id = None
    self.radii_scene_id = None
    #self.scene_id = None
    #self.rotation_mx = matrix.identity(3)
    self.rot_recip_zvec = None
    self.rot_zvec = None
    self.meanradius = -1
    self.past = time.time()
    self.orientmessage = None
    self.high_quality = True
    if 'high_quality' in kwds:
      self.high_quality = kwds['high_quality']
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
    self.sleeptime = 0.025
    self.colstraliases = ""
    self.binvals = []
    self.binvalsboundaries = []
    self.proc_arrays = []
    self.HKLscenes = []
    self.HKLscenesdict = {}
    self.HKLscenesMaxdata = []
    self.HKLscenesMindata = []
    self.HKLscenesMaxsigmas = []
    self.HKLscenesMinsigmas = []
    self.bindata = None
    self.reciproc_scale = 1.0
    self.realspace_scale = 1.0
    self.sceneisdirty = True
    self.hkl_scenes_info = []
    self.match_valarrays = []
    self.array_infostrs = []
    self.array_infotpls = []
    self.binstrs = []
    self.nuniqueval = 0
    self.bin_infotpls = []
    self.mapcoef_fom_dict = {}
    self.sceneid_from_arrayid = []
    self.was_disconnected = False
    self.parent = None
    if 'parent' in kwds:
      self.parent = kwds['parent']
    self.verbose = 0
    if 'verbose' in kwds:
      self.verbose = eval(kwds['verbose'])
    self.mprint = sys.stdout.write
    if 'mprint' in kwds:
      self.mprint = kwds['mprint']
    self.nbinvalsboundaries = 0
    tempdir = tempfile.gettempdir()
    self.hklfname = os.path.join(tempdir, "hkl.htm" )
    if os.path.isfile(self.hklfname):
      os.remove(self.hklfname)
    if 'htmlfname' in kwds and kwds['htmlfname']:
      self.hklfname = kwds['htmlfname']
    self.hklfname = os.path.abspath( self.hklfname )
    self.jscriptfname = os.path.join(tempdir, "hkljstr.js")
    if os.path.isfile(self.jscriptfname):
      os.remove(self.jscriptfname)
    if 'jscriptfname' in kwds and kwds['jscriptfname'] != "":
      self.jscriptfname = kwds['jscriptfname']
    self.websockport = 7894
    if 'websockport' in kwds:
      self.websockport = kwds['websockport']
    self.send_info_to_gui = None
    if 'send_info_to_gui' in kwds:
      self.send_info_to_gui = kwds['send_info_to_gui']
    self.mprint('Output will be written to \"%s\"\n' \
      'rendering WebGL with JavaScript in \"%s\"' %(self.hklfname, self.jscriptfname))
    self.hklhtml = r"""
<!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 4.0 Transitional//EN">

<html>
<head>
<meta charset="utf-8" />
</head>

<body>
<script src="%s" type="text/javascript"></script>

<script src="%s" type="text/javascript"></script>
    """
    self.htmldiv = """
<div id="viewport" style="width:100%; height:100%;"></div>

</body></html>

    """
    self.colourgradientvalues = []
    self.isinjected = False
    self.UseOSBrowser = True
    if 'UseOSBrowser' in kwds:
      self.UseOSBrowser = eval(kwds['UseOSBrowser'])
    self.viewmtrx = None
    self.lastviewmtrx = None
    self.currentRotmx = matrix.identity(3)
    self.HKLscenesKey = ( 0, False,
                          self.viewerparams.expand_anomalous, self.viewerparams.expand_to_p1  )
    self.msgqueue = []
    self.websockclient = None
    self.handshakewait = 5
    if 'handshakewait' in kwds:
      self.handshakewait = eval(kwds['handshakewait'])
    self.lastmsg = "" # "Ready"
    self.browserisopen = False
    self.msgdelim = ":\n"
    self.msgqueuethrd = None
    self.StartWebsocket()
    self.isterminating = False
    self.javascriptcleaned = False


  def __exit__(self, exc_type, exc_value, traceback):
    # not called unless instantiated with a "with hklview_3d ... " statement
    self.JavaScriptCleanUp()
    nwait = 0
    while not self.isterminating and nwait < 5:
      sleep(self.sleeptime)
      nwait += self.sleeptime
    if os.path.isfile(self.hklfname):
      os.remove(self.hklfname)
    self.mprint("Destroying hklview_3d", 1)


  def SendInfoToGUI(self, mydict):
    if self.send_info_to_gui:
      self.send_info_to_gui( mydict )


  def update_settings(self, diff_phil, curphilparam) :
    self.ngl_settings = curphilparam.viewer.NGL
    self.viewerparams = curphilparam.viewer
    self.params = curphilparam
    self.diff_phil = diff_phil

    if has_phil_path(diff_phil,
                     "openfilename",
                     "spacegroup_choice",
                     "using_space_subgroup",
                     "merge_data",
                     "camera_type",
                     "miller_array_operations",
                     ) \
     or has_phil_path(diff_phil, "viewer") \
     and has_phil_path(diff_phil,
                       "show_data_over_sigma",
                       "show_missing",
                       "show_only_missing",
                       "show_systematic_absences",
                       "slice_axis",
                       "slice_mode",
                       "slice_index",
                       "scene_id",
                       "scale",
                       "nth_power_scale_radii"
                       ) \
     or self.viewerparams.inbrowser==False and \
      ( has_phil_path(diff_phil,
                     "expand_anomalous",
                     "expand_to_p1",
                     "show_anomalous_pairs")
       ):
        if curphilparam.viewer.slice_mode and self.viewerparams.inbrowser:
          self.viewerparams.inbrowser = False
        self.sceneisdirty = True
        self.ConstructReciprocalSpace(curphilparam, merge=self.merge )
    msg = ""
    if self.viewerparams.scene_id is not None and \
      ( has_phil_path(diff_phil,
                      "show_missing",
                      "show_only_missing",
                      "show_systematic_absences",
                      "scene_bin_thresholds",
                      "bin_scene_label",
                      "nbins"
                      )
       ):
      self.binvals, self.nuniqueval = self.calc_bin_thresholds(curphilparam.bin_scene_label, curphilparam.nbins)
      self.sceneisdirty = True

    if has_phil_path(diff_phil, "camera_type"):
      self.set_camera_type()

    if has_phil_path(diff_phil, "miller_array_operations"):
      self.viewerparams.scene_id = len(self.HKLscenes)-1
      self.set_scene(self.viewerparams.scene_id)

    if self.viewerparams.scene_id is not None:
      if not self.isinjected:
        self.scene = self.HKLscenes[self.viewerparams.scene_id]
      self.DrawNGLJavaScript()
      msg = "Rendered %d reflections\n" % self.scene.points.size()
      #if not has_phil_path(diff_phil, "scene_id"):
# set_volatile_params() is already called when we receive the AutoViewSet message when loading a new scene
      msg += self.set_volatile_params()
    return msg, curphilparam


  def set_volatile_params(self):
    msg = ""
    if self.viewerparams.scene_id is not None:
      if has_phil_path(self.diff_phil, "angle_around_vector"): # no need to redraw any clip plane
        return msg
      self.fix_orientation(self.viewerparams.NGL.fixorientation)
      self.SetMouseSpeed(self.viewerparams.NGL.mouse_sensitivity)
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
        #self.clip_plane_hkl_vector(hkl[0], hkl[1], hkl[2], clipwidth=200,
        #                 fixorientation = self.viewerparams.NGL.fixorientation)
      if self.viewerparams.inbrowser and not self.viewerparams.slice_mode:
        msg += self.ExpandInBrowser(P1= self.viewerparams.expand_to_p1,
                              friedel_mate= self.viewerparams.expand_anomalous)
      if self.params.clip_plane.clipwidth:
        clipwidth = self.params.clip_plane.clipwidth
        hkldist = self.params.clip_plane.hkldist
        R = flex.vec3_double( [(self.params.clip_plane.h, self.params.clip_plane.k, self.params.clip_plane.l)])
        if self.params.clip_plane.fractional_vector == "realspace" or self.params.clip_plane.fractional_vector == "tncs":
          isreciprocal = False

      self.clip_plane_vector(R[0][0], R[0][1], R[0][2], hkldist,
          clipwidth, self.viewerparams.NGL.fixorientation, self.params.clip_plane.is_parallel,
          isreciprocal)
      msg += self.SetOpacities(self.viewerparams.NGL.bin_opacities )
      if self.params.real_space_unit_cell_scale_fraction is None:
        scale = None
      else:
        scale = (self.realspace_scale - 1.0)*self.params.real_space_unit_cell_scale_fraction + 1.0
      msg += self.DrawUnitCell(scale )
      if self.params.reciprocal_unit_cell_scale_fraction is None:
        scale = None
      else:
        scale = (self.reciproc_scale - 1.0)*self.params.reciprocal_unit_cell_scale_fraction + 1.0
      msg += self.DrawReciprocalUnitCell(scale )
      self.set_tooltip_opacity()
      self.set_show_tooltips()
      #self.SetAutoView()
    return msg


  def set_scene(self, scene_id):
    self.binvals = []
    self.isinjected = False
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
      self.isinjected = False
    if self.viewerparams and self.viewerparams.scene_id is not None and self.viewerparams.scene_id >= 0 and self.HKLscenes:
      self.miller_array = self.HKLscenes[self.viewerparams.scene_id].miller_array
      self.scene = self.HKLscenes[self.viewerparams.scene_id]
    self.merge = merge
    if (self.miller_array is None):
      return
    self.identify_suitable_fomsarrays()
    self.GetUnitcellScales()
    self.d_min = self.miller_array.d_min()
    array_info = self.miller_array.info()
    self.sg = self.miller_array.space_group()
    self.symops = self.sg.all_ops()
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
      datval = None
      if hkl in hklscene.work_array.indices():
        datval = hklscene.work_array.data_at_first_index(hkl)
      else:
        if id >= hklscene.data.size():
          continue
        datval = hklscene.data[id]
      if datval is not None and (not (math.isnan( abs(datval) ) or datval == display.inanval)):
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
        elif hklscene.work_array.sigmas() is not None:
          sigma = hklscene.work_array.sigma_at_first_index(hkl)
          spbufttip += str(roundoff(datval, 2)) + ", " + str(roundoff(sigma, 2))
        else:
          spbufttip += str(roundoff(datval, 2))
    spbufttip += '\\n\\n%d,%d,%d' %(id, sym_id, anomalous) # compared by the javascript
    spbufttip += '\''
    return spbufttip


  def get_col_fomcol(self, idx):
    if len(self.hkl_scenes_info) == 0:
      return -1, -1
    return self.hkl_scenes_info[idx][6], self.hkl_scenes_info[idx][7]


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


  def ConstructReciprocalSpace(self, curphilparam, merge=None, scene_id=None):
    self.HKLscenesKey = (curphilparam.spacegroup_choice,
                         curphilparam.using_space_subgroup,
                         curphilparam.merge_data,
                         self.viewerparams.expand_anomalous,
                         self.viewerparams.expand_to_p1,
                         self.viewerparams.inbrowser,
                         self.viewerparams.slice_axis,
                         self.viewerparams.slice_mode,
                         self.viewerparams.slice_index,
                         self.viewerparams.show_missing,
                         self.viewerparams.show_only_missing,
                         self.viewerparams.show_systematic_absences,
                         self.viewerparams.scale,
                         self.viewerparams.nth_power_scale_radii
                         )
    if self.HKLscenesKey in self.HKLscenesdict and not self.has_new_miller_array:
      (
        self.HKLscenes,
        self.tooltipstringsdict,
        self.HKLscenesMaxdata,
        self.HKLscenesMindata,
        self.HKLscenesMaxsigmas,
        self.HKLscenesMinsigmas,
        self.hkl_scenes_info
      ) =  self.HKLscenesdict[self.HKLscenesKey]
      self.mprint("Scene key is already present", verbose=1)
      #self.sceneisdirty = False
      return True
    self.mprint("Constructing HKL scenes", verbose=0)

    HKLscenes = []
    HKLscenesMaxdata = []
    HKLscenesMindata = []
    HKLscenesMaxsigmas = []
    HKLscenesMinsigmas = []
    hkl_scenes_info = []
    tooltipstringsdict = {}
    i = 0
    assert(self.proc_arrays)
    self.mprint("\nReflection data scenes:", verbose=0)
    arrayid = None
    if scene_id is not None:
      arrayid = self.scene_id_to_array_id(scene_id)
    for (idx,e) in enumerate(self.proc_arrays):
      if arrayid is not None and arrayid != idx:
        continue
      (hklscenes, scenemaxdata,
        scenemindata, scenemaxsigmas,
         sceneminsigmas, scenearrayinfos
         ) = MakeHKLscene( e.deep_copy(), idx, copy.deepcopy(self.viewerparams), self.mapcoef_fom_dict, merge, self.mprint )

      HKLscenesMaxdata.extend(scenemaxdata)
      HKLscenesMindata.extend(scenemindata)
      HKLscenesMaxsigmas.extend(scenemaxsigmas)
      HKLscenesMinsigmas.extend(sceneminsigmas)
      hkl_scenes_info.extend(scenearrayinfos)
      HKLscenes.extend(hklscenes)
      for i,inf in enumerate(scenearrayinfos):
        self.mprint("%d, %s" %(idx+i+1, inf[0]), verbose=0)

    tooltipstringsdict = {}
    self.colstraliases = "var hk = \'H,K,L: \';\n"
    self.HKLscenesdict[self.HKLscenesKey] = (
                HKLscenes,
                tooltipstringsdict,
                HKLscenesMaxdata,
                HKLscenesMindata,
                HKLscenesMaxsigmas,
                HKLscenesMinsigmas,
                hkl_scenes_info
                )
    (
      self.HKLscenes,
      self.tooltipstringsdict,
      self.HKLscenesMaxdata,
      self.HKLscenesMindata,
      self.HKLscenesMaxsigmas,
      self.HKLscenesMinsigmas,
      self.hkl_scenes_info
    ) =  self.HKLscenesdict[self.HKLscenesKey]
    #import code, traceback; code.interact(local=locals(), banner="".join( traceback.format_stack(limit=10) ) )
    self.sceneisdirty = True
    self.SendInfoToGUI({ "hklscenes_arrays": self.hkl_scenes_info, "NewHKLscenes" : True })
    self.has_new_miller_array = False
    return True


  def identify_suitable_fomsarrays(self):
    self.mprint("Matching complex arrays to suitable FOM arrays")
    self.mapcoef_fom_dict = {}
    self.sceneid_from_arrayid = []
    for k,proc_array in enumerate(self.proc_arrays):
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


  def calc_bin_thresholds(self, bin_scene_label, nbins):
    self.binscenelabel = bin_scene_label
    if self.binscenelabel=="Resolution":
      warray = self.HKLscenes[int(self.viewerparams.scene_id)].work_array
      dres = self.HKLscenes[int(self.viewerparams.scene_id)].dres
      uc = warray.unit_cell()
      indices = self.HKLscenes[int(self.viewerparams.scene_id)].indices
      binning = miller.binning( uc, nbins, indices, flex.max(dres), flex.min(dres) )
      binvals = [ binning.bin_d_range(n)[0] for n in binning.range_all() ]
      binvals = [ e for e in binvals if e != -1.0] # delete dummy limit
      binvals = list( 1.0/flex.double(binvals) )
      nuniquevalues = len(set(list(dres)))
    else:
      bindata = self.HKLscenes[int(self.binscenelabel)].data.deep_copy()
      if isinstance(bindata, flex.complex_double):
        raise Sorry("Cannot order complex data values for binning.")
      selection = flex.sort_permutation( bindata )
      bindata_sorted = bindata.select(selection)
      # get binvals by dividing bindata_sorted with nbins
      binvals = [bindata_sorted[0]] * nbins #
      for i,e in enumerate(bindata_sorted):
        idiv = int(nbins*float(i)/len(bindata_sorted))
        binvals[idiv] = e
      nuniquevalues = len(set(list(bindata)))
    binvals.sort()
    return binvals, nuniquevalues


  def UpdateBinValues(self, binvals = [], nuniquevalues = 0):
    if binvals:
      binvals.sort()
      self.binvals = binvals
    else: # ensure default resolution interval includes all data by avoiding rounding errors
      self.binvals = [ 1.0/(self.miller_array.d_max_min()[0]*1.001),
                       1.0/(self.miller_array.d_max_min()[1]*0.999) ]
    self.nuniqueval = nuniquevalues


  def MatchBinArrayToSceneArray(self, ibinarray):
    # match bindata with data(scene_id)
    if self.binscenelabel=="Resolution":
      return 1.0/self.scene.dres
    # get the array id that is mapped through an HKLscene id
    binarraydata = self.HKLscenes[ibinarray].data
    scenearraydata = self.HKLscenes[self.viewerparams.scene_id].data
    matchindices = miller.match_indices(self.HKLscenes[self.viewerparams.scene_id].indices, self.HKLscenes[ibinarray].indices )
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
    data = millarr.data()
    sigmas = millarr.sigmas()
    dres = millarr.unit_cell().d( millarr.indices() )
    newarray = millarr.deep_copy()
    self.mprint("Creating new miller array through the operation: %s" %operation)
    try:
      newdata = None
      newsigmas = None
      exec(operation)
      newarray._data = newdata
      newarray._sigmas = newsigmas
      return newarray
    except Exception as e:
      self.mprint( str(e), verbose=0)
      return None


  def OperateOn2MillerArrays(self, millarr1, millarr2, operation):
    # lets user specify a one line python expression operating on data1 and data2
    matchindices = miller.match_indices(millarr1.indices(), millarr2.indices() )
    matcharr1 = millarr1.select( matchindices.pairs().column(0) )
    matcharr2 = millarr2.select( matchindices.pairs().column(1) )
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
      exec(operation)
      newarray._data = newdata
      newarray._sigmas = newsigmas
      return newarray
    except Exception as e:
      self.mprint( str(e), verbose=0)
      return None


  def DrawNGLJavaScript(self, blankscene=False):
    if not self.scene or not self.sceneisdirty:
      return
    if self.miller_array is None :
      self.mprint( "Select an HKL scene to display reflections" )
      return
    self.mprint("Composing JavaScript...")

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
    # make arrow font size roughly proportional to radius of highest resolution shell
    #fontsize = str(1.0 + roundoff(math.pow( max(self.miller_array.index_span().max()), 1.0/3.0)))
    if not self.miller_array:
      fontsize = str(1.0)
    else:
      fontsize = str(1.0 + roundoff(math.pow( max(self.miller_array.index_span().max()), 1.0/2.0)))

    if blankscene:
      axisfuncstr = "\nvar MakeHKL_Axis = function() { };\n"
    else:
      axisfuncstr = """
var fontsize = %s;
function MakeHKL_Axis(mshape)
{
  // xyz arrows
  // mshape.addSphere( [0,0,0] , [ 1, 1, 1 ], 0.3, 'Origin');
  //blue-x
  mshape.addArrow( %s, %s , [ 0, 0, 1 ], 0.1);
  //green-y
  mshape.addArrow( %s, %s , [ 0, 1, 0 ], 0.1);
  //red-z
  mshape.addArrow( %s, %s , [ 1, 0, 0 ], 0.1);

  mshape.addText( %s, [ 0, 0, 1 ], fontsize, 'h');
  mshape.addText( %s, [ 0, 1, 0 ], fontsize, 'k');
  mshape.addText( %s, [ 1, 0, 0 ], fontsize, 'l');
};
    """ %(fontsize, str(Hstararrowstart), str(Hstararrowend), str(Kstararrowstart),
          str(Kstararrowend), str(Lstararrowstart), str(Lstararrowend), Hstararrowtxt,
          Kstararrowtxt, Lstararrowtxt)

    # Make colour gradient array used for drawing a bar of colours next to associated values on the rendered html
    mincolourscalar = self.HKLscenesMindata[self.colour_scene_id]
    maxcolourscalar = self.HKLscenesMaxdata[self.colour_scene_id]
    if self.viewerparams.sigma_color:
      mincolourscalar = self.HKLscenesMinsigmas[self.colour_scene_id]
      maxcolourscalar = self.HKLscenesMaxsigmas[self.colour_scene_id]
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
    if self.HKLscenes[self.colour_scene_id].miller_array.is_complex_array():
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
      if self.HKLscenes[self.colour_scene_id].isUsingFOMs():
        fomln = 50
        fom = 1.0
        fomdecr = 1.0/(fomln-1.0)
      # make fomln fom arrays of size len(colourscalararray) when calling colour_by_phi_FOM
        for j in range(fomln):
          fomarrays.append( flex.double(len(colourscalararray), fom) )
          fom -= fomdecr
        for j in range(fomln):
          colourgradarrays.append( graphics_utils.colour_by_phi_FOM( colourscalararray*(math.pi/180.0), fomarrays[j] ) * 255.0)
      else:
        fomln =1
        fomarrays = [1.0]
        colourgradarrays.append( graphics_utils.colour_by_phi_FOM( colourscalararray*(math.pi/180.0) ) * 255.0)
    else:
      fomln = 1
      fomarrays = [1.0]
      colourgradarrays.append(graphics_utils.color_by_property(
        properties= flex.double(colourscalararray),
        selection=flex.bool( len(colourscalararray), True),
        color_all=False,
        gradient_type= self.viewerparams.color_scheme) * 255.0)

    colors = self.HKLscenes[self.colour_scene_id].colors
    radii = self.HKLscenes[self.radii_scene_id].radii
    self.meanradius = flex.mean(radii)

    if blankscene:
      points = flex.vec3_double( [ ] )
      colors = flex.vec3_double( [ ] )
      radii = flex.double( [ ] )
      self.binscenelabel = "Resolution"
    else:
      points = self.scene.points

    nrefls = points.size()
    hkls = self.scene.indices
    dres = self.scene.dres
    if self.binscenelabel=="Resolution":
      colstr = "dres"
    else:
      colstr = self.HKLscenes[ int(self.binscenelabel) ].work_array.info().label_string()
    data = self.scene.data
    colourlabel = self.HKLscenes[self.colour_scene_id].colourlabel
    fomlabel = self.HKLscenes[self.colour_scene_id].fomlabel
    #import code, traceback; code.interact(local=locals(), banner="".join( traceback.format_stack(limit=10) ) )
    assert (colors.size() == radii.size() == nrefls)
    colours = []
    positions = []
    radii2 = []
    spbufttips = []

    self.binvalsboundaries = []
    if self.binscenelabel=="Resolution":
      self.binvalsboundaries = self.binvals
      self.bindata = 1.0/self.scene.dres
    else:
      ibinarray= int(self.binscenelabel)
      self.binvalsboundaries = [ self.HKLscenesMindata[ibinarray] - 0.1 , self.HKLscenesMaxdata[ibinarray] + 0.1 ]
      self.binvalsboundaries.extend( self.binvals )
      self.binvalsboundaries.sort()
      if self.binvalsboundaries[0] < 0.0:
        self.binvalsboundaries.append(0.0)
        self.binvalsboundaries.sort()
      #self.bindata = self.HKLscenes[ibinarray].data
      self.bindata = self.MatchBinArrayToSceneArray(ibinarray)
      if self.HKLscenes[ibinarray].work_array.is_complex_array():
        self.bindata = self.HKLscenes[ibinarray].ampl

    self.nbinvalsboundaries = len(self.binvalsboundaries)
    # Un-binnable data is scene data values where there's no matching reflection in the bin data
    # Put these in a separate bin and be diligent with the book keeping!
    for ibin in range(self.nbinvalsboundaries+1): # adding the extra bin for un-binnable data
      colours.append([]) # colours and positions are 3 x size of data()
      positions.append([])
      radii2.append([])
      spbufttips.append([])

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

    flexbinvalsboundaries = flex.double(self.binvalsboundaries)
    start_time = time.time()
    for i, hklstars in enumerate(points):
      # bin currently displayed data according to the values of another miller array
      ibin = data2bin( self.bindata[i], self.binvalsboundaries, self.nbinvalsboundaries )
      positions[ibin].extend( graphics_utils.flt_roundoffvec3(hklstars, 2) )
      colours[ibin].extend( graphics_utils.flt_roundoffvec3(colors[i], 2) )
      radii2[ibin].append( graphics_utils.flt_roundoff(radii[i], 2) )
      spbufttips[ibin].append( i )

    elapsed_time = time.time() - start_time
    self.mprint("elapsed time: %s" %elapsed_time, verbose=2)

    spherebufferstr = self.colstraliases
    negativeradiistr = ""
    cntbin = 0
    self.binstrs = []
    self.bin_infotpls = []
    if self.nuniqueval < self.params.nbins:
      self.mprint("%d bins was requested but %s data has only %d unique value(s)!" %(self.params.nbins, colstr, self.nuniqueval), 0)
    for ibin in range(self.nbinvalsboundaries+1):
      mstr =""
      nreflsinbin = len(radii2[ibin])
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

      spherebufferstr += "\n// %s\n" %mstr
      #spherebufferstr += "  ttips.push( [ ] );"
      ttlst = [-1]
      ttlst.extend(spbufttips[ibin])
      ttipsobj = "{ ids: " + str( ttlst ) + """,
       getPosition: function() { return { x:0, y:0 }; } // dummy function to avoid crash
  } """
      # + ",\n      getPosition: function() { return stage.mouseObserver.canvasPosition; } }"
      spherebufferstr += "  ttips.push( %s );" %str( ttipsobj )
      spherebufferstr += """
  positions.push( new Float32Array( %s ) );
  colours.push( new Float32Array( %s ) );
  radii.push( new Float32Array( %s ) );
  shapebufs.push( new NGL.%s({
    position: positions[%d],
    color: colours[%d], """ %(str(positions[ibin]), str(colours[ibin]), \
        str(radii2[ibin]), self.primitivetype, cntbin, \
        cntbin)
      if self.primitivetype == "SphereBuffer":
        spherebufferstr += "\n    radius: radii[%d]," %cntbin
      spherebufferstr += "\n    picking: ttips[%d]," %cntbin
      if self.primitivetype == "PointBuffer":
        spherebufferstr += "\n  }, {pointSize: %1.2f})\n" %self.viewerparams.scale
      else:
        if self.high_quality:
          spherebufferstr += """
    })
  );
  """
        else:
          spherebufferstr += """
    }, { disableImpostor: true
   ,    sphereDetail: 0 }) // rather than default value of 2 icosahedral subdivisions
  );
  """
      spherebufferstr += "shape.addBuffer(shapebufs[%d]);\n  alphas.push(1.0);\n" %cntbin

      if ibin <self.nbinvalsboundaries and self.binvalsboundaries[ibin] < 0.0:
        negativeradiistr += "shapebufs[%d].setParameters({metalness: 1});\n" %cntbin
      cntbin += 1

    #self.ngl_settings.bin_opacities = str([ "1.0, %d"%e for e in range(cntbin) ])
    self.ngl_settings.bin_opacities = str([ (1.0, e) for e in range(cntbin) ])
    self.SendInfoToGUI( { "bin_opacities": self.ngl_settings.bin_opacities,
                          "bin_infotpls": self.bin_infotpls,
                          "bin_data_label": colstr,
                          "tooltip_opacity": self.ngl_settings.tooltip_alpha
                         } )

    spherebufferstr += """
// create tooltip element and add to the viewer canvas
  stage.viewer.container.appendChild(tooltip);


  stage.signals.clicked.add(
    PickingProxyfunc
  );


  stage.mouseObserver.signals.dragged.add(
    function ( deltaX, deltaY)
    {
      if (clipFixToCamPosZ === true)
      {
        stage.viewer.parameters.clipNear = origclipnear + (origcameraZpos -stage.viewer.camera.position.z);
        stage.viewer.parameters.clipFar = origclipfar + (origcameraZpos -stage.viewer.camera.position.z);
        stage.viewer.requestRender();
      }
      msg = getOrientMsg();
      rightnow = timefunc();
      if (rightnow - timenow > 250)
      { // only post every 250 milli second as not to overwhelm python
        postrotmxflag = true;
        WebsockSendMsg('CurrentViewOrientation:\\n' + msg );
        timenow = timefunc();
      }
    }
  );


  stage.mouseObserver.signals.clicked.add(
    function (x, y)
    {
      msg = getOrientMsg();
      WebsockSendMsg('CurrentViewOrientation:\\n' + msg );
    }
  );


  stage.mouseObserver.signals.scrolled.add(
    function (delta)
    {
      if (clipFixToCamPosZ === true)
      {
        stage.viewer.parameters.clipNear = origclipnear + (origcameraZpos -stage.viewer.camera.position.z);
        stage.viewer.parameters.clipFar = origclipfar + (origcameraZpos -stage.viewer.camera.position.z);
        stage.viewer.requestRender();
      }
      msg = getOrientMsg();
      rightnow = timefunc();
      if (rightnow - timenow > 250)
      { // only post every 250 milli second as not to overwhelm python
        WebsockSendMsg('CurrentViewOrientation:\\n' + msg );
        timenow = timefunc();
      }
    }
  );


  stage.viewer.signals.rendered.add(
    function()
    {
      if (postrotmxflag === true)
      {
        msg = getOrientMsg();
        WebsockSendMsg('CurrentViewOrientation:\\n' + msg );
        postrotmxflag = false;
      }
    }
  );


  stage.viewerControls.signals.changed.add(
    function()
    {
      msg = getOrientMsg();
      rightnow = timefunc();
      if (rightnow - timenow > 250)
      { // only post every 250 milli second as not to overwhelm python
        WebsockSendMsg('CurrentViewOrientation:\\n' + msg );
        //ReturnClipPlaneDistances();
        sleep(250).then(()=> {
            ReturnClipPlaneDistances();
          }
        );
        timenow = timefunc();
      }
    }
  );

    """

    colourgradstrs = "colourgradvalarray = new Array(%s);\n" %fomln
    # if displaying phases from map coefficients together with fom values then
    for g,colourgradarray in enumerate(colourgradarrays):
      self.colourgradientvalues = []
      for j,e in enumerate(colourgradarray):
        self.colourgradientvalues.append( [colourscalararray[j], e] )
      self.colourgradientvalues = roundoff( self.colourgradientvalues )
      fom = fomarrays[g]
      colourgradstr = []
      for j,val in enumerate(self.colourgradientvalues):
        vstr = ""
        alpha = 1.0
        rgb = (int(val[1][0]), int(val[1][1]), int(val[1][2]) )
        gradval = "rgba(%s, %s, %s, %s)" %(rgb[0], rgb[1], rgb[2], alpha)
        if j%10 == 0 or j==len(self.colourgradientvalues)-1 :
          vstr = str( roundoff(val[0], 2, as_string=True) )
        colourgradstr.append([vstr , gradval])
      colourgradstrs += "  colourgradvalarray[%s] = %s;\n" %(g, str(colourgradstr) )
    if blankscene:
      colourscriptstr = ""
    else:
      colourscriptstr = """

  //colourgradvalarrays
  %s

  ColourChart("%s", "%s");

    """ % (colourgradstrs, colourlabel, fomlabel)


    #negativeradiistr = ""
    #for ibin in range(self.nbinvalsboundaries):
    #  if self.binvalsboundaries[ibin] < 0.0:
    #    negativeradiistr += "shapebufs[%d].setParameters({metalness: 1})\n" %ibin
    qualitystr = """ , { disableImpostor: true
                  , sphereDetail: 0 } // rather than default value of 2 icosahedral subdivisions
            """
    if self.high_quality:
      qualitystr = ""

    self.NGLscriptstr = """


function createElement(name, properties, style)
{
// utility function used in for loop over colourgradvalarray
  var el = document.createElement(name);
  Object.assign(el, properties);
  Object.assign(el.style, style);
  Object.assign(el.style,
  {
      display: "block",
      position: "absolute",
      fontFamily: "sans-serif",
      //fontSize: "smaller",
      fontSize: "12px",
  }
  );
  return el;
}


function addElement(el)
{
// utility function used in for loop over colourgradvalarray
  Object.assign(el.style,
  {
    position: "absolute",
    zIndex: 10
  }
  );
  stage.viewer.container.appendChild(el);
}


function addDivBox(txt, t, l, w, h, bgcolour="rgba(255, 255, 255, 0.0)")
{
  divbox = createElement("div",
  {
    innerText: txt
  },
  {
    backgroundColor: bgcolour,
    color:  "rgba(0, 0, 0, 1.0)",
    top: t.toString() + "px",
    left: l.toString() + "px",
    width: w.toString() + "px",
    height: h.toString() + "px",
  }
  );
  addElement(divbox);
}

var dbgmsg = "";
// debug message window
var debugmessage = document.createElement("div");
Object.assign(debugmessage.style, {
  position: "absolute",
  zIndex: 10,
  pointerEvents: "none",
  backgroundColor: "rgba(255, 255, 255, 0.8 )",
  color: "black",
  padding: "0.1em",
  fontFamily: "sans-serif",
  bottom: "10px",
  left: "10px",
  fontSize: "smaller",
  display: "block"
});



// Microsoft Edge users follow instructions on
// https://stackoverflow.com/questions/31772564/websocket-to-localhost-not-working-on-microsoft-edge
// to enable websocket connection

var pagename = location.pathname.substring(1);
var mysocket;

try
{
  mysocket = new WebSocket('ws://127.0.0.1:%s/');
}
catch(err)
{
  alert('JavaScriptError: ' + err.stack );
  addDivBox("Error!", window.innerHeight - 50, 20, 40, 20, rgba(100, 100, 100, 0.0));
}



function WebsockSendMsg(msg)
{
  try
  {
    // Avoid "WebSocket is already in CLOSING or CLOSED state" errors when using QWebEngineView
    // See https://stackoverflow.com/questions/48472977/how-to-catch-and-deal-with-websocket-is-already-in-closing-or-closed-state-in
    if (mysocket.readyState === mysocket.OPEN)
    {
      mysocket.send(msg);
      mysocket.send( 'Ready ' + pagename + '\\n' );
    }
  }
  catch(err)
  {
    alert('JavaScriptError: ' + err.stack );
    addDivBox("Error!", window.innerHeight - 50, 20, 40, 20, rgba(100, 100, 100, 0.0));
  }
}


mysocket.onopen = function(e)
{
  msg = '%s now connected via websocket to ' + pagename + '\\n';
  WebsockSendMsg(msg);
  dbgmsg =msg;
  rerendered = false;
  //ReRender();
};


mysocket.onclose = function(e)
{
  msg = '%s now disconnecting from websocket ' + pagename + '\\n';
  WebsockSendMsg(msg);
  dbgmsg =msg;
};


function sleep(ms) {
  return new Promise(resolve => setTimeout(resolve, ms));
}

async function ReRender()
{
  await sleep(500);
  if (shapeComp != null && rerendered==false) // workaround for QTWebEngine bug sometimes failing to render scene
  {
    shapeComp.autoView();
    rerendered = true;
    WebsockSendMsg( 'AutoViewSet ' + pagename );
  }
}


async function RenderRequest()
{
  await sleep(100);
  stage.viewer.requestRender();
  WebsockSendMsg( 'RenderRequest ' + pagename );
}

// Log errors to debugger of your browser
mysocket.onerror = function(error)
{
  msg = 'WebSocket Error ' + error;
  console.log(msg);
  dbgmsg =msg;
};


var stage;
var shape;
var shapeComp;
var vectorshape = null;
var repr;
var AA = String.fromCharCode(197); // short for angstrom
var DGR = String.fromCharCode(176); // short for degree symbol
var current_ttip = "";
var ttips = [];
var vectorreprs = [];
var vectorshapeComps = [];
var positions = [];
var br_positions = [];
var br_colours = [];
var br_radii = [];
var br_ttips = [];
var colours = [];
var alphas = [];
var radii = [];
var shapebufs = [];
var br_shapebufs = [];
var nrots = 0;
var postrotmxflag = false;
var cvorient = new NGL.Matrix4();
var oldmsg = "";
var clipFixToCamPosZ = false;
var origclipnear;
var origclipfar;
var origcameraZpos;
var nbins = %s;
var rerendered = false;
var expstate = "";
var current_ttip_ids;
var isdebug = %s;
var tdelay = 100;
var displaytooltips = true;

function timefunc() {
  var d = new Date();
  var now = d.getTime();
  return now
}

var timenow = timefunc();
var rightnow = timefunc();


window.addEventListener( 'resize',
  function( event ){
      stage.handleResize();
  },
  false
);


if (isdebug)
{
  var script=document.createElement('script');
  script.src='https://rawgit.com/paulirish/memory-stats.js/master/bookmarklet.js';
  document.head.appendChild(script);
}


// define tooltip element
var tooltip = document.createElement("div");
Object.assign(tooltip.style, {
  display: "none",
  position: "absolute",
  zIndex: 10,
  pointerEvents: "none",
  backgroundColor: "rgba(255, 255, 255, %s )",
  color: "black",
  padding: "0.1em",
  fontFamily: "sans-serif"
});


%s


function getOrientMsg()
{
  cvorientmx = stage.viewerControls.getOrientation();
  if (cvorientmx.determinant() == 0)
      return oldmsg; // don't return invalid matrix

  cvorient = cvorientmx.elements;
  for (j=0; j<16; j++)
  {
    if (Number.isNaN( cvorient[j]) )
      return oldmsg; // don't return invalid matrix
  }

  if (stage.viewer.cDist != 0
        && stage.viewer.parameters.clipFar > stage.viewer.cDist
        && stage.viewer.cDist > stage.viewer.parameters.clipNear)
    cameradist = stage.viewer.cDist;
  else if (stage.viewer.camera.position.z != 0
        && stage.viewer.parameters.clipFar > -stage.viewer.camera.position.z
        && -stage.viewer.camera.position.z > stage.viewer.parameters.clipNear)
    cameradist = -stage.viewer.camera.position.z;
  else
    cameradist = cvorient[14]; // fall back if stage.viewer.camera.position.z is corrupted
  cvorient.push( cameradist );
  msg = String(cvorient);
  oldmsg = msg;
  return msg;
}


  // listen to `hovered` signal to move tooltip around and change its text
PickingProxyfunc = function(pickingProxy)
{
  if (pickingProxy
        && (Object.prototype.toString.call(pickingProxy.picker["ids"]) === '[object Array]' )
        && displaytooltips )
  {
    var cp = pickingProxy.canvasPosition;
    var sym_id = -1;
    var hkl_id = -1;
    var ttipid = "";
    if (pickingProxy.picker["ids"].length > 0)
    { // get stored id number of symmetry operator applied to this hkl
      sym_id = pickingProxy.picker["ids"][0];
      var ids = pickingProxy.picker["ids"].slice(1);
      var is_friedel_mate = 0;
      hkl_id = ids[ pickingProxy.pid %% ids.length ];
      if (pickingProxy.pid >= ids.length)
        is_friedel_mate = 1;
    }
    // tell python the id of the hkl and id number of the symmetry operator
    rightnow = timefunc();
    if (rightnow - timenow > tdelay)
    { // only post every 50 milli second as not to overwhelm python
      ttipid = String([hkl_id, sym_id, is_friedel_mate]);
      WebsockSendMsg( 'tooltip_id: [' + ttipid + ']' );
      timenow = timefunc();
    }

    if (isdebug)
      console.log( "current_ttip_ids: " + String(current_ttip_ids) + ", ttipid: " + String(ttipid) );
    if (current_ttip !== "" && current_ttip_ids == ttipid )
    {
      tooltip.innerText = current_ttip;
      tooltip.style.bottom = cp.y + 7 + "px";
      tooltip.style.left = cp.x + 8 + "px";
      tooltip.style.fontSize = "smaller";
      tooltip.style.display = "block";
    }
  }
  else
  {
    tooltip.style.display = "none";
    current_ttip = "";
  }
};


function ColourChart(millerlabel, fomlabel)
{
  var ih = 3,
  topr = 35,
  topr2 = 10,
  lp = 10,
  wp = 60,
  lp2 = lp + wp,
  gl = 3,
  wp2 = gl,
  fomlabelheight = 25;
  if (colourgradvalarray.length === 1)
  {
    wp2 = 15;
    fomlabelheight = 0;
  }

  var wp3 = wp + colourgradvalarray.length * wp2 + 2;

  totalheight = ih*colourgradvalarray[0].length + 35 + fomlabelheight;
  // make a white box on top of which boxes with transparent background are placed
  // containing the colour values at regular intervals as well as label legend of
  // the displayed miller array
  addDivBox("", topr2, lp, wp3, totalheight, 'rgba(255, 255, 255, 1.0)');

  // print label of the miller array used for colouring
  addDivBox(millerlabel, topr2, lp, wp, 20);

  if (colourgradvalarray.length > 1)
  {
    // print FOM label, 1, 0.5 and 0.0 values below colour chart
    fomtop = topr2 + totalheight - 18;
    fomlp = lp + wp;
    fomwp = wp3;
    fomtop2 = fomtop - 13;
    // print the 1 number
    addDivBox("1", fomtop2, fomlp, fomwp, 20);
    // print the 0.5 number
    leftp = fomlp + 0.48 * gl * colourgradvalarray.length;
    addDivBox("0.5", fomtop2, leftp, fomwp, 20);
    // print the FOM label
    addDivBox(fomlabel, fomtop, fomlp, fomwp, 20);
    // print the 0 number
    leftp = fomlp + 0.96 * gl * colourgradvalarray.length;
    addDivBox("0", fomtop2, leftp, fomwp, 20);
  }

  for (j = 0; j < colourgradvalarray[0].length; j++)
  {
    rgbcol = colourgradvalarray[0][j][1];
    val = colourgradvalarray[0][j][0];
    topv = j*ih + topr;
    toptxt = topv - 5;
    // print value of miller array if present in colourgradvalarray[0][j][0]
    addDivBox(val, toptxt, lp, wp, ih);
  }

  // draw the colour gradient
  for (g = 0; g < colourgradvalarray.length; g++)
  {
    leftp = g*gl + lp + wp;
    // if FOM values are supplied draw colour gradients with decreasing
    // saturation values as stored in the colourgradvalarray[g] arrays
    for (j = 0; j < colourgradvalarray[g].length; j++)
    {
      rgbcol = colourgradvalarray[g][j][1];
      val = colourgradvalarray[g][j][0];
      topv = j*ih + topr;
      addDivBox("", topv, leftp, wp2, ih, rgbcol);
    }
  }
}


function HKLscene()
{
  shape = new NGL.Shape('shape');
  //vectorshape = new NGL.Shape('vectorshape');
  stage = new NGL.Stage('viewport', {  backgroundColor: "rgb(128, 128, 128)",
                                      tooltip:false, // create our own tooltip from a div element
                                      fogNear: 100, fogFar: 100 });
  stage.setParameters( { cameraType: "%s" } );

  /*
  canvas = stage.viewer.renderer.domElement;
  const ctx = canvas.getContext('webgl', {
    desynchronized: true,
    preserveDrawingBuffer: true
  });
  */

  MakeHKL_Axis(shape);

  %s

  shapeComp = stage.addComponentFromObject(shape);
  repr = shapeComp.addRepresentation('buffer');
  shapeComp.autoView();
  repr.update();

  // if some radii are negative draw them with wireframe
  %s
  %s

  if (isdebug)
    stage.viewer.container.appendChild(debugmessage);

  // avoid NGL zoomFocus messing up clipplanes positions. So reassign those signals to zoomDrag
  stage.mouseControls.remove("drag-shift-right");
  stage.mouseControls.add("drag-shift-right", NGL.MouseActions.zoomDrag);
  stage.mouseControls.remove("drag-middle");
  stage.mouseControls.add("drag-middle", NGL.MouseActions.zoomDrag);
  stage.mouseControls.remove('clickPick-left'); // avoid undefined move-pick when clicking on a sphere

  stage.viewer.requestRender();
  if (isdebug)
    debugmessage.innerText = dbgmsg;
}


function OnUpdateOrientation()
{
  msg = getOrientMsg();
  WebsockSendMsg('MouseMovedOrientation:\\n' + msg );
}


try
{
  document.addEventListener('DOMContentLoaded', function() { HKLscene() }, false );
  document.addEventListener('mouseup', function() { OnUpdateOrientation() }, false );
  document.addEventListener('wheel', function(e) { OnUpdateOrientation() }, false );
  document.addEventListener('scroll', function(e) { OnUpdateOrientation() }, false );
  // mitigate flickering on some PCs when resizing
  document.addEventListener('resize', function() { RenderRequest() }, false );
}
catch(err)
{
  WebsockSendMsg('JavaScriptError: ' + err.stack );
}

    """ % (self.websockport, self.__module__, self.__module__, cntbin, str(self.verbose>=2).lower(), \
           self.ngl_settings.tooltip_alpha, axisfuncstr, self.camera_type, spherebufferstr, \
           negativeradiistr, colourscriptstr)

    WebsockMsgHandlestr = """

function ReturnClipPlaneDistances()
{
  if (stage.viewer.parameters.clipScale == 'relative')
    cameradist = stage.viewer.cDist;
  if (stage.viewer.parameters.clipScale == 'absolute')
    if (stage.viewer.cDist != 0
         && stage.viewer.parameters.clipFar > stage.viewer.cDist
         && stage.viewer.cDist > stage.viewer.parameters.clipNear)
      cameradist = stage.viewer.cDist;
    else if (stage.viewer.camera.position.z != 0
         && stage.viewer.parameters.clipFar > -stage.viewer.camera.position.z
         && -stage.viewer.camera.position.z > stage.viewer.parameters.clipNear)
      cameradist = stage.viewer.camera.position.z;
    else if (stage.viewer.camera.position.z == -stage.viewer.cDist)
      cameradist = stage.viewer.cDist;
    else
      return;

  msg = String( [stage.viewer.parameters.clipNear,
                  stage.viewer.parameters.clipFar,
                  cameradist ] )
  WebsockSendMsg('ReturnClipPlaneDistances:\\n' + msg );
}


mysocket.onmessage = function(e)
{
  var c,
  si;
  WebsockSendMsg('\\n    Browser: Got ' + e.data ); // tell server what it sent us
  try
  {
    var datval = e.data.split(":\\n");
    var msgtype = datval[0];
    var val = datval[1].split(",");

    if (msgtype === "alpha")
    {
      bin = parseInt(val[0]);
      alphas[bin] = parseFloat(val[1]);
      shapebufs[bin].setParameters({opacity: alphas[bin]});
      for (var g=0; g < nrots; g++ )
        br_shapebufs[bin][g].setParameters({opacity: alphas[bin]});
      //stage.viewer.requestRender();
      RenderRequest();
    }

    if (msgtype === "colour")
    {
      bin = parseInt(val[0]);
      si =  parseInt(val[1]);
      colours[bin][3*si] = parseFloat(val[2]);
      colours[bin][3*si+1] = parseFloat(val[3]);
      colours[bin][3*si+2] = parseFloat(val[4]);
      shapebufs[bin].setAttributes({ color: colours[bin] });

      for (var g=0; g < nrots; g++ )
      {
        br_colours[bin][3*si] = parseFloat(val[2]);
        br_colours[bin][3*si+1] = parseFloat(val[3]);
        br_colours[bin][3*si+2] = parseFloat(val[4]);
        br_shapebufs[bin][g].setAttributes({ color: br_colours[bin] });
      }
      //stage.viewer.requestRender();
      RenderRequest();
    }

    if (msgtype === "DisplayTooltips")
    {
      displaytooltips = val[0];
      stage.signals.hovered.removeAll();
      stage.signals.clicked.removeAll();
      if (displaytooltips == "click")
        stage.signals.clicked.add( PickingProxyfunc );
      if (displaytooltips == "hover")
        stage.signals.hovered.add( PickingProxyfunc );
    }

    if (msgtype === "ShowThisTooltip")
    {
      current_ttip = eval(datval[1]).split("\\n\\n")[0];
      current_ttip_ids = eval(datval[1]).split("\\n\\n")[1];
    }

    if (msgtype === "TooltipOpacity")
    {
      Object.assign(tooltip.style, {
        backgroundColor: "rgba(255, 255, 255, " + val[0] + " )",
      });
    }

    if (msgtype === "Redraw")
    {
      //stage.viewer.requestRender();
      RenderRequest();
    }

    if (msgtype === "ReOrient")
    {
      WebsockSendMsg( 'Reorienting ' + pagename );
      sm = new Float32Array(16);
      for (j=0; j<16; j++)
      {
        sm[j] = parseFloat(val[j]);
        if (isNaN( sm[j] ))
          return; // do nothing just in case
      }

      var m = new NGL.Matrix4();
      m.fromArray(sm);
      stage.viewerControls.orient(m);
      //stage.viewer.renderer.setClearColor( 0xffffff, 0.01);
      //stage.viewer.requestRender();
      RenderRequest();
      msg = getOrientMsg();
      WebsockSendMsg('CurrentViewOrientation:\\n' + msg );
    }

    if (msgtype === "Reload")
    {
    // refresh browser with the javascript file
      msg = getOrientMsg();
      WebsockSendMsg('OrientationBeforeReload:\\n' + msg );
      WebsockSendMsg( 'Refreshing ' + pagename );
      window.location.reload(true);
      // Now we are gone. A new javascript file has been loaded in the browser
    }

    if (msgtype.includes("Expand") )
    {
      WebsockSendMsg( 'Expanding data...' );

      if (msgtype == "Expand" && expstate == "")
        return;

      if (msgtype == "ExpandP1" && expstate == "isP1Expanded")
        return;

      if (msgtype == "ExpandFriedel" && expstate == "isFriedelExpanded")
        return;

      if (msgtype == "ExpandP1Friedel" && expstate == "isP1FriedelExpanded")
        return;

      // delete the shapebufs[] that holds the positions[] arrays
      shapeComp.removeRepresentation(repr);
      // remove shapecomp from stage first
      stage.removeComponent(shapeComp);

      br_positions = [];
      br_colours = [];
      br_radii = [];
      br_ttips = [];
      br_shapebufs = [];
      var nexpandrefls = 0;

      //alert('rotations:\\n' + val);
      // Rotation matrices are concatenated to a string of floats
      // separated by line breaks between each roation matrix
      rotationstrs = datval[1].split("\\n");
      var Rotmats = [];
      var r = new NGL.Vector3();

      for (var rotmxidx=0; rotmxidx < rotationstrs.length; rotmxidx++ )
      {
        Rotmats.push( new NGL.Matrix3() );
        // convert string of rotation matrix elements into a Matrix3
        var elmstrs = rotationstrs[rotmxidx].split(",");
        for (j=0; j<9; j++)
          Rotmats[rotmxidx].elements[j] = parseFloat(elmstrs[j]);
      }

      var Imx = new NGL.Matrix3();
      Imx.identity(); // for testing
      if ( !(msgtype.includes("P1")) && rotationstrs.length == 1 && Rotmats[0].equals(Imx) )
        throw "Only the identity matrix is provided. That means no P1 expansion of reflections!";

      for (var bin=0; bin<nbins; bin++)
      {
        var nsize = positions[bin].length/3; // number of reflections in each bin
        var csize = nsize*3;
        var nsize3 = nsize*3;
        var anoexp = false;

        if (msgtype.includes("Friedel") )
        {
          anoexp = true;
          csize = nsize*6;
        }
        br_positions.push( [] );
        br_shapebufs.push( [] );
        br_colours.push( [] );
        br_radii.push( [] );
        br_ttips.push( [] );

        br_colours[bin] = colours[bin];
        br_radii[bin] = radii[bin];
        if (anoexp)
        {
          var colarr = [];
          var cl = colours[bin].length;
          for (var i=0; i<cl; i++)
          {
            colarr[i] = colours[bin][i];
            colarr[i+cl] = colours[bin][i];
          }
          br_colours[bin] = new Float32Array(colarr);

          var radiiarr = [];
          var rl = radii[bin].length;
          for (var i=0; i<rl; i++)
          {
            radiiarr[i] = radii[bin][i];
            radiiarr[i+rl] = radii[bin][i];
          }
          br_radii[bin] = new Float32Array(radiiarr);
        }

        nrots = 0;
        nexpandrefls = 0;
        for (var rotmxidx=0; rotmxidx < rotationstrs.length; rotmxidx++ )
        {
          if (rotationstrs[rotmxidx].length < 1 )
            continue;
          nrots++;

          br_positions[bin].push( [] );
          br_shapebufs[bin].push( [] );
          br_ttips[bin].push( [] );
          Object.assign(br_ttips[bin][rotmxidx], ttips[bin]); // deep copy the ttips[bin] object
          br_ttips[bin][rotmxidx][0] = rotmxidx;
          br_positions[bin][rotmxidx] = new Float32Array( csize );
          nexpandrefls += csize;

          for (var i=0; i<nsize; i++)
          {
            idx= i*3;
            r.x = positions[bin][idx];
            r.y = positions[bin][idx+1];
            r.z = positions[bin][idx+2];

            r.applyMatrix3(Rotmats[rotmxidx])

            br_positions[bin][rotmxidx][idx] = r.x;
            br_positions[bin][rotmxidx][idx + 1] = r.y;
            br_positions[bin][rotmxidx][idx + 2] = r.z;

            if (anoexp)
            {
              r.negate(); // inversion for anomalous pair
              br_positions[bin][rotmxidx][nsize3 + idx] = r.x;
              br_positions[bin][rotmxidx][nsize3 + idx + 1] = r.y;
              br_positions[bin][rotmxidx][nsize3 + idx + 2] = r.z;
            }
          }

          br_shapebufs[bin][rotmxidx] = new NGL.SphereBuffer({
              position: br_positions[bin][rotmxidx],
              color: br_colours[bin],
              radius: br_radii[bin],
              // rotmxidx works as the id of the rotation of applied symmetry operator when creating tooltip for an hkl
              picking: br_ttips[bin][rotmxidx],
              } %s  );
          shape.addBuffer(br_shapebufs[bin][rotmxidx]);
          //WebsockSendMsg( 'Memory usage: ' + String(window.performance.memory.totalJSHeapSize) +
          //        ', ' + String(window.performance.memory.totalJSHeapSize) );
        }
        if (nexpandrefls == nsize*3)
          expstate = "";
        if (nexpandrefls == nsize*6)
          expstate = "isFriedelExpanded";
        if (nexpandrefls == nsize*3*nrots && nrots > 1)
          expstate = "isP1Expanded";
        if (nexpandrefls == nsize*6*nrots && nrots > 1)
          expstate = "isP1FriedelExpanded";
      }
      MakeHKL_Axis(shape);

      shapeComp = stage.addComponentFromObject(shape);
      repr = shapeComp.addRepresentation('buffer');

      for (var bin=0; bin<nbins; bin++)
      {
        for (var rotmxidx=0; rotmxidx < nrots; rotmxidx++ )
        {
          br_shapebufs[bin][rotmxidx].setParameters({opacity: alphas[bin]});
        }
      }

      //stage.viewer.requestRender();
      RenderRequest();
      WebsockSendMsg( 'Done ' + msgtype );
    }

    if (msgtype === "DisableMouseRotation")
    {
      WebsockSendMsg( 'Fix mouse rotation' + pagename );
      stage.mouseControls.remove("drag-left");
      stage.mouseControls.remove("scroll-ctrl");
      stage.mouseControls.remove("scroll-shift");
    }

    if (msgtype === "EnableMouseRotation")
    {
      WebsockSendMsg( 'Can mouse rotate ' + pagename );
      stage.mouseControls.add("drag-left", NGL.MouseActions.rotateDrag);
      stage.mouseControls.add("scroll-ctrl", NGL.MouseActions.scrollCtrl);
      stage.mouseControls.add("scroll-shift", NGL.MouseActions.scrollShift);
    }

    if (msgtype === "RotateStage")
    {
      WebsockSendMsg( 'Rotating stage ' + pagename );

      var sm = new Float32Array(9);
      var m4 = new NGL.Matrix4();

      for (j=0; j<9; j++)
        sm[j] = parseFloat(val[j]);

      // GL matrices are the transpose of conventional rotation matrices
      m4.set( sm[0], sm[3], sm[6], 0.0,
              sm[1], sm[4], sm[7], 0.0,
              sm[2], sm[5], sm[8], 0.0,
              0.0,   0.0,   0.0,   1.0
      );
      stage.viewerControls.orient(m4);
      if (val[9]=="verbose")
        postrotmxflag = true;
      ReturnClipPlaneDistances();
      //ReRender();
      //stage.viewer.requestRender();
      RenderRequest();
      sleep(100).then(()=> {
          msg = getOrientMsg();
          WebsockSendMsg('CurrentViewOrientation:\\n' + msg );
        }
      );
    }

    if (msgtype === "SpinAnimate")
    {
      WebsockSendMsg( 'SpinAnimating ' + pagename );
      //strs = datval[1].split("\\n");
      var r = new Float32Array(3);
      //var elmstrs = strs[0].split(",");
      for (j=0; j<3; j++)
        r[j] = parseFloat(val[j]);
      if (r[0] == 0.0 && r[1] == 0.0 && r[2] == 0.0)
      {
        // default bindings as per ngl\src\controls\mouse-actions.ts
        stage.mouseControls.add("drag-ctrl-left", NGL.MouseActions.panDrag);
        stage.mouseControls.add("drag-ctrl-right", NGL.MouseActions.focusScroll);
        stage.mouseControls.add("drag-shift-left", NGL.MouseActions.zoomDrag);
        stage.mouseControls.add("drag-shift-right", NGL.MouseActions.zoomDrag);
        stage.mouseControls.add("drag-middle", NGL.MouseActions.zoomDrag);
        stage.mouseControls.add("drag-right", NGL.MouseActions.panDrag);
        stage.mouseControls.add("drag-left", NGL.MouseActions.rotateDrag);
        stage.mouseControls.add("scroll-ctrl", NGL.MouseActions.scrollCtrl);
        stage.mouseControls.add("scroll-shift", NGL.MouseActions.scrollShift);
        stage.setSpin(false);
      }
      else
      {
        stage.spinAnimation.axis.set(r[0], r[1], r[2]);

        stage.mouseControls.remove("drag-ctrl-left");
        stage.mouseControls.remove("drag-ctrl-right");
        stage.mouseControls.remove("drag-shift-left");
        stage.mouseControls.remove("drag-shift-right");
        stage.mouseControls.remove("drag-middle");
        stage.mouseControls.remove("drag-right");
        stage.mouseControls.remove("drag-left");
        stage.mouseControls.remove("scroll-ctrl");
        stage.mouseControls.remove("scroll-shift");
        stage.setSpin(true);
      }
    }

    if (msgtype === "TranslateHKLpoints")
    {
      WebsockSendMsg( 'Translating HKLs ' + pagename );
      strs = datval[1].split("\\n");
      var sm = new Float32Array(3);
      var elmstrs = strs[0].split(",");
      for (j=0; j<3; j++)
        sm[j] = parseFloat(elmstrs[j]);
      shapeComp.setPosition([ sm[0], sm[1], sm[2] ]);
      //stage.viewer.requestRender();
      RenderRequest();
      sleep(100).then(()=> {
          msg = getOrientMsg();
          WebsockSendMsg('CurrentViewOrientation:\\n' + msg );
        }
      );
    }

    function DeleteVectors(reprname)
    {
      thisrepr = stage.getRepresentationsByName(reprname);
      var wasremoved = false;
      for (i=0; i<stage.compList.length; i++)
        if (stage.compList[i].reprList[0].name == reprname)
        {
          thiscomp = stage.compList[i];
          thiscomp.removeRepresentation(thisrepr);
          stage.removeComponent(thiscomp);
          wasremoved = true;
        }
      return wasremoved;
    };

    if (msgtype === "AddVector")
    {
      var r1 = new Float32Array(3);
      var r2 = new Float32Array(3);
      var rgb = new Float32Array(3);
      for (j=0; j<3; j++)
      {
        r1[j] = parseFloat(val[j]);
        r2[j] = parseFloat(val[j+3]);
        rgb[j]= parseFloat(val[j+6]);
      }
      radius = parseFloat(val[11]);

      if (vectorshape == null)
        vectorshape = new NGL.Shape('vectorshape');

      vectorshape.addArrow( r1, r2 , [rgb[0], rgb[1], rgb[2]], radius);
      if (val[6] !== "") {
        var txtR = [ (r1[0] + r2[0])/2.0, (r1[1] + r2[1])/2.0, (r1[2] + r2[2])/2.0 ];
        vectorshape.addText( txtR, [rgb[0], rgb[1], rgb[2]], fontsize, val[9] );
      }
      // if reprname is supplied with a vector then make a representation named reprname
      // of this and all pending vectors stored in vectorshape and render them.
      // Otherwise just accummulate the new vector
      var reprname = val[10].trim();
      if (reprname != "")
      {
        DeleteVectors(reprname); // delete any existing vectors with the same name
        vectorshapeComps.push( stage.addComponentFromObject(vectorshape) );
        vectorreprs.push(
          vectorshapeComps[vectorshapeComps.length-1].addRepresentation('vecbuf',
                                                                      { name: reprname} )
        );
        vectorshape = null;
        //stage.viewer.requestRender();
        RenderRequest();
      }
    }

    if (msgtype === "RemoveVectors")
    {
      var reprname = val[0].trim(); // elmstrs[0].trim();
      // if reprname is supplied only remove vectors with that name
      var reprnamegone = false;
      var clipvecgone = false;
      var unitcellgone = false;
      var reciprocunitcellgone = false;
      if (reprname != "")
        reprnamegone = DeleteVectors(reprname);
      else // otherwise remove all vectors
      {
        clipvecgone = DeleteVectors("clip_vector");
        unitcellgone = DeleteVectors("unitcell");
        reciprocunitcellgone = DeleteVectors("reciprocal_unitcell");
      }
      if (reprnamegone || clipvecgone || unitcellgone || reciprocunitcellgone)
        //stage.viewer.requestRender();
        RenderRequest();
    }

    if (msgtype === "SetMouseSpeed")
    {
      stage.trackballControls.rotateSpeed = parseFloat(val[0]);
    }

    if (msgtype === "GetMouseSpeed")
    {
      msg = String( [stage.trackballControls.rotateSpeed] )
      WebsockSendMsg('ReturnMouseSpeed:\\n' + msg );
    }

    if (msgtype === "SetClipPlaneDistances")
    {
      var near = parseFloat(val[0]);
      var far = parseFloat(val[1]);
      origcameraZpos = parseFloat(val[2]);
      stage.viewer.parameters.clipMode =  'camera';
      // clipScale = 'absolute' means clip planes are using scene dimensions
      stage.viewer.parameters.clipScale = 'absolute';
      clipFixToCamPosZ = true;

      if (near >= far )
      { // default to no clipping if near >= far
        stage.viewer.parameters.clipMode = 'scene';
      // clipScale = 'relative' means clip planes are in percentage
        stage.viewer.parameters.clipScale = 'relative';
        clipFixToCamPosZ = false;
        near = 0;
        far = 100;
      }
      else
        stage.viewer.camera.position.z = origcameraZpos;
      stage.viewer.parameters.clipNear = near;
      stage.viewer.parameters.clipFar = far;
      origclipnear = near;
      origclipfar = far;
      //stage.viewer.requestRender();
      RenderRequest();
    }

    if (msgtype === "GetClipPlaneDistances")
      ReturnClipPlaneDistances();

    if (msgtype === "GetBoundingBox")
    {
      msg = String( [stage.viewer.boundingBoxSize.x,
                     stage.viewer.boundingBoxSize.y,
                     stage.viewer.boundingBoxSize.z]
                  )
      WebsockSendMsg('ReturnBoundingBox:\\n' + msg );
    }

    if (msgtype ==="JavaScriptCleanUp")
    {
      stage.removeAllComponents();
      stage.mouseObserver.dispose();
      ttips = [];
      vectorreprs = [];
      vectorshapeComps = [];
      positions = [];
      br_positions = [];
      br_colours = [];
      br_radii = [];
      br_ttips = [];
      colours = [];
      alphas = [];
      radii = [];
      shapebufs = [];
      br_shapebufs = [];
      shape = null;
      shapeComp = null;
      vectorshape = null;
      repr = null;
      stage.dispose();
      stage = null;
      WebsockSendMsg('CleanUp:\\nDestroying JavaScript objects');
      document = null;
    }

    if (msgtype === "InjectNewReflections")
    {
      WebsockSendMsg( 'Rendering new reflections ' + pagename );
      var nrefl = parseInt(val.length/7);
      if (nrefl !== val.length/7)
      {
        alert("Mismatch in array of reflections, colours and radii!")
        return;
      }

      // delete the shapebufs[] that holds the positions[] arrays
      shapeComp.removeRepresentation(repr);
      // remove shapecomp from stage first
      stage.removeComponent(shapeComp);

      positions = [];
      colours = [];
      radii = [];
      alphas = [];
      shapebufs = [];
      ttips = [];
      shapebufs = [];
      nbins = 1; // currently no binning when injecting reflections

      positions_ = []; // dummy variables for conforming to binning scheme above
      colours_ = [];   // as used when expanding reflections
      radii_ = [];
      ttips_ = [-1]

      for (j=0; j<nrefl; j++)
      {
        positions_.push( parseFloat(val[7*j]) );
        positions_.push( parseFloat(val[7*j+1]) );
        positions_.push( parseFloat(val[7*j+2]) );
        colours_.push( parseFloat(val[7*j+3]) );
        colours_.push( parseFloat(val[7*j+4]) );
        colours_.push( parseFloat(val[7*j+5]) );
        radii_.push( parseFloat(val[7*j+6]) );
        ttips_.push(j)
      }

      positions.push( new Float32Array( positions_ ));
      colours.push( new Float32Array( colours_ ));
      radii.push( new Float32Array( radii_ ));
      ttips.push(ttips_);

      shapebufs.push( new NGL.SphereBuffer({
        position: positions[0],
        color: colours[0],
        radius: radii[0],
        picking: ttips[0],
        })
      );
      shape.addBuffer(shapebufs[0]);
      alphas.push(1.0);

      MakeHKL_Axis(shape);
      shapeComp = stage.addComponentFromObject(shape);
      repr = shapeComp.addRepresentation('buffer');
      //stage.viewer.requestRender();
      RenderRequest();
      WebsockSendMsg('Injected new reflections');
    }

    if (msgtype === "SetAutoView")
    {
      if (shapeComp != null) // workaround for QTWebEngine bug sometimes failing to render scene
        shapeComp.autoView();
      WebsockSendMsg('AutoViewSet ' + pagename);
    }

    if (msgtype === "Testing")
    {
      // test something new
      /*
      var newradii = radii[0].map(function(element) {
        return element*1.5;
      });
      shapebufs[0].setAttributes({
          radius: newradii
      })
      repr = shapeComp.addRepresentation('buffer');
      //stage.viewer.requestRender();
      RenderRequest();
      */
    }
    WebsockSendMsg('Received message: ' + msgtype );
    if (isdebug)
      debugmessage.innerText = dbgmsg;
  }

  catch(err)
  {
    WebsockSendMsg('JavaScriptError: ' + err.stack );
  }

};

    """ %qualitystr

    self.NGLscriptstr += WebsockMsgHandlestr
    if self.jscriptfname:
      with open( self.jscriptfname, "w") as f:
        f.write( self.NGLscriptstr )
    self.ReloadNGL()
    if not blankscene:
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
    self.sceneisdirty = False
    self.lastscene_id = self.viewerparams.scene_id


  def OnWebsocketClientMessage(self, client, server, message):
    if self.viewerparams.scene_id is None or self.miller_array is None:
      return
    try:
      if message != "":
        if "Orientation" in message:
          self.ProcessOrientationMessage(message)
        elif 'Received message:' in message:
          self.mprint( message, verbose=2)
        elif "websocket" in message:
          self.mprint( message, verbose=1)
        elif "AutoViewSet" in message:
          self.set_volatile_params()
        elif "JavaScriptCleanUp:" in message:
          self.mprint( message, verbose=1)
          self.StopThreads()
        elif "JavaScriptError:" in message:
          self.mprint( message, verbose=0)
        elif "Expand" in message:
          self.mprint( message, verbose=0)
          #raise Sorry(message)
        elif "ReturnClipPlaneDistances:" in message:
          datastr = message[ message.find("\n") + 1: ]
          lst = datastr.split(",")
          flst = [float(e) for e in lst]
          self.clipNear = flst[0]
          self.clipFar = flst[1]
          self.cameraPosZ = flst[2]
        elif "ReturnBoundingBox:" in message:
          datastr = message[ message.find("\n") + 1: ]
          lst = datastr.split(",")
          flst = [float(e) for e in lst]
          self.boundingX = flst[0]
          self.boundingY = flst[1]
          self.boundingZ = flst[2]
        elif "ReturnMouseSpeed" in message:
          datastr = message[ message.find("\n") + 1: ]
          lst = datastr.split(",")
          flst = [float(e) for e in lst]
          if flst[0] is not None and not cmath.isnan(flst[0]):
            self.ngl_settings.mouse_sensitivity = flst[0]
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
          #self.send_msg_to_browser("ShowThisTooltip", ttip)
          self.AddToBrowserMsgQueue("ShowThisTooltip", ttip)
        else:
          if "Ready " in message:
            self.mprint( message, verbose=5)
        self.lastmsg = message
    except Exception as e:
      self.mprint( to_str(e) + "\n" + traceback.format_exc(limit=10), verbose=0)


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
      return
    else:
      cameradist = math.pow(rotdet, 1.0/3.0)
    self.mprint("Scale distance: %s" %roundoff(cameradist), verbose=3)
    currentRotmx = matrix.identity(3)
    if cameradist > 0.0:
      currentRotmx = ScaleRotMx/cameradist
      cameraPosZ = cameradist
    return cameraPosZ, currentRotmx, cameratranslation


  def ProcessOrientationMessage(self, message):
    if message.find("NaN")>=0 or message.find("undefined")>=0:
      return
    if "OrientationBeforeReload:" in message:
      #sleep(0.2)
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
    #while not self.websockclient :
    while not self.browserisopen :
      time.sleep(self.sleeptime)
      nwait += self.sleeptime
      if nwait > sec:
        return False
    return True


  def AddToBrowserMsgQueue(self, msgtype, msg=""):
    self.msgqueue.append( (msgtype, msg) )


  def WebBrowserMsgQueue(self):
    try:
      while True:
        nwait = 0.0
        sleep(self.sleeptime)
        if self.javascriptcleaned:
          self.mprint("Shutting down WebBrowser message queue", verbose=1)
          return
        if len(self.msgqueue):
          pendingmessagetype, pendingmessage = self.msgqueue[0]
          self.send_msg_to_browser(pendingmessagetype, pendingmessage)
          while not self.browserisopen:  #self.websockclient:
            sleep(self.sleeptime)
            nwait += self.sleeptime
            if nwait > self.handshakewait or self.javascriptcleaned or not self.viewerparams.scene_id is not None:
              return
          self.msgqueue.remove( self.msgqueue[0] )
          #if self.was_disconnected:
          #  nwait2 = 0.0
          #  while nwait2 < self.handshakewait:
          #    nwait2 += self.sleeptime
          #  self.ReloadNGL()
# if the html content is huge the browser will be unresponsive until it has finished
# reading the html content. This may crash this thread. So try restarting this thread until
# browser is ready
    except Exception as e:
      self.mprint( str(e) + ", Restarting WebBrowserMsgQueue\n" \
                         + traceback.format_exc(limit=10), verbose=2)
      self.websockclient = None
      self.WebBrowserMsgQueue()


  def OnConnectWebsocketClient(self, client, server):
    self.websockclient = client
    self.mprint( "Browser connected:" + str( self.websockclient ), verbose=1 )
    if self.was_disconnected:
      self.was_disconnected = False
    if self.lastviewmtrx and self.viewerparams.scene_id is not None:
      self.set_volatile_params()
      self.mprint( "Reorienting client after refresh:" + str( self.websockclient ), verbose=2 )
      self.AddToBrowserMsgQueue("ReOrient", self.lastviewmtrx)
    else:
      self.SetAutoView()


  def OnDisconnectWebsocketClient(self, client, server):
    self.mprint( "Browser disconnected:" + str( client ), verbose=1 )
    self.was_disconnected = True


  def send_msg_to_browser(self, msgtype, msg=""):
    message = u"" + msgtype + self.msgdelim + str(msg)
    if self.websockclient:
      nwait = 0.0
      while not ("Ready" in self.lastmsg or "tooltip_id" in self.lastmsg \
        or "CurrentViewOrientation" in self.lastmsg or "AutoViewSet" in self.lastmsg \
        or "ReOrient" in self.lastmsg):
        sleep(self.sleeptime)
        nwait += self.sleeptime
        if nwait > self.handshakewait and self.browserisopen:
          self.mprint("ERROR: No handshake from browser!", verbose=0 )
          self.mprint("failed sending " + msgtype, verbose=1)
          break
      self.server.send_message(self.websockclient, message )
    else:
      self.OpenBrowser()


  def OpenBrowser(self):
    if self.viewerparams.scene_id is not None and not (self.websockclient or self.browserisopen):
      NGLlibpath = libtbx.env.under_root(os.path.join("modules","cctbx_project","crys3d","hklview","ngl.js") )
      htmlstr = self.hklhtml %(NGLlibpath, os.path.abspath( self.jscriptfname))
      htmlstr += self.htmldiv
      with open(self.hklfname, "w") as f:
        f.write( htmlstr )
      self.url = "file:///" + os.path.abspath( self.hklfname )
      self.url = self.url.replace("\\", "/")
      self.mprint( "Writing %s and connecting to its websocket client" %self.hklfname, verbose=1)
      if self.UseOSBrowser:
        webbrowser.open(self.url, new=1)
      self.SendInfoToGUI({ "html_url": self.url } )
      self.browserisopen = True
      self.isnewfile = False


  def StartWebsocket(self):
    self.server = WebsocketServer(self.websockport, host='127.0.0.1')
    if not self.server:
      raise Sorry("Could not connect to web browser")
    self.server.set_fn_new_client(self.OnConnectWebsocketClient)
    self.server.set_fn_client_left(self.OnDisconnectWebsocketClient)
    self.server.set_fn_message_received(self.OnWebsocketClientMessage)
    self.wst = threading.Thread(target=self.server.run_forever)
    self.wst.daemon = True
    self.wst.start()
    self.msgqueuethrd = threading.Thread(target = self.WebBrowserMsgQueue )
    self.msgqueuethrd.daemon = True
    self.msgqueuethrd.start()


  def StopThreads(self):
    try:
      if self.websockclient: # might not have been created if program is closed before a data set is shown
        self.websockclient['handler'].send_text(u"", opcode=0x8)
    except Exception as e:
      self.mprint( to_str(e) + "\n" + traceback.format_exc(limit=10), verbose=0)
    self.mprint("Shutting down Websocket listening thread", verbose=1)
    self.server.shutdown()
    self.javascriptcleaned = True
    self.msgqueuethrd.join()
    self.mprint("Shutting down WebsocketServer", verbose=1)
    self.wst.join()
    self.isterminating = True


  def set_camera_type(self):
    self.camera_type = self.ngl_settings.camera_type


  def set_show_tooltips(self):
    msg = "%s" %self.ngl_settings.show_tooltips
    self.send_msg_to_browser("DisplayTooltips", msg)


  def set_tooltip_opacity(self):
    msg = "%f" %self.ngl_settings.tooltip_alpha
    self.send_msg_to_browser("TooltipOpacity", msg)


  def SetOpacities(self, bin_opacities_str):
    retstr = ""
    if self.miller_array and bin_opacities_str and not self.isinjected:
      self.ngl_settings.bin_opacities = bin_opacities_str
      bin_opacitieslst = eval(self.ngl_settings.bin_opacities)
      for binopacity in bin_opacitieslst:
        alpha = binopacity[0] # float(binopacity.split(",")[0])
        bin = binopacity[1] # int(binopacity.split(",")[1])
        retstr += self.set_opacity(bin, alpha)
      self.SendInfoToGUI( { "bin_opacities": self.ngl_settings.bin_opacities } )
    return retstr


  def set_opacity(self, bin, alpha):
    if bin > self.nbinvalsboundaries-1:
      return "There are only %d bins present\n" %self.nbinvalsboundaries
    msg = "%d, %f" %(bin, alpha)
    self.send_msg_to_browser("alpha", msg)
    return "Opacity %s set on bin[%s]\n" %(alpha, bin)


  def RedrawNGL(self):
    self.AddToBrowserMsgQueue("Redraw")


  def ReloadNGL(self): # expensive as javascript may be several Mbytes large
    self.mprint("Rendering JavaScript...", verbose=1)
    self.AddToBrowserMsgQueue("Reload")


  def JavaScriptCleanUp(self):
    self.AddToBrowserMsgQueue("JavaScriptCleanUp")
    if self.viewerparams.scene_id is None: # nothing is showing in the browser
      self.StopThreads()


  def ExpandInBrowser(self, P1=True, friedel_mate=True):
    retmsg = "Not expanding in browser\n"
    if self.sceneisdirty:
      return retmsg
    uc = self.miller_array.unit_cell()
    OrtMx = matrix.sqr( uc.orthogonalization_matrix())
    InvMx = OrtMx.inverse()
    msgtype = "Expand"
    msg = ""
    unique_rot_ops = []
    if P1:
      msgtype += "P1"
      unique_rot_ops = self.symops[ 0 : self.sg.order_p() ] # avoid duplicate rotation matrices
      retmsg = "Expanding to P1 in browser\n"
      if not self.miller_array.is_unique_set_under_symmetry():
        retmsg += "Not all reflections are in the same asymmetric unit in reciprocal space.\n"
        retmsg += "Some reflections might be displayed on top of one another.\n"
    else:
      unique_rot_ops = [ self.symops[0] ] # No P1 expansion. So only submit the identity matrix
    if friedel_mate and not self.miller_array.anomalous_flag():
      msgtype += "Friedel"
      retmsg = "Expanding Friedel mates in browser\n"
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
    return retmsg


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
      return "Removing real space unit cell\n"
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
    return "Adding real space unit cell\n"


  def DrawReciprocalUnitCell(self, scale=1):
    if scale is None:
      self.RemoveVectors("reciprocal_unitcell")
      return "Removing reciprocal unit cell\n"
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
    return "Adding reciprocal unit cell\n"


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


  def fix_orientation(self, val):
    if val:
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


  def SetMouseSpeed(self, trackspeed):
    msg = str(trackspeed)
    self.AddToBrowserMsgQueue("SetMouseSpeed", msg)
    #self.GetMouseSpeed() # TODO: fix wait time


  def GetMouseSpeed(self):
    self.ngl_settings.mouse_sensitivity = None
    self.AddToBrowserMsgQueue("GetMouseSpeed", "")
    if self.WaitforHandshake():
      nwait = 0
      while self.ngl_settings.mouse_sensitivity is None and nwait < 5:
        time.sleep(self.sleeptime)
        nwait += self.sleeptime


  def SetClipPlaneDistances(self, near, far, cameraPosZ=None):
    if cameraPosZ is None:
      cameraPosZ = self.cameraPosZ
    msg = str(near) + ", " + str(far) + ", " + str(cameraPosZ)
    self.AddToBrowserMsgQueue("SetClipPlaneDistances", msg)


  def GetClipPlaneDistances(self):
    self.clipNear = None
    self.clipFar = None
    self.cameraPosZ = None
    self.AddToBrowserMsgQueue("GetClipPlaneDistances", "")
    if self.WaitforHandshake():
      nwait = 0
      while self.clipFar is None and nwait < self.handshakewait:
        time.sleep(self.sleeptime)
        nwait += self.sleeptime
      self.mprint("clipnear, clipfar, cameraPosZ: %s, %s %s" \
                 %(self.clipNear, self.clipFar, self.cameraPosZ), 2)
    return (self.clipNear, self.clipFar, self.cameraPosZ)


  def GetBoundingBox(self):
    self.boundingX = 0.0
    self.boundingY = 0.0
    self.boundingZ = 0.0
    self.AddToBrowserMsgQueue("GetBoundingBox", "")
    if self.WaitforHandshake():
      nwait = 0
      while self.boundingX is None and nwait < self.handshakewait:
        time.sleep(self.sleeptime)
        nwait += self.sleeptime
      self.mprint("boundingXYZ: %s, %s %s" \
         %(self.boundingX, self.boundingY, self.boundingZ), verbose=2)
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
    #self.parent.update_settings()


  def DisableMouseRotation(self): # disable rotating with the mouse
    self.send_msg_to_browser("DisableMouseRotation")


  def EnableMouseRotation(self): # enable rotating with the mouse
    self.send_msg_to_browser("EnableMouseRotation")


  def ReOrientStage(self):
    if self.viewmtrx:
      #self.send_msg_to_browser("ReOrient", self.viewmtrx)
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


  def InjectNewReflections(self, proc_array):
    (hklscenes, scenemaxdata,
      scenemindata, scenemaxsigmas,
        sceneminsigmas, scenearrayinfos
     ) = MakeHKLscene(proc_array, 0, copy.deepcopy(self.viewerparams), { } , None)

    strdata = ""
    hklscene = hklscenes[0]
    self.scene = hklscene

    for i,radius in enumerate(hklscene.radii):
      ftuple = (hklscene.points[i][0], hklscene.points[i][1], hklscene.points[i][2],
               hklscene.colors[i][0], hklscene.colors[i][1], hklscene.colors[i][2], radius )
      strdata += "%s,%s,%s,%s,%s,%s,%s," % roundoff(ftuple, 2)
    strdata = strdata[:-1] # avoid the last comma
    self.isinjected = True
    self.AddToBrowserMsgQueue("InjectNewReflections", strdata)



ngl_philstr = """
  mouse_sensitivity = 0.2
    .type = float
  bin_opacities = ""
    .type = str
  tooltip_alpha = 0.85
    .type = float
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

websocket.enableTrace(True)
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


import asyncio
import math
import websockets

async def time(websocket, path):
  x = 0
  for i in range(1000):
    x += 0.2
    alpha =  (math.cos(x) +1.0 )/2.0
    msg = u"alpha, 2, %f" %alpha
    await websocket.send( msg )
    r = (math.cos(x) +1.0 )/2.0
    g = (math.cos(x+1) +1.0 )/2.0
    b = (math.cos(x+2) +1.0 )/2.0
    msg = u"colour, 1, %d, %f, %f, %f" %(i,r,g,b)
    await websocket.send( msg )
    message = await websocket.recv()
    print( message)
    await asyncio.sleep(0.2)

start_server = websockets.serve(time, '127.0.0.1', 7894)
asyncio.get_event_loop().run_until_complete(start_server)
asyncio.get_event_loop().run_forever()

"""
