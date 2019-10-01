
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
from libtbx import easy_mp
import webbrowser, tempfile
from six.moves import range



def has_phil_path(philobj, path):
  return [ e.path for e in philobj.all_definitions() if path in e.path ]


class ArrayInfo:
  def __init__(self, millarr, mprint=sys.stdout.write, fomlabel=None):
    from iotbx.gui_tools.reflections import get_array_description
    data = millarr.data()
    if (isinstance(data, flex.int)):
      data = [e for e in data if e!= display.inanval]
    if millarr.is_complex_array():
      data = flex.abs(millarr.data())
    data = [e for e in data if not math.isnan(e)]
    self.maxdata =max( data )
    self.mindata =min( data )
    self.maxsigmas = self.minsigmas = None
    if millarr.sigmas() is not None:
      data = millarr.sigmas()
      data = [e for e in data if not math.isnan(e)]
      self.maxsigmas =max( data )
      self.minsigmas =min( data )
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
      b = flex.bool([bool(math.isnan(e[0]) + math.isnan(e[1]) + math.isnan(e[2])) for e in hklscene.colors])
      hklscene.colors = hklscene.colors.set_selected(b, (1.0, 1.0, 1.0))
      b = flex.bool([bool(math.isnan(e)) for e in hklscene.radii])
      hklscene.radii = hklscene.radii.set_selected(b, 0.2)
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
    b = flex.bool([bool(math.isnan(e)) for e in phases])
    # replace the nan values with an arbitrary float value
    phases = phases.set_selected(b, 42.4242)
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
      if not tooltipstringsdict.has_key(hkl):
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
    self.viewerparams = None
    self.params = None
    self.miller_array = None
    self.symops = []
    self.sg = None
    self.tooltipstrings = []
    self.tooltipstringsdict = {}
    self.d_min = None
    self.scene = None
    self.merge = False
    self.NGLscriptstr = ""
    self.camera_type = "orthographic"
    self.primitivetype = "SphereBuffer"
    self.script_has_tooltips = False
    self.url = ""
    self.binscenelabel = "Resolution"
    self.colour_scene_id = None
    self.radii_scene_id = None
    #self.scene_id = None
    self.rotation_mx = matrix.identity(3)
    self.rot_recip_zvec = None
    self.rot_zvec = None
    self.meanradius = -1
    self.past = time.time()
    self.orientmessage = None
    self.high_quality = True
    if kwds.has_key('high_quality'):
      self.high_quality = kwds['high_quality']
    self.cameradist = 0.0
    self.clipNear = None
    self.clipFar = None
    self.cameraPosZ = None
    self.boundingX = None
    self.boundingY = None
    self.boundingZ = None
    self.OrigClipNear = None
    self.OrigClipFar = None
    self.cameratranslation = ( 0,0,0 )
    self.angle_x_svec = 0.0
    self.angle_y_svec = 0.0
    self.angle_z_svec = 0.0
    self.angle_z_yzvec = 0.0
    self.angle_y_yzvec = 0.0
    self.angle_y_xyvec = 0.0
    self.angle_x_xyvec = 0.0
    self.unit_h_axis = None
    self.unit_k_axis = None
    self.unit_l_axis = None
    self.normal_hk = None
    self.normal_kl = None
    self.normal_lh = None
    self.isnewfile = False
    self.has_new_miller_array = False
    self.sleeptime = 0.1
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
    self.sceneisdirty = True
    self.hkl_scenes_info = []
    self.match_valarrays = []
    self.array_infostrs = []
    self.array_infotpls = []
    self.binstrs = []
    self.bin_infotpls = []
    self.mapcoef_fom_dict = {}
    self.verbose = 0
    if kwds.has_key('verbose'):
      self.verbose = kwds['verbose']
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
    self.guisocket = None
    if 'guisocket' in kwds:
      self.guisocket = kwds['guisocket']
    self.mprint('Output will be written to \"%s\"\n' \
      'including reference to NGL JavaScript \"%s\"' %(self.hklfname, self.jscriptfname))
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
      self.UseOSBrowser = kwds['UseOSBrowser']
    self.viewmtrx = None
    self.HKLscenesKey = ( 0, False,
                          self.settings.expand_anomalous, self.settings.expand_to_p1  )
    self.msgqueue = []
    self.websockclient = None
    self.handshakewait = 5
    if 'handshakewait' in kwds:
      self.handshakewait = kwds['handshakewait']
    self.lastmsg = "" # "Ready"
    self.browserisopen = False
    self.msgdelim = ":\n"
    self.msgqueuethrd = None
    self.StartWebsocket()


  def __exit__(self, exc_type, exc_value, traceback):
    # not called unless instantiated with a "with hklview_3d ... " statement
    self.server.shutdown()
    if os.path.isfile(self.hklfname):
      os.remove(self.hklfname)


  def SendInfoToGUI(self, mydict):
    if self.guisocket:
      self.guisocket.send( str(mydict).encode("utf-8") )


  def update_settings(self, diff_phil, curphilparam) :
    self.ngl_settings = curphilparam.viewer.NGL
    self.viewerparams = curphilparam.viewer
    self.params = curphilparam
    if has_phil_path(diff_phil, "filename") \
     or has_phil_path(diff_phil, "spacegroup_choice") \
     or has_phil_path(diff_phil, "merge_data") \
     or has_phil_path(diff_phil, "miller_array_operation") \
     or has_phil_path(diff_phil, "scene_id")  \
     or has_phil_path(diff_phil, "nbins") \
     or has_phil_path(diff_phil, "camera_type") \
     or has_phil_path(diff_phil, "bin_scene_label") \
     or has_phil_path(diff_phil, "scene_bin_thresholds") \
     or has_phil_path(diff_phil, "spacegroup_choice") \
     or has_phil_path(diff_phil, "using_space_subgroup") \
     or has_phil_path(diff_phil, "viewer") \
     and ( \
      has_phil_path(diff_phil, "show_data_over_sigma") \
     or has_phil_path(diff_phil, "show_missing") \
     or has_phil_path(diff_phil, "show_only_missing") \
     or has_phil_path(diff_phil, "show_systematic_absences") \
     or has_phil_path(diff_phil, "slice_axis") \
     or has_phil_path(diff_phil, "slice_mode") \
     or has_phil_path(diff_phil, "slice_index") \
     or has_phil_path(diff_phil, "scale") \
     or has_phil_path(diff_phil, "nth_power_scale_radii") \
     or self.settings.inbrowser==False and \
               ( has_phil_path(diff_phil, "expand_anomalous") or \
                has_phil_path(diff_phil, "expand_to_p1") )\
     or has_phil_path(diff_phil, "show_anomalous_pairs") \
      ):
        if curphilparam.viewer.slice_mode and self.settings.inbrowser:
          self.settings.inbrowser = False
        self.sceneisdirty = True
        self.ConstructReciprocalSpace(curphilparam, merge=self.merge)
    msg = ""
    if has_phil_path(diff_phil, "show_missing") \
     or has_phil_path(diff_phil, "show_only_missing") \
     or has_phil_path(diff_phil, "show_systematic_absences"):
      self.binvals = self.calc_bin_thresholds(curphilparam.bin_scene_label, curphilparam.nbins)

    if has_phil_path(diff_phil, "camera_type"):
      self.set_camera_type()

    if self.viewerparams.scene_id >=0:
      if not self.isinjected:
        self.scene = self.HKLscenes[self.viewerparams.scene_id]
      self.DrawNGLJavaScript()
      msg = "Rendered %d reflections\n" % self.scene.points.size()
      if has_phil_path(diff_phil, "fixorientation"):
        self.fix_orientation(curphilparam.viewer.NGL.fixorientation)
      if has_phil_path(diff_phil, "mouse_sensitivity"):
        self.SetTrackBallRotateSpeed(curphilparam.viewer.NGL.mouse_sensitivity)

      if curphilparam.viewer.slice_mode: # explicit slicing
        if curphilparam.viewer.slice_axis=="h": hkl = [1,0,0]
        if curphilparam.viewer.slice_axis=="k": hkl = [0,1,0]
        if curphilparam.viewer.slice_axis=="l": hkl = [0,0,1]
        self.clip_plane_hkl_vector(hkl[0], hkl[1], hkl[2], clipwidth=200,
                         fixorientation = curphilparam.viewer.NGL.fixorientation)
    if self.settings.inbrowser and not curphilparam.viewer.slice_mode:
      msg += self.ExpandInBrowser(P1= self.settings.expand_to_p1,
                            friedel_mate= self.settings.expand_anomalous)

    if curphilparam.clip_plane.clipwidth:
      if  curphilparam.clip_plane.is_real_space_frac_vec:
        self.clip_plane_abc_vector(curphilparam.clip_plane.h, curphilparam.clip_plane.k,
          curphilparam.clip_plane.l, curphilparam.clip_plane.hkldist,
          curphilparam.clip_plane.clipwidth, curphilparam.viewer.NGL.fixorientation,
          curphilparam.clip_plane.is_parallel)
      else:
        self.clip_plane_hkl_vector(curphilparam.clip_plane.h, curphilparam.clip_plane.k,
          curphilparam.clip_plane.l, curphilparam.clip_plane.hkldist,
          curphilparam.clip_plane.clipwidth, curphilparam.viewer.NGL.fixorientation,
          curphilparam.clip_plane.is_parallel)

    msg += self.SetOpacities(curphilparam.viewer.NGL.bin_opacities )
    if has_phil_path(diff_phil, "tooltip_alpha"):
      self.set_tooltip_opacity()
    return msg, curphilparam


  def set_miller_array(self, scene_id=None, merge=None, details=""):
    if scene_id is not None:
      self.viewerparams.scene_id = scene_id
      self.isinjected = False
    if self.viewerparams and self.viewerparams.scene_id >= 0 and self.HKLscenes:
      self.miller_array = self.HKLscenes[self.viewerparams.scene_id].miller_array
      self.scene = self.HKLscenes[self.viewerparams.scene_id]
    self.merge = merge
    if (self.miller_array is None):
      return
    self.identify_suitable_fomsarrays()
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


  def MakeToolTips(self, HKLscenes):
    #import code, traceback; code.interact(local=locals(), banner="".join( traceback.format_stack(limit=10) ) )
    allcolstraliases = "var hk = \'H,K,L: \';\n"
    alltooltipstringsdict = {}
    if self.script_has_tooltips:
      # large data sets will make javascript file very large with risk of crashing browser
      self.mprint( "making tooltips")
      tooltipstringsdict = {}
      for j,hklscene in enumerate(HKLscenes):
        #tooltipstringsdict, colstraliases = MakeTtips(hklscene, j)
        #"""
        if hklscene.isUsingFOMs():
          continue # already have tooltips for the scene without the associated fom
        colstraliases = "\n  var st%d = '\\n%s: ';" %(j, hklscene.work_array.info().label_string() )
        ocolstr = hklscene.work_array.info().label_string()
        if hklscene.work_array.is_complex_array():
          ampl = flex.abs(hklscene.data)
          phases = flex.arg(hklscene.data) * 180.0/math.pi
          # purge nan values from array to avoid crash in fmod_positive()
          b = flex.bool([bool(math.isnan(e)) for e in phases])
          # replace the nan values with an arbitrary float value
          phases = phases.set_selected(b, 42.4242)
          # Cast negative degrees to equivalent positive degrees
          phases = flex.fmod_positive(phases, 360.0)
        sigmas = hklscene.sigmas
        for i,datval in enumerate(hklscene.data):
          hkl = hklscene.indices[i]
          if not tooltipstringsdict.has_key(hkl):
            spbufttip = '\'+hk+\'%s, %s, %s' %(hkl[0], hkl[1], hkl[2])
            spbufttip += '\ndres: %s ' %str(roundoff(hklscene.dres[i], 2) )
            spbufttip += '\'+AA+\'' # javascript alias for angstrom
            tooltipstringsdict[hkl] = spbufttip
          od =""
          if hklscene.work_array.is_complex_array():
            od = str(roundoff(ampl[i], 2)) + ", " + str(roundoff(phases[i], 1)) + \
              "\'+DGR+\'"
          elif sigmas is not None:
            od = str(roundoff(datval, 2)) + ", " + str(roundoff(sigmas[i], 2))
          else:
            od = str(roundoff(datval, 2))
          if not (math.isnan( abs(datval) ) or datval == display.inanval):
            # st1, st2,... are javascript aliases for miller array labelstrings as declared in self.colstraliases
            tooltipstringsdict[hkl] += '\'+st%d+\'%s' %(j, od)
        #"""
        alltooltipstringsdict.update( tooltipstringsdict )
        allcolstraliases += colstraliases
      allcolstraliases += "\n"

    return alltooltipstringsdict, allcolstraliases


  def GetTooltipOnTheFly(self, id, rotmx=None, anomalous=False):
    hkl = self.scene.indices[id]
    hklvec = flex.vec3_double( [(hkl[0], hkl[1], hkl[2])])
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
      if datval and (not (math.isnan( abs(datval) ) or datval == display.inanval)):
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
    spbufttip += '\''
    return spbufttip


  def get_col_fomcol(self, idx):
    if len(self.hkl_scenes_info) == 0:
      return -1, -1
    return self.hkl_scenes_info[idx][6], self.hkl_scenes_info[idx][7]


  def ConstructReciprocalSpace(self, curphilparam, merge=None):
    #self.miller_array = self.match_valarrays[self.scene_id]
    #self.miller_array = self.proc_arrays[self.scene_id]
    self.HKLscenesKey = (curphilparam.filename,
                         curphilparam.spacegroup_choice,
                         curphilparam.using_space_subgroup,
                         curphilparam.merge_data,
                         self.settings.expand_anomalous,
                         self.settings.expand_to_p1,
                         self.settings.inbrowser,
                         self.settings.slice_axis,
                         self.settings.slice_mode,
                         self.settings.slice_index,
                         self.settings.show_missing,
                         self.settings.show_only_missing,
                         self.settings.show_systematic_absences,
                         self.settings.scale,
                         self.settings.nth_power_scale_radii
                         )

    if self.HKLscenesdict.has_key(self.HKLscenesKey) and not self.has_new_miller_array:
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
    # arguments tuple for multi_core_run
    assert(self.proc_arrays)
    argstuples = [ (e.deep_copy(), idx, copy.deepcopy(self.settings), self.mapcoef_fom_dict, merge, self.mprint) \
                     for (idx,e) in enumerate(self.proc_arrays)]
    #import code, traceback; code.interact(local=locals(), banner="".join( traceback.format_stack(limit=10) ) )
    """
    for (i, (args, res, errstr)) in enumerate( easy_mp.multi_core_run( MakeHKLscene, argstuples, 8)):
      if errstr:
        self.mprint(errstr)
      (hkl_scenes, scenemaxdata,
        scenemindata, scenemaxsigmas,
        sceneminsigmas, scenearrayinfos
      ) = res
      HKLscenesMaxdata.extend(scenemaxdata)
      HKLscenesMindata.extend(scenemindata)
      HKLscenesMaxsigmas.extend(scenemaxsigmas)
      HKLscenesMinsigmas.extend(sceneminsigmas)
      hkl_scenes_info.extend(scenearrayinfos)
      HKLscenes.extend(hkl_scenes)
      for inf in scenearrayinfos:
        self.mprint("%d, %s" %(i, inf) )
        i += 1

    """
    for j,proc_array in enumerate(self.proc_arrays):
      (hklscenes, scenemaxdata,
        scenemindata, scenemaxsigmas,
         sceneminsigmas, scenearrayinfos
         ) = MakeHKLscene(argstuples[j][0], argstuples[j][1], argstuples[j][2], argstuples[j][3], argstuples[j][4], argstuples[j][5] )
         #) = MakeHKLscene(proc_array, j, copy.deepcopy(self.settings), self.mapcoef_fom_dict, None)

      HKLscenesMaxdata.extend(scenemaxdata)
      HKLscenesMindata.extend(scenemindata)
      HKLscenesMaxsigmas.extend(scenemaxsigmas)
      HKLscenesMinsigmas.extend(sceneminsigmas)
      hkl_scenes_info.extend(scenearrayinfos)
      HKLscenes.extend(hklscenes)
      #for inf in scenearrayinfos:
      #  self.mprint("%d, %s" %(i, inf) )
      #  i += 1

    tooltipstringsdict, self.colstraliases = self.MakeToolTips(HKLscenes)
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
    self.mprint("\nReflection data scenes:", verbose=0)
    for j,inf in enumerate(hkl_scenes_info):
      self.mprint("%d, %s" %(j, inf[0]), verbose=0)
    self.sceneisdirty = True
    self.SendInfoToGUI({ "hklscenes_arrays": self.hkl_scenes_info, "NewHKLscenes" : True })
    self.has_new_miller_array = False
    return True


  def identify_suitable_fomsarrays(self):
    self.mprint("Matching complex arrays to suitable FOM arrays")
    self.mapcoef_fom_dict = {}
    for proc_array in self.proc_arrays:
      fom_arrays_idx = []
      for i,foms_array in enumerate(self.proc_arrays):
        if not proc_array.is_complex_array() or not foms_array.is_real_array():
          continue
        if proc_array.size() != foms_array.size():
          continue
        if  min(foms_array.data()) < 0.0 or max(foms_array.data()) > 1.0:
          continue
        fom_arrays_idx.append( (foms_array, i) )
      self.mapcoef_fom_dict[proc_array.info().label_string()] = fom_arrays_idx


  def calc_bin_thresholds(self, bin_scene_label, nbins):
    self.binscenelabel = bin_scene_label
    if self.binscenelabel=="Resolution":
      warray = self.HKLscenes[int(self.viewerparams.scene_id)].work_array
      dres = self.HKLscenes[int(self.viewerparams.scene_id)].dres
      uc = warray.unit_cell()
      indices = self.HKLscenes[int(self.viewerparams.scene_id)].indices
      binning = miller.binning( uc, nbins, indices, max(dres), min(dres) )
      binvals = [ binning.bin_d_range(n)[0]  for n in binning.range_all() ]
      binvals = [ e for e in binvals if e != -1.0] # delete dummy limit
      binvals = list( 1.0/flex.double(binvals) )
    else:
      bindata = self.HKLscenes[int(self.binscenelabel)].data.deep_copy()
      selection = flex.sort_permutation( bindata )
      bindata_sorted = bindata.select(selection)
      # get binvals by dividing bindata_sorted with nbins
      binvals = [bindata_sorted[0]] * nbins #
      for i,e in enumerate(bindata_sorted):
        idiv = int(nbins*float(i)/len(bindata_sorted))
        binvals[idiv] = e
    binvals.sort()
    return binvals


  def UpdateBinValues(self, binvals = [] ):
    if binvals:
      binvals.sort()
      self.binvals = binvals
    else: # ensure default resolution interval includes all data by avoiding rounding errors
      self.binvals = [ 1.0/(self.miller_array.d_max_min()[0]*1.001),
                       1.0/(self.miller_array.d_max_min()[1]*0.999) ]


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
var MakeHKL_Axis = function()
{
  // xyz arrows
  shape.addSphere( [0,0,0] , [ 1, 1, 1 ], 0.3, 'Origo');
  //blue-x
  shape.addArrow( %s, %s , [ 0, 0, 1 ], 0.1);
  //green-y
  shape.addArrow( %s, %s , [ 0, 1, 0 ], 0.1);
  //red-z
  shape.addArrow( %s, %s , [ 1, 0, 0 ], 0.1);

  shape.addText( %s, [ 0, 0, 1 ], fontsize, 'h');
  shape.addText( %s, [ 0, 1, 0 ], fontsize, 'k');
  shape.addText( %s, [ 1, 0, 0 ], fontsize, 'l');
};
    """ %(fontsize, str(Hstararrowstart), str(Hstararrowend), str(Kstararrowstart),
          str(Kstararrowend), str(Lstararrowstart), str(Lstararrowend), Hstararrowtxt,
          Kstararrowtxt, Lstararrowtxt)

    # Make colour gradient array used for drawing a bar of colours next to associated values on the rendered html
    mincolourscalar = self.HKLscenesMindata[self.colour_scene_id]
    maxcolourscalar = self.HKLscenesMaxdata[self.colour_scene_id]
    if self.settings.sigma_color:
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
        gradient_type= self.settings.color_scheme) * 255.0)

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
    # Un-binnable data is scene data values where the bin array has no corresponding miller index
    # Just put these in a separate bin and pay attention to the book keeping!
    for ibin in range(self.nbinvalsboundaries+1): # adding the extra bin for un-binnable data
      colours.append([]) # colours and positions are 3 x size of data()
      positions.append([])
      radii2.append([])
      spbufttips.append([])

    def data2bin(d):
      for ibin, binval in enumerate(self.binvalsboundaries):
        if math.isnan(d): # NaN values are un-binnable. Tag them for an additional last bin
          return self.nbinvalsboundaries
        if (ibin+1) == self.nbinvalsboundaries:
          return ibin
        if d > binval and d <= self.binvalsboundaries[ibin+1]:
          return ibin
      #import code, traceback; code.interact(local=locals(), banner="".join( traceback.format_stack(limit=10) ) )
      raise Sorry("Should never get here")

    if nrefls > 0 and self.bindata.size() != points.size():
      raise Sorry("Not the same number of reflections in bin-data and displayed data")

    for i, hklstars in enumerate(points):
      # bin currently displayed data according to the values of another miller array
      ibin = data2bin( self.bindata[i] )
      positions[ibin].extend( roundoff(list(hklstars), 2) )
      colours[ibin].extend( roundoff(list( colors[i] ), 2) )
      radii2[ibin].append( roundoff(radii[i], 2) )
      #spbufttips[ibin].append(self.tooltipstrings[i] )
      if self.script_has_tooltips:
        spbufttips[ibin].append(self.tooltipstringsdict[hkls[i]])
      else:
        spbufttips[ibin].append( i )

    #spherebufferstr = ""
    spherebufferstr = self.colstraliases
    negativeradiistr = ""
    cntbin = 0
    self.binstrs = []
    self.bin_infotpls = []
    for ibin in range(self.nbinvalsboundaries+1):
      mstr =""
      nreflsinbin = len(radii2[ibin])
      if nreflsinbin == 0:
        continue
      bin2 = float("nan"); bin1= float("nan") # indicates un-binned data
      if ibin == self.nbinvalsboundaries:
        mstr= "bin[%d] has %d un-matching reflections with %s in ]%2.3f; %2.3f]" %(cntbin, nreflsinbin, \
                colstr, bin1, bin2)
      if ibin < (self.nbinvalsboundaries-1):
        bin1= self.binvalsboundaries[ibin]
        bin2= self.binvalsboundaries[ibin+1]
        if colstr=="dres":
          bin1= 1.0/self.binvalsboundaries[ibin]
          bin2= 1.0/self.binvalsboundaries[ibin+1]
        mstr= "bin[%d] has %d reflections with %s in ]%2.3f; %2.3f]" %(cntbin, nreflsinbin, \
                colstr, bin1, bin2)
      self.bin_infotpls.append( roundoff((nreflsinbin, bin1, bin2 )) )
      self.binstrs.append(mstr)
      self.mprint(mstr, verbose=0)

      spherebufferstr += "\n// %s\n" %mstr
      if self.script_has_tooltips:
        uncrustttips = str(spbufttips[ibin]).replace('\"', '\'')
        uncrustttips = uncrustttips.replace("\'\'+", "")
        spherebufferstr += "  ttips.push( %s );" %uncrustttips
      else:
        #spherebufferstr += "  ttips.push( [ ] );"
        ttlst = [-1]
        ttlst.extend(spbufttips[ibin])
        spherebufferstr += "  ttips.push( %s );" %str( ttlst )
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
        spherebufferstr += "\n  }, {pointSize: %1.2f})\n" %self.settings.scale
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
  // listen to `hovered` signal to move tooltip around and change its text
  stage.signals.hovered.add(
    function (pickingProxy)
    {
      //tooltip.style.display = "none";
      if (pickingProxy && (Object.prototype.toString.call(pickingProxy.picker) === '[object Array]'  ))
      {
        var cp = pickingProxy.canvasPosition;
    """
    if self.script_has_tooltips:
      spherebufferstr += """
        tooltip.innerText = pickingProxy.picker[pickingProxy.pid];
    """
    else:
      spherebufferstr += """
        var sym_id = -1;
        var hkl_id = -1
        if (pickingProxy.picker.length > 0)
        { // get stored id number of symmetry operator applied to this hkl
          sym_id = pickingProxy.picker[0];
          var ids = pickingProxy.picker.slice(1);
          var is_friedel_mate = 0;
          hkl_id = ids[ pickingProxy.pid % ids.length ];
          if (pickingProxy.pid >= ids.length)
            is_friedel_mate = 1;
        }
        // tell python the id of the hkl and id number of the symmetry operator
        rightnow = timefunc();
        if (rightnow - timenow > 250)
        { // only post every 250 milli second as not to overwhelm python
          WebsockSendMsg( 'tooltip_id: [' + String([sym_id, hkl_id, is_friedel_mate]) + ']' );
          timenow = timefunc();
        }

        if (current_ttip !== "" )
        {
          tooltip.innerText = current_ttip;
    """
    spherebufferstr += """      tooltip.style.bottom = cp.y + 7 + "px";
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
    }
  );


  stage.mouseObserver.signals.dragged.add(
    function ( deltaX, deltaY)
    {
      if (clipFixToCamPosZ === true)
      {
        stage.viewer.parameters.clipNear = origclipnear + (origcameraZpos - stage.viewer.camera.position.z);
        stage.viewer.parameters.clipFar = origclipfar + (origcameraZpos - stage.viewer.camera.position.z);
        stage.viewer.requestRender();
      }
      cvorient = stage.viewerControls.getOrientation().elements;
      msg = String(cvorient);
      WebsockSendMsg('CurrentViewOrientation:\\n' + msg );
    }
  );


  stage.mouseObserver.signals.scrolled.add(
    function (delta)
    {
      if (clipFixToCamPosZ === true)
      {
        stage.viewer.parameters.clipNear = origclipnear + (origcameraZpos - stage.viewer.camera.position.z);
        stage.viewer.parameters.clipFar = origclipfar + (origcameraZpos - stage.viewer.camera.position.z);
        stage.viewer.requestRender();
      }
      cvorient = stage.viewerControls.getOrientation().elements;
      msg = String(cvorient);
      WebsockSendMsg('CurrentViewOrientation:\\n' + msg );
    }
  );


    """

    spherebufferstr += """
  stage.signals.clicked.add(
    function (pickingProxy)
    {
      if (pickingProxy && (Object.prototype.toString.call(pickingProxy.picker) === '[object Array]'  ))
      {
"""
    if self.script_has_tooltips:
      spherebufferstr += "   var innerText = pickingProxy.picker[pickingProxy.pid];\n"
    else:
      spherebufferstr += "        var innerText = pickingProxy.pid;"
    spherebufferstr += """
        WebsockSendMsg( innerText);
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
        rgb = roundoff(val[1], 1)
        gradval = "rgba(%s, %s, %s, %s)" %(rgb[0], rgb[1], rgb[2], alpha)
        if j%10 == 0 or j==len(self.colourgradientvalues)-1 :
          vstr = str( roundoff(val[0], 2) )
        colourgradstr.append([vstr , gradval])

      colourgradstrs += "  colourgradvalarray[%s] = %s;\n" %(g, str(colourgradstr) )
    if blankscene:
      colourscriptstr = ""
    else:
      colourscriptstr = """

  //colourgradvalarrays
  %s

  var ih = 3,
  topr = 35,
  topr2 = 10,
  lp = 10,
  wp = 40,
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
  addDivBox("", topr2, lp, wp3, totalheight, 'rgba(255.0, 255.0, 255.0, 1.0)');

  // print label of the miller array used for colouring
  addDivBox("%s", topr2, lp, wp, 20);

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
    addDivBox("%s", fomtop, fomlp, fomwp, 20);
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
      color: "black",
      fontFamily: "sans-serif",
      fontSize: "smaller",
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


function addDivBox(txt, t, l, w, h, bgcolour='rgba(255.0, 255.0, 255.0, 0.0)')
{
  divbox = createElement("div",
  {
    innerText: txt
  },
  {
    backgroundColor: bgcolour,
    color:  'rgba(0.0, 0.0, 0.0, 1.0)',
    top: t.toString() + "px",
    left: l.toString() + "px",
    width: w.toString() + "px",
    height: h.toString() + "px",
  }
  );
  addElement(divbox);
}

// Microsoft Edge users follow instructions on
// https://stackoverflow.com/questions/31772564/websocket-to-localhost-not-working-on-microsoft-edge
// to enable websocket connection

var pagename = location.pathname.substring(1);
var mysocket = new WebSocket('ws://127.0.0.1:%s/');


function WebsockSendMsg(msg)
{
  try
  {
    mysocket.send(msg);
    mysocket.send( 'Ready ' + pagename + '\\n' );
  }

  catch(err)
  {
    alert('JavaScriptError: ' + err.stack );
    addDivBox("Error!", window.innerHeight - 50, 20, 40, 20, rgba(100.0, 100.0, 100.0, 0.0));
  }
}


mysocket.onopen = function(e)
{
  WebsockSendMsg('%s now connected via websocket to ' + pagename + '\\n');
};

mysocket.onclose = function(e)
{
  WebsockSendMsg('%s now disconnecting from websocket ' + pagename + '\\n');
};

// Log errors to debugger of your browser
mysocket.onerror = function(error)
{
  console.log('WebSocket Error ' + error);
};


window.addEventListener( 'resize',
  function( event ){
      stage.handleResize();
  },
  false
);




var stage;
var shape;
var shapeComp;
var vectorreprs = [];
var vectorshape;
var vectorshapeComps = [];
var repr;
var AA = String.fromCharCode(197); // short for angstrom
var DGR = String.fromCharCode(176); // short for degree symbol
var ttips = [];
var current_ttip = "";
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
var cvorient = new NGL.Matrix4();
var clipFixToCamPosZ = false;
var origclipnear;
var origclipfar;
var origcameraZpos;
var nbins = %s;

function timefunc() {
  var d = new Date();
  var now = d.getTime();
  return now
}

var timenow = timefunc();
var rightnow = timefunc();



///var script=document.createElement('script');
//script.src='https://rawgit.com/paulirish/memory-stats.js/master/bookmarklet.js';
//document.head.appendChild(script);


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


var hklscene = function()
{
  shape = new NGL.Shape('shape');
  vectorshape = new NGL.Shape('vectorshape');
  stage = new NGL.Stage('viewport', { backgroundColor: "grey", tooltip:false,
                                      fogNear: 100, fogFar: 100 });
  stage.setParameters( { cameraType: "%s" } );

  MakeHKL_Axis();

  %s

  shapeComp = stage.addComponentFromObject(shape);
  repr = shapeComp.addRepresentation('buffer');
  shapeComp.autoView();
  repr.update();

  // if some radii are negative draw them with wireframe
  %s

  %s
  stage.viewer.requestRender();
}


try
{
  document.addEventListener('DOMContentLoaded', function() { hklscene() }, false );
}
catch(err)
{
  WebsockSendMsg('JavaScriptError: ' + err.stack );
}

    """ % (self.websockport, self.__module__, self.__module__, cntbin, self.ngl_settings.tooltip_alpha, \
            axisfuncstr, self.camera_type, spherebufferstr, negativeradiistr, colourscriptstr)



    WebsockMsgHandlestr = """
mysocket.onmessage = function (e)
{
  var c,
  si;
  WebsockSendMsg('\\n    Browser: Got ' + e.data ); // tell server what it sent us
  try
  {
    var datval = e.data.split(":\\n");
    //alert('received2:\\n' + datval);
    var msgtype = datval[0];
    var val = datval[1].split(",");

    if (msgtype === "alpha")
    {
      bin = parseInt(val[0]);
      alphas[bin] = parseFloat(val[1]);
      shapebufs[bin].setParameters({opacity: alphas[bin]});
      for (var g=0; g < nrots; g++ )
        br_shapebufs[bin][g].setParameters({opacity: alphas[bin]});
      stage.viewer.requestRender();
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
      stage.viewer.requestRender();
    }

    if (msgtype === "ShowTooltip")
    {
      current_ttip = eval( String(val));
    }

    if (msgtype === "Redraw")
    {
      stage.viewer.requestRender();
    }

    if (msgtype === "ReOrient")
    {
      WebsockSendMsg( 'Reorienting ' + pagename );
      sm = new Float32Array(16);
      //alert('ReOrienting: ' + val)
      for (j=0; j<16; j++)
        sm[j] = parseFloat(val[j]);

      var m = new NGL.Matrix4();
      m.fromArray(sm);
      stage.viewerControls.orient(m);
      stage.viewer.renderer.setClearColor( 0xffffff, 0.01);
      stage.viewer.requestRender();
    }

    if (msgtype === "Reload")
    {
    // refresh browser with the javascript file
      cvorient = stage.viewerControls.getOrientation().elements;
      msg = String(cvorient);
      WebsockSendMsg('OrientationBeforeReload:\\n' + msg );
      WebsockSendMsg( 'Refreshing ' + pagename );
      window.location.reload(true);
    }

    if (msgtype.includes("Expand") )
    {
      WebsockSendMsg( 'Expanding data...' );
      // delete the shapebufs[] that holds the positions[] arrays
      shapeComp.removeRepresentation(repr);
      // remove shapecomp from stage first
      stage.removeComponent(shapeComp);

      br_positions = [];
      br_colours = [];
      br_radii = [];
      br_ttips = [];
      br_shapebufs = [];

      //alert('rotations:\\n' + val);
      // Rotation matrices are concatenated to a string of floats
      // separated by line breaks between each roation matrix
      rotationstrs = datval[1].split("\\n");
      var Rotmat = new NGL.Matrix3();
      var sm = new Float32Array(9);
      var r = new NGL.Vector3();

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
        // if there is only the identity matrix that means no P1 expansion
        for (var rotmxidx=0; rotmxidx < rotationstrs.length; rotmxidx++ )
        {
          if (rotationstrs[rotmxidx] < 1 )
            continue;
          nrots++;

          br_positions[bin].push( [] );
          br_shapebufs[bin].push( [] );
          br_ttips[bin].push( [] );
          br_ttips[bin][rotmxidx] = ttips[bin].slice(); // deep copy with slice()
          br_ttips[bin][rotmxidx][0] = rotmxidx;
          br_positions[bin][rotmxidx] = new Float32Array( csize );

          // convert string of rotation matrix elements into a Matrix3
          var elmstrs = rotationstrs[rotmxidx].split(",");
          //alert('rot' + rotmxidx + ': ' + elmstrs);
          for (j=0; j<9; j++)
            sm[j] = parseFloat(elmstrs[j]);
          Rotmat.fromArray(sm);

          for (var i=0; i<nsize; i++)
          {
            idx= i*3;
            r.x = positions[bin][idx];
            r.y = positions[bin][idx+1];
            r.z = positions[bin][idx+2];

            r.applyMatrix3(Rotmat)

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
      }
      MakeHKL_Axis();

      shapeComp = stage.addComponentFromObject(shape);
      repr = shapeComp.addRepresentation('buffer');

      for (var bin=0; bin<nbins; bin++)
      {
        for (var rotmxidx=0; rotmxidx < nrots; rotmxidx++ )
        {
          br_shapebufs[bin][rotmxidx].setParameters({opacity: alphas[bin]});
        }
      }

      stage.viewer.requestRender();
      WebsockSendMsg( 'Expanded data' );
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

      strs = datval[1].split("\\n");
      var sm = new Float32Array(9);
      var m4 = new NGL.Matrix4();
      var elmstrs = strs[0].split(",");
      //alert('rot: ' + elmstrs);

      for (j=0; j<9; j++)
        sm[j] = parseFloat(elmstrs[j]);

      /* GL matrices are the transpose of the conventional rotation matrices
      m4.set( sm[0], sm[1], sm[2], 0.0,
              sm[3], sm[4], sm[5], 0.0,
              sm[6], sm[7], sm[8], 0.0,
              0.0,   0.0,   0.0,   1.0
      );
      */
      m4.set( sm[0], sm[3], sm[6], 0.0,
              sm[1], sm[4], sm[7], 0.0,
              sm[2], sm[5], sm[8], 0.0,
              0.0,   0.0,   0.0,   1.0
      );
      stage.viewerControls.orient(m4);
      stage.viewer.requestRender();

    }

    if (msgtype === "TranslateHKLpoints")
    {
      WebsockSendMsg( 'Translating HKLs ' + pagename );
      strs = datval[1].split("\\n");
      var sm = new Float32Array(3);
      var elmstrs = strs[0].split(",");
      //alert('trans: ' + elmstrs);
      for (j=0; j<3; j++)
        sm[j] = parseFloat(elmstrs[j]);
      shapeComp.setPosition([ sm[0], sm[1], sm[2] ])
      stage.viewer.requestRender();
    }

    if (msgtype === "AddVector")
    {
      strs = datval[1].split("\\n");
      var r1 = new Float32Array(3);
      var r2 = new Float32Array(3);
      var rgb = new Float32Array(3);
      var elmstrs = strs[0].split(",");
      for (j=0; j<3; j++)
      {
        r1[j] = parseFloat(elmstrs[j]);
        r2[j] = parseFloat(elmstrs[j+3]);
        rgb[j]= parseFloat(elmstrs[j+6]);
      }

      vectorshape.addArrow( r1, r2 , [rgb[0], rgb[1], rgb[2]], 0.15);
      if (elmstrs[6] !== "") {
        var txtR = [ r1[0] + r2[0], r1[1] + r2[1], r1[2] + r2[2] ];
        vectorshape.addText( txtR, [rgb[0], rgb[1], rgb[2]], fontsize/2.0, elmstrs[9] );
      }

      vectorshapeComps.push( stage.addComponentFromObject(vectorshape) );
      vectorreprs.push( vectorshapeComps[vectorshapeComps.length-1].addRepresentation('vectorbuffer') );
      stage.viewer.requestRender();
    }

    if (msgtype === "RemoveAllVectors")
    {
      for (i=0; i<vectorshapeComps.length; i++)
      {
        vectorshapeComps[i].removeRepresentation(vectorreprs[i]);
        stage.removeComponent(vectorshapeComps[i]);
      }
      vectorshapeComps = [];
      vectorreprs = [];
      clipFixToCamPosZ = false;
      stage.viewer.requestRender();
    }

    if (msgtype === "TooltipOpacity")
    {
      strs = datval[1].split("\\n");
      var elmstrs = strs[0].split(",");
      Object.assign(tooltip.style, {
        backgroundColor: "rgba(255, 255, 255, " + elmstrs[0] + " )",
      });
    }

    if (msgtype === "SetTrackBallRotateSpeed")
    {
      strs = datval[1].split("\\n");
      var elmstrs = strs[0].split(",");
      stage.trackballControls.rotateSpeed = parseFloat(elmstrs[0]);
    }

    if (msgtype === "GetTrackBallRotateSpeed")
    {
      msg = String( [stage.trackballControls.rotateSpeed] )
      WebsockSendMsg('ReturnTrackBallRotateSpeed:\\n' + msg );
    }

    if (msgtype === "SetClipPlaneDistances")
    {
      strs = datval[1].split("\\n");
      var elmstrs = strs[0].split(",");
      var near = parseFloat(elmstrs[0]);
      var far = parseFloat(elmstrs[1]);
      origcameraZpos = parseFloat(elmstrs[2]);
      stage.viewer.parameters.clipMode = 'camera';
      stage.viewer.parameters.clipScale = 'absolute';

      if (near >= far )
      { // default to no clipping if near >= far
        stage.viewer.parameters.clipMode = 'scene';
        stage.viewer.parameters.clipScale = 'relative';
        near = -1000;
        far = 1000;
      }
      stage.viewer.parameters.clipNear = near;
      stage.viewer.parameters.clipFar = far;
      origclipnear = near;
      origclipfar = far;
      clipFixToCamPosZ = true;
      stage.viewer.camera.position.z = origcameraZpos;
      stage.viewer.requestRender();
    }

    if (msgtype === "GetClipPlaneDistances")
    {
      msg = String( [stage.viewer.parameters.clipNear,
                     stage.viewer.parameters.clipFar,
                     stage.viewer.camera.position.z] )
      WebsockSendMsg('ReturnClipPlaneDistances:\\n' + msg );
    }

    if (msgtype === "GetBoundingBox")
    {
      msg = String( [stage.viewer.boundingBoxSize.x,
                     stage.viewer.boundingBoxSize.y,
                     stage.viewer.boundingBoxSize.z]
                  )
      WebsockSendMsg('ReturnBoundingBox:\\n' + msg );
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

      MakeHKL_Axis();
      shapeComp = stage.addComponentFromObject(shape);
      repr = shapeComp.addRepresentation('buffer');
      stage.viewer.requestRender();
      WebsockSendMsg('Injected new reflections');
    }

    if (msgtype === "Testing")
    {
      // test something new
      WebsockSendMsg( 'Testing something new ' + pagename );
      /*
      var newradii = radii[0].map(function(element) {
        return element*1.5;
      });
      shapebufs[0].setAttributes({
          radius: newradii
      })
      repr = shapeComp.addRepresentation('buffer');
      stage.viewer.requestRender();
      */
    }

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
      self.GetClipPlaneDistances()
      self.GetBoundingBox()
      self.OrigClipFar = self.clipFar
      self.OrigClipNear = self.clipNear
      self.SetTrackBallRotateSpeed( self.ngl_settings.mouse_sensitivity )
    self.sceneisdirty = False


  def OnConnectWebsocketClient(self, client, server):
    #if not self.websockclient:
    self.websockclient = client
    self.mprint( "Browser connected:" + str( self.websockclient ), verbose=1 )
    #else:
    #  self.mprint( "Unexpected browser connection was rejected" )


  def OnWebsocketClientMessage(self, client, server, message):
    if self.viewerparams.scene_id is None or self.miller_array is None:
      return
    try:
      if message != "":
        if "Orientation" in message:
          self.orientmessage = message
          self.ProcessOrientationMessage()
        else:
          self.mprint( message, verbose=4)
        self.lastmsg = message
      if "JavaScriptError:" in message:
        self.mprint( message, verbose=0)
        #raise Sorry(message)
      if "OrientationBeforeReload:" in message:
        #sleep(0.2)
        self.mprint( "Reorienting client after refresh:" + str( self.websockclient ), verbose=2 )
        if not self.isnewfile:
          self.viewmtrx = self.orientmessage[ self.orientmessage.find("\n") + 1: ]
          #import code, traceback; code.interact(local=locals(), banner="".join( traceback.format_stack(limit=10) ) )
          self.msgqueue.append( ("ReOrient", self.viewmtrx) )
        self.isnewfile = False
      if "ReturnClipPlaneDistances:" in message:
        datastr = message[ message.find("\n") + 1: ]
        lst = datastr.split(",")
        flst = [float(e) for e in lst]
        self.clipNear = flst[0]
        self.clipFar = flst[1]
        self.cameraPosZ = flst[2]
        self.params.clip_plane.clipwidth = None
      if "ReturnBoundingBox:" in message:
        datastr = message[ message.find("\n") + 1: ]
        lst = datastr.split(",")
        flst = [float(e) for e in lst]
        self.boundingX = flst[0]
        self.boundingY = flst[1]
        self.boundingZ = flst[2]
      if "ReturnTrackBallRotateSpeed" in message:
        datastr = message[ message.find("\n") + 1: ]
        lst = datastr.split(",")
        flst = [float(e) for e in lst]
        self.ngl_settings.mouse_sensitivity = flst[0]
      if "tooltip_id:" in message:
        #import code, traceback; code.interact(local=locals(), banner="".join( traceback.format_stack(limit=10) ) )
        sym_id = eval(message.split("tooltip_id:")[1])[0]
        id = eval(message.split("tooltip_id:")[1])[1]
        is_friedel_mate = eval(message.split("tooltip_id:")[1])[2]
        rotmx = None
        if sym_id >= 0 and sym_id < len(self.symops):
          rotmx = self.symops[sym_id].r()
        hkls = self.scene.indices
        if not is_friedel_mate:
          ttip = self.GetTooltipOnTheFly(id, rotmx)
        else:
          # if id > len(hkls) then these hkls are added as the friedel mates during the
          # "if (anoexp)" condition in the javascript code
          id = id % len(hkls)
          ttip = "id: %d" %id
          #ttip = self.GetTooltipOnTheFly(hkls[id], rotmx, anomalous=True)
          ttip = self.GetTooltipOnTheFly(id, rotmx, anomalous=True)
        self.SendMsgToBrowser("ShowTooltip", ttip)
        #import code, traceback; code.interact(local=locals(), banner="".join( traceback.format_stack(limit=10) ) )
    except Exception as e:
      self.mprint( to_str(e) + "\n" + traceback.format_exc(limit=10), verbose=0)


  def ProcessOrientationMessage(self):
    if self.orientmessage is None:
      return

    if self.orientmessage.find("NaN")>=0:
      return

    self.viewmtrx = self.orientmessage[ self.orientmessage.find("\n") + 1: ]
    lst = self.viewmtrx.split(",")
    flst = [float(e) for e in lst]
    ScaleRotMx = matrix.sqr( (flst[0], flst[4], flst[8],
                          flst[1], flst[5], flst[9],
                          flst[2], flst[6], flst[10]
                          )
    )
    self.cameratranslation = (flst[12], flst[13], flst[14])
    self.mprint("translation: %s" %str(roundoff(self.cameratranslation)), verbose=3)
    self.cameradist = math.pow(ScaleRotMx.determinant(), 1.0/3.0)
    self.mprint("distance: %s" %roundoff(self.cameradist), verbose=3)
    self.rotation_mx = ScaleRotMx/self.cameradist
    rotlst = roundoff(self.rotation_mx.elems)
    self.mprint("""Rotation matrix:
  %s,  %s,  %s
  %s,  %s,  %s
  %s,  %s,  %s
    """ %rotlst, verbose=3)
    alllst = roundoff(flst)
    self.mprint("""OrientationMatrix matrix:
  %s,  %s,  %s,  %s
  %s,  %s,  %s,  %s
  %s,  %s,  %s,  %s
  %s,  %s,  %s,  %s
    """ %tuple(alllst), verbose=4)
    self.params.mouse_moved = True
    if self.rotation_mx.is_r3_rotation_matrix():
      angles = self.rotation_mx.r3_rotation_matrix_as_x_y_z_angles(deg=True)
      self.mprint("angles: %s" %str(roundoff(angles)), verbose=3)
      z_vec = flex.vec3_double( [(0,0,1)])
      self.rot_zvec = z_vec * self.rotation_mx
      self.mprint("Rotated cartesian Z direction : %s" %str(roundoff(self.rot_zvec[0])), verbose=3)
      rfracmx = matrix.sqr( self.miller_array.unit_cell().reciprocal().fractionalization_matrix() )
      self.rot_recip_zvec = self.rot_zvec * rfracmx
      self.rot_recip_zvec = (1.0/self.rot_recip_zvec.norm()) * self.rot_recip_zvec
      self.mprint("Rotated reciprocal L direction : %s" %str(roundoff(self.rot_recip_zvec[0])), verbose=3)


  def WaitforHandshake(self, sec):
    nwait = 0
    while not self.websockclient:
      time.sleep(self.sleeptime)
      nwait += self.sleeptime
      if nwait > sec:
        return False
    return True


  def WebBrowserMsgQueue(self):
    try:
      while True:
        nwait = 0.0
        sleep(self.sleeptime)
        #if time.time() - self.past > 1.0:
        #  self.ProcessOrientationMessage() # avoid doing too often as this is expensive
        #self.past = time.time()
        if len(self.msgqueue):
          #print("self.msgqueue: " + str(self.msgqueue))
          pendingmessagetype, pendingmessage = self.msgqueue[0]
          self.SendMsgToBrowser(pendingmessagetype, pendingmessage)
          while not (self.browserisopen and self.websockclient):
            sleep(self.sleeptime)
            nwait += self.sleeptime
            if nwait > self.handshakewait and self.browserisopen:
              self.mprint("ERROR: No handshake from browser! Security settings may have to be adapted", verbose=0 )
              return
              #break
          self.msgqueue.remove( self.msgqueue[0] )
# if the html content is huge the browser will be unresponsive until it has finished
# reading the html content. This may crash this thread. So try restarting this thread until
# browser is ready
    except Exception as e:
      self.mprint( str(e) + ", Restarting WebBrowserMsgQueue\n" \
                         + traceback.format_exc(limit=10), verbose=2)
      self.WebBrowserMsgQueue()


  def SendMsgToBrowser(self, msgtype, msg=""):
    #import code, traceback; code.interact(local=locals(), banner="".join( traceback.format_stack(limit=10) ) )
    #print "self.server.clients: ", self.server.clients
    #print "self.websockclient: ",
    message = u"" + msgtype + self.msgdelim + str(msg)
    if self.websockclient:
      while not ("Ready" in self.lastmsg or "tooltip_id" in self.lastmsg \
        or "CurrentViewOrientation" in self.lastmsg):
        sleep(self.sleeptime)
      self.server.send_message(self.websockclient, message )
    else:
      self.OpenBrowser()


  def StartWebsocket(self):
    self.server = WebsocketServer(self.websockport, host='127.0.0.1')
    if not self.server:
      raise Sorry("Could not connect to web browser")
    self.server.set_fn_new_client(self.OnConnectWebsocketClient)
    self.server.set_fn_message_received(self.OnWebsocketClientMessage)
    self.wst = threading.Thread(target=self.server.run_forever)
    self.wst.daemon = True
    self.wst.start()
    self.msgqueuethrd = threading.Thread(target = self.WebBrowserMsgQueue )
    self.msgqueuethrd.daemon = True
    self.msgqueuethrd.start()


  def set_camera_type(self):
    self.camera_type = self.ngl_settings.camera_type


  def set_tooltip_opacity(self):
    msg = "%f" %self.ngl_settings.tooltip_alpha
    self.SendMsgToBrowser("TooltipOpacity", msg)


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
    self.SendMsgToBrowser("alpha", msg)
    return "Opacity %s set on bin[%s]\n" %(alpha, bin)


  def RedrawNGL(self):
    #self.SendMsgToBrowser("Redraw")
    self.msgqueue.append( ("Redraw", "") )


  def ReloadNGL(self): # expensive as javascript may be several Mbytes large
    self.mprint("Rendering JavaScript...", verbose=1)
    #self.SendMsgToBrowser("Reload")
    self.msgqueue.append( ("Reload", "") )


  def OpenBrowser(self):
    if not self.browserisopen:
      #NGLlibpath = libtbx.env.under_root(os.path.join("modules","cctbx_project","crys3d","hklview","ngl.js") )
      NGLlibpath = libtbx.env.under_root(os.path.join("modules","cctbx_project","crys3d","hklview","ngl.dev.js") )
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
      self.isnewfile = False
      self.browserisopen = True


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
      msg += str_rot + "\n"
    self.msgqueue.append( (msgtype, msg) )
    self.GetBoundingBox() # bounding box changes when the extent of the displayed lattice changes
    return retmsg


  def AddVector(self, s1, s2, s3, t1, t2, t3, isreciprocal=True, label="", r=0, g=0, b=0):
    """
    place vector from {s1, s2, s3] to [t1, t2, t3] with colour r,g,b and label
    """
    uc = self.miller_array.unit_cell()
    vec1 = (s1*self.scene.renderscale, s2*self.scene.renderscale, s3*self.scene.renderscale)
    vec2 = (t1*self.scene.renderscale, t2*self.scene.renderscale, t3*self.scene.renderscale)
    #svec = list(vec)
    if isreciprocal:
      # uc.reciprocal_space_vector() only takes integer miller indices so compute
      # the cartesian coordinates for real valued miller indices with the transpose of the fractionalization matrix
      vec1 = list( vec1 * matrix.sqr(uc.fractionalization_matrix()).transpose() )
      vec2 = list( vec2 * matrix.sqr(uc.fractionalization_matrix()).transpose() )
      svec1 = [ vec1[0], vec1[1], vec1[2] ]
      svec2 = [ vec2[0], vec2[1], vec2[2] ]
    else:
      vec1 = list( vec1 * matrix.sqr(uc.orthogonalization_matrix()) )
      vec2 = list( vec2 * matrix.sqr(uc.orthogonalization_matrix()) )
      vscale = 200.0/uc.volume()
      # TODO: find suitable scale factor for displaying real space vector together with reciprocal vectors
      svec1 = [ vscale*vec1[0], vscale*vec1[1], vscale*vec1[2] ]
      svec2 = [ vscale*vec2[0], vscale*vec2[1], vscale*vec2[2] ]
    self.mprint("cartesian vector is: %s to %s" %(str(roundoff(svec1)), str(roundoff(svec2))), verbose=1)
    svec = [svec2[0]-svec1[0], svec2[1]-svec1[1], svec2[2]-svec1[2] ]
    xyvec = svec[:] # deep copying
    xyvec[2] = 0.0 # projection vector of svec in the xy plane
    xyvecnorm = math.sqrt( xyvec[0]*xyvec[0] + xyvec[1]*xyvec[1] )
    if xyvecnorm > 0.0:
      self.angle_x_xyvec = math.acos( xyvec[0]/xyvecnorm )*180.0/math.pi
      self.angle_y_xyvec = math.acos( xyvec[1]/xyvecnorm )*180.0/math.pi
    else:
      self.angle_x_xyvec = 90.0
      self.angle_y_xyvec = 90.0
    yzvec = svec[:]
    yzvec[0] = 0.0 # projection vector of svec in the yz plane
    yzvecnorm = math.sqrt( yzvec[1]*yzvec[1] + yzvec[2]*yzvec[2] )
    if yzvecnorm > 0.0:
      self.angle_y_yzvec = math.acos( yzvec[1]/yzvecnorm )*180.0/math.pi
      self.angle_z_yzvec = math.acos( yzvec[2]/yzvecnorm )*180.0/math.pi
    else:
      self.angle_y_yzvec = 90.0
      self.angle_z_yzvec = 90.0
    svecnorm = math.sqrt( svec[0]*svec[0] + svec[1]*svec[1] + svec[2]*svec[2] )
    self.angle_x_svec = math.acos( svec[0]/svecnorm )*180.0/math.pi
    self.angle_y_svec = math.acos( svec[1]/svecnorm )*180.0/math.pi
    self.angle_z_svec = math.acos( svec[2]/svecnorm )*180.0/math.pi
    if self.angle_y_svec > 90.0:
      self.angle_x_xyvec = -self.angle_x_xyvec

    self.mprint("angles in xy plane to x,y axis are: %s, %s" %(self.angle_x_xyvec, self.angle_y_xyvec), verbose=2)
    self.mprint("angles in yz plane to y,z axis are: %s, %s" %(self.angle_y_yzvec, self.angle_z_yzvec), verbose=2)
    self.mprint("angles to x,y,z axis are: %s, %s, %s" %(self.angle_x_svec, self.angle_y_svec, self.angle_z_svec ), verbose=2)

    self.msgqueue.append( ("AddVector", "%s, %s, %s, %s, %s, %s, %s, %s, %s, %s" \
         %tuple(svec1 + svec2 + [r, g, b, label]) ))


  def PointVectorPerpendicularToClipPlane(self):
    self.RotateStage(( self.angle_x_xyvec, self.angle_z_svec, 0.0 ))


  def PointVectorParallelToClipPlane(self):
    self.RotateStage(( self.angle_x_xyvec, self.angle_z_svec+90.0, 90.0 ))


  def DrawUnitCell(self):
    self.AddVector(0,0,0, 1,0,0, False, label="200a/V", r=0.5, g=0.8, b=0.8)
    self.AddVector(0,0,0, 0,1,0, False, label="200b/V", r=0.8, g=0.5, b=0.8)
    self.AddVector(0,0,0, 0,0,1, False, label="200c/V", r=0.8, g=0.8, b=0.5)
    self.AddVector(1,0,0, 1,1,0, False, r=0.8, g=0.5, b=0.8)
    self.AddVector(0,1,0, 1,1,0, False, r=0.5, g=0.8, b=0.8)
    self.AddVector(0,0,1, 1,0,1, False, r=0.5, g=0.8, b=0.8)
    self.AddVector(0,0,1, 0,1,1, False, r=0.8, g=0.5, b=0.8)
    self.AddVector(0,1,1, 1,1,1, False, r=0.5, g=0.8, b=0.8)
    self.AddVector(1,0,1, 1,1,1, False, r=0.8, g=0.5, b=0.8)
    self.AddVector(1,0,0, 1,0,1, False, r=0.8, g=0.8, b=0.5)
    self.AddVector(0,1,0, 0,1,1, False, r=0.8, g=0.8, b=0.5)
    self.AddVector(1,1,0, 1,1,1, False, r=0.8, g=0.8, b=0.5)


  def DrawReciprocalUnitCell(self):
    n=2
    self.AddVector(0,0,0, n,0,0, label="2a*", r=0.5, g=0.3, b=0.3)
    self.AddVector(0,0,0, 0,n,0, label="2b*", r=0.3, g=0.5, b=0.3)
    self.AddVector(0,0,0, 0,0,n, label="2c*", r=0.3, g=0.3, b=0.5)
    self.AddVector(n,0,0, n,n,0, r=0.3, g=0.5, b=0.3)
    self.AddVector(0,n,0, n,n,0, r=0.5, g=0.3, b=0.3)
    self.AddVector(0,0,n, n,0,n, r=0.5, g=0.3, b=0.3)
    self.AddVector(0,0,n, 0,n,n, r=0.3, g=0.5, b=0.3)
    self.AddVector(0,n,n, n,n,n, r=0.5, g=0.3, b=0.3)
    self.AddVector(n,0,n, n,n,n, r=0.3, g=0.5, b=0.3)
    self.AddVector(n,0,0, n,0,n, r=0.3, g=0.3, b=0.5)
    self.AddVector(0,n,0, 0,n,n, r=0.3, g=0.3, b=0.5)
    self.AddVector(n,n,0, n,n,n, r=0.3, g=0.3, b=0.5)


  def fix_orientation(self, val):
    if val:
      self.DisableMouseRotation()
    else:
      self.EnableMouseRotation()


  def clip_plane_hkl_vector(self, h, k, l, hkldist=0.0,
             clipwidth=None, fixorientation=True, is_parallel=False):
    # create clip plane that is normal to the reciprocal hkl vector
    if h==0.0 and k==0.0 and l==0.0 or clipwidth <= 0.0:
      self.RemoveNormalVectorToClipPlane()
      return
    self.RemoveAllReciprocalVectors()
    R = -l * self.normal_hk + h * self.normal_kl + k * self.normal_lh
    self.AddVector(0, 0, 0, R[0][0], R[0][1], R[0][2], isreciprocal=False)
    if fixorientation:
      self.DisableMouseRotation()
    else:
      self.EnableMouseRotation()
    if is_parallel:
      self.PointVectorParallelToClipPlane()
    else:
      self.PointVectorPerpendicularToClipPlane()
    halfdist = -self.cameraPosZ  - hkldist # self.viewer.boundingZ*0.5
    if clipwidth is None:
      clipwidth = self.meanradius
    clipNear = halfdist - clipwidth # 50/self.viewer.boundingZ
    clipFar = halfdist + clipwidth  #50/self.viewer.boundingZ
    self.SetClipPlaneDistances(clipNear, clipFar, self.cameraPosZ)
    self.TranslateHKLpoints(R[0][0], R[0][1], R[0][2], hkldist)
    self.DrawUnitCell()
    self.DrawReciprocalUnitCell()


  def clip_plane_abc_vector(self, a, b, c, hkldist=0.0,
             clipwidth=None, fixorientation=True, is_parallel=False):
    # create clip plane that is normal to the realspace fractional abc vector
    if a==0.0 and b==0.0 and c==0.0 or clipwidth <= 0.0:
      self.RemoveNormalVectorToClipPlane()
      return
    self.RemoveAllReciprocalVectors()
    self.AddVector(0, 0, 0, a, b, c, isreciprocal=False)
    if fixorientation:
      self.DisableMouseRotation()
    else:
      self.EnableMouseRotation()
    if is_parallel:
      self.PointVectorParallelToClipPlane()
    else:
      self.PointVectorPerpendicularToClipPlane()
    halfdist = -self.cameraPosZ  - hkldist # self.viewer.boundingZ*0.5
    if clipwidth is None:
      clipwidth = self.meanradius
    clipNear = halfdist - clipwidth # 50/self.viewer.boundingZ
    clipFar = halfdist + clipwidth  #50/self.viewer.boundingZ
    self.SetClipPlaneDistances(clipNear, clipFar, self.cameraPosZ)
    self.TranslateHKLpoints(a,b,c, hkldist)
    self.DrawUnitCell()
    self.DrawReciprocalUnitCell()


  def clip_plane_to_HKL_vector(self, h, k, l, hkldist=0.0,
             clipwidth=None, fixorientation=True):
    if h==0.0 and k==0.0 and l==0.0 or clipwidth==None:
      self.RemoveNormalVectorToClipPlane()
      return
    self.RemoveAllReciprocalVectors()
    self.AddVector(0, 0, 0, h, k, l, isreciprocal=False)
    if fixorientation:
      self.DisableMouseRotation()
    else:
      self.EnableMouseRotation()
    self.PointVectorPerpendicularToClipPlane()
    halfdist = -self.cameraPosZ  - hkldist # self.viewer.boundingZ*0.5
    if clipwidth is None:
      clipwidth = self.meanradius
    clipNear = halfdist - clipwidth # 50/self.viewer.boundingZ
    clipFar = halfdist + clipwidth  #50/self.viewer.boundingZ
    self.SetClipPlaneDistances(clipNear, clipFar, self.cameraPosZ)
    self.TranslateHKLpoints(h,k,l, hkldist)


  def RemoveNormalVectorToClipPlane(self):
    self.EnableMouseRotation()
    self.RemoveAllReciprocalVectors()
    self.SetClipPlaneDistances(0, 0)
    self.TranslateHKLpoints(0, 0, 0, 0.0)


  def SetTrackBallRotateSpeed(self, trackspeed):
    msg = str(trackspeed)
    self.msgqueue.append( ("SetTrackBallRotateSpeed", msg) )
    self.GetTrackBallRotateSpeed()


  def GetTrackBallRotateSpeed(self):
    self.ngl_settings.mouse_sensitivity = None
    self.msgqueue.append( ("GetTrackBallRotateSpeed", "") )
    if self.WaitforHandshake(5):
      nwait = 0
      while self.ngl_settings.mouse_sensitivity is None and nwait < 5:
        time.sleep(self.sleeptime)
        nwait += self.sleeptime


  def SetClipPlaneDistances(self, near, far, cameraPosZ=None):
    if cameraPosZ is None:
      cameraPosZ = self.cameraPosZ
    msg = str(near) + ", " + str(far) + ", " + str(cameraPosZ)
    self.msgqueue.append( ("SetClipPlaneDistances", msg) )


  def GetClipPlaneDistances(self):
    self.clipNear = None
    self.clipFar = None
    self.cameraPosZ = None
    self.msgqueue.append( ("GetClipPlaneDistances", "") )
    if self.WaitforHandshake(5):
      nwait = 0
      while self.clipFar is None and nwait < 5:
        time.sleep(self.sleeptime)
        nwait += self.sleeptime
      self.mprint("clipnear, clipfar, cameraPosZ: %s, %s %s" \
                 %(self.clipNear, self.clipFar, self.cameraPosZ), 2)
    return (self.clipNear, self.clipFar, self.cameraPosZ)


  def GetBoundingBox(self):
    self.boundingX = 0.0
    self.boundingY = 0.0
    self.boundingZ = 0.0
    self.msgqueue.append( ("GetBoundingBox", "") )
    if self.WaitforHandshake(5):
      nwait = 0
      while self.boundingX is None and nwait < 5:
        time.sleep(self.sleeptime)
        nwait += self.sleeptime
      self.mprint("boundingXYZ: %s, %s %s" \
         %(self.boundingX, self.boundingY, self.boundingZ), verbose=2)
    return (self.boundingX, self.boundingY, self.boundingZ)


  def RemoveAllReciprocalVectors(self):
    self.msgqueue.append( ("RemoveAllVectors", "" ))


  def TestNewFunction(self):
    self.SendMsgToBrowser("Testing")


  def DisableMouseRotation(self): # disable rotating with the mouse
    self.SendMsgToBrowser("DisableMouseRotation")


  def EnableMouseRotation(self): # enable rotating with the mouse
    self.SendMsgToBrowser("EnableMouseRotation")


  def RotateStage(self, eulerangles):
    eulerangles1 = eulerangles
    #uc = self.miller_array.unit_cell()
    #OrtMx = matrix.sqr( uc.orthogonalization_matrix())
    #InvMx = OrtMx.inverse()
    #RotMx = matrix.sqr( (1.0, 0.0, 0.0, 0.0, -1.0, 0.0, 0.0, 0.0, -1.0 ) )
    #ortrot = (OrtMx * RotMx * InvMx).as_mat3()
    radangles = [e*math.pi/180.0 for e in eulerangles1]
    RotMx = scitbx.math.euler_angles_as_matrix(radangles)
    scaleRot = RotMx * self.cameradist
    ortrot = scaleRot.as_mat3()
    if scaleRot.determinant() < 1.0:
      self.mprint("Waiting for scaleRot determinant > 0")
      return
      sleep(1)
      self.RotateStage(eulerangles)
      return
    str_rot = str(ortrot)
    str_rot = str_rot.replace("(", "")
    str_rot = str_rot.replace(")", "")
    msg = str_rot + "\n"
    self.msgqueue.append( ("RotateStage", msg) )


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
    self.msgqueue.append( ("TranslateHKLpoints", msg) )


  def InjectNewReflections(self, proc_array):
    (hklscenes, scenemaxdata,
      scenemindata, scenemaxsigmas,
        sceneminsigmas, scenearrayinfos
     ) = MakeHKLscene(proc_array, 0, copy.deepcopy(self.settings), { } , None)

    strdata = ""
    hklscene = hklscenes[0]
    self.scene = hklscene

    for i,radius in enumerate(hklscene.radii):
      ftuple = (hklscene.points[i][0], hklscene.points[i][1], hklscene.points[i][2],
               hklscene.colors[i][0], hklscene.colors[i][1], hklscene.colors[i][2], radius )
      strdata += "%s,%s,%s,%s,%s,%s,%s," % roundoff(ftuple, 2)
    strdata = strdata[:-1] # avoid the last comma
    self.isinjected = True
    self.msgqueue.append( ("InjectNewReflections", strdata) )



ngl_philstr = """
  mouse_sensitivity = 0.2
    .type = float
  bin_opacities = ""
    .type = str
  tooltip_alpha = 0.85
    .type = float
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
