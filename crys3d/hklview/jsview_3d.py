
from __future__ import division
from libtbx.math_utils import roundoff
from cctbx.miller import display2 as display
from cctbx.array_family import flex
from cctbx import miller
from scitbx import graphics_utils
from scitbx import matrix
from libtbx.utils import Sorry, to_str
from websocket_server import WebsocketServer
import threading, math, sys
from time import sleep
import os.path, time, copy
import libtbx
from libtbx import easy_mp
import webbrowser, tempfile



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
    dmin = 0.0
    dmax = 0.0
    try:
      self.span = ( millarr.index_span().min(), millarr.index_span().max())
      dmin = millarr.d_max_min()[1]
      dmax = millarr.d_max_min()[0]
    except Exception, e:
      mprint(to_str(e))
    issymunique = millarr.is_unique_set_under_symmetry()
    self.infotpl = (self.labels, self.desc, millarr.indices().size(), self.span,
     self.minmaxdata, self.minmaxsigs, (roundoff(dmin), roundoff(dmax)), issymunique )
    self.infostr = "%s (%s), %s HKLs: %s, MinMax: %s, MinMaxSigs: %s, d_minmax: %s, SymUnique: %d" %self.infotpl




def MakeHKLscene( proc_array, pidx, setts, mapcoef_fom_dict, merge, mprint=sys.stdout.write):
  scenemaxdata =[]
  scenemindata =[]
  scenemaxsigmas = []
  sceneminsigmas = []
  scenearrayinfos = []
  hklscenes = []
  fomsarrays_idx = [(None, [])]
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
      infostr = ainf.infostr
      scenemaxdata.append( ainf.maxdata )
      scenemindata.append( ainf.mindata )
      scenemaxsigmas.append(ainf.maxsigmas)
      sceneminsigmas.append(ainf.minsigmas)
      scenearrayinfos.append((infostr, pidx, fidx))
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
    self.miller_array = None
    self.tooltipstrings = []
    self.tooltipstringsdict = {}
    self.d_min = None
    self.scene = None
    self.merge = False
    self.NGLscriptstr = ""
    self.cameratype = "orthographic"
    self.primitivetype = "SphereBuffer"
    self.usingtooltips = True
    self.url = ""
    self.binarray = "Resolution"
    self.icolourcol = None
    self.iradiicol = None
    self.iarray = None
    self.isnewfile = False
    self.colstraliases = ""
    self.binvals = []
    self.workingbinvals = []
    self.proc_arrays = []
    self.HKLscenes = []
    self.HKLscenesdict = {}
    self.HKLscenesMaxdata = []
    self.HKLscenesMindata = []
    self.HKLscenesMaxsigmas = []
    self.HKLscenesMinsigmas = []
    self.sceneisdirty = False
    self.hkl_scenes_info = []
    self.match_valarrays = []
    self.binstrs = []
    self.mapcoef_fom_dict = {}
    self.verbose = True
    if kwds.has_key('verbose'):
      self.verbose = kwds['verbose']
    self.mprint = sys.stdout.write
    if kwds.has_key('mprint'):
      self.mprint = kwds['mprint']
    self.nbin = 0
    tempdir = tempfile.gettempdir()
    self.hklfname = os.path.join(tempdir, "hkl.htm" )
    if os.path.isfile(self.hklfname):
      os.remove(self.hklfname)
    if kwds.has_key('htmlfname'):
      self.hklfname = kwds['htmlfname']
    self.hklfname = os.path.abspath( self.hklfname )
    self.jscriptfname = os.path.join(tempdir, "hkljstr.js")
    if os.path.isfile(self.jscriptfname):
      os.remove(self.jscriptfname)
    if kwds.has_key('jscriptfname'):
      self.jscriptfname = kwds['jscriptfname']
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
    self.UseOSBrowser = True
    if kwds.has_key('UseOSBrowser'):
      self.UseOSBrowser = kwds['UseOSBrowser']
    self.viewmtrxelms = None
    self.HKLscenesKey = ( 0, False,
                          self.settings.expand_anomalous, self.settings.expand_to_p1  )
    self.pendingmessage = None
    self.pendingmessagetype = None
    self.websockclient = None
    self.lastmsg = ""
    self.msgdelim = ":\n"
    self.msgqueuethrd = None
    self.StartWebsocket()


  def __exit__(self, exc_type, exc_value, traceback):
    # not called unless instantiated with a "with hklview_3d ... " statement
    self.server.shutdown()
    if os.path.isfile(self.hklfname):
      os.remove(self.hklfname)


  def update_settings(self, diffphil, currentphil) :
    if hasattr(diffphil, "filename") \
      or hasattr(diffphil, "spacegroup_choice") \
      or hasattr(diffphil, "merge_data") \
      or hasattr(diffphil, "column")  \
      or hasattr(diffphil, "spacegroup_choice") \
      or hasattr(diffphil, "using_space_subgroup") \
      or hasattr(diffphil, "viewer") \
      and ( \
       hasattr(diffphil.viewer, "show_data_over_sigma") \
      or hasattr(diffphil.viewer, "show_missing") \
      or hasattr(diffphil.viewer, "show_only_missing") \
      or hasattr(diffphil.viewer, "show_systematic_absences") \
      or hasattr(diffphil.viewer, "slice_axis") \
      or hasattr(diffphil.viewer, "slice_mode") \
      or hasattr(diffphil.viewer, "slice_index") \
      or hasattr(diffphil.viewer, "scale") \
      or hasattr(diffphil.viewer, "nth_power_scale_radii") \
      or hasattr(diffphil.viewer, "expand_anomalous") \
      or hasattr(diffphil.viewer, "expand_to_p1")
      or hasattr(diffphil.viewer, "show_anomalous_pairs") \
      ):
        self.sceneisdirty = True
        if not self.ConstructReciprocalSpace(currentphil, merge=self.merge) or \
         self.miller_array is None or self.iarray < 0:
          return
    #import code, traceback; code.interact(local=locals(), banner="".join( traceback.format_stack(limit=10) ) )
    msg = ""
    if self.iarray >=0:
      self.scene = self.HKLscenes[self.iarray]
      self.DrawNGLJavaScript()
      msg = "Rendered %d reflections\n" % self.scene.points.size()
    return msg


  def set_miller_array(self, col=None, merge=None, details=""):
    if col is not None:
      self.iarray = col
    if self.iarray >= 0:
      self.miller_array = self.HKLscenes[self.iarray].miller_array
      self.scene = self.HKLscenes[self.iarray]
    self.merge = merge
    if (self.miller_array is None):
      return
    self.identify_suitable_fomsarrays()
    self.d_min = self.miller_array.d_min()
    array_info = self.miller_array.info()
    self.binvals = [ 1.0/self.miller_array.d_max_min()[0], 1.0/self.miller_array.d_max_min()[1]  ]
    #import code, traceback; code.interact(local=locals(), banner="".join( traceback.format_stack(limit=10) ) )
    uc = "a=%g b=%g c=%g angles=%g,%g,%g" % self.miller_array.unit_cell().parameters()
    self.mprint( "Data: %s %s, %d reflections in space group: %s, unit Cell: %s" \
      % (array_info.label_string(), details, self.miller_array.indices().size(), \
          self.miller_array.space_group_info(), uc) )


  def MakeToolTips(self, HKLscenes):
    #import code, traceback; code.interact(local=locals(), banner="".join( traceback.format_stack(limit=10) ) )
    allcolstraliases = "\n  var hk = \'H,K,L: \'"
    alltooltipstringsdict = {}
    if self.usingtooltips:
      self.mprint( "making tooltips")
      tooltipstringsdict = {}
      for j,hklscene in enumerate(HKLscenes):
        #tooltipstringsdict, colstraliases = MakeTtips(hklscene, j)
        #"""
        if hklscene.isUsingFOMs():
          continue # already have tooltips for the scene without the associated fom
        colstraliases = "\n  var st%d = '\\n%s: '" %(j, hklscene.work_array.info().label_string() )
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


  def get_col_fomcol(self, idx):
    if len(self.hkl_scenes_info) == 0:
      return -1, -1
    return self.hkl_scenes_info[idx][6], self.hkl_scenes_info[idx][7]


  def ConstructReciprocalSpace(self, currentphil, merge=None):
    self.mprint("Constructing HKL scenes")
    #self.miller_array = self.match_valarrays[self.iarray]
    #self.miller_array = self.proc_arrays[self.iarray]
    if not self.sceneisdirty:
      self.mprint("Scene is clean", verbose=True)
      return

    self.HKLscenesKey = (currentphil.filename,
                         currentphil.spacegroup_choice,
                         currentphil.using_space_subgroup,
                         currentphil.merge_data,
                         self.settings.expand_anomalous,
                         self.settings.expand_to_p1,
                         self.settings.slice_axis,
                         self.settings.slice_mode,
                         self.settings.slice_index,
                         self.settings.show_missing,
                         self.settings.show_only_missing,
                         self.settings.show_systematic_absences,
                         self.settings.scale,
                         self.settings.nth_power_scale_radii
                         )
    if self.HKLscenesdict.has_key(self.HKLscenesKey):
      (
        self.HKLscenes,
        self.tooltipstringsdict,
        self.HKLscenesMaxdata,
        self.HKLscenesMindata,
        self.HKLscenesMaxsigmas,
        self.HKLscenesMinsigmas,
        self.hkl_scenes_info
      ) =  self.HKLscenesdict[self.HKLscenesKey]
      self.mprint("Scene key is already present", verbose=True)
      return True

    HKLscenes = []
    HKLscenesMaxdata = []
    HKLscenesMindata = []
    HKLscenesMaxsigmas = []
    HKLscenesMinsigmas = []
    hkl_scenes_info = []
    tooltipstringsdict = {}
    #import code, traceback; code.interact(local=locals(), banner="".join( traceback.format_stack(limit=10) ) )
    #self.scene = display.scene(miller_array=self.miller_array, merge=merge,
    # settings=self.settings, foms_array=foms_array)
    # compute scenes for each of the miller arrays
    #for i,match_valarray in enumerate(self.match_valarrays):
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
         #) = MakeHKLscene(proc_array, copy.deepcopy(self.settings), self.mapcoef_fom_dict, merge)

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
    self.sceneisdirty = False
    for j,inf in enumerate(hkl_scenes_info):
      self.mprint("%d, %s" %(j, inf[0]) )

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


  def UpdateBinValues(self, binvals = [] ):
    if binvals:
      self.binvals = binvals
    else:
      self.binvals = [ 1.0/self.miller_array.d_max_min()[0], 1.0/self.miller_array.d_max_min()[1] ]


  def DrawNGLJavaScript(self):
    if not self.scene:
      return
    if self.miller_array is None :
      self.mprint( "A miller array must be selected for rendering the reflections" )
      return
    self.mprint("Composing NGL JavaScript...")
    h_axis = self.scene.axes[0]
    k_axis = self.scene.axes[1]
    l_axis = self.scene.axes[2]
    nrefls = self.scene.points.size()

    Hstararrowstart = roundoff( [-h_axis[0]*100, -h_axis[1]*100, -h_axis[2]*100] )
    Hstararrowend = roundoff( [h_axis[0]*100, h_axis[1]*100, h_axis[2]*100] )
    Hstararrowtxt  = roundoff( [h_axis[0]*102, h_axis[1]*102, h_axis[2]*102] )
    Kstararrowstart = roundoff( [-k_axis[0]*100, -k_axis[1]*100, -k_axis[2]*100] )
    Kstararrowend = roundoff( [k_axis[0]*100, k_axis[1]*100, k_axis[2]*100] )
    Kstararrowtxt  = roundoff( [k_axis[0]*102, k_axis[1]*102, k_axis[2]*102] )
    Lstararrowstart = roundoff( [-l_axis[0]*100, -l_axis[1]*100, -l_axis[2]*100] )
    Lstararrowend = roundoff( [l_axis[0]*100, l_axis[1]*100, l_axis[2]*100] )
    Lstararrowtxt  = roundoff( [l_axis[0]*102, l_axis[1]*102, l_axis[2]*102] )
    # make arrow font size roughly proportional to radius of highest resolution shell
    fontsize = str(1.0 + roundoff(math.pow( max(self.miller_array.index_span().max()), 1.0/3.0)))

    arrowstr = """
    // xyz arrows
    shape.addSphere( [0,0,0] , [ 1, 1, 1 ], 0.3, 'Origo');
    //blue-x
    shape.addArrow( %s, %s , [ 0, 0, 1 ], 0.1);
    //green-y
    shape.addArrow( %s, %s , [ 0, 1, 0 ], 0.1);
    //red-z
    shape.addArrow( %s, %s , [ 1, 0, 0 ], 0.1);

    shape.addText( %s, [ 0, 0, 1 ], %s, 'h');
    shape.addText( %s, [ 0, 1, 0 ], %s, 'k');
    shape.addText( %s, [ 1, 0, 0 ], %s, 'l');
    """ %(str(Hstararrowstart), str(Hstararrowend), str(Kstararrowstart), str(Kstararrowend),
          str(Lstararrowstart), str(Lstararrowend), Hstararrowtxt, fontsize,
          Kstararrowtxt, fontsize, Lstararrowtxt, fontsize)

    # Make colour gradient array used for drawing a bar of colours next to associated values on the rendered html
    mincolourscalar = self.HKLscenesMindata[self.icolourcol]
    maxcolourscalar = self.HKLscenesMaxdata[self.icolourcol]
    if self.settings.sigma_color:
      mincolourscalar = self.HKLscenesMinsigmas[self.icolourcol]
      maxcolourscalar = self.HKLscenesMaxsigmas[self.icolourcol]
    span = maxcolourscalar - mincolourscalar
    ln = 60
    incr = span/ln
    colourgradarrays = []
    val = mincolourscalar
    colourscalararray =flex.double()
    colourscalararray.append( val )
    for j,sc in enumerate(range(ln)):
      val += incr
      colourscalararray.append( val )
    if self.HKLscenes[self.icolourcol].miller_array.is_complex_array():
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
      if self.HKLscenes[self.icolourcol].isUsingFOMs():
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

    colors = self.HKLscenes[self.icolourcol].colors
    radii = self.HKLscenes[self.iradiicol].radii
    points = self.scene.points
    hkls = self.scene.indices
    dres = self.scene.dres
    colstr = self.scene.miller_array.info().label_string()
    data = self.scene.data
    colourlabel = self.HKLscenes[self.icolourcol].colourlabel
    fomlabel = self.HKLscenes[self.icolourcol].fomlabel
    #import code, traceback; code.interact(local=locals(), banner="".join( traceback.format_stack(limit=10) ) )
    assert (colors.size() == radii.size() == nrefls)
    colours = []
    positions = []
    radii2 = []
    spbufttips = []

    self.workingbinvals = []
    if not self.binarray=="Resolution":
      ibinarray= int(self.binarray)
      self.workingbinvals = [ self.HKLscenesMindata[ibinarray] - 0.1 , self.HKLscenesMaxdata[ibinarray] + 0.1 ]
      self.workingbinvals.extend( self.binvals )
      self.workingbinvals.sort()
      if self.workingbinvals[0] < 0.0:
         self.workingbinvals.append(0.0)
         self.workingbinvals.sort()
      bindata = self.HKLscenes[ibinarray].data
      if self.HKLscenes[ibinarray].work_array.is_complex_array():
        bindata = self.HKLscenes[ibinarray].ampl
    else:
      self.workingbinvals = self.binvals
      colstr = "dres"
      bindata = 1.0/dres
    self.nbin = len(self.workingbinvals)

    for ibin in range(self.nbin):
      colours.append([]) # colours and positions are 3 x size of data()
      positions.append([])
      radii2.append([])
      spbufttips.append([])

    def data2bin(d):
      for ibin, binval in enumerate(self.workingbinvals):
        if (ibin+1) == self.nbin:
          return ibin
        if d > binval and d <= self.workingbinvals[ibin+1]:
          return ibin
      #import code, traceback; code.interact(local=locals(), banner="".join( traceback.format_stack(limit=10) ) )
      raise Sorry("Should never get here")

    for i, hklstars in enumerate(points):
      # bin currently displayed data according to the values of another miller array
      ibin = data2bin( bindata[i] )
      positions[ibin].extend( roundoff(list(hklstars), 2) )
      colours[ibin].extend( roundoff(list( colors[i] ), 2) )
      radii2[ibin].append( roundoff(radii[i], 2) )
      #spbufttips[ibin].append(self.tooltipstrings[i] )
      spherebufferstr = ""
      if self.usingtooltips:
        spbufttips[ibin].append(self.tooltipstringsdict[hkls[i]])

    spherebufferstr += self.colstraliases
    negativeradiistr = ""
    cntbin = 0
    self.binstrs = []
    for ibin in range(self.nbin):
      mstr =""
      nreflsinbin = len(radii2[ibin])
      if (ibin+1) < self.nbin and nreflsinbin > 0:
        bin1= self.workingbinvals[ibin]
        bin2= self.workingbinvals[ibin+1]
        if colstr=="dres":
          bin1= 1.0/self.workingbinvals[ibin]
          bin2= 1.0/self.workingbinvals[ibin+1]
        mstr= "bin:%d, %d reflections with %s in ]%2.2f; %2.2f]" %(cntbin, nreflsinbin, \
                colstr, bin1, bin2)
        self.binstrs.append(mstr)
        self.mprint(mstr, verbose=True)

        spherebufferstr += "// %s\n" %mstr
        if self.usingtooltips:
          uncrustttips = str(spbufttips[ibin]).replace('\"', '\'')
          uncrustttips = uncrustttips.replace("\'\'+", "")
          spherebufferstr += "  ttips.push( %s );" %uncrustttips
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
        if self.usingtooltips:
          spherebufferstr += "\n    picking: ttips[%d]," %cntbin
        if self.primitivetype == "PointBuffer":
          spherebufferstr += "\n  }, {pointSize: %1.2f})\n" %self.settings.scale
        else:
          spherebufferstr += """\n  })
  //}, { disableImpostor: true // to enable changing sphereDetail
  //, sphereDetail: 0 }) // rather than default value of 2 icosahedral subdivisions
  //}, { disableImpostor: true }) // if true allows wireframe spheres but does not allow resizing spheres
  );
  """
        spherebufferstr += "  shape.addBuffer(shapebufs[%d])" %cntbin

        if self.workingbinvals[ibin] < 0.0:
          negativeradiistr += "shapebufs[%d].setParameters({metalness: 1})\n" %cntbin
        cntbin += 1

    if self.usingtooltips:
      spherebufferstr += """

// create tooltip element and add to the viewer canvas
  tooltip = document.createElement("div");
  Object.assign(tooltip.style, {
    display: "none",
    position: "absolute",
    zIndex: 10,
    pointerEvents: "none",
    backgroundColor: "rgba(255, 255, 255, 0.75)",
    color: "black",
    padding: "0.1em",
    fontFamily: "sans-serif"
  });

  stage.viewer.container.appendChild(tooltip);
  // listen to `hovered` signal to move tooltip around and change its text
  stage.signals.hovered.add(function (pickingProxy) {
    if (pickingProxy && (Object.prototype.toString.call(pickingProxy.picker) === '[object Array]'  )){
      var sphere = pickingProxy.sphere;
      var cp = pickingProxy.canvasPosition;
      tooltip.innerText = pickingProxy.picker[pickingProxy.pid];
      tooltip.style.bottom = cp.y + 7 + "px";
      tooltip.style.left = cp.x + 8 + "px";
      tooltip.style.fontSize = "smaller";
      tooltip.style.display = "block";
    }else{
      tooltip.style.display = "none";
    }
  });

    """

    spherebufferstr += """
  stage.signals.clicked.add(function (pickingProxy) {
  if (pickingProxy && (Object.prototype.toString.call(pickingProxy.picker) === '[object Array]'  )){
    var innerText = pickingProxy.picker[pickingProxy.pid];
    mysocket.send( innerText);
  }
  });

    """
    colourgradstrs = "colourgradvalarray = new Array(%s)\n" %fomln
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

      colourgradstrs += "  colourgradvalarray[%s] = %s\n" %(g, str(colourgradstr) )

    #negativeradiistr = ""
    #for ibin in range(self.nbin):
    #  if self.workingbinvals[ibin] < 0.0:
    #    negativeradiistr += "shapebufs[%d].setParameters({metalness: 1})\n" %ibin



    self.NGLscriptstr = """
// Microsoft Edge users follow instructions on
// https://stackoverflow.com/questions/31772564/websocket-to-localhost-not-working-on-microsoft-edge
// to enable websocket connection

var pagename = location.pathname.substring(1);
var mysocket = new WebSocket('ws://127.0.0.1:7894/');

mysocket.onopen = function (e) {
  mysocket.send('%s now connected via websocket to ' + pagename + '\\n');
};

mysocket.onclose = function (e) {
  mysocket.send('%s now disconnecting from websocket ' + pagename + '\\n');
};

// Log errors to debugger of your browser
mysocket.onerror = function (error) {
  console.log('WebSocket Error ' + error);
};


window.addEventListener( 'resize', function( event ){
    stage.handleResize();
}, false );




var stage;
var shape;
var shapeComp;
var repr;
var AA = String.fromCharCode(197); // short for angstrom
var DGR = String.fromCharCode(176); // short for degree symbol
var ttips = [];
var positions = [];
var br_positions = [];
var br_colours = [];
var br_radii = [];
var colours = [];
var radii = [];
var shapebufs = [];
var br_shapebufs = [];


function createElement (name, properties, style) {
// utility function used in for loop over colourgradvalarray
  var el = document.createElement(name)
  Object.assign(el, properties)
  Object.assign(el.style, style)
  Object.assign(el.style, {
      display: "block",
      position: "absolute",
      color: "black",
      fontFamily: "sans-serif",
      fontSize: "smaller",
    }
  )
  return el
}


function addElement (el) {
// utility function used in for loop over colourgradvalarray
  Object.assign(el.style, {
    position: "absolute",
    zIndex: 10
  })
  stage.viewer.container.appendChild(el)
}


function addDivBox(txt, t, l, w, h, bgcolour='rgba(255.0, 255.0, 255.0, 0.0)') {
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
  addElement(divbox)
}


var hklscene = function () {
  shape = new NGL.Shape('shape');
  stage = new NGL.Stage('viewport', { backgroundColor: "grey", tooltip:false,
                                      fogNear: 100, fogFar: 100 });
  stage.setParameters( { cameraType: "%s" } );

  %s

  %s

  shapeComp = stage.addComponentFromObject(shape);
  repr = shapeComp.addRepresentation('buffer');
  shapeComp.autoView();
  repr.update()

  // if some radii are negative draw them with wireframe
  %s

  //colourgradvalarrays
  %s

  var ih = 3;
  var topr = 35
  var topr2 = 10
  var lp = 10
  var wp = 40
  var lp2 = lp + wp
  var gl = 3
  var wp2 = gl
  var fomlabelheight = 25
  if (colourgradvalarray.length === 1) {
    wp2 = 15
    fomlabelheight = 0
  }

  var wp3 = wp + colourgradvalarray.length * wp2 + 2

  totalheight = ih*colourgradvalarray[0].length + 35 + fomlabelheight
  // make a white box on top of which boxes with transparent background are placed
  // containing the colour values at regular intervals as well as label legend of
  // the displayed miller array
  addDivBox("", topr2, lp, wp3, totalheight, 'rgba(255.0, 255.0, 255.0, 1.0)');

  // print label of the miller array used for colouring
  addDivBox("%s", topr2, lp, wp, 20);

  if (colourgradvalarray.length > 1) {
    // print FOM label, 1, 0.5 and 0.0 values below colour chart
    fomtop = topr2 + totalheight - 18
    fomlp = lp + wp
    fomwp = wp3
    fomtop2 = fomtop - 13
    // print the 1 number
    addDivBox("1", fomtop2, fomlp, fomwp, 20)
    // print the 0.5 number
    leftp = fomlp + 0.48 * gl * colourgradvalarray.length
    addDivBox("0.5", fomtop2, leftp, fomwp, 20)
    // print the FOM label
    addDivBox("%s", fomtop, fomlp, fomwp, 20);
    // print the 0 number
    leftp = fomlp + 0.96 * gl * colourgradvalarray.length
    addDivBox("0", fomtop2, leftp, fomwp, 20)
  }

  for (j = 0; j < colourgradvalarray[0].length; j++) {
    rgbcol = colourgradvalarray[0][j][1];
    val = colourgradvalarray[0][j][0]
    topv = j*ih + topr
    toptxt = topv - 5
    // print value of miller array if present in colourgradvalarray[0][j][0]
    addDivBox(val, toptxt, lp, wp, ih);
  }

  // draw the colour gradient
  for (g = 0; g < colourgradvalarray.length; g++) {
    leftp = g*gl + lp + wp

    // if FOM values are supplied draw colour gradients with decreasing
    // saturation values as stored in the colourgradvalarray[g] arrays
    for (j = 0; j < colourgradvalarray[g].length; j++) {
      rgbcol = colourgradvalarray[g][j][1];
      val = colourgradvalarray[g][j][0]
      topv = j*ih + topr
      addDivBox("", topv, leftp, wp2, ih, rgbcol);
    }
  }

}

document.addEventListener('DOMContentLoaded', function() { hklscene() }, false );


mysocket.onmessage = function (e) {
  //alert('received:\\n' + e.data);
  var c
  var alpha
  var si
  mysocket.send('got ' + e.data ); // tell server what it sent us
  try {
    var datval = e.data.split(":\\n");
    //alert('received2:\\n' + datval);
    var msgtype = datval[0];
    //alert('received3:\\n' + msgtype);
    var val = datval[1].split(",");

    if (msgtype === "alpha")
    {
      ibin = parseInt(val[0])
      alpha = parseFloat(val[1])
      shapebufs[ibin].setParameters({opacity: alpha})
      stage.viewer.requestRender()
    }

    if (msgtype === "colour")
    {
      ibin = parseInt(val[0])
      si =  parseInt(val[1])
      colours[ibin][3*si] = parseFloat(val[2])
      colours[ibin][3*si+1] = parseFloat(val[3])
      colours[ibin][3*si+2] = parseFloat(val[4])
      shapebufs[ibin].setAttributes({ color: colours[ibin] })
      stage.viewer.requestRender()
    }

    if (msgtype === "Redraw")
    {
      stage.viewer.requestRender()
    }

    if (msgtype === "ReOrient")
    {
      mysocket.send( 'Reorienting ' + pagename );
      sm = new Float32Array(16);
      //alert('ReOrienting: ' + val)
      for (j=0; j<16; j++)
        sm[j] = parseFloat(val[j])

      var m = new NGL.Matrix4();
      m.fromArray(sm);
      stage.viewerControls.orient(m);
      stage.viewer.requestRender();
    }

    if (msgtype === "Reload")
    {
    // refresh browser with the javascript file
      cvorient = stage.viewerControls.getOrientation().elements
      msg = String(cvorient)
      mysocket.send('Current vieworientation:\\n' + msg );

      mysocket.send( 'Refreshing ' + pagename );
      window.location.reload(true);
    }


    if (msgtype.includes("Expand") !== -1  )
    {
      // test something new
      mysocket.send( 'Expanding data...' );

      if (br_positions.length > 0)
      {
        mysocket.send('Data has already been expanded' );
        return;
      }

      var nsize = positions[0].length/3;

      shapeComp.removeRepresentation(repr);

      //alert('rotations:\\n' + val)
      strs = datval[1].split("\\n");

      nrots = strs.length;
      bins = %d -1;

      var Rotmat = new NGL.Matrix3();
      var sm = new Float32Array(9);
      var r = new NGL.Vector3();

      var csize = nsize*3;
      var nsize3 = nsize*3;
      var anoexp = false;

      if (msgtype.includes("Friedel") !== -1  )
      {
        anoexp = true;
        csize = nsize*6;
      }

      for (var h=0; h<bins; h++)
      {
        br_positions.push( [] )
        br_shapebufs.push( [] )
        br_colours.push( [] )
        br_radii.push( [] )

        br_colours[h] = colours[h];
        br_radii[h] = radii[h]
        if (anoexp)
        {
          var colarr = [];
          var cl = colours[h].length;
          for (var i=0; i<cl; i++)
          {
            colarr[i] = colours[h][i];
            colarr[i+cl] = colours[h][i];
          }
          br_colours[h] = new Float32Array(colarr);

          var radiiarr = [];
          var rl = radii[h].length;
          for (var i=0; i<rl; i++)
          {
            radiiarr[i] = radii[h][i];
            radiiarr[i+rl] = radii[h][i];
          }
          br_radii[h] = new Float32Array(radiiarr);
        }


        for (var g=0; g < nrots; g++ )
        {
          if (strs[g] < 1 )
            continue
          br_positions[h].push( [] )
          br_shapebufs[h].push( [] )
          var elmstrs = strs[g].split(",");

          for (j=0; j<9; j++)
            sm[j] = parseFloat(elmstrs[j])
          Rotmat.fromArray(sm);
          //alert('rot' + g + ': ' + elmstrs)

          br_positions[h][g] = new Float32Array( csize );

          for (var i=0; i<nsize; i++)
          {
            idx= i*3;
            r.x = positions[h][idx];
            r.y = positions[h][idx+1];
            r.z = positions[h][idx+2];

            r.applyMatrix3(Rotmat)

            br_positions[h][g][idx] = r.x
            br_positions[h][g][idx + 1] = r.y
            br_positions[h][g][idx + 2] = r.z

            if (anoexp)
            {
              r.negate(); // inversion for anomalous pair
              br_positions[h][g][nsize3 + idx] = r.x
              br_positions[h][g][nsize3 + idx + 1] = r.y
              br_positions[h][g][nsize3 + idx + 2] = r.z
            }

          }

          br_shapebufs[h][g] = new NGL.SphereBuffer({
              position: br_positions[h][g],
              color: br_colours[h],
              radius: br_radii[h],
            });

          shape.addBuffer(br_shapebufs[h][g]);

        }
      }

      repr = shapeComp.addRepresentation('buffer');

      stage.viewer.requestRender()
    }






    if (msgtype === "Testing") {
      // test something new
      mysocket.send( 'Testing something new ' + pagename );

      var newradii = radii[0].map(function(element) {
        return element*1.5;
      });
      shapebufs[0].setAttributes({
          radius: newradii
      })

      repr = shapeComp.addRepresentation('buffer');

      stage.viewer.requestRender()
    }







  }  catch(err) {
    mysocket.send('error: ' + err );
  }
};


    """ % (self.__module__, self.__module__, self.cameratype, arrowstr, spherebufferstr, \
            negativeradiistr, colourgradstrs, colourlabel, fomlabel, self.nbin)
    if self.jscriptfname:
      with open( self.jscriptfname, "w") as f:
        f.write( self.NGLscriptstr )
    self.ReloadNGL()


  def OnConnectWebsocketClient(self, client, server):
    #if not self.websockclient:
    self.websockclient = client
    self.mprint( "Browser connected:" + str( self.websockclient ), self.verbose )
    #else:
    #  self.mprint( "Unexpected browser connection was rejected" )



  def OnWebsocketClientMessage(self, client, server, message):
    if message != "":
      self.mprint( message, verbose=False)
    self.lastmsg = message
    if "Current vieworientation:" in message:
      # The NGL.Matrix4 with the orientation is a list of floats.
      self.viewmtrxelms = message[ message.find("\n") + 1: ]
      sleep(0.2)
      self.mprint( "Reorienting client after refresh:" + str( self.websockclient ) )
      if not self.isnewfile:
        #import code, traceback; code.interact(local=locals(), banner="".join( traceback.format_stack(limit=10) ) )
        self.pendingmessagetype = "ReOrient"
        self.pendingmessage = self.viewmtrxelms
      self.isnewfile = False


  def WebBrowserMsgQueue(self):
    try:
      while True:
        sleep(1)
        if self.pendingmessagetype:
          self.SendWebSockMsg(self.pendingmessagetype, self.pendingmessage)
          self.pendingmessage = None
          self.pendingmessagetype = None
# if the html content is huge the browser will be unresponsive until it has finished
# reading the html content. This may crash this thread. So try restarting this thread until
# browser is ready
    except Exception, e:
      self.mprint( str(e) + ", Restarting WebBrowserMsgQueue", self.verbose)
      self.WebBrowserMsgQueue()


  def StartWebsocket(self):
    self.server = WebsocketServer(7894, host='127.0.0.1')
    if not self.server:
      raise Sorry("Could not connect socket to web browser")
    self.server.set_fn_new_client(self.OnConnectWebsocketClient)
    self.server.set_fn_message_received(self.OnWebsocketClientMessage)
    self.wst = threading.Thread(target=self.server.run_forever)
    self.wst.daemon = True
    self.wst.start()
    self.msgqueuethrd = threading.Thread(target = self.WebBrowserMsgQueue )
    self.msgqueuethrd.daemon = True
    self.msgqueuethrd.start()


  def SendWebSockMsg(self, msgtype, msg=""):
    #import code, traceback; code.interact(local=locals(), banner="".join( traceback.format_stack(limit=10) ) )
    #print "self.server.clients: ", self.server.clients
    #print "self.websockclient: ",
    message = u""
    message = msgtype + self.msgdelim + msg
    if self.websockclient:
      while "Refreshing" in self.lastmsg:
        sleep(0.5)
        self.lastmsg = ""
      self.server.send_message(self.websockclient, message )
    else:
      self.OpenBrowser()


  def SetOpacity(self, bin, alpha):
    if bin > self.nbin:
      self.mprint( "There are only %d bins of data present" %self.nbin )
      return
    msg = u"alpha, %d, %f" %(bin, alpha)
    self.SendWebSockMsg(msg)


  def RedrawNGL(self):
    self.SendWebSockMsg("Redraw")


  def ReloadNGL(self): # expensive as javascript may be several Mbytes large
    self.mprint("Rendering JavaScript...")
    self.SendWebSockMsg("Reload")


  def OpenBrowser(self):
    NGLlibpath = libtbx.env.under_root(os.path.join("modules","cctbx_project","crys3d","hklview","ngl.js") )
    htmlstr = self.hklhtml %(NGLlibpath, os.path.abspath( self.jscriptfname))
    htmlstr += self.htmldiv
    with open(self.hklfname, "w") as f:
      f.write( htmlstr )
    self.url = "file://" + os.path.abspath( self.hklfname )
    self.mprint( "Writing %s and connecting to its websocket client" %self.hklfname, self.verbose )
    if self.UseOSBrowser:
      webbrowser.open(self.url, new=1)
    self.isnewfile = False


  def ExpandInBrowser(self, P1=True, friedel_mate=True):
    sg = self.miller_array.space_group()
    symops = sg.all_ops()
    uc = self.miller_array.unit_cell()
    OrtMx = matrix.sqr( uc.orthogonalization_matrix())
    InvMx = OrtMx.inverse()

    msgtype = "Expand"
    msg = ""
    unique_rot_ops = []
    if P1:
      unique_rot_ops = symops[ 0 : sg.order_p() ]
    if friedel_mate:
      msgtype += "Friedel"

    for symop in unique_rot_ops:
      RotMx = matrix.sqr( symop.r().as_double())

      ortrot = (OrtMx * RotMx * InvMx).as_mat3()
      str_rot = str(ortrot)
      str_rot = str_rot.replace("(", "")
      str_rot = str_rot.replace(")", "")
      msg += str_rot + "\n"

    self.SendWebSockMsg(msgtype, msg)


  def TestNewFunction(self):
    self.SendWebSockMsg("Testing")





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
